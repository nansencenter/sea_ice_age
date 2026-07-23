from __future__ import annotations

import numpy as np
from matplotlib.tri import Triangulation
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from .io import velocity_km_day
from .vectors import east_north_to_projected


def signed_area(x: np.ndarray, y: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    points = np.column_stack((x, y))[triangles]
    return 0.5 * (
        (points[:, 1, 0] - points[:, 0, 0]) *
        (points[:, 2, 1] - points[:, 0, 1]) -
        (points[:, 2, 0] - points[:, 0, 0]) *
        (points[:, 1, 1] - points[:, 0, 1])
    )


def orient_counterclockwise(
    x: np.ndarray, y: np.ndarray, triangles: np.ndarray
) -> np.ndarray:
    triangles = np.asarray(triangles, dtype=np.int32).copy()
    reverse = signed_area(x, y, triangles) < 0
    triangles[reverse] = triangles[reverse][:, [0, 2, 1]]
    return triangles


def boundary_nodes(triangles: np.ndarray, layers: int = 0) -> np.ndarray:
    """Return boundary nodes, optionally expanded by adjacent node layers."""
    edges = np.sort(np.vstack((
        triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]
    )), axis=1)
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    selected = set(np.unique(unique_edges[counts == 1]).tolist())
    for _ in range(layers):
        selected_array = np.fromiter(selected, dtype=np.int32)
        touching = np.isin(unique_edges, selected_array).any(axis=1)
        selected.update(np.unique(unique_edges[touching]).tolist())
    return np.asarray(sorted(selected), dtype=np.int32)


def projected_velocity_grid(
    longitude: np.ndarray,
    latitude: np.ndarray,
    eastward: np.ndarray,
    northward: np.ndarray,
    target_crs: str,
) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray]:
    u_xy, v_xy = east_north_to_projected(
        longitude, latitude, eastward, northward, target_crs
    )
    return velocity_km_day(u_xy), velocity_km_day(v_xy)


def sample_velocity_nodes(
    x_grid_km: np.ndarray,
    y_grid_km: np.ndarray,
    u_grid_km_day: np.ndarray,
    v_grid_km_day: np.ndarray,
    node_x_km: np.ndarray,
    node_y_km: np.ndarray,
    interpolation_mode: str = "legacy_zero_fill",
    missing_policy: str = "zero",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample drift on nodes with selectable model/legacy missing-data rules.

    ``model_linear`` renormalizes bilinear weights over finite source ocean
    cells.  It therefore does not blend model drift with artificial zeroes at
    land. ``legacy_zero_fill`` preserves the earlier implementation exactly.
    The strict ``error`` policy leaves missing values as NaN so the movement
    stage can stop after exempting deliberately fixed boundary nodes.
    """
    if interpolation_mode not in {"model_linear", "legacy_zero_fill"}:
        raise ValueError(f"Unsupported drift interpolation mode: {interpolation_mode}")
    if missing_policy not in {"error", "error_on_ice", "nearest_on_ice", "zero"}:
        raise ValueError(f"Unsupported missing drift policy: {missing_policy}")
    points = np.column_stack((node_y_km, node_x_km))
    u_masked = np.ma.asarray(u_grid_km_day, dtype=np.float64)
    v_masked = np.ma.asarray(v_grid_km_day, dtype=np.float64)
    u_values = np.ma.filled(u_masked, np.nan)
    v_values = np.ma.filled(v_masked, np.nan)
    valid_grid = ~(
        np.ma.getmaskarray(u_masked) |
        np.ma.getmaskarray(v_masked) |
        ~np.isfinite(u_values) |
        ~np.isfinite(v_values)
    )
    u_filled = np.where(valid_grid, u_values, 0.0)
    v_filled = np.where(valid_grid, v_values, 0.0)
    u = RegularGridInterpolator(
        (y_grid_km, x_grid_km), u_filled,
        method="linear", bounds_error=False, fill_value=0.0,
    )(points)
    v = RegularGridInterpolator(
        (y_grid_km, x_grid_km), v_filled,
        method="linear", bounds_error=False, fill_value=0.0,
    )(points)
    coverage = RegularGridInterpolator(
        (y_grid_km, x_grid_km), valid_grid.astype(float),
        method="linear", bounds_error=False, fill_value=0.0,
    )(points)
    if interpolation_mode == "model_linear":
        good = coverage > np.finfo(np.float64).eps
        u[good] /= coverage[good]
        v[good] /= coverage[good]
        missing = ~good
    else:
        missing = coverage < 0.5
    if missing_policy == "zero":
        u[missing] = 0.0
        v[missing] = 0.0
    else:
        u[missing] = np.nan
        v[missing] = np.nan
    return u, v, missing


def fill_missing_velocity_nearest(
    x_grid_km: np.ndarray,
    y_grid_km: np.ndarray,
    u_grid_km_day: np.ndarray,
    v_grid_km_day: np.ndarray,
    node_x_km: np.ndarray,
    node_y_km: np.ndarray,
    u_nodes: np.ndarray,
    v_nodes: np.ndarray,
    fill_mask: np.ndarray,
    maximum_distance_km: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fill selected nodes from the nearest finite remapped drift cell.

    This is a bounded repair for small target-grid coverage holes, applied only
    after the caller has established that the node touches model sea ice.
    Unresolved nodes are returned for a mandatory error.
    """
    u = np.asarray(u_nodes, dtype=np.float64).copy()
    v = np.asarray(v_nodes, dtype=np.float64).copy()
    requested = np.asarray(fill_mask, dtype=bool)
    unresolved = requested.copy()
    if not requested.any():
        return u, v, unresolved, np.empty(0, dtype=np.float64)

    u_grid = np.ma.asarray(u_grid_km_day, dtype=np.float64)
    v_grid = np.ma.asarray(v_grid_km_day, dtype=np.float64)
    u_values = np.ma.filled(u_grid, np.nan)
    v_values = np.ma.filled(v_grid, np.nan)
    valid = ~(
        np.ma.getmaskarray(u_grid) |
        np.ma.getmaskarray(v_grid) |
        ~np.isfinite(u_values) |
        ~np.isfinite(v_values)
    )
    if not valid.any():
        return u, v, unresolved, np.empty(0, dtype=np.float64)

    grid_x, grid_y = np.meshgrid(x_grid_km, y_grid_km)
    tree = cKDTree(np.column_stack((grid_x[valid], grid_y[valid])))
    node_indices = np.flatnonzero(requested)
    distances, nearest = tree.query(
        np.column_stack((node_x_km[requested], node_y_km[requested]))
    )
    accepted = np.isfinite(distances) & (distances <= maximum_distance_km)
    accepted_nodes = node_indices[accepted]
    valid_u = u_values[valid]
    valid_v = v_values[valid]
    u[accepted_nodes] = valid_u[nearest[accepted]]
    v[accepted_nodes] = valid_v[nearest[accepted]]
    unresolved[accepted_nodes] = False
    return u, v, unresolved, distances[accepted]


def advect_nodes_carefully(
    x: np.ndarray,
    y: np.ndarray,
    u_km_day: np.ndarray,
    v_km_day: np.ndarray,
    triangles: np.ndarray,
    maximum_iterations: int = 100,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Apply one daily displacement while preventing inverted elements.

    The correction follows the reference routine: speeds at nodes belonging to
    inverted elements are repeatedly halved until all signed areas are positive.
    """
    u = np.asarray(u_km_day, dtype=np.float64).copy()
    v = np.asarray(v_km_day, dtype=np.float64).copy()
    affected: set[int] = set()
    for iteration in range(maximum_iterations + 1):
        x_new = x + u
        y_new = y + v
        inverted = signed_area(x_new, y_new, triangles) <= 0
        if not inverted.any():
            return x_new, y_new, iteration, len(affected)
        nodes = np.unique(triangles[inverted])
        affected.update(nodes.tolist())
        u[nodes] *= 0.5
        v[nodes] *= 0.5
    raise RuntimeError("Could not prevent inverted elements during daily advection")


def sample_mapping_points(
    x_grid_km: np.ndarray,
    y_grid_km: np.ndarray,
    domain_mask: np.ndarray,
    subdivisions: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Create a stable sub-cell point cloud for overlap-weight estimation."""
    dx = float(np.median(np.abs(np.diff(x_grid_km))))
    dy = float(np.median(np.abs(np.diff(y_grid_km))))
    offsets_x = (np.arange(subdivisions) + 0.5) * dx / subdivisions - dx / 2
    offsets_y = (np.arange(subdivisions) + 0.5) * dy / subdivisions - dy / 2
    x_centres, y_centres = np.meshgrid(x_grid_km, y_grid_km)
    centres = np.column_stack((x_centres[domain_mask], y_centres[domain_mask]))
    offsets = np.array(np.meshgrid(offsets_x, offsets_y)).reshape(2, -1).T
    points = (centres[:, None, :] + offsets[None, :, :]).reshape(-1, 2)
    return points[:, 0], points[:, 1]


def compute_element_mapping(
    source_x: np.ndarray,
    source_y: np.ndarray,
    source_triangles: np.ndarray,
    destination_x: np.ndarray,
    destination_y: np.ndarray,
    destination_triangles: np.ndarray,
    sample_x: np.ndarray,
    sample_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate normalized source-to-destination overlap weights."""
    source_finder = Triangulation(
        source_x, source_y, source_triangles
    ).get_trifinder()
    destination_finder = Triangulation(
        destination_x, destination_y, destination_triangles
    ).get_trifinder()
    source_index = source_finder(sample_x, sample_y)
    destination_index = destination_finder(sample_x, sample_y)
    valid = (source_index >= 0) & (destination_index >= 0)
    pairs = np.column_stack((source_index[valid], destination_index[valid])).astype(np.int32)
    pairs, counts = np.unique(pairs, axis=0, return_counts=True)

    totals = np.zeros(len(destination_triangles), dtype=np.float64)
    np.add.at(totals, pairs[:, 1], counts)
    weights = counts.astype(np.float64) / totals[pairs[:, 1]]

    missing = np.where(totals == 0)[0]
    if len(missing):
        centroids = np.column_stack((destination_x, destination_y))[
            destination_triangles[missing]
        ].mean(axis=1)
        fallback_source = source_finder(centroids[:, 0], centroids[:, 1])
        valid_fallback = fallback_source >= 0
        fallback_pairs = np.column_stack(
            (fallback_source[valid_fallback], missing[valid_fallback])
        ).astype(np.int32)
        pairs = np.vstack((pairs, fallback_pairs))
        weights = np.concatenate((weights, np.ones(len(fallback_pairs))))
    return pairs, weights


def mapped_area_ratio(
    previous_x: np.ndarray,
    previous_y: np.ndarray,
    advected_x: np.ndarray,
    advected_y: np.ndarray,
    source_triangles: np.ndarray,
    destination_count: int,
    mapping: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    initial_area = np.abs(signed_area(previous_x, previous_y, source_triangles))
    advected_area = np.abs(signed_area(advected_x, advected_y, source_triangles))
    source_ratio = advected_area / initial_area
    result = np.ones(destination_count, dtype=np.float64)
    accumulated = np.zeros(destination_count, dtype=np.float64)
    np.add.at(accumulated, mapping[:, 1], source_ratio[mapping[:, 0]] * weights)
    covered = np.zeros(destination_count, dtype=np.float64)
    np.add.at(covered, mapping[:, 1], weights)
    good = covered > 0
    result[good] = accumulated[good] / covered[good]
    return result

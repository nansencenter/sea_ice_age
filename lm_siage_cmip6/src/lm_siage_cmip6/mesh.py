from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation
from scipy.spatial import Delaunay


@dataclass(frozen=True)
class MeshQuality:
    vertices: int
    triangles: int
    minimum_area_km2: float
    median_area_km2: float
    maximum_area_km2: float
    minimum_edge_km: float
    median_edge_km: float
    maximum_edge_km: float
    minimum_quality: float
    median_quality: float


def build_domain_mask(
    maximum_concentration: np.ndarray,
    land_mask: np.ndarray,
    latitude: np.ndarray,
    concentration_threshold: float,
    buffer_cells: int,
    minimum_latitude: float,
) -> np.ndarray:
    """Build an ocean domain around all ice observed in the input period."""
    ice = np.asarray(maximum_concentration) >= concentration_threshold
    buffered = binary_dilation(ice, iterations=buffer_cells)
    return buffered & ~np.asarray(land_mask, dtype=bool) & (latitude >= minimum_latitude)


def triangulate_cell_centres(
    x_km: np.ndarray,
    y_km: np.ndarray,
    domain_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Triangulate selected regular-grid cell centres without crossing gaps."""
    x_grid, y_grid = np.meshgrid(x_km, y_km)
    x = x_grid[domain_mask]
    y = y_grid[domain_mask]
    points = np.column_stack((x, y))
    triangles = Delaunay(points, qhull_options="QJ").simplices.astype(np.int32)

    dx = float(np.median(np.abs(np.diff(x_km))))
    dy = float(np.median(np.abs(np.diff(y_km))))
    maximum_edge = 1.1 * np.hypot(dx, dy)
    triangle_points = points[triangles]
    edges = np.roll(triangle_points, -1, axis=1) - triangle_points
    short_enough = np.linalg.norm(edges, axis=2).max(axis=1) <= maximum_edge
    twice_area = np.abs(
        (triangle_points[:, 1, 0] - triangle_points[:, 0, 0]) *
        (triangle_points[:, 2, 1] - triangle_points[:, 0, 1]) -
        (triangle_points[:, 2, 0] - triangle_points[:, 0, 0]) *
        (triangle_points[:, 1, 1] - triangle_points[:, 0, 1])
    )
    nondegenerate = twice_area > np.finfo(np.float64).eps

    centroids = triangle_points.mean(axis=1)
    ix = np.abs(x_km[None, :] - centroids[:, 0, None]).argmin(axis=1)
    iy = np.abs(y_km[None, :] - centroids[:, 1, None]).argmin(axis=1)
    inside = domain_mask[iy, ix]
    triangles = triangles[short_enough & nondegenerate & inside]

    used = np.unique(triangles)
    remap = np.full(len(points), -1, dtype=np.int32)
    remap[used] = np.arange(len(used), dtype=np.int32)
    return x[used], y[used], remap[triangles]


def mesh_quality(x: np.ndarray, y: np.ndarray, triangles: np.ndarray) -> MeshQuality:
    points = np.column_stack((x, y))[triangles]
    vectors = np.roll(points, -1, axis=1) - points
    edges = np.linalg.norm(vectors, axis=2)
    twice_area = np.abs(
        (points[:, 1, 0] - points[:, 0, 0]) *
        (points[:, 2, 1] - points[:, 0, 1]) -
        (points[:, 2, 0] - points[:, 0, 0]) *
        (points[:, 1, 1] - points[:, 0, 1])
    )
    area = 0.5 * twice_area
    quality = 4.0 * np.sqrt(3.0) * area / np.sum(edges * edges, axis=1)
    return MeshQuality(
        vertices=len(x),
        triangles=len(triangles),
        minimum_area_km2=float(area.min()),
        median_area_km2=float(np.median(area)),
        maximum_area_km2=float(area.max()),
        minimum_edge_km=float(edges.min()),
        median_edge_km=float(np.median(edges)),
        maximum_edge_km=float(edges.max()),
        minimum_quality=float(quality.min()),
        median_quality=float(np.median(quality)),
    )

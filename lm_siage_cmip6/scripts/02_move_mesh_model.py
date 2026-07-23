#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from bamg_remesher import BamgOpts, BamgRemesher
from lm_siage_cmip6.advection import (
    advect_nodes_carefully,
    boundary_nodes,
    compute_element_mapping,
    fill_missing_velocity_nearest,
    mapped_area_ratio,
    orient_counterclockwise,
    projected_velocity_grid,
    sample_mapping_points,
    sample_velocity_nodes,
)
from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.interpolation import sample_sic_to_elements
from lm_siage_cmip6.io import concentration_fraction
from lm_siage_cmip6.mesh import mesh_quality


def mesh_path(root: Path, date_key: str) -> Path:
    return root / date_key[:4] / f"mesh_{date_key}.npz"


def load_mesh(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path) as data:
        triangles_name = "t" if "t" in data else "triangles"
        return data["x"], data["y"], data[triangles_name]


def main() -> None:
    parser = argparse.ArgumentParser(description="Advect and remesh the CMCC ice mesh")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--start", help="First source date, YYYYMMDD")
    parser.add_argument("--stop", help="Last destination date, YYYYMMDD")
    parser.add_argument("--days", type=int, help="Limit the number of daily intervals")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    mesh_config = config["mesh"]
    catalog = InputCatalog.from_config(config)
    records = catalog.select(args.start, args.stop)
    if args.days is not None:
        records = records[: args.days + 1]
    if len(records) < 2:
        raise ValueError("At least two consecutive input days are required")

    initial_path = Path(mesh_config["initial_mesh_file"])
    if not initial_path.is_absolute():
        initial_path = config["repo_root"] / initial_path
    daily_root = Path(mesh_config["daily_mesh_dir"])
    if not daily_root.is_absolute():
        daily_root = config["repo_root"] / daily_root
    with np.load(initial_path) as initial:
        initial_x = initial["x"]
        initial_y = initial["y"]
        initial_t = initial["triangles"]
        x_grid_km = initial["x_grid_km"]
        y_grid_km = initial["y_grid_km"]
        domain_mask = initial["domain_mask"]
    initial_t = orient_counterclockwise(initial_x, initial_y, initial_t)

    first_path = mesh_path(daily_root, records[0].key)
    first_path.parent.mkdir(parents=True, exist_ok=True)
    if args.force or not first_path.exists():
        np.savez_compressed(first_path, x=initial_x, y=initial_y, t=initial_t)

    sample_x, sample_y = sample_mapping_points(
        x_grid_km, y_grid_km, domain_mask,
        int(mesh_config["mapping_subdivisions"]),
    )
    options = BamgOpts(
        hmin=float(mesh_config["bamg_hmin_km"]),
        hmax=float(mesh_config["bamg_hmax_km"]),
    )
    remesher = BamgRemesher()
    variables = config["variables"]
    drift_interpolation = mesh_config.get(
        "drift_interpolation_mode", "legacy_zero_fill"
    )
    missing_drift_policy = mesh_config.get("missing_drift_policy", "zero")

    for source_record, destination_record in zip(records[:-1], records[1:]):
        source_path = mesh_path(daily_root, source_record.key)
        destination_path = mesh_path(daily_root, destination_record.key)
        if destination_path.exists() and not args.force:
            print(f"Skipping existing {destination_path}")
            continue
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source mesh {source_path}")
        x, y, triangles = load_mesh(source_path)
        triangles = orient_counterclockwise(x, y, triangles)

        with Dataset(source_record.path) as dataset:
            longitude = dataset.variables[variables["longitude"]][:]
            latitude = dataset.variables[variables["latitude"]][:]
            eastward = dataset.variables[variables["eastward_velocity"]][source_record.index]
            northward = dataset.variables[variables["northward_velocity"]][source_record.index]
            concentration = concentration_fraction(
                dataset.variables[variables["concentration"]][source_record.index]
            ) * 100.0
            land = np.asarray(
                dataset.variables[variables["status_flag"]][source_record.index]
            ) == 1
        u_grid, v_grid = projected_velocity_grid(
            longitude, latitude, eastward, northward, config["projection"]["proj4"]
        )
        u, v, missing = sample_velocity_nodes(
            x_grid_km, y_grid_km, u_grid, v_grid, x, y,
            interpolation_mode=drift_interpolation,
            missing_policy=missing_drift_policy,
        )
        fixed_nodes = boundary_nodes(
            triangles, int(mesh_config["fixed_boundary_layers"])
        )
        u[fixed_nodes] = 0.0
        v[fixed_nodes] = 0.0
        active_missing = missing.copy()
        active_missing[fixed_nodes] = False
        blocking_missing = active_missing.copy()
        repaired_distances = np.empty(0, dtype=np.float64)
        if missing_drift_policy in {"error_on_ice", "nearest_on_ice"}:
            element_points = np.column_stack((
                y[triangles].mean(axis=1), x[triangles].mean(axis=1)
            ))
            element_sic, _ = sample_sic_to_elements(
                x_grid_km, y_grid_km, concentration, land, element_points,
                mode="model_linear",
            )
            ice_elements = np.isfinite(element_sic) & (
                element_sic >= 100.0 * float(config["ice"]["concentration_threshold"])
            )
            ice_nodes = np.zeros(len(x), dtype=bool)
            if ice_elements.any():
                ice_nodes[np.unique(triangles[ice_elements])] = True
            blocking_missing &= ice_nodes
            inactive_missing = active_missing & ~ice_nodes
            u[inactive_missing] = 0.0
            v[inactive_missing] = 0.0
        if missing_drift_policy == "nearest_on_ice" and blocking_missing.any():
            u, v, blocking_missing, repaired_distances = fill_missing_velocity_nearest(
                x_grid_km, y_grid_km, u_grid, v_grid, x, y, u, v,
                blocking_missing,
                float(mesh_config.get("nearest_drift_max_distance_km", 50.0)),
            )
        if missing_drift_policy in {"error", "error_on_ice", "nearest_on_ice"} and blocking_missing.any():
            raise RuntimeError(
                f"{source_record.key}: drift is missing at "
                f"{int(blocking_missing.sum())} protected mesh nodes. "
                "The configured model-data policy will not replace it with zero motion."
            )
        advected_x, advected_y, correction_iterations, affected_nodes = (
            advect_nodes_carefully(
                x, y, u, v, triangles,
                int(mesh_config["maximum_advection_iterations"]),
            )
        )
        destination_x, destination_y, destination_t = remesher.remesh(
            advected_x, advected_y, triangles, options
        )
        destination_t = orient_counterclockwise(
            destination_x, destination_y, destination_t
        )
        mapping, weights = compute_element_mapping(
            advected_x, advected_y, triangles,
            destination_x, destination_y, destination_t,
            sample_x, sample_y,
        )
        area_ratio = mapped_area_ratio(
            x, y, advected_x, advected_y, triangles,
            len(destination_t), mapping, weights,
        )
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            destination_path,
            x=destination_x,
            y=destination_y,
            t=destination_t,
            src2dst=mapping,
            weights=weights.astype(np.float32),
            ar=area_ratio.astype(np.float32),
        )
        quality = mesh_quality(destination_x, destination_y, destination_t)
        speed = np.hypot(u, v)
        print(
            f"{source_record.key}->{destination_record.key}: "
            f"nodes={len(destination_x)}, triangles={len(destination_t)}, "
            f"missing_drift={active_missing.mean():.1%}, "
            f"missing_protected={int(blocking_missing.sum())}, "
            f"nearest_drift_repairs={len(repaired_distances)}, "
            f"maximum_repair_distance={repaired_distances.max() if len(repaired_distances) else 0:.2f} km, "
            f"fixed_boundary_nodes={len(fixed_nodes)}, "
            f"median_speed={np.median(speed[~missing]) if (~missing).any() else 0:.3f} km/day, "
            f"correction_iterations={correction_iterations}, "
            f"affected_nodes={affected_nodes}, "
            f"area_ratio={area_ratio.min():.3f}..{area_ratio.max():.3f}, "
            f"minimum_quality={quality.minimum_quality:.3f}"
        )


if __name__ == "__main__":
    main()

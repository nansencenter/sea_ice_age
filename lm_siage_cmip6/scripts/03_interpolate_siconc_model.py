#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.interpolation import sample_sic_to_elements
from lm_siage_cmip6.io import concentration_fraction
from lm_siage_cmip6.pipeline import dated_path, save_npz


def main() -> None:
    parser = argparse.ArgumentParser(description="Interpolate CMCC SIC onto mesh elements")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--start")
    parser.add_argument("--stop")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    records = InputCatalog.from_config(config).select(args.start, args.stop)
    mesh_root = config["repo_root"] / config["mesh"]["daily_mesh_dir"]
    sic_root = config["repo_root"] / config["ice"]["daily_sic_dir"]
    variables = config["variables"]
    interpolation_mode = config["ice"].get(
        "sic_interpolation_mode", "legacy_zero_fill"
    )
    missing_policy = config["ice"].get("missing_sic_policy", "zero")
    if missing_policy not in {"error", "error_on_ice_domain", "zero"}:
        raise ValueError(f"Unsupported missing SIC policy: {missing_policy}")

    with Dataset(records[0].path) as dataset:
        x_grid = np.asarray(dataset.variables["x"][:]) / 1000.0
        y_grid = np.asarray(dataset.variables["y"][:]) / 1000.0
        land = np.asarray(dataset.variables[variables["status_flag"]][0]) == 1
    initial_mesh = config["repo_root"] / config["mesh"]["initial_mesh_file"]
    with np.load(initial_mesh) as initial:
        maximum_concentration = initial["maximum_concentration"] * 100.0
    potential_ice = RegularGridInterpolator(
        (y_grid, x_grid), maximum_concentration,
        method="linear", bounds_error=False, fill_value=0.0,
    )
    threshold_percent = 100.0 * float(config["ice"]["concentration_threshold"])
    for record in records:
        mesh_file = dated_path(mesh_root, "mesh", record.key)
        output = dated_path(sic_root, "sic", record.key)
        if output.exists() and not args.force:
            continue
        if not mesh_file.exists():
            raise FileNotFoundError(mesh_file)
        with np.load(mesh_file) as mesh:
            x, y, triangles = mesh["x"], mesh["y"], mesh["t"]
        centroids_x = x[triangles].mean(axis=1)
        centroids_y = y[triangles].mean(axis=1)
        points = np.column_stack((centroids_y, centroids_x))
        with Dataset(record.path) as dataset:
            concentration = concentration_fraction(
                dataset.variables[variables["concentration"]][record.index]
            )
        concentration_percent = concentration * 100.0
        sic, missing = sample_sic_to_elements(
            x_grid, y_grid, concentration_percent, land, points,
            mode=interpolation_mode,
        )
        blocking_missing = missing.copy()
        if missing_policy == "error_on_ice_domain":
            blocking_missing &= potential_ice(points) >= threshold_percent
            sic[missing & ~blocking_missing] = 0.0
        if missing_policy in {"error", "error_on_ice_domain"} and blocking_missing.any():
            raise RuntimeError(
                f"{record.key}: SIC is missing at {int(blocking_missing.sum())} "
                "protected ocean mesh elements. The configured model-data "
                "policy will not fill them."
            )
        if missing_policy == "zero":
            sic[missing] = 0.0
        sic = np.clip(sic, 0.0, 100.0).astype(np.float32)
        save_npz(output, {"sic": sic})
        print(
            f"{record.key}: elements={len(sic)}, "
            f"ice_elements={(sic >= 15).sum()}, missing={int(missing.sum())}, "
            f"missing_protected={int(blocking_missing.sum())}, "
            f"mode={interpolation_mode}, range={sic.min():.2f}..{sic.max():.2f}%"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from matplotlib.tri import Triangulation
import numpy as np
from netCDF4 import Dataset

from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.pipeline import dated_path, load_npz, save_npz


def main() -> None:
    parser = argparse.ArgumentParser(description="Grid element sea-ice age fields")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--start")
    parser.add_argument("--stop")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    mesh_root = config["repo_root"] / config["mesh"]["daily_mesh_dir"]
    age_root = config["repo_root"] / config["age"]["daily_age_dir"]
    grid_root = config["repo_root"] / config["age"]["gridded_age_dir"]
    variables = config["variables"]

    with Dataset(config["input_file"]) as dataset:
        x_km = np.asarray(dataset.variables["x"][:]) / 1000.0
        y_km = np.asarray(dataset.variables["y"][:]) / 1000.0
        land = np.asarray(dataset.variables[variables["status_flag"]][0]) == 1
    x_grid, y_grid = np.meshgrid(x_km, y_km)

    for age_file in sorted(age_root.glob("*/age_*.npz")):
        date_key = age_file.stem.split("_")[-1]
        if args.start and date_key < args.start:
            continue
        if args.stop and date_key > args.stop:
            continue
        age_data = load_npz(age_file)
        if "age" not in age_data or "f" not in age_data:
            continue
        output = dated_path(grid_root, "grid", date_key)
        if output.exists() and not args.force:
            continue
        mesh = load_npz(dated_path(mesh_root, "mesh", date_key))
        finder = Triangulation(mesh["x"], mesh["y"], mesh["t"]).get_trifinder()
        element_index = finder(x_grid, y_grid)
        invalid = (element_index < 0) | land
        safe_index = element_index.copy()
        safe_index[invalid] = 0

        output_data: dict[str, np.ndarray] = {}
        age_grid = age_data["age"][safe_index].astype(np.float32)
        age_grid[invalid] = np.nan
        output_data["age"] = age_grid
        fractions = age_data["f"]
        for index, fraction in enumerate(fractions):
            name = f"sic_{len(fractions) - index}yi"
            fraction_grid = fraction[safe_index].astype(np.float32)
            fraction_grid[invalid] = np.nan
            output_data[name] = fraction_grid
        save_npz(output, output_data)
        print(
            f"{date_key}: gridded age and {len(fractions)} fractions; "
            f"valid cells={np.isfinite(age_grid).sum()}"
        )


if __name__ == "__main__":
    main()

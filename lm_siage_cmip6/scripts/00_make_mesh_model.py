#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.io import concentration_fraction
from lm_siage_cmip6.mesh import build_domain_mask, mesh_quality, triangulate_cell_centres


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the initial model sea-ice mesh")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument(
        "--use-bamg", action="store_true",
        help="Remesh the initial triangulation with BAMG after constructing it",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    mesh_config = config["mesh"]
    output = Path(mesh_config["initial_mesh_file"])
    if not output.is_absolute():
        output = config["repo_root"] / output
    if output.exists() and not args.force:
        print(f"Keeping existing initial mesh: {output}")
        return

    catalog = InputCatalog.from_config(config)
    with Dataset(catalog.records[0].path) as ds:
        x_km = np.asarray(ds.variables["x"][:], dtype=np.float64) / 1000.0
        y_km = np.asarray(ds.variables["y"][:], dtype=np.float64) / 1000.0
        latitude = np.asarray(ds.variables[config["variables"]["latitude"]][:])
        status = ds.variables[config["variables"]["status_flag"]][0]
        land_mask = np.asarray(status) == 1

    maximum_concentration = np.zeros(land_mask.shape, dtype=np.float32)
    concentration_name = config["variables"]["concentration"]
    records_by_path = {}
    for record in catalog.records:
        records_by_path.setdefault(record.path, []).append(record.index)
    for path, indices in records_by_path.items():
        with Dataset(path) as ds:
            for index in indices:
                concentration = concentration_fraction(ds.variables[concentration_name][index])
                np.maximum(
                    maximum_concentration,
                    np.ma.filled(concentration, 0.0),
                    out=maximum_concentration,
                )

    domain_mask = build_domain_mask(
        maximum_concentration,
        land_mask,
        latitude,
        config["ice"]["concentration_threshold"],
        int(mesh_config["ice_domain_buffer_cells"]),
        float(mesh_config["minimum_latitude"]),
    )
    x, y, triangles = triangulate_cell_centres(x_km, y_km, domain_mask)

    if args.use_bamg:
        from bamg_remesher import BamgOpts, BamgRemesher

        options = BamgOpts(
            hmin=float(mesh_config["bamg_hmin_km"]),
            hmax=float(mesh_config["bamg_hmax_km"]),
        )
        x, y, triangles = BamgRemesher().remesh(x, y, triangles, options)

    quality = mesh_quality(x, y, triangles)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        x=x,
        y=y,
        triangles=triangles,
        x_grid_km=x_km,
        y_grid_km=y_km,
        domain_mask=domain_mask,
        land_mask=land_mask,
        maximum_concentration=maximum_concentration,
    )
    print(f"Initial mesh written: {output}")
    print(f"BAMG applied: {args.use_bamg}")
    print(f"domain cells: {int(domain_mask.sum())}")
    for key, value in asdict(quality).items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the non-interactive production pipeline.

Plotting and drift snapshots are intentionally excluded. Use
01_plot_drift_snapshot.py explicitly during testing.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.pipeline import (
    current_experiment_start,
    latest_complete_output,
    require_checkpoint_state,
)


STAGES = [
    ("remap", "00_remap_model_data.py"),
    ("mesh", "00_make_mesh_model.py"),
    ("move", "02_move_mesh_model.py"),
    ("sic", "03_interpolate_siconc_model.py"),
    ("minimum", "04_find_autumn_minima_model.py"),
    ("age", "05_compute_age_model.py"),
    ("grid", "11_grid_age_model.py"),
    ("netcdf", "12_netcdf_age_model.py"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CMIP6 sea-ice-age pipeline")
    parser.add_argument(
        "--config", default="config/noresm2_mm_historical_r1i1p1f1.yaml"
    )
    parser.add_argument("--from-stage", choices=[name for name, _ in STAGES])
    parser.add_argument("--to-stage", choices=[name for name, _ in STAGES], default="netcdf")
    parser.add_argument("--start", help="Optional YYYYMMDD range start")
    parser.add_argument("--stop", help="Optional YYYYMMDD range stop")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--no-auto-restart", action="store_true",
        help="Do not start from the latest complete annual output",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    script_dir = Path(__file__).resolve().parent
    names = [name for name, _ in STAGES]
    first_stage = args.from_stage or ("remap" if "preprocess" in config else "mesh")
    selected = STAGES[names.index(first_stage): names.index(args.to_stage) + 1]

    effective_start = args.start
    output_start = args.start
    parent_handoff = False
    catalog = None
    for name, filename in selected:
        if name != "remap" and catalog is None:
            catalog = InputCatalog.from_config(config)
            if effective_start is None and not args.no_auto_restart and not args.force:
                checkpoint = latest_complete_output(config, catalog)
                if checkpoint is not None:
                    require_checkpoint_state(config, checkpoint)
                    effective_start = checkpoint.date_key
                    parent_handoff = (
                        checkpoint.experiment != config["dataset"]["experiment"]
                    )
                    output_start = (
                        current_experiment_start(config, catalog)
                        if parent_handoff else effective_start
                    )
                    print(
                        f"Restarting from complete {checkpoint.year} output "
                        f"({checkpoint.path.name}); state date={effective_start}",
                        flush=True,
                    )
                    if effective_start >= catalog.records[-1].key:
                        print("All available input years already have complete outputs.")
                        return
        command = [sys.executable, str(script_dir / filename), "--config", str(config["config_path"])]
        if name == "remap":
            if args.start:
                command.extend(("--start-year", args.start[:4]))
            if args.stop:
                command.extend(("--stop-year", args.stop[:4]))
        stage_start = effective_start
        if parent_handoff and name in {"sic", "minimum", "grid", "netcdf"}:
            stage_start = output_start
        if name in {"move", "sic", "age", "grid", "netcdf"}:
            if stage_start:
                command.extend(("--start", stage_start))
            if args.stop:
                command.extend(("--stop", args.stop))
        if name == "minimum" and (stage_start or args.stop):
            years = sorted({
                record.date.year
                for record in (catalog or InputCatalog.from_config(config)).select(stage_start, args.stop)
            })
            if years:
                command.append("--year")
                command.extend(str(year) for year in years)
        if name == "mesh":
            command.append("--use-bamg")
        if args.force:
            command.append("--force")
        print(f"\n=== Stage {name}: {filename} ===", flush=True)
        subprocess.run(command, cwd=config["repo_root"], check=True)


if __name__ == "__main__":
    main()

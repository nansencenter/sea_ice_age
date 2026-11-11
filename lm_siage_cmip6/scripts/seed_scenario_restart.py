#!/usr/bin/env python3
"""Seed an emissions scenario from the final historical Lagrangian state."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil

import numpy as np
from netCDF4 import Dataset

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.pipeline import (
    configured_path,
    dated_path,
    latest_complete_output,
    require_checkpoint_state,
)


def atomic_copy(source: Path, destination: Path, force: bool) -> None:
    if destination.exists() and not force:
        print(f"Keeping existing restart state: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    shutil.copy2(source, temporary)
    os.replace(temporary, destination)
    print(f"Copied restart state: {source} -> {destination}")


def copy_time_record(
    source: Path,
    destination: Path,
    index: int,
    *,
    parent_experiment: str,
    state_date: str,
    force: bool,
) -> None:
    """Copy one remapped forcing record without altering its date."""
    if destination.exists() and not force:
        print(f"Keeping existing restart bridge: {destination}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".part.nc")
    if temporary.exists():
        temporary.unlink()
    with Dataset(source) as src, Dataset(temporary, "w", format="NETCDF4") as dst:
        for name, dimension in src.dimensions.items():
            if name == "time":
                dst.createDimension(name, None if dimension.isunlimited() else 1)
            else:
                dst.createDimension(name, None if dimension.isunlimited() else len(dimension))
        dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
        for name, variable in src.variables.items():
            options = {}
            if "_FillValue" in variable.ncattrs():
                options["fill_value"] = variable.getncattr("_FillValue")
            if variable.ndim >= 2:
                options.update(zlib=True, complevel=1)
            output = dst.createVariable(name, variable.datatype, variable.dimensions, **options)
            output.setncatts({
                attribute: variable.getncattr(attribute)
                for attribute in variable.ncattrs()
                if attribute != "_FillValue"
            })
            if "time" in variable.dimensions:
                axis = variable.dimensions.index("time")
                source_selection = [slice(None)] * variable.ndim
                destination_selection = [slice(None)] * variable.ndim
                source_selection[axis] = slice(index, index + 1)
                destination_selection[axis] = slice(0, 1)
                output[tuple(destination_selection)] = variable[tuple(source_selection)]
            else:
                output[...] = variable[...]
        dst.restart_bridge_record = np.int8(1)
        dst.parent_experiment_id = parent_experiment
        dst.parent_state_date = state_date
        dst.comment = (
            "Single final historical forcing record used only to advance the "
            "historical Lagrangian state into the first scenario day."
        )
    os.replace(temporary, destination)
    print(f"Wrote restart forcing bridge: {destination}")


def validate_state(mesh_path: Path, age_path: Path) -> None:
    with np.load(mesh_path) as mesh:
        if not {"x", "y", "t"}.issubset(mesh.files):
            raise ValueError(f"Historical mesh checkpoint is incomplete: {mesh_path}")
    with np.load(age_path) as age:
        required = {"age", "f"}
        if not required.issubset(age.files) or not any(
            name.startswith("sic") and len(name) == 11 for name in age.files
        ):
            raise ValueError(f"Historical age checkpoint is incomplete: {age_path}")


def validate_scientific_switches(parent: dict, scenario: dict) -> None:
    """Require one algorithm definition across the historical/scenario seam."""
    settings = [
        ("preprocess", "target_mask_mode"),
        ("preprocess", "drift_extent_mode"),
        ("preprocess", "remap_neighbours"),
        ("preprocess", "remap_power"),
        ("preprocess", "remap_max_distance_m"),
        ("preprocess", "minimum_valid_drift_weight"),
        ("preprocess", "fix_index_block"),
        ("ice", "sic_interpolation_mode"),
        ("ice", "missing_sic_policy"),
        ("ice", "autumn_search_start"),
        ("ice", "autumn_search_stop"),
        ("age", "clip_autumn_area_ratio"),
        ("age", "autumn_area_ratio_clip"),
        ("age", "maximum_age_classes"),
        ("age", "year_length_days"),
        ("mesh", "drift_interpolation_mode"),
        ("mesh", "missing_drift_policy"),
        ("mesh", "nearest_drift_max_distance_km"),
        ("mesh", "fixed_boundary_layers"),
    ]
    mismatches = []
    for section, key in settings:
        parent_value = parent.get(section, {}).get(key)
        scenario_value = scenario.get(section, {}).get(key)
        if parent_value != scenario_value:
            mismatches.append(
                f"  {section}.{key}: parent={parent_value!r}, "
                f"scenario={scenario_value!r}"
            )
    if parent.get("projection", {}).get("proj4") != scenario.get("projection", {}).get("proj4"):
        mismatches.append("  projection.proj4 differs")
    if mismatches:
        raise ValueError(
            "Historical and scenario scientific switches differ:\n"
            + "\n".join(mismatches)
            + "\nUse the same algorithm settings on both sides of the handoff."
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Copy the final historical mesh/cohort checkpoint and one forcing "
            "record into an emissions-scenario work tree"
        )
    )
    parser.add_argument("--config", required=True, help="Scenario configuration")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    scenario = load_config(args.config)
    restart = scenario.get("restart") or {}
    required_restart = {
        "parent_config", "parent_experiment", "parent_state_date", "bridge_input_file"
    }
    missing = sorted(required_restart - set(restart))
    if missing:
        raise ValueError(f"Scenario restart config is missing: {', '.join(missing)}")

    parent_config_path = configured_path(scenario, restart["parent_config"])
    parent = load_config(parent_config_path)
    scenario_dataset = scenario["dataset"]
    parent_dataset = parent["dataset"]
    for field in ("model", "variant", "grid_label"):
        if scenario_dataset.get(field) != parent_dataset.get(field):
            raise ValueError(
                f"Scenario {field}={scenario_dataset.get(field)!r} does not match "
                f"parent {field}={parent_dataset.get(field)!r}"
            )
    if restart["parent_experiment"] != parent_dataset["experiment"]:
        raise ValueError("restart.parent_experiment does not match the parent config")
    if scenario_dataset["experiment"] == parent_dataset["experiment"]:
        raise ValueError("Scenario and parent experiments must be different")
    validate_scientific_switches(parent, scenario)

    state_date = str(restart["parent_state_date"])
    if len(state_date) != 8 or not state_date.isdigit():
        raise ValueError("restart.parent_state_date must be quoted YYYYMMDD")

    parent_catalog = InputCatalog.from_config(parent)
    checkpoint = latest_complete_output(parent, parent_catalog)
    if checkpoint is None:
        raise FileNotFoundError("No complete annual output exists for the parent run")
    if checkpoint.date_key != state_date:
        raise ValueError(
            f"Configured parent state is {state_date}, but the latest complete "
            f"parent output ends on {checkpoint.date_key}: {checkpoint.path}"
        )
    require_checkpoint_state(parent, checkpoint)

    parent_mesh_root = configured_path(parent, parent["mesh"]["daily_mesh_dir"])
    parent_age_root = configured_path(parent, parent["age"]["daily_age_dir"])
    scenario_mesh_root = configured_path(scenario, scenario["mesh"]["daily_mesh_dir"])
    scenario_age_root = configured_path(scenario, scenario["age"]["daily_age_dir"])
    source_mesh = dated_path(parent_mesh_root, "mesh", state_date)
    source_age = dated_path(parent_age_root, "age", state_date)
    validate_state(source_mesh, source_age)

    atomic_copy(
        configured_path(parent, parent["mesh"]["initial_mesh_file"]),
        configured_path(scenario, scenario["mesh"]["initial_mesh_file"]),
        args.force,
    )
    atomic_copy(source_mesh, dated_path(scenario_mesh_root, "mesh", state_date), args.force)
    atomic_copy(source_age, dated_path(scenario_age_root, "age", state_date), args.force)

    try:
        source_record = next(
            record for record in parent_catalog.records if record.key == state_date
        )
    except StopIteration as error:
        raise ValueError(f"Parent input does not contain {state_date}") from error
    copy_time_record(
        source_record.path,
        configured_path(scenario, restart["bridge_input_file"]),
        source_record.index,
        parent_experiment=parent_dataset["experiment"],
        state_date=state_date,
        force=args.force,
    )
    print(
        f"Scenario {scenario_dataset['experiment']} is seeded from "
        f"{parent_dataset['experiment']} state {state_date}."
    )


if __name__ == "__main__":
    main()

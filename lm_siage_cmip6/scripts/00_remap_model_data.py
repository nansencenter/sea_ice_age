#!/usr/bin/env python3
"""Select CMIP6 native-grid inputs and remap them to annual LAEA files.

The scientific remapping is delegated unchanged to the established
perturbation/Jacobian remapper.  This script adds only selection, annual
slicing, restart-safe output handling, and the small verification subset.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import shlex
import subprocess
import sys
import tempfile

import numpy as np
from netCDF4 import Dataset

from lm_siage_cmip6.cmip6 import AnnualTriplet, discover_annual_triplets
from lm_siage_cmip6.config import load_config


def resolve_path(config: dict, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (config["repo_root"] / path).resolve()


def copy_time_slice(source: Path, destination: Path, start: int, stop: int) -> None:
    """Copy a contiguous time slice while retaining CMIP6 metadata."""
    with Dataset(source) as src, Dataset(destination, "w", format="NETCDF4") as dst:
        for name, dimension in src.dimensions.items():
            if name == "time":
                dst.createDimension(name, None if dimension.isunlimited() else stop - start)
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
                attr: variable.getncattr(attr)
                for attr in variable.ncattrs()
                if attr != "_FillValue"
            })
            if "time" in variable.dimensions:
                axis = variable.dimensions.index("time")
                source_selection = [slice(None)] * variable.ndim
                destination_selection = [slice(None)] * variable.ndim
                source_selection[axis] = slice(start, stop)
                destination_selection[axis] = slice(0, stop - start)
                output[tuple(destination_selection)] = variable[tuple(source_selection)]
            else:
                output[...] = variable[...]


def output_is_complete(
    path: Path,
    expected_count: int,
    require_drift_completion: bool = False,
) -> bool:
    if not path.is_file():
        return False
    try:
        with Dataset(path) as dataset:
            complete = (
                int(getattr(dataset, "remapping_complete", 0)) == 1
                and int(getattr(dataset, "static_status_repeated", 0)) == 1
                and len(dataset.dimensions["time"]) == expected_count
                and {"siu", "siv", "siconc"}.issubset(dataset.variables)
            )
            if require_drift_completion:
                unresolved = int(
                    getattr(
                        dataset,
                        "drift_completion_unresolved_cell_times",
                        -1,
                    )
                )
                complete = (
                    complete
                    and int(getattr(dataset, "drift_completion_enabled", 0)) == 1
                    and int(getattr(dataset, "drift_completion_complete", 0)) == 1
                    and unresolved >= 0
                    and (
                        unresolved == 0
                        or int(
                            getattr(dataset, "drift_exclusion_complete", 0)
                        ) == 1
                    )
                    and {"drift_supported", "siconc_for_age"}.issubset(
                        dataset.variables
                    )
                )
            return complete
    except (OSError, KeyError, ValueError):
        return False


def repeat_static_status(path: Path, variable_name: str = "status_flag") -> None:
    """Repeat the target-grid status mask across the annual time dimension.

    The established remapper copies the static target status into record zero
    of its time-dependent output variable. NetCDF fill values otherwise remain
    in later records even though the target land/invalid mask is constant.
    """
    with Dataset(path, "a") as dataset:
        if variable_name not in dataset.variables:
            dataset.static_status_repeated = np.int8(0)
            return
        variable = dataset.variables[variable_name]
        if "time" not in variable.dimensions:
            dataset.static_status_repeated = np.int8(1)
            return
        time_axis = variable.dimensions.index("time")
        if time_axis != 0:
            raise ValueError(
                f"Expected time to be the first {variable_name} dimension; "
                f"found {variable.dimensions}"
            )
        status = variable[0]
        if np.ma.count(status) == 0:
            raise ValueError(f"{variable_name} record zero contains no valid values")
        for index in range(1, len(dataset.dimensions["time"])):
            variable[index] = status
        variable.comment = (
            "Static target-grid status repeated for every output record by "
            "lm_siage_cmip6/scripts/00_remap_model_data.py."
        )
        dataset.static_status_repeated = np.int8(1)


def annotate_output(path: Path, selection: dict, triplet: AnnualTriplet) -> None:
    with Dataset(path, "a") as dataset:
        dataset.source_id = selection["model"]
        dataset.experiment_id = selection["experiment"]
        dataset.variant_label = selection["variant"]
        dataset.table_id = selection.get("table", "SIday")
        dataset.native_grid_label = selection.get("grid_label", "gn")
        dataset.input_start_date = triplet.date_keys[0]
        dataset.input_end_date = triplet.date_keys[-1]
        dataset.remapping_complete = np.int8(1)
        dataset.remapping_driver = "lm_siage_cmip6/scripts/00_remap_model_data.py"


def write_verification(
    source: Path,
    destination: Path,
    count: int,
    data_variables: set[str] | None = None,
) -> None:
    """Write a compact NetCDF containing selected first remapped records."""
    data_variables = data_variables or {
        "siu", "siv", "siconc", "siconc_for_age", "drift_supported",
    }
    keep = {
        "time", "time_bnds", "x", "y", "longitude", "latitude",
        "status_flag", "Lambert_Azimuthal_Equal_Area",
    } | data_variables
    temporary = destination.with_suffix(".part.nc")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if temporary.exists():
        temporary.unlink()
    with Dataset(source) as src, Dataset(temporary, "w", format="NETCDF4") as dst:
        actual = min(count, len(src.dimensions["time"]))
        used_dimensions = {
            dimension
            for name, variable in src.variables.items()
            if name in keep
            for dimension in variable.dimensions
        }
        for name in used_dimensions:
            dimension = src.dimensions[name]
            size = actual if name == "time" else len(dimension)
            dst.createDimension(name, size)
        dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
        for name, variable in src.variables.items():
            if name not in keep:
                continue
            options = {}
            if "_FillValue" in variable.ncattrs():
                options["fill_value"] = variable.getncattr("_FillValue")
            if variable.ndim >= 2:
                options.update(zlib=True, complevel=4)
            output = dst.createVariable(name, variable.datatype, variable.dimensions, **options)
            output.setncatts({
                attr: variable.getncattr(attr)
                for attr in variable.ncattrs()
                if attr != "_FillValue"
            })
            selection = tuple(
                slice(0, actual) if dim == "time" else slice(None)
                for dim in variable.dimensions
            )
            output[...] = variable[selection] if selection else variable[...]
        names = ", ".join(sorted(data_variables))
        dst.title = f"Verification subset of remapped CMIP6 {names} input"
        dst.purpose = (
            f"Manual inspection of the first remapped {names} records; not used "
            "by the production sea-ice-age calculation."
        )
        dst.verification_time_steps = np.int32(actual)
        dst.verification_complete = np.int8(1)
    os.replace(temporary, destination)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select and annually remap native CMIP6 sea-ice input to LAEA"
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--start-year", type=int)
    parser.add_argument("--stop-year", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    selection = config["dataset"]
    preprocess = config["preprocess"]
    source_dir = resolve_path(config, preprocess["source_dir"])
    output_dir = resolve_path(config, preprocess["remapped_dir"])
    target_grid = resolve_path(config, preprocess["target_grid_file"])
    remapper = resolve_path(config, preprocess["remapper_script"])
    weights = resolve_path(config, preprocess["weights_file"])
    verification = resolve_path(config, preprocess["verification_file"])
    sic_verification_value = preprocess.get("sic_verification_file")
    sic_verification = (
        resolve_path(config, sic_verification_value)
        if sic_verification_value else None
    )

    triplets = discover_annual_triplets(
        source_dir,
        model=selection["model"],
        experiment=selection["experiment"],
        variant=selection["variant"],
        table=selection.get("table", "SIday"),
        grid_label=selection.get("grid_label"),
        recursive=bool(preprocess.get("recursive", False)),
    )
    triplets = [
        item for item in triplets
        if (args.start_year is None or item.year >= args.start_year)
        and (args.stop_year is None or item.year <= args.stop_year)
    ]
    if not triplets:
        raise ValueError("No selected years remain after applying the requested range")

    print(
        f"Selected {selection['model']} {selection['experiment']} "
        f"{selection['variant']}: {triplets[0].year}..{triplets[-1].year}"
    )
    temp_root = resolve_path(config, preprocess.get("temporary_dir", "work/remapping/tmp"))
    complete_drift_on_ice = bool(
        preprocess.get("complete_drift_on_ice", True)
    )
    ice_threshold_percent = 100.0 * float(
        preprocess.get(
            "drift_completion_concentration_threshold",
            config["ice"]["concentration_threshold"],
        )
    )
    drift_completion_max_distance_km = float(
        preprocess.get(
            "drift_completion_max_distance_km",
            config.get("mesh", {}).get(
                "nearest_drift_max_distance_km",
                300.0,
            ),
        )
    )
    unresolved_ice_drift_policy = preprocess.get(
        "unresolved_ice_drift_policy",
        "error",
    )
    if unresolved_ice_drift_policy not in {"error", "exclude"}:
        raise ValueError(
            "preprocess.unresolved_ice_drift_policy must be error or exclude"
        )
    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        weights.parent.mkdir(parents=True, exist_ok=True)
        temp_root.mkdir(parents=True, exist_ok=True)

    first_output: Path | None = None
    for triplet in triplets:
        filename = (
            f"RemappedInput_{selection['model']}_{selection['experiment']}_"
            f"{selection['variant']}_{triplet.year}_LAEA.nc"
        )
        output = output_dir / filename
        first_output = first_output or output
        if output_is_complete(
            output,
            len(triplet.date_keys),
            require_drift_completion=complete_drift_on_ice,
        ) and not args.force:
            print(f"Skipping complete remapped year {triplet.year}: {output}")
            continue
        print(
            f"Year {triplet.year}: {triplet.date_keys[0]}..{triplet.date_keys[-1]} "
            f"({len(triplet.date_keys)} records) -> {output}"
        )
        if args.dry_run:
            continue
        with tempfile.TemporaryDirectory(
            prefix=f"remap_{triplet.year}_", dir=temp_root
        ) as temporary_directory:
            temporary_directory = Path(temporary_directory)
            annual_inputs = {}
            for variable, source_slice in triplet.sources.items():
                annual_path = temporary_directory / source_slice.path.name
                copy_time_slice(
                    source_slice.path,
                    annual_path,
                    source_slice.start_index,
                    source_slice.stop_index,
                )
                annual_inputs[variable] = annual_path

            partial = output.with_suffix(".part.nc")
            if partial.exists():
                partial.unlink()
            command = [
                sys.executable, str(remapper),
                "--siu", str(annual_inputs["siu"]),
                "--siv", str(annual_inputs["siv"]),
                "--siconc", str(annual_inputs["siconc"]),
                "--grid", str(target_grid),
                "--out", str(partial),
                "--weights", str(weights),
                "--time_chunk", str(preprocess.get("time_chunk", 50)),
                "--target_mask_mode", preprocess.get("target_mask_mode", "land_only"),
                "--drift_extent_mode", preprocess.get("drift_extent_mode", "nearest"),
                "--k", str(preprocess.get("remap_neighbours", 4)),
                "--power", str(preprocess.get("remap_power", 1.0)),
                "--max_dist", str(preprocess.get("remap_max_distance_m", 500000.0)),
                "--min_valid_weight", str(preprocess.get("minimum_valid_drift_weight", 0.5)),
            ]
            if preprocess.get("fix_index_block", True):
                command.extend([
                    "--fix_index_block",
                    "--index_block_start", str(preprocess.get("index_block_start", 214)),
                    "--index_block_end", str(preprocess.get("index_block_end", 218)),
                    "--index_block_pad", str(preprocess.get("index_block_pad", 3)),
                    "--block_method", preprocess.get("block_method", "idw"),
                    "--block_idw_power", str(preprocess.get("block_idw_power", 2.0)),
                ])
            if preprocess.get("auto_clip_siconc", True):
                command.append("--auto_clip_siconc")
            if complete_drift_on_ice:
                command.extend([
                    "--complete_drift_on_ice",
                    "--ice_threshold_percent", str(ice_threshold_percent),
                    "--drift_fill_max_distance",
                    str(1000.0 * drift_completion_max_distance_km),
                    "--drift_fill_donor_domain",
                    preprocess.get(
                        "drift_completion_donor_domain",
                        "non_land",
                    ),
                    "--unresolved_ice_drift_policy",
                    unresolved_ice_drift_policy,
                ])
            print(shlex.join(command), flush=True)
            subprocess.run(command, cwd=remapper.parent, check=True)
            repeat_static_status(partial, preprocess.get("target_mask_var", "status_flag"))
            annotate_output(partial, selection, triplet)
            os.replace(partial, output)
            print(f"Completed remapped year {triplet.year}: {output}")

    if not args.dry_run and first_output is not None:
        if args.force or not verification.is_file():
            write_verification(
                first_output, verification, int(preprocess.get("verification_timesteps", 10))
            )
            print(f"Wrote remapping verification data: {verification}")
        else:
            print(f"Keeping existing verification data: {verification}")
        if sic_verification is not None:
            if args.force or not sic_verification.is_file():
                write_verification(
                    first_output,
                    sic_verification,
                    int(preprocess.get("verification_timesteps", 10)),
                    data_variables={"siconc"},
                )
                print(f"Wrote SIC verification data: {sic_verification}")
            else:
                print(f"Keeping existing SIC verification data: {sic_verification}")


if __name__ == "__main__":
    main()

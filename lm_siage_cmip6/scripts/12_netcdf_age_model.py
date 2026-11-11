#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.pipeline import (
    annual_age_filename,
    configured_path,
    load_npz,
)


def set_attributes(variable, attributes: dict) -> None:
    variable.setncatts(attributes)


def output_is_complete(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        with Dataset(path) as dataset:
            return (
                int(getattr(dataset, "processing_complete", 0)) == 1
                and len(dataset.dimensions["time"]) > 0
                and {"sea_ice_age", "status_flag"}.issubset(dataset.variables)
            )
    except (OSError, KeyError, TypeError, ValueError):
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Export annual model sea-ice-age NetCDF files")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--start")
    parser.add_argument("--stop")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    grid_root = configured_path(config, config["age"]["gridded_age_dir"])
    output_config = config["output"]
    catalog = InputCatalog.from_config(config)
    record_by_key = {record.key: record for record in catalog.records}
    variables = config["variables"]
    selection = config.get("dataset", {})
    clip_enabled = bool(config["age"].get("clip_autumn_area_ratio", True))
    clip_values = config["age"].get("autumn_area_ratio_clip")
    clip_description = (
        f"enabled: {clip_values}" if clip_enabled and clip_values is not None
        else "disabled"
    )

    grid_by_key = {
        path.stem.split("_")[-1]: path
        for path in sorted(grid_root.glob("*/grid_*.npz"))
    }
    if not grid_by_key:
        raise FileNotFoundError(f"No gridded age fields below {grid_root}")
    first_age_key = min(grid_by_key)
    requested_keys = [
        key for key in grid_by_key
        if (args.start is None or key >= args.start)
        and (args.stop is None or key <= args.stop)
    ]
    years = sorted({int(key[:4]) for key in requested_keys})
    if not years:
        raise ValueError("No gridded age years lie in the requested range")

    template_path = catalog.records[0].path
    with Dataset(template_path) as template:
        x = np.asarray(template.variables["x"][:])
        y = np.asarray(template.variables["y"][:])
        longitude = np.asarray(template.variables[variables["longitude"]][:])
        latitude = np.asarray(template.variables[variables["latitude"]][:])
        source_status = np.asarray(template.variables[variables["status_flag"]][0])
        source_global = {name: template.getncattr(name) for name in template.ncattrs()}

    model = selection.get("model", source_global.get("source_id", "unknown-model"))
    experiment = selection.get(
        "experiment", source_global.get("experiment_id", "unknown-experiment")
    )
    variant = selection.get(
        "variant", source_global.get("variant_label", "unknown-variant")
    )
    output_root = configured_path(config, config["output_dir"]) / model
    output_root.mkdir(parents=True, exist_ok=True)

    all_catalog_keys = [record.key for record in catalog.records]
    for year in years:
        expected_keys = [
            key for key in all_catalog_keys
            if int(key[:4]) == year and key >= first_age_key
        ]
        actual_keys = sorted(key for key in grid_by_key if int(key[:4]) == year)
        if actual_keys != expected_keys:
            missing = sorted(set(expected_keys) - set(actual_keys))
            raise RuntimeError(
                f"Cannot mark {year} complete: {len(missing)} gridded daily fields "
                f"are missing (first: {missing[:3]})"
            )

        output = output_root / annual_age_filename(config, year)
        if output_is_complete(output) and not args.force:
            print(f"Skipping complete annual output: {output}")
            continue
        partial = output.with_suffix(".part.nc")
        if partial.exists():
            partial.unlink()

        fraction_names: set[str] = set()
        for key in actual_keys:
            fraction_names.update(
                name for name in load_npz(grid_by_key[key])
                if name.startswith("sic_") and name.endswith("yi")
            )
        ordered_fractions = sorted(
            fraction_names, key=lambda name: int(name.split("_")[1].replace("yi", ""))
        )

        time_values = []
        time_bounds = []
        time_units = None
        calendar = None
        open_sources: dict[Path, Dataset] = {}
        try:
            for key in actual_keys:
                record = record_by_key[key]
                if record.path not in open_sources:
                    open_sources[record.path] = Dataset(record.path)
                source = open_sources[record.path]
                time = source.variables["time"]
                units = time.units
                file_calendar = getattr(time, "calendar", "standard")
                if time_units is None:
                    time_units, calendar = units, file_calendar
                elif units != time_units or file_calendar != calendar:
                    raise ValueError("Time units or calendar change within an annual output")
                value = float(time[record.index])
                time_values.append(value)
                if "time_bnds" in source.variables:
                    time_bounds.append(
                        np.asarray(source.variables["time_bnds"][record.index], dtype=float)
                    )
                else:
                    time_bounds.append(np.array([value - 0.5, value + 0.5]))
        finally:
            for source in open_sources.values():
                source.close()

        with Dataset(partial, "w", format="NETCDF4") as dataset:
            dataset.createDimension("time", len(actual_keys))
            dataset.createDimension("y", len(y))
            dataset.createDimension("x", len(x))
            dataset.createDimension("nv", 2)

            projection = dataset.createVariable("Lambert_Azimuthal_Equal_Area", "i1")
            set_attributes(projection, {
                "grid_mapping_name": "lambert_azimuthal_equal_area",
                "longitude_of_projection_origin": 0.0,
                "latitude_of_projection_origin": 90.0,
                "false_easting": 0.0,
                "false_northing": 0.0,
                "semi_major_axis": 6378137.0,
                "inverse_flattening": 298.257223563,
                "proj4_string": config["projection"]["proj4"],
            })
            x_variable = dataset.createVariable("x", "f8", ("x",))
            y_variable = dataset.createVariable("y", "f8", ("y",))
            x_variable[:] = x
            y_variable[:] = y
            set_attributes(x_variable, {
                "standard_name": "projection_x_coordinate", "units": "m", "axis": "X"
            })
            set_attributes(y_variable, {
                "standard_name": "projection_y_coordinate", "units": "m", "axis": "Y"
            })
            time_variable = dataset.createVariable("time", "f8", ("time",))
            time_variable[:] = time_values
            set_attributes(time_variable, {
                "standard_name": "time", "long_name": "time", "units": time_units,
                "calendar": calendar, "axis": "T", "bounds": "time_bnds",
            })
            bounds_variable = dataset.createVariable("time_bnds", "f8", ("time", "nv"))
            bounds_variable[:] = np.asarray(time_bounds)

            age_variable = dataset.createVariable(
                output_config["sea_ice_age_name"], "f4", ("time", "y", "x"),
                zlib=True, complevel=4, fill_value=np.float32(1.0e20),
                chunksizes=(1, len(y), len(x)),
            )
            set_attributes(age_variable, {
                "standard_name": "age_of_sea_ice",
                "long_name": "Weighted Average of Sea Ice Age",
                "name": output_config["sea_ice_age_name"],
                "comment": "The weighted average is computed over all available fractions.",
                "units": "years", "valid_min": np.float32(0),
                "valid_max": np.float32(20), "ancillary_variables": "status_flag",
                "coverage_content_type": "modelResult",
                "grid_mapping": "Lambert_Azimuthal_Equal_Area",
            })
            status_variable = dataset.createVariable(
                "status_flag", "i1", ("time", "y", "x"), zlib=True,
                chunksizes=(1, len(y), len(x)),
            )
            set_attributes(status_variable, {
                "standard_name": "age_of_sea_ice status_flag",
                "long_name": "status flag array for sea ice age",
                "valid_min": np.int8(1), "valid_max": np.int8(3),
                "flag_masks": np.array([1, 2, 3], dtype=np.int8),
                "flag_meanings": "nominal land invalid",
                "flag_descriptions": "flag = 1: Nominal retrieval by the SIA algorithm; flag = 2: Position is over land; flag = 3: Pixel is invalid",
                "coverage_content_type": "qualityInformation",
                "grid_mapping": "Lambert_Azimuthal_Equal_Area",
            })

            if output_config["export_longitude_latitude"]:
                lon_variable = dataset.createVariable("longitude", "f4", ("y", "x"), zlib=True)
                lat_variable = dataset.createVariable("latitude", "f4", ("y", "x"), zlib=True)
                lon_variable[:] = longitude
                lat_variable[:] = latitude
                set_attributes(lon_variable, {
                    "standard_name": "longitude", "long_name": "longitude",
                    "units": "degrees_east",
                })
                set_attributes(lat_variable, {
                    "standard_name": "latitude", "long_name": "latitude",
                    "units": "degrees_north",
                })

            numerals = {
                1: "First", 2: "Second", 3: "Third", 4: "Fourth",
                5: "Fifth", 6: "Sixth", 7: "Seventh", 8: "Eighth",
            }
            fraction_variables = {}
            if output_config["export_age_fractions"]:
                for source_name in ordered_fractions:
                    age_class = int(source_name.split("_")[1].replace("yi", ""))
                    variable_name = f"conc_{age_class}yi"
                    variable = dataset.createVariable(
                        variable_name, "i2", ("time", "y", "x"), zlib=True,
                        complevel=4, fill_value=np.int16(-32768),
                        chunksizes=(1, len(y), len(x)),
                    )
                    variable.scale_factor = np.float32(0.001)
                    variable.add_offset = np.float32(0.0)
                    variable.set_auto_maskandscale(False)
                    set_attributes(variable, {
                        "long_name": f"Concentration of {numerals[age_class]} Year Sea Ice",
                        "name": variable_name, "units": "1",
                        "valid_min": np.int16(0), "valid_max": np.int16(1100),
                        "ancillary_variables": "status_flag",
                        "coverage_content_type": "modelResult",
                        "standard_name": "sea_ice_area_fraction",
                        "grid_mapping": "Lambert_Azimuthal_Equal_Area",
                    })
                    fraction_variables[source_name] = variable

            for index, key in enumerate(actual_keys):
                grid = load_npz(grid_by_key[key])
                age = grid["age"]
                status = np.ones(age.shape, dtype=np.int8)
                status[source_status == 1] = 2
                status[(source_status != 1) & ~np.isfinite(age)] = 3
                age_variable[index] = np.ma.masked_invalid(age)
                status_variable[index] = status
                for source_name, variable in fraction_variables.items():
                    packed = np.full(age.shape, -32768, dtype=np.int16)
                    if source_name in grid:
                        fraction = np.asarray(grid[source_name], dtype=float) / 100.0
                        valid = np.isfinite(fraction)
                        packed[valid] = np.rint(
                            np.clip(fraction[valid], 0.0, 1.0) * 1000.0
                        ).astype(np.int16)
                    variable[index] = packed

            created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            dataset.setncatts({
                "title": f"{model} Lagrangian Sea Ice Age",
                "summary": "Sea-ice age computed from CMIP6 sea-ice concentration and drift using Lagrangian advection of a triangular mesh.",
                "source": source_global.get("source", model),
                "source_id": model,
                "experiment_id": experiment,
                "variant_label": variant,
                "institution": source_global.get("institution", "unknown"),
                "Conventions": "CF-1.7, ACDD-1.3",
                "algorithm": "lagrangian_sea_ice_age_v2p1 adapted for CMIP6",
                "configuration_label": config.get("configuration_label", "default"),
                "remap_neighbours": int(
                    config["preprocess"].get("remap_neighbours", 4)
                ),
                "remap_power": float(
                    config["preprocess"].get("remap_power", 1.0)
                ),
                "remap_max_distance_m": float(
                    config["preprocess"].get("remap_max_distance_m", 500000.0)
                ),
                "drift_extent_mode": config["preprocess"].get(
                    "drift_extent_mode", "nearest"
                ),
                "sic_interpolation_mode": config["ice"].get(
                    "sic_interpolation_mode", "legacy_zero_fill"
                ),
                "missing_sic_policy": config["ice"].get(
                    "missing_sic_policy", "zero"
                ),
                "drift_interpolation_mode": config["mesh"].get(
                    "drift_interpolation_mode", "legacy_zero_fill"
                ),
                "missing_drift_policy": config["mesh"].get(
                    "missing_drift_policy", "zero"
                ),
                "nearest_drift_max_distance_km": float(
                    config["mesh"].get("nearest_drift_max_distance_km", 0.0)
                ),
                "autumn_minimum_search_window": (
                    f"{config['ice']['autumn_search_start']} to "
                    f"{config['ice']['autumn_search_stop']}"
                ),
                "autumn_area_ratio_clipping": clip_description,
                "history": f"Created {created} by lm_siage_cmip6",
                "date_created": created,
                "processing_level": "Level 4",
                "processing_complete": np.int8(1),
                "last_complete_date": actual_keys[-1],
                "time_coverage_start": actual_keys[0],
                "time_coverage_end": actual_keys[-1],
                "comment": "Uncertainty variables are omitted because the configured model input provides no concentration or drift uncertainties.",
            })
        os.replace(partial, output)
        print(
            f"Wrote complete annual output: {output} "
            f"({len(actual_keys)} daily age fields)"
        )


if __name__ == "__main__":
    main()

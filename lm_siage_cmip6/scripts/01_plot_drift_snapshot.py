#!/usr/bin/env python3
"""Create testing-only drift-vector snapshots.

This script is deliberately not part of the production pipeline. It must be
invoked explicitly and writes only to the ignored diagnostics directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset, num2date

from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.io import concentration_fraction, velocity_km_day
from lm_siage_cmip6.vectors import east_north_to_projected


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot testing-only model drift vectors in projected km/day"
    )
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument(
        "--time-index", type=int, nargs="+", default=[0],
        help="One or more zero-based time indices",
    )
    parser.add_argument("--stride", type=int, default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    variables = config["variables"]
    model = config.get("dataset", {}).get("model", "model")
    stride = args.stride or int(config["diagnostics"]["vector_stride"])
    output_dir = Path(config["diagnostics"]["drift_snapshot_dir"])
    if not output_dir.is_absolute():
        output_dir = config["repo_root"] / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    with Dataset(config["input_file"]) as ds:
        time = ds.variables["time"]
        dates = num2date(time[:], time.units, calendar=time.calendar)
        x_km = np.asarray(ds.variables["x"][:]) / 1000.0
        y_km = np.asarray(ds.variables["y"][:]) / 1000.0
        longitude = np.asarray(ds.variables[variables["longitude"]][:])
        latitude = np.asarray(ds.variables[variables["latitude"]][:])
        land = np.asarray(ds.variables[variables["status_flag"]][0]) == 1
        projection = config["projection"]["proj4"]

        for index in args.time_index:
            if index < 0 or index >= len(dates):
                raise IndexError(f"time index {index} is outside 0..{len(dates)-1}")
            concentration = concentration_fraction(
                ds.variables[variables["concentration"]][index]
            )
            eastward = ds.variables[variables["eastward_velocity"]][index]
            northward = ds.variables[variables["northward_velocity"]][index]
            u_xy, v_xy = east_north_to_projected(
                longitude, latitude, eastward, northward, projection
            )
            u_day = velocity_km_day(u_xy)
            v_day = velocity_km_day(v_xy)
            ice_present = np.ma.filled(
                concentration >= config["ice"]["concentration_threshold"], False
            )
            valid = ~(
                np.ma.getmaskarray(u_day) | np.ma.getmaskarray(v_day) | land
            ) & ice_present

            sample = np.s_[::stride, ::stride]
            x_sample, y_sample = np.meshgrid(x_km[::stride], y_km[::stride])
            u_sample = np.ma.array(u_day[sample], mask=~valid[sample])
            v_sample = np.ma.array(v_day[sample], mask=~valid[sample])

            figure, axis = plt.subplots(figsize=(9, 9), constrained_layout=True)
            image = axis.pcolormesh(
                x_km, y_km, np.ma.masked_where(land, concentration),
                shading="auto", cmap="Blues", vmin=0.0, vmax=1.0,
            )
            axis.contour(x_km, y_km, land.astype(float), levels=[0.5], colors="0.3", linewidths=0.4)
            arrows = axis.quiver(
                x_sample, y_sample, u_sample, v_sample,
                color="tab:red", scale=350.0, width=0.0022,
                headwidth=3.5, headlength=4.5,
            )
            axis.quiverkey(arrows, 0.82, 0.04, 10.0, "10 km/day", coordinates="axes")
            figure.colorbar(image, ax=axis, shrink=0.72, label="Sea-ice concentration")
            date_label = f"{dates[index].year:04d}-{dates[index].month:02d}-{dates[index].day:02d}"
            axis.set(
                title=f"{model} projected sea-ice drift — {date_label}",
                xlabel="LAEA x (km)", ylabel="LAEA y (km)", aspect="equal",
            )
            axis.text(
                0.01, 0.01,
                "Red arrows: east/north drift rotated to LAEA x/y and converted to km/day",
                transform=axis.transAxes, fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            )
            output = output_dir / f"drift_{dates[index].year:04d}{dates[index].month:02d}{dates[index].day:02d}.png"
            figure.savefig(output, dpi=160)
            plt.close(figure)
            speed = np.ma.array(
                np.ma.sqrt(u_day * u_day + v_day * v_day), mask=~valid
            )
            print(f"Wrote {output}")
            print(
                f"  valid vectors={speed.count()}, "
                f"median={float(np.ma.median(speed)):.3f} km/day, "
                f"p99={float(np.percentile(speed.compressed(), 99)):.3f} km/day, "
                f"max={float(speed.max()):.3f} km/day"
            )


if __name__ == "__main__":
    main()

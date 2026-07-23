#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a sea-ice-age NetCDF")
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    required = {
        "time", "time_bnds", "x", "y", "Lambert_Azimuthal_Equal_Area",
        "sea_ice_age", "status_flag",
    }
    with Dataset(args.path) as dataset:
        missing = sorted(required.difference(dataset.variables))
        if missing:
            raise ValueError(f"Missing variables: {', '.join(missing)}")
        age = dataset.variables["sea_ice_age"][:]
        if age.count() and (float(age.min()) < 0 or float(age.max()) > 20):
            raise ValueError("sea_ice_age lies outside 0..20 years")
        status = np.asarray(dataset.variables["status_flag"][:])
        if not set(np.unique(status)).issubset({1, 2, 3}):
            raise ValueError("status_flag contains values outside 1, 2, 3")
        fractions = sorted(
            name for name in dataset.variables
            if name.startswith("conc_") and name.endswith("yi")
        )
        fraction_sum = np.ma.zeros(age.shape, dtype=np.float64)
        for name in fractions:
            values = dataset.variables[name][:]
            if values.count() and (float(values.min()) < 0 or float(values.max()) > 1.001):
                raise ValueError(f"{name} lies outside 0..1")
            fraction_sum += values
        # Fractions are independently packed at a resolution of 0.001. Their
        # decoded sum may therefore exceed one by at most half a packing unit
        # per exported class even though the unpacked scientific sum is <= 1.
        packing_tolerance = 0.0005 * len(fractions) + 1.0e-6
        if fraction_sum.count() and float(fraction_sum.max()) > 1.0 + packing_tolerance:
            raise ValueError(
                "Decoded age fractions exceed one beyond packing tolerance: "
                f"{float(fraction_sum.max()):.6f} > "
                f"{1.0 + packing_tolerance:.6f}"
            )
        print(f"Output validation passed: {args.path}")
        print(f"variables: {', '.join(dataset.variables)}")
        print(f"valid age cells: {age.count()}")
        print(f"age range: {float(age.min()):.4f}..{float(age.max()):.4f} years")
        print(f"age fractions: {', '.join(fractions)}")


if __name__ == "__main__":
    main()

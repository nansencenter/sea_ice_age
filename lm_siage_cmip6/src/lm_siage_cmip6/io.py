from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date


EXPECTED_DIMS = ("time", "y", "x")


@dataclass(frozen=True)
class InputSummary:
    path: Path
    shape: tuple[int, int, int]
    first_date: str
    last_date: str
    calendar: str
    concentration_min: float
    concentration_max: float
    max_speed_m_s: float
    x_direction: str
    y_direction: str


def concentration_fraction(values: np.ndarray) -> np.ma.MaskedArray:
    """Convert valid CMIP concentration percentages to fractions."""
    values = np.ma.masked_invalid(np.ma.asarray(values, dtype=np.float64))
    result = values / 100.0
    return np.ma.masked_outside(result, 0.0, 1.0)


def velocity_km_day(values: np.ndarray) -> np.ma.MaskedArray:
    """Convert velocity from metres per second to kilometres per day."""
    return np.ma.masked_invalid(np.ma.asarray(values, dtype=np.float64)) * 86.4


def validate_input(path: str | Path) -> InputSummary:
    """Validate the structural and unit assumptions used by the CMCC workflow."""
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)

    required = {
        "time", "x", "y", "longitude", "latitude", "status_flag",
        "siu", "siv", "siconc", "Lambert_Azimuthal_Equal_Area",
    }
    with Dataset(path) as ds:
        missing = sorted(required.difference(ds.variables))
        if missing:
            raise ValueError(f"Missing required variables: {', '.join(missing)}")

        for name in ("siu", "siv", "siconc", "status_flag"):
            if ds.variables[name].dimensions != EXPECTED_DIMS:
                raise ValueError(
                    f"{name} dimensions are {ds.variables[name].dimensions}, "
                    f"expected {EXPECTED_DIMS}"
                )

        if ds.variables["siconc"].units != "%":
            raise ValueError("siconc must be stored as percent (%)")
        for name in ("siu", "siv"):
            if ds.variables[name].units != "m s-1":
                raise ValueError(f"{name} must use m s-1")
        for name in ("x", "y"):
            if ds.variables[name].units != "m":
                raise ValueError(f"{name} must use metres")

        time = ds.variables["time"]
        calendar = getattr(time, "calendar", "standard")
        dates = num2date(time[:], time.units, calendar=calendar)
        if len(dates) < 2:
            raise ValueError("At least two daily records are required")
        day_steps = np.diff(np.asarray(time[:], dtype=np.float64))
        if not np.allclose(day_steps, 1.0):
            raise ValueError("Input time steps must be consecutive daily records")

        x = np.asarray(ds.variables["x"][:])
        y = np.asarray(ds.variables["y"][:])
        def coordinate_direction(values: np.ndarray, name: str) -> str:
            differences = np.diff(values)
            if np.all(differences > 0):
                return "increasing"
            if np.all(differences < 0):
                return "decreasing"
            raise ValueError(f"{name} coordinate must be strictly monotonic")

        x_direction = coordinate_direction(x, "x")
        y_direction = coordinate_direction(y, "y")

        cmin, cmax = np.inf, -np.inf
        max_speed = 0.0
        for index in range(len(dates)):
            concentration = concentration_fraction(ds.variables["siconc"][index])
            if concentration.count():
                cmin = min(cmin, float(concentration.min()))
                cmax = max(cmax, float(concentration.max()))
            u = np.ma.filled(ds.variables["siu"][index], np.nan)
            v = np.ma.filled(ds.variables["siv"][index], np.nan)
            speed = np.ma.masked_invalid(np.hypot(u, v))
            if speed.count():
                max_speed = max(max_speed, float(speed.max()))

        if not np.isfinite(cmin) or not np.isfinite(cmax):
            raise ValueError("siconc contains no valid concentration values")

        return InputSummary(
            path=path,
            shape=(len(dates), len(y), len(x)),
            first_date=str(dates[0]),
            last_date=str(dates[-1]),
            calendar=calendar,
            concentration_min=cmin,
            concentration_max=cmax,
            max_speed_m_s=max_speed,
            x_direction=x_direction,
            y_direction=y_direction,
        )

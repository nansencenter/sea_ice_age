from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import glob

import cftime
import numpy as np
from netCDF4 import Dataset, num2date


@dataclass(frozen=True)
class TimeRecord:
    path: Path
    index: int
    date: object
    ordinal: float

    @property
    def key(self) -> str:
        return f"{self.date.year:04d}{self.date.month:02d}{self.date.day:02d}"


class InputCatalog:
    def __init__(self, records: list[TimeRecord], calendar: str):
        self.records = records
        self.calendar = calendar

    @classmethod
    def from_config(cls, config: dict) -> "InputCatalog":
        pattern = Path(config.get("input_pattern", config["input_file"]))
        if not pattern.is_absolute():
            pattern = config["repo_root"] / pattern
        paths = [Path(path).resolve() for path in sorted(glob.glob(str(pattern)))]
        if not paths:
            fallback = Path(config["input_file"]).resolve()
            if fallback.is_file():
                paths = [fallback]
        if not paths:
            raise FileNotFoundError(f"No input files match {pattern}")

        restart = config.get("restart") or {}
        bridge_value = restart.get("bridge_input_file")
        if bridge_value:
            bridge = Path(bridge_value).expanduser()
            if not bridge.is_absolute():
                bridge = config["repo_root"] / bridge
            bridge = bridge.resolve()
            if not bridge.is_file():
                raise FileNotFoundError(
                    f"Scenario restart bridge is missing: {bridge}\n"
                    "Run scripts/seed_scenario_restart.py for this config first."
                )
            if bridge not in paths:
                paths.append(bridge)

        records: list[TimeRecord] = []
        calendar: str | None = None
        for path in paths:
            with Dataset(path) as dataset:
                time = dataset.variables["time"]
                file_calendar = getattr(time, "calendar", "standard")
                if calendar is None:
                    calendar = file_calendar
                elif file_calendar != calendar:
                    raise ValueError(
                        f"Calendar changed from {calendar} to {file_calendar} in {path}"
                    )
                dates = num2date(time[:], time.units, calendar=file_calendar)
                ordinals = cftime.date2num(
                    dates, "days since 0001-01-01", calendar=file_calendar
                )
                records.extend(
                    TimeRecord(path, index, date, float(ordinal))
                    for index, (date, ordinal) in enumerate(zip(dates, ordinals))
                )

        records.sort(key=lambda record: record.ordinal)
        ordinals = np.asarray([record.ordinal for record in records])
        if len(np.unique(ordinals)) != len(ordinals):
            raise ValueError("Input series contains duplicate dates")
        if len(records) > 1 and not np.allclose(np.diff(ordinals), 1.0):
            gaps = np.where(~np.isclose(np.diff(ordinals), 1.0))[0]
            first = int(gaps[0])
            raise ValueError(
                f"Input series is not daily-contiguous between "
                f"{records[first].key} and {records[first + 1].key}"
            )
        return cls(records, calendar or "standard")

    def select(self, start: str | None = None, stop: str | None = None) -> list[TimeRecord]:
        return [
            record for record in self.records
            if (start is None or record.key >= start) and
            (stop is None or record.key <= stop)
        ]

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from pathlib import Path
import re

import cftime
import numpy as np
from netCDF4 import Dataset, num2date


REQUIRED_VARIABLES = ("siu", "siv", "siconc")
_DATE_RANGE = re.compile(r"^(?P<start>\d{8})-(?P<stop>\d{8})$")
_VARIANT = re.compile(r"^r\d+i\d+p\d+f\d+$")


@dataclass(frozen=True)
class CMIP6Name:
    path: Path
    variable: str
    table: str
    model: str
    experiment: str
    variant: str
    grid_label: str
    start: date
    stop: date


@dataclass(frozen=True)
class SourceSlice:
    path: Path
    start_index: int
    stop_index: int

    @property
    def count(self) -> int:
        return self.stop_index - self.start_index


@dataclass(frozen=True)
class AnnualTriplet:
    year: int
    date_keys: tuple[str, ...]
    sources: dict[str, SourceSlice]


def parse_cmip6_filename(path: str | Path) -> CMIP6Name:
    """Parse a CMIP6 filename from the right-hand, fixed fields.

    Parsing from the right is important because model names commonly contain
    hyphens and may contain additional underscores.
    """
    path = Path(path)
    if path.suffix != ".nc":
        raise ValueError(f"Not a NetCDF filename: {path.name}")
    parts = path.stem.split("_")
    if len(parts) < 7:
        raise ValueError(f"Not a CMIP6 daily filename: {path.name}")
    match = _DATE_RANGE.fullmatch(parts[-1])
    if match is None or _VARIANT.fullmatch(parts[-3]) is None:
        raise ValueError(f"Not a supported CMIP6 daily filename: {path.name}")
    model = "_".join(parts[2:-4])
    if not model:
        raise ValueError(f"Model field is empty in {path.name}")
    return CMIP6Name(
        path=path.resolve(),
        variable=parts[0],
        table=parts[1],
        model=model,
        experiment=parts[-4],
        variant=parts[-3],
        grid_label=parts[-2],
        start=_parse_date(match.group("start")),
        stop=_parse_date(match.group("stop")),
    )


def _parse_date(value: str) -> date:
    return date(int(value[:4]), int(value[4:6]), int(value[6:8]))


def discover_annual_triplets(
    source_dir: str | Path,
    *,
    model: str,
    experiment: str,
    variant: str,
    table: str = "SIday",
    grid_label: str | None = None,
    recursive: bool = False,
) -> list[AnnualTriplet]:
    """Select matching siu/siv/siconc files and split their records by year."""
    source_dir = Path(source_dir).expanduser().resolve()
    iterator = source_dir.rglob("*.nc") if recursive else source_dir.glob("*.nc")
    selected: dict[str, list[CMIP6Name]] = {name: [] for name in REQUIRED_VARIABLES}
    for path in sorted(iterator):
        try:
            item = parse_cmip6_filename(path)
        except ValueError:
            continue
        if (
            item.variable in selected
            and item.model == model
            and item.experiment == experiment
            and item.variant == variant
            and item.table == table
            and (grid_label is None or item.grid_label == grid_label)
        ):
            selected[item.variable].append(item)

    missing = [name for name, items in selected.items() if not items]
    if missing:
        raise FileNotFoundError(
            f"No matching {', '.join(missing)} input for "
            f"{model}/{experiment}/{variant} in {source_dir}"
        )

    records_by_variable = {
        variable: _read_source_records(items, variable, model, experiment, variant)
        for variable, items in selected.items()
    }
    reference_keys = set(records_by_variable["siconc"])
    for variable in ("siu", "siv"):
        keys = set(records_by_variable[variable])
        if keys != reference_keys:
            missing_dates = sorted(reference_keys - keys)
            extra_dates = sorted(keys - reference_keys)
            raise ValueError(
                f"{variable} dates do not match siconc; "
                f"missing={missing_dates[:3]}, extra={extra_dates[:3]}"
            )

    ordered_keys = sorted(reference_keys)
    reference_records = records_by_variable["siconc"]
    ordinals = np.asarray([reference_records[key][2] for key in ordered_keys])
    if len(ordinals) > 1 and not np.allclose(np.diff(ordinals), 1.0):
        index = int(np.where(~np.isclose(np.diff(ordinals), 1.0))[0][0])
        raise ValueError(
            f"Selected input is not daily-contiguous between "
            f"{ordered_keys[index]} and {ordered_keys[index + 1]}"
        )

    years = sorted({int(key[:4]) for key in ordered_keys})
    result: list[AnnualTriplet] = []
    for year in years:
        date_keys = tuple(key for key in ordered_keys if int(key[:4]) == year)
        sources: dict[str, SourceSlice] = {}
        for variable in REQUIRED_VARIABLES:
            entries = [records_by_variable[variable][key] for key in date_keys]
            paths = {entry[0] for entry in entries}
            if len(paths) != 1:
                raise ValueError(
                    f"{variable} records for {year} span multiple files; "
                    "prepare files so that a year is contained in one file"
                )
            indices = [entry[1] for entry in entries]
            if indices != list(range(indices[0], indices[0] + len(indices))):
                raise ValueError(f"{variable} indices for {year} are not contiguous")
            sources[variable] = SourceSlice(
                path=entries[0][0], start_index=indices[0], stop_index=indices[-1] + 1
            )
        result.append(AnnualTriplet(year=year, date_keys=date_keys, sources=sources))
    return result


def _read_source_records(
    items: list[CMIP6Name],
    variable: str,
    model: str,
    experiment: str,
    variant: str,
) -> dict[str, tuple[Path, int, float]]:
    result: dict[str, tuple[Path, int, float]] = {}
    calendar: str | None = None
    for item in items:
        with Dataset(item.path) as dataset:
            if variable not in dataset.variables:
                raise KeyError(f"{item.path} does not contain {variable}")
            metadata = {
                "model": getattr(dataset, "source_id", getattr(dataset, "model_id", None)),
                "experiment": getattr(dataset, "experiment_id", None),
                "variant": getattr(dataset, "variant_label", None),
            }
            expected = {"model": model, "experiment": experiment, "variant": variant}
            for field, value in expected.items():
                if metadata[field] is not None and metadata[field] != value:
                    raise ValueError(
                        f"{field} metadata in {item.path.name} is {metadata[field]!r}, "
                        f"expected {value!r}"
                    )
            time = dataset.variables["time"]
            file_calendar = getattr(time, "calendar", "standard")
            if calendar is None:
                calendar = file_calendar
            elif calendar != file_calendar:
                raise ValueError(f"Calendar changes from {calendar} to {file_calendar}")
            dates = num2date(time[:], time.units, calendar=file_calendar)
            ordinals = cftime.date2num(
                dates, "days since 0001-01-01", calendar=file_calendar
            )
            for index, (decoded, ordinal) in enumerate(zip(dates, ordinals)):
                key = f"{decoded.year:04d}{decoded.month:02d}{decoded.day:02d}"
                if key in result:
                    raise ValueError(f"Duplicate {variable} date {key}")
                result[key] = (item.path, index, float(ordinal))
    return result

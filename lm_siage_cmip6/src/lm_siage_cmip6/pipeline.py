from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset, num2date


def dated_path(root: Path, prefix: str, date_key: str) -> Path:
    return root / date_key[:4] / f"{prefix}_{date_key}.npz"


@dataclass(frozen=True)
class OutputCheckpoint:
    path: Path
    date_key: str
    year: int
    experiment: str


def configured_path(config: dict, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (config["repo_root"] / path).resolve()


def annual_age_filename(config: dict, year: int, experiment: str | None = None) -> str:
    dataset = config["dataset"]
    values = {
        "model": dataset["model"],
        "experiment": experiment or dataset["experiment"],
        "variant": dataset["variant"],
        "year": year,
    }
    template = config["output"].get(
        "filename_template", "SeaIceAge_{model}_{experiment}_{variant}_{year}.nc"
    )
    return template.format(**values)


def latest_complete_output(config: dict, catalog) -> OutputCheckpoint | None:
    """Find a safe restart output for the selected experiment.

    A current-experiment output always takes precedence.  A different
    experiment is eligible only when it is named explicitly as the parent in
    ``restart.parent_experiment``.  This prevents one emissions scenario from
    accidentally continuing from another scenario for the same model/variant.
    """
    if "dataset" not in config:
        return None
    dataset = config["dataset"]
    output_dir = configured_path(config, config["output_dir"]) / dataset["model"]
    if not output_dir.is_dir():
        return None
    pattern = re.compile(
        rf"^SeaIceAge_{re.escape(dataset['model'])}_(?P<experiment>.+)_"
        rf"{re.escape(dataset['variant'])}_(?P<year>\d{{4}})\.nc$"
    )
    catalog_keys = {record.key for record in catalog.records}
    candidates: list[OutputCheckpoint] = []
    for path in output_dir.glob("SeaIceAge_*.nc"):
        match = pattern.fullmatch(path.name)
        if match is None:
            continue
        try:
            with Dataset(path) as dataset_file:
                if int(getattr(dataset_file, "processing_complete", 0)) != 1:
                    continue
                if not {"time", "sea_ice_age", "status_flag"}.issubset(
                    dataset_file.variables
                ) or len(dataset_file.dimensions["time"]) == 0:
                    continue
                time = dataset_file.variables["time"]
                decoded = num2date(
                    time[-1], time.units, calendar=getattr(time, "calendar", "standard")
                )
                key = f"{decoded.year:04d}{decoded.month:02d}{decoded.day:02d}"
        except (OSError, KeyError, TypeError, ValueError):
            continue
        if key not in catalog_keys or int(match.group("year")) != int(key[:4]):
            continue
        candidates.append(OutputCheckpoint(
            path=path,
            date_key=key,
            year=int(match.group("year")),
            experiment=match.group("experiment"),
        ))
    current_experiment = dataset["experiment"]
    current = [
        item for item in candidates if item.experiment == current_experiment
    ]
    if current:
        return max(current, key=lambda item: item.date_key)

    parent_experiment = (config.get("restart") or {}).get("parent_experiment")
    if parent_experiment:
        parent = [
            item for item in candidates if item.experiment == parent_experiment
        ]
        if parent:
            return max(parent, key=lambda item: item.date_key)
    return None


def current_experiment_start(config: dict, catalog) -> str:
    """Return the first non-bridge input date in a scenario catalog."""
    bridge_value = (config.get("restart") or {}).get("bridge_input_file")
    if not bridge_value:
        return catalog.records[0].key
    bridge = configured_path(config, bridge_value)
    records = [record for record in catalog.records if record.path != bridge]
    if not records:
        raise ValueError("The input catalog contains only the restart bridge record")
    return records[0].key


def require_checkpoint_state(config: dict, checkpoint: OutputCheckpoint) -> None:
    """Ensure the Lagrangian state needed to continue is still available.

    The annual NetCDF contains gridded age, but not the moving mesh and cohort
    arrays.  Restarting scientifically therefore requires the matching internal
    checkpoint files as well as the public annual product.
    """
    mesh_root = configured_path(config, config["mesh"]["daily_mesh_dir"])
    age_root = configured_path(config, config["age"]["daily_age_dir"])
    required = [
        dated_path(mesh_root, "mesh", checkpoint.date_key),
        dated_path(age_root, "age", checkpoint.date_key),
    ]
    missing = [path for path in required if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(
            f"The latest complete output is {checkpoint.path}, but its restart "
            f"state is missing:\n{formatted}\n"
            "The gridded annual NetCDF alone cannot reconstruct mesh/cohort state."
        )


def load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {name: data[name] for name in data.files}


def save_npz(path: Path, data: dict[str, np.ndarray], merge: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    output = load_npz(path) if merge and path.exists() else {}
    output.update(data)
    temporary = path.with_suffix(path.suffix + ".tmp.npz")
    np.savez_compressed(temporary, **output)
    temporary.replace(path)


def propagate_element_field(
    source: np.ndarray,
    mapping: np.ndarray,
    weights: np.ndarray,
    area_ratio: np.ndarray,
    observed_limit: np.ndarray | None = None,
    area_ratio_clip: tuple[float, float] | None = None,
) -> np.ndarray:
    """Map an extensive element field and correct for element-area change.

    Dividing by ``area_ratio`` is required to express conserved transported
    ice as concentration after convergence or divergence. ``area_ratio_clip``
    is a separate optional numerical clamp inherited by the autumn-minimum
    stage; disabling the clamp must not disable the physical ratio correction.
    """
    ratio = np.asarray(area_ratio, dtype=np.float64)
    if area_ratio_clip is not None:
        ratio = np.clip(ratio, *area_ratio_clip)
    ratio = np.where(ratio > 0, ratio, 0.01)
    destination = np.zeros(len(ratio), dtype=np.float64)
    np.add.at(
        destination,
        mapping[:, 1],
        np.asarray(source)[mapping[:, 0]] * weights,
    )
    destination /= ratio
    destination = np.clip(destination, 0.0, 100.0)
    if observed_limit is not None:
        destination = np.minimum(destination, observed_limit)
    return destination


def enforce_cumulative_nesting(
    cumulative: list[np.ndarray],
    observed: np.ndarray,
) -> list[np.ndarray]:
    """Ensure older cumulative ice cohorts are subsets of newer cohorts.

    Age fractions are differences between cumulative cohort fields.  Those
    differences are physical only when oldest <= ... <= newest <= observed.
    Propagated older ice is therefore capped by each successively newer annual
    minimum, working backwards from the observed concentration ceiling.
    """
    ceiling = np.clip(np.asarray(observed, dtype=np.float64), 0.0, 100.0)
    nested: list[np.ndarray] = [np.empty(0)] * len(cumulative)
    for index in range(len(cumulative) - 1, -1, -1):
        values = np.clip(np.asarray(cumulative[index], dtype=np.float64), 0.0, ceiling)
        nested[index] = values
        ceiling = values
    return nested


def select_recent_cohort_names(
    data: dict[str, np.ndarray], maximum_classes: int
) -> list[str]:
    """Return at most the newest annual cumulative cohort names, oldest first.

    The original age notebook transports six annual cumulative cohorts.  On
    the autumn reference date the newly calculated minimum is stored beside
    that day's age result, then becomes one of the six transported cohorts on
    the following day.  Keeping this selection explicit prevents a transient
    seventh cumulative cohort from entering the age fractions at year end.
    """
    names = sorted(
        name for name in data if name.startswith("sic") and len(name) == 11
    )
    return names[-maximum_classes:]

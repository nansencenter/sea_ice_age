#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.pipeline import (
    dated_path,
    enforce_cumulative_nesting,
    load_npz,
    propagate_element_field,
    resolve_age_start,
    save_npz,
    select_recent_cohort_names,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Propagate ice cohorts and compute age")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--start")
    parser.add_argument("--stop")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    catalog = InputCatalog.from_config(config)
    mesh_root = config["repo_root"] / config["mesh"]["daily_mesh_dir"]
    sic_root = config["repo_root"] / config["ice"]["daily_sic_dir"]
    age_root = config["repo_root"] / config["age"]["daily_age_dir"]
    maximum_classes = int(config["age"]["maximum_age_classes"])
    year_days = float(config["age"]["year_length_days"])
    reference_month, reference_day = (
        int(value) for value in config["ice"]["autumn_search_stop"].split("-")
    )

    # Match only finalized YYYYMMDD checkpoints, never interrupted .tmp.npz files.
    initial_files = sorted(age_root.glob("*/age_????????.npz"))
    initial_files = [path for path in initial_files if any(
        name.startswith("sic") for name in load_npz(path)
    )]
    if not initial_files:
        raise FileNotFoundError("Run 04_find_autumn_minima_model.py first")
    initial_key = initial_files[0].stem.split("_")[-1]
    start = resolve_age_start(args.start, initial_key)
    if args.start is not None and start != args.start:
        print(
            f"Requested age start {args.start} precedes the first autumn "
            f"checkpoint; starting from {start} instead"
        )
    records = catalog.select(start, args.stop)
    if not records or records[0].key != start:
        raise ValueError(f"Start date {start} is not in the input catalog")

    for source_record, destination_record in zip(records[:-1], records[1:]):
        source_path = dated_path(age_root, "age", source_record.key)
        destination_path = dated_path(age_root, "age", destination_record.key)
        if destination_path.exists() and "age" in load_npz(destination_path) and not args.force:
            continue
        source_age = load_npz(source_path)
        cohort_names = select_recent_cohort_names(source_age, maximum_classes)
        mesh = load_npz(dated_path(mesh_root, "mesh", destination_record.key))
        observed = load_npz(dated_path(sic_root, "sic", destination_record.key))["sic"]
        destination_data = {}
        for name in cohort_names:
            destination_data[name] = propagate_element_field(
                source_age[name], mesh["src2dst"], mesh["weights"], mesh["ar"],
                observed_limit=observed,
            ).astype(np.float32)

        # Match the original notebook: compute today's age from the cohorts
        # transported into today.  A new autumn minimum stored on this date is
        # preserved below, but joins the transported set only tomorrow.  This
        # keeps the configured six cumulative cohorts and seven age fractions.
        cohort_names = sorted(destination_data)
        cumulative = enforce_cumulative_nesting(
            [destination_data[name] for name in cohort_names], observed
        )
        for name, values in zip(cohort_names, cumulative):
            destination_data[name] = values.astype(np.float32)
        fractions = np.diff(
            np.asarray([np.zeros_like(observed)] + cumulative + [observed]), axis=0
        )
        negative_count = int((fractions < -1.0e-3).sum())
        month_day = (destination_record.date.month, destination_record.date.day)
        reference_month_day = (reference_month, reference_day)
        autumn_year = (
            destination_record.date.year
            if month_day >= reference_month_day
            else destination_record.date.year - 1
        )
        autumn_type = destination_record.date.__class__
        autumn_date = autumn_type(autumn_year, reference_month, reference_day)
        offset = (destination_record.date - autumn_date).days / year_days
        class_ages = np.arange(len(fractions), 0, -1) - 1 + offset
        age = np.sum(fractions * class_ages[:, None], axis=0) / 100.0
        destination_data["age"] = age.astype(np.float32)
        destination_data["f"] = fractions.astype(np.float32)
        annual_name = f"sic{destination_record.key}"
        if destination_path.exists():
            existing = load_npz(destination_path)
            if annual_name in existing:
                destination_data[annual_name] = existing[annual_name]
        # Replace rather than merge so a forced rerun cannot retain stale
        # cohorts from an earlier configuration or interrupted calculation.
        save_npz(destination_path, destination_data, merge=False)
        print(
            f"{destination_record.key}: cohorts={len(cohort_names)}, "
            f"fractions={len(fractions)}, negative_fraction_values={negative_count}, "
            f"age={np.nanmin(age):.3f}..{np.nanmax(age):.3f} years"
        )


if __name__ == "__main__":
    main()

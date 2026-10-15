#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.pipeline import dated_path, load_npz, propagate_element_field, save_npz


def main() -> None:
    parser = argparse.ArgumentParser(description="Find annual autumn ice minima")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    parser.add_argument("--year", type=int, nargs="+")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    catalog = InputCatalog.from_config(config)
    mesh_root = config["repo_root"] / config["mesh"]["daily_mesh_dir"]
    sic_root = config["repo_root"] / config["ice"]["daily_sic_dir"]
    age_root = config["repo_root"] / config["age"]["daily_age_dir"]
    start_mmdd = config["ice"]["autumn_search_start"].replace("-", "")
    stop_mmdd = config["ice"]["autumn_search_stop"].replace("-", "")
    # Keep the area-ratio correction itself: it carries concentration through
    # physical mesh convergence/divergence. This switch controls only the
    # inherited autumn-only clamp that limited the ratio to roughly +/-2.4%
    # as an observational-noise safeguard. Model tests normally disable that
    # clamp so autumn and daily cohort propagation use the same deformation.
    clip_enabled = bool(config["age"].get("clip_autumn_area_ratio", True))
    clip_values = config["age"].get("autumn_area_ratio_clip")
    clip = tuple(clip_values) if clip_enabled and clip_values is not None else None
    if clip is not None and (len(clip) != 2 or clip[0] <= 0 or clip[0] > clip[1]):
        raise ValueError("autumn_area_ratio_clip must contain two ordered positive values")
    print(
        "Autumn area-ratio clipping: "
        + (f"enabled at {clip}" if clip is not None else "disabled")
    )
    years = args.year or sorted({record.date.year for record in catalog.records})

    by_key = {record.key: record for record in catalog.records}
    for year in years:
        start_key = f"{year:04d}{start_mmdd}"
        stop_key = f"{year:04d}{stop_mmdd}"
        candidates = [
            record for record in catalog.records
            if start_key <= record.key <= stop_key
        ]
        if not candidates or candidates[-1].key != stop_key:
            print(f"Skipping {year}: incomplete autumn search window")
            continue
        output = dated_path(age_root, "age", stop_key)
        cohort_name = f"sic{stop_key}"
        if output.exists() and cohort_name in load_npz(output) and not args.force:
            continue
        propagated_candidates = []
        for candidate in candidates:
            concentration = load_npz(dated_path(sic_root, "sic", candidate.key))["sic"]
            current_key = candidate.key
            while current_key < stop_key:
                current_index = next(
                    index for index, record in enumerate(catalog.records)
                    if record.key == current_key
                )
                destination = catalog.records[current_index + 1]
                mesh = load_npz(dated_path(mesh_root, "mesh", destination.key))
                observed = load_npz(dated_path(sic_root, "sic", destination.key))["sic"]
                concentration = propagate_element_field(
                    concentration, mesh["src2dst"], mesh["weights"], mesh["ar"],
                    observed_limit=observed, area_ratio_clip=clip,
                )
                current_key = destination.key
            propagated_candidates.append(concentration)
        minimum = np.min(np.asarray(propagated_candidates), axis=0).astype(np.float32)
        save_npz(output, {cohort_name: minimum}, merge=True)
        print(
            f"{year}: {len(candidates)} candidates, "
            f"{cohort_name} range={minimum.min():.2f}..{minimum.max():.2f}%"
        )


if __name__ == "__main__":
    main()

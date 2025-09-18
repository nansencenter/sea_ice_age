#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stitch two nextsim restart files with different meshes.
The first file is remapped to the mesh of the second file.
The output is saved in a npz file in mesh/YYYY/mesh_YYYYMMDD.npz

Usage:
    python stitch_restarts.py <input_file1> <input_file2> <output_dir> [--search-dist 15000] [--cores 1] [--force 0]

Created on Thu Jun 20 10:00:00 2024
"""

import argparse
import os
import datetime as dt
from joblib import Parallel, delayed

from stitch_restarts import(
        get_dst_file,
        stitch_restart_pair,
        DEFAULT_SEARCH_DIST,
        DEFAULT_CORES as DEFAULT_CORES_1PAIR,
        )


DEFAULT_CORES = 10

DEFAULT_PATTERN = "%Y%m%d/inputs/field_%Y%m%dT000000Z.bin"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process two input files and specify an output directory.")
    parser.add_argument('input_root_dir', type=str, help='Path to the root directory containing directories organised by dates with format yyyymmdd')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')

    valid_date = lambda s: dt.datetime.strptime(s, "%Y%m%d")
    parser.add_argument('start_date', type=valid_date, help="Don't process dates before this date (given in yyyymmdd format)")
    parser.add_argument('end_date', type=valid_date, help="Don't process dates after this date (given in yyyymmdd format)")
    default_pattern = DEFAULT_PATTERN.replace("%", "%%")
    parser.add_argument('--pattern', type=str,
            default=DEFAULT_PATTERN,
            help=f'pattern for finding input restart file (relative to input_root_dir) from date (default: {default_pattern})')
    parser.add_argument('--search-dist', type=float,
            default=DEFAULT_SEARCH_DIST,
            help=f'Search distance for compute_mapping (default: {DEFAULT_SEARCH_DIST})')
    parser.add_argument('--cores-one-pair', type=int,
            default=DEFAULT_CORES_1PAIR,
            help=f'Number of cores to use for stitching one pair of restarts (default: {DEFAULT_CORES_1PAIR})')
    parser.add_argument('--cores', type=int,
            default=DEFAULT_CORES,
            help=f'Number of pairs to process at once (default: {DEFAULT_CORES})')
    parser.add_argument('--force', action="store_true",
            help='Force overwrite of output files')
    return parser.parse_args()


def get_jobs(input_root_dir, start_date, end_date, pattern,
        output_dir, cores_one_pair, force, search_dist, **kwargs):

    job_opts_common = {
            "output_dir": output_dir,
            "cores": cores_one_pair,
            "force": force,
            "search_dist": search_dist,
            }

    def get_filename(date):
        rel_path = date.strftime(pattern)
        return os.path.join(input_root_dir, rel_path)

    date = start_date
    while date <= end_date:
        next_date = date + dt.timedelta(1)
        f1 = get_filename(date)
        f2 = get_filename(next_date)
        date = next_date
        output = get_dst_file(f2, output_dir)
        if os.path.exists(output) and not force:
            print(f"{output} exists: stopping. Use --force to overwrite.")
            continue
        if not (os.path.exists(f1) and os.path.exists(f2)):
            print(f"WARNING: input_files {f1} and {f2} not present.")
            continue
        yield job_opts_common | {
                "input_file1": f1,
                "input_file2": f2,
                }


def main():
    args = vars(parse_arguments())
    Parallel(n_jobs=args["cores"], verbose=10)(
            delayed(stitch_restart_pair)(**job_opts)
            for job_opts in get_jobs(**args)
            )


if __name__ == "__main__":
    main()

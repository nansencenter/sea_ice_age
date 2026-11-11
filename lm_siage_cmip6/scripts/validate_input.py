#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict

from lm_siage_cmip6.config import load_config
from lm_siage_cmip6.io import validate_input


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the CMCC age input")
    parser.add_argument("--config", default="config/cmcc_1950.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    summary = validate_input(config["input_file"])
    print("Model input validation passed")
    for key, value in asdict(summary).items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()

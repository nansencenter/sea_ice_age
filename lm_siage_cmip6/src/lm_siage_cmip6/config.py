from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a workflow configuration and resolve paths relative to the repo."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)

    if not isinstance(config, dict):
        raise ValueError(f"Configuration must contain a mapping: {config_path}")

    repo_root = config_path.parent.parent
    config["config_path"] = config_path
    config["repo_root"] = repo_root
    for key in ("input_file", "work_dir", "output_dir"):
        value = Path(config[key]).expanduser()
        config[key] = value if value.is_absolute() else (repo_root / value).resolve()
    return config

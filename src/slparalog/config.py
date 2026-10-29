"""Locate the project root and read pre-registered config files."""

from __future__ import annotations

from pathlib import Path

import yaml

# src/slparalog/config.py -> project root is three levels up (editable/src layout).
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = PROJECT_ROOT / "config"


def load_yaml(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_success(path: str | Path | None = None) -> dict:
    """Return the pre-registered success-criterion constants (config/success.yaml)."""
    return load_yaml(path or CONFIG_DIR / "success.yaml")

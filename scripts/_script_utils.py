#!/usr/bin/env python3
"""Shared helpers for script path and argument handling."""

from __future__ import annotations

from argparse import ArgumentTypeError
from pathlib import Path


def get_project_root() -> Path:
    """Return repository root (one level above scripts/)."""
    return Path(__file__).resolve().parents[1]


def resolve_path(root: Path, path: Path | None) -> Path | None:
    """Resolve path relative to repository root when not absolute."""
    if path is None:
        return None
    return path if path.is_absolute() else (root / path)


def resolve_data_dir(root: Path, data_dir: Path | None = None) -> Path:
    """
    Resolve input data directory.

    Priority:
    1) User provided path
    2) data/raw/202511 (legacy default)
    3) Latest subdirectory under data/raw/
    4) data/raw/
    """
    candidate = resolve_path(root, data_dir)
    if candidate is not None:
        return candidate

    raw_root = root / "data" / "raw"
    preferred = raw_root / "202511"
    if preferred.exists():
        return preferred

    if raw_root.exists():
        subdirs = sorted([p for p in raw_root.iterdir() if p.is_dir()])
        if subdirs:
            return subdirs[-1]

    return raw_root


def ensure_dir(path: Path) -> Path:
    """Create directory if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def positive_float(value: str) -> float:
    """Argparse type: strictly positive float."""
    parsed = float(value)
    if parsed <= 0:
        raise ArgumentTypeError("must be > 0")
    return parsed


def unit_interval(value: str) -> float:
    """Argparse type: float in [0, 1]."""
    parsed = float(value)
    if parsed < 0 or parsed > 1:
        raise ArgumentTypeError("must be between 0 and 1")
    return parsed


def positive_int(value: str) -> int:
    """Argparse type: strictly positive integer."""
    parsed = int(value)
    if parsed <= 0:
        raise ArgumentTypeError("must be a positive integer")
    return parsed

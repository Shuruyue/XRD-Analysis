"""Tests for scripts/_script_utils.py helpers."""

import sys
from argparse import ArgumentTypeError
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from _script_utils import (
    positive_float,
    positive_int,
    resolve_data_dir,
    resolve_path,
    unit_interval,
)


def test_resolve_path_relative(tmp_path):
    resolved = resolve_path(tmp_path, Path("outputs/plots"))
    assert resolved == tmp_path / "outputs/plots"


def test_resolve_data_dir_prefers_202511(tmp_path):
    (tmp_path / "data" / "raw" / "202511").mkdir(parents=True)
    (tmp_path / "data" / "raw" / "old").mkdir(parents=True)
    data_dir = resolve_data_dir(tmp_path, None)
    assert data_dir == tmp_path / "data" / "raw" / "202511"


def test_resolve_data_dir_uses_latest_subdir(tmp_path):
    (tmp_path / "data" / "raw" / "202401").mkdir(parents=True)
    (tmp_path / "data" / "raw" / "202512").mkdir(parents=True)
    data_dir = resolve_data_dir(tmp_path, None)
    assert data_dir == tmp_path / "data" / "raw" / "202512"


def test_resolve_data_dir_honors_user_path(tmp_path):
    explicit = tmp_path / "custom_raw"
    explicit.mkdir(parents=True)
    data_dir = resolve_data_dir(tmp_path, explicit)
    assert data_dir == explicit


def test_positive_float_validator():
    assert positive_float("0.1") == 0.1
    with pytest.raises(ArgumentTypeError):
        positive_float("0")


def test_unit_interval_validator():
    assert unit_interval("0") == 0.0
    assert unit_interval("1") == 1.0
    with pytest.raises(ArgumentTypeError):
        unit_interval("1.1")


def test_positive_int_validator():
    assert positive_int("1") == 1
    with pytest.raises(ArgumentTypeError):
        positive_int("0")

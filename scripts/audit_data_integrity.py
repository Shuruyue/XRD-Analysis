#!/usr/bin/env python3
"""
Data integrity audit for XRD workflows.

Checks:
1) Raw file presence and filename pattern consistency
2) Empty/corrupted file detection
3) Standard scan availability for instrument calibration
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


FILENAME_PATTERN = re.compile(
    r"^\d{8}_(\d+(\.\d+)?)ml_(\d+)h(_(\d+)min)?$"
)


def audit_raw(raw_dir: Path) -> dict:
    files = sorted(raw_dir.glob("*.txt"))
    bad_names = []
    empty_files = []

    for fp in files:
        stem = fp.stem
        if not FILENAME_PATTERN.match(stem):
            bad_names.append(fp.name)
        if fp.stat().st_size == 0:
            empty_files.append(fp.name)

    return {
        "raw_count": len(files),
        "bad_names": bad_names,
        "empty_files": empty_files,
    }


def audit_standards(standards_dir: Path) -> dict:
    if not standards_dir.exists():
        return {
            "exists": False,
            "standard_files": [],
            "has_calibration_standard": False,
        }

    files = [p.name for p in standards_dir.iterdir() if p.is_file() and p.name != ".gitkeep"]
    has_standard = any(
        key in name.lower()
        for name in files
        for key in ("lab6", "la_b6", "srm660", "si_standard", "silicon_standard")
    )
    return {
        "exists": True,
        "standard_files": files,
        "has_calibration_standard": has_standard,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit XRD data integrity and calibration readiness.")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/202511"))
    parser.add_argument("--standards-dir", type=Path, default=Path("data/standards"))
    args = parser.parse_args()

    raw_stats = audit_raw(args.raw_dir)
    std_stats = audit_standards(args.standards_dir)

    print("Data Integrity Audit")
    print(f"  Raw files: {raw_stats['raw_count']}")
    print(f"  Filename violations: {len(raw_stats['bad_names'])}")
    print(f"  Empty files: {len(raw_stats['empty_files'])}")
    if raw_stats["bad_names"]:
        print("  Bad filename samples:")
        for name in raw_stats["bad_names"][:10]:
            print(f"    - {name}")

    print(f"  Standards dir exists: {std_stats['exists']}")
    print(f"  Standard files found: {len(std_stats['standard_files'])}")
    print(f"  Calibration standard ready: {std_stats['has_calibration_standard']}")
    if std_stats["standard_files"]:
        for name in std_stats["standard_files"]:
            print(f"    - {name}")

    if raw_stats["empty_files"]:
        return 2
    if not std_stats["has_calibration_standard"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate per-sample fitting diagnosis plots."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _script_utils import ensure_dir, get_project_root, resolve_data_dir, resolve_path

from xrd_analysis.visualization.generate_fitting_diagnosis import (
    generate_sample_fitting_plot,
)


def generate_diagnosis_plots(data_dir: Path, output_dir: Path, clean: bool = True) -> int:
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    individual_dir = ensure_dir(output_dir / "individual")
    if clean:
        for file in individual_dir.glob("*.png"):
            file.unlink()

    data_files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(data_files)} data files in {data_dir}")
    if not data_files:
        return 1

    print(f"Generating {len(data_files)} diagnosis plots...")
    for idx, filepath in enumerate(data_files, 1):
        print(f"  {idx}/{len(data_files)}: {filepath.stem}")
        try:
            generate_sample_fitting_plot(filepath, output_dir=individual_dir)
        except Exception as exc:
            print(f"    [ERROR] {filepath.stem}: {exc}")

    print(f"Done! Plots saved to: {individual_dir}")
    return 0


def main() -> int:
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Generate fitting diagnosis plots.")
    parser.add_argument("--data-dir", type=Path, default=None, help="Directory containing raw .txt scans.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/plots/diagnostics/fitting"),
        help="Output diagnosis directory.",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Do not remove old png files before generation.",
    )
    args = parser.parse_args()

    data_dir = resolve_data_dir(root, args.data_dir)
    output_dir = resolve_path(root, args.output_dir)
    assert output_dir is not None
    ensure_dir(output_dir)
    return generate_diagnosis_plots(data_dir, output_dir, clean=not args.no_clean)


if __name__ == "__main__":
    raise SystemExit(main())

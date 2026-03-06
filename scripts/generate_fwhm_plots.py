#!/usr/bin/env python3
"""Generate FWHM evolution plots from pipeline results."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _script_utils import ensure_dir, get_project_root, resolve_data_dir, resolve_path

from xrd_analysis.analysis import batch_analyze
from xrd_analysis.visualization import plot_fwhm_by_concentration, plot_fwhm_evolution


def generate_fwhm_plots(data_dir: Path, output_dir: Path) -> int:
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    ensure_dir(output_dir)
    for file in output_dir.glob("*.png"):
        file.unlink()

    data_files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(data_files)} data files in {data_dir}")
    if not data_files:
        return 1

    results = batch_analyze([str(f) for f in data_files])
    plot_data = []
    for result in results:
        sample = {
            "name": result.sample_name,
            "concentration": result.leveler_concentration or 0.0,
            "time": result.sample_age_hours or 0.0,
            "peaks": [],
        }
        for peak in result.peaks:
            hkl_str = f"({peak.hkl[0]}{peak.hkl[1]}{peak.hkl[2]})"
            sample["peaks"].append(
                {
                    "hkl": hkl_str,
                    "fwhm": peak.fwhm,
                    "fwhm_error": 0.01,
                }
            )
        plot_data.append(sample)

    plot_fwhm_evolution(
        plot_data,
        x_param="time",
        output_path=str(output_dir / "fwhm_evolution.png"),
        show=False,
    )
    plot_fwhm_by_concentration(
        plot_data,
        output_path=str(output_dir / "fwhm_compare_grid.png"),
        dpi=300,
        show=False,
    )
    print(f"Done! Plots saved to: {output_dir}")
    return 0


def main() -> int:
    root = get_project_root()
    parser = argparse.ArgumentParser(
        description="Generate FWHM plots from batch scans."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing raw .txt scans.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/plots/analysis/fwhm"),
        help="Output plot directory.",
    )
    args = parser.parse_args()

    data_dir = resolve_data_dir(root, args.data_dir)
    output_dir = resolve_path(root, args.output_dir)
    assert output_dir is not None
    return generate_fwhm_plots(data_dir, output_dir)


if __name__ == "__main__":
    raise SystemExit(main())

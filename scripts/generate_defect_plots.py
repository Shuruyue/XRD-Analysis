#!/usr/bin/env python3
"""Generate defect-related plots (stacking fault evolution)."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _script_utils import (
    ensure_dir,
    get_project_root,
    positive_float,
    positive_int,
    resolve_data_dir,
    resolve_path,
    unit_interval,
)

from xrd_analysis.analysis.pipeline import load_bruker_txt, parse_filename
from xrd_analysis.fitting.peak_fitter import fit_peak_with_diagnosis
from xrd_analysis.methods.defect_analysis import analyze_stacking_faults
from xrd_analysis.visualization.defect_plots import (
    plot_stacking_fault_evolution,
)

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

TARGET_PEAKS = {
    (1, 1, 1): 43.316,
    (2, 0, 0): 50.448,
    (2, 2, 0): 74.124,
}


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def process_raw_files(
    data_dir: Path,
    window: float,
    min_r2: float,
    doublet_max_iterations: int,
) -> list[dict]:
    """Process all raw files and build plotting records."""
    samples: list[dict] = []

    if not data_dir.exists():
        logger.error("Data directory not found: %s", data_dir)
        return []

    files = sorted(data_dir.glob("*.txt"))
    logger.info("Found %d raw files in %s", len(files), data_dir)

    for idx, filepath in enumerate(files, 1):
        try:
            file_info = parse_filename(str(filepath))
            two_theta, intensity = load_bruker_txt(str(filepath))

            fitted_peaks: dict[tuple[int, int, int], dict] = {}
            for hkl, center in TARGET_PEAKS.items():
                fit_res = fit_peak_with_diagnosis(
                    two_theta,
                    intensity,
                    center,
                    window=window,
                    use_doublet=True,
                    doublet_max_iterations=doublet_max_iterations,
                )
                if fit_res.get("success", False) and fit_res.get("r_squared", 0.0) >= min_r2:
                    fitted_peaks[hkl] = fit_res

            if not fitted_peaks:
                continue

            sample_entry = {
                "concentration": file_info.get("concentration_ml", 0.0),
                "time": file_info.get("time_hours", 0.0),
                "sample_name": filepath.stem,
                "high_quality": True,
                "stacking_fault": None,
            }

            if (1, 1, 1) in fitted_peaks and (2, 0, 0) in fitted_peaks:
                sf_res = analyze_stacking_faults(
                    fitted_peaks[(1, 1, 1)]["center"],
                    fitted_peaks[(2, 0, 0)]["center"],
                )
                sample_entry["stacking_fault"] = {
                    "alpha_probability": sf_res.alpha_probability,
                    "alpha_percent": sf_res.alpha_percent,
                    "severity": sf_res.severity.value,
                    "peak_separation": sf_res.peak_separation_deg,
                }

            samples.append(sample_entry)
            if idx % 10 == 0:
                logger.info("  Processed %d/%d", idx, len(files))
        except Exception as exc:
            logger.warning("  Skip %s: %s", filepath.name, exc)

    logger.info("Successfully processed %d samples.", len(samples))
    return samples


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate stacking-fault evolution plots.")
    parser.add_argument("--data-dir", type=Path, default=None, help="Directory containing raw .txt scans.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Plot output directory.")
    parser.add_argument(
        "--window",
        type=positive_float,
        default=2.5,
        help="Peak fitting window (deg, >0).",
    )
    parser.add_argument(
        "--min-r2",
        type=unit_interval,
        default=0.8,
        help="Minimum accepted fit R^2 (0 to 1).",
    )
    parser.add_argument(
        "--doublet-max-iterations",
        type=positive_int,
        default=20000,
        help="Maximum iterations for doublet fitting (>0).",
    )
    args = parser.parse_args()

    root_dir = get_project_root()
    config = load_config(root_dir / "config.yaml")

    data_dir = resolve_data_dir(root_dir, args.data_dir)
    default_output = (
        Path(config.get("output", {}).get("paths", {}).get("root", "outputs"))
        / "plots"
        / "microstructure"
        / "defect"
    )
    plots_dir = resolve_path(root_dir, args.output_dir or default_output)
    assert plots_dir is not None
    ensure_dir(plots_dir)

    for old in plots_dir.glob("*.png"):
        old.unlink()

    print("=" * 60)
    print("Generating Defect Analysis Plots")
    print("=" * 60)

    samples = process_raw_files(
        data_dir=data_dir,
        window=args.window,
        min_r2=args.min_r2,
        doublet_max_iterations=args.doublet_max_iterations,
    )
    if not samples:
        logger.error("No valid samples processed.")
        return 1

    plot_stacking_fault_evolution(
        samples,
        output_path=str(plots_dir / "stacking_fault_evolution.png"),
        dpi=300,
        show=False,
    )
    print(f"Done! Plots saved to: {plots_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

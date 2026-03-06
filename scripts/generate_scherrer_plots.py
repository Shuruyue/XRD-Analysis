#!/usr/bin/env python3
"""Generate Scherrer size evolution plots from batch scans."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

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
from xrd_analysis.core.config import ParameterConfig
from xrd_analysis.fitting.peak_fitter import fit_peak_with_diagnosis
from xrd_analysis.methods.scherrer import ScherrerCalculator
from xrd_analysis.visualization.scherrer_plots import (
    plot_scherrer_by_concentration,
    plot_scherrer_evolution_by_peak,
)

warnings.filterwarnings("ignore")

TARGET_PEAKS = {
    (1, 1, 1): 43.316,
    (2, 0, 0): 50.448,
    (2, 2, 0): 74.124,
}


def generate_scherrer_plots(
    data_dir: Path,
    output_dir: Path,
    window: float,
    min_r2: float,
    doublet_max_iterations: int,
) -> int:
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    ensure_dir(output_dir)
    for file in output_dir.glob("*.png"):
        file.unlink()

    cfg = ParameterConfig()
    calculator = ScherrerCalculator(
        use_cubic_habit=True,
        caglioti_params=(
            cfg.instrument.caglioti_u,
            cfg.instrument.caglioti_v,
            cfg.instrument.caglioti_w,
        ),
    )

    files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(files)} samples in {data_dir}")
    if not files:
        return 1

    results: list[dict] = []
    for idx, filepath in enumerate(files, 1):
        file_info = parse_filename(str(filepath))
        two_theta, intensity = load_bruker_txt(str(filepath))

        sample_data = {
            "filename": filepath.name,
            "concentration": file_info.get("concentration_ml", 0.0),
            "time": file_info.get("time_hours", 0.0),
            "peaks": [],
        }

        for hkl, center in TARGET_PEAKS.items():
            fit_res = fit_peak_with_diagnosis(
                two_theta,
                intensity,
                center,
                window=window,
                use_doublet=True,
                doublet_max_iterations=doublet_max_iterations,
            )
            if not fit_res.get("success", False):
                continue
            if fit_res.get("r_squared", 0.0) < min_r2:
                continue

            scherrer_res = calculator.calculate(
                two_theta=fit_res["center"],
                fwhm_observed=fit_res["fwhm"],
                hkl=hkl,
                eta_observed=fit_res.get("eta"),
                eta_instrumental=0.0,
            )
            if scherrer_res.size_nm <= 0:
                continue

            sample_data["peaks"].append(
                {
                    "hkl": hkl,
                    "two_theta": fit_res["center"],
                    "size_nm": scherrer_res.size_nm,
                    "size_err": 0.0,
                    "k_factor": scherrer_res.k_factor,
                }
            )

        if sample_data["peaks"]:
            results.append(sample_data)
        if idx % 10 == 0:
            print(f"  Processed {idx}/{len(files)}")

    print(f"Successfully processed {len(results)} samples.")
    if not results:
        return 2

    plot_scherrer_evolution_by_peak(
        results,
        output_path=str(output_dir / "size_evolution_by_peak.png"),
        show=False,
    )
    plot_scherrer_by_concentration(
        results,
        output_path=str(output_dir / "size_evolution_by_concentration.png"),
        show=False,
    )
    print(f"Done! Plots saved to: {output_dir}")
    return 0


def main() -> int:
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Generate Scherrer size plots.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing raw .txt scans.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/plots/microstructure/scherrer"),
        help="Output plot directory.",
    )
    parser.add_argument(
        "--window",
        type=positive_float,
        default=2.0,
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

    data_dir = resolve_data_dir(root, args.data_dir)
    output_dir = resolve_path(root, args.output_dir)
    assert output_dir is not None
    return generate_scherrer_plots(
        data_dir=data_dir,
        output_dir=output_dir,
        window=args.window,
        min_r2=args.min_r2,
        doublet_max_iterations=args.doublet_max_iterations,
    )


if __name__ == "__main__":
    raise SystemExit(main())

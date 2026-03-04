#!/usr/bin/env python3
"""Generate Williamson-Hall plots from batch scans."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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
from xrd_analysis.methods.scherrer import ScherrerCalculator
from xrd_analysis.methods.williamson_hall import analyze_williamson_hall
from xrd_analysis.visualization.wh_plots import (
    plot_strain_evolution,
    plot_williamson_hall,
)

TARGET_PEAKS = {
    (1, 1, 1): 43.316,
    (2, 0, 0): 50.448,
    (2, 2, 0): 74.124,
    (3, 1, 1): 89.930,
    (2, 2, 2): 95.140,
}


def _load_caglioti(root: Path) -> tuple[float, float, float]:
    cfg_path = root / "config.yaml"
    if not cfg_path.exists():
        return 0.0, 0.0, 0.003
    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    caglioti = cfg.get("instrument", {}).get("caglioti", {})
    u = float(caglioti.get("U") or 0.0)
    v = float(caglioti.get("V") or 0.0)
    w = float(caglioti.get("W") or 0.003)
    return u, v, w


def _hkl_label(hkl: tuple[int, int, int]) -> str:
    return f"({hkl[0]}{hkl[1]}{hkl[2]})"


def generate_wh_plots(
    root: Path,
    data_dir: Path,
    output_dir: Path,
    window: float,
    min_r2: float,
    doublet_max_iterations: int,
) -> int:
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    ensure_dir(output_dir)
    dir_individual = ensure_dir(output_dir / "individual_plots")
    for file in output_dir.glob("*.png"):
        file.unlink()
    for file in dir_individual.glob("*.png"):
        file.unlink()

    files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(files)} samples in {data_dir}")
    if not files:
        return 1

    scherrer = ScherrerCalculator(caglioti_params=_load_caglioti(root))
    results_summary: list[dict] = []

    for idx, filepath in enumerate(files, 1):
        try:
            two_theta, intensity = load_bruker_txt(str(filepath))
            file_info = parse_filename(str(filepath))

            centers: list[float] = []
            fwhm_sample: list[float] = []
            hkl_labels: list[str] = []
            hkl_list: list[tuple[int, int, int]] = []

            for hkl, center_guess in TARGET_PEAKS.items():
                fit_res = fit_peak_with_diagnosis(
                    two_theta,
                    intensity,
                    center_guess,
                    window=window,
                    use_doublet=True,
                    doublet_max_iterations=doublet_max_iterations,
                )
                if not fit_res.get("success", False):
                    continue
                if fit_res.get("r_squared", 0.0) < min_r2:
                    continue

                sch_res = scherrer.calculate(
                    two_theta=fit_res["center"],
                    fwhm_observed=fit_res["fwhm"],
                    hkl=hkl,
                    eta_observed=fit_res.get("eta"),
                    eta_instrumental=0.0,
                )
                if sch_res.fwhm_sample <= 0:
                    continue

                centers.append(float(fit_res["center"]))
                fwhm_sample.append(float(sch_res.fwhm_sample))
                hkl_labels.append(_hkl_label(hkl))
                hkl_list.append(hkl)

            if len(centers) < 3:
                print(f"  [SKIP] {filepath.stem}: not enough valid peaks ({len(centers)}/5)")
                continue

            wh_result = analyze_williamson_hall(
                np.array(centers),
                np.array(fwhm_sample),
                hkl_list=hkl_list,
            )
            results_summary.append(
                {
                    "name": filepath.stem,
                    "concentration": file_info.get("concentration_ml", 0.0),
                    "time": file_info.get("time_hours", 0.0),
                    "microstrain": wh_result.microstrain,
                    "strain_error": wh_result.strain_error,
                    "fitted_peaks": len(centers),
                    "r_squared": wh_result.r_squared,
                }
            )

            plot_williamson_hall(
                np.array(centers),
                np.array(fwhm_sample),
                fit_result={
                    "slope": wh_result.slope,
                    "intercept": wh_result.intercept,
                    "r_squared": wh_result.r_squared,
                },
                hkl_labels=hkl_labels,
                output_path=str(dir_individual / f"{filepath.stem}_WH.png"),
                sample_name=filepath.stem,
                dpi=300,
                show=False,
            )
            plt.close("all")

            if idx % 10 == 0:
                print(f"  Processed {idx}/{len(files)}")
        except Exception as exc:
            print(f"  [ERROR] {filepath.stem}: {exc}")

    print(f"Successfully processed {len(results_summary)} samples.")
    if not results_summary:
        return 2

    plot_strain_evolution(
        results_summary,
        output_path=str(output_dir / "microstrain_evolution.png"),
        show=False,
    )
    print(f"Done! Plots saved to: {output_dir}")
    return 0


def main() -> int:
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Generate Williamson-Hall plots.")
    parser.add_argument("--data-dir", type=Path, default=None, help="Directory containing raw .txt scans.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/plots/microstructure/williamson_hall"),
        help="Output plot directory.",
    )
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

    data_dir = resolve_data_dir(root, args.data_dir)
    output_dir = resolve_path(root, args.output_dir)
    assert output_dir is not None
    return generate_wh_plots(
        root=root,
        data_dir=data_dir,
        output_dir=output_dir,
        window=args.window,
        min_r2=args.min_r2,
        doublet_max_iterations=args.doublet_max_iterations,
    )


if __name__ == "__main__":
    raise SystemExit(main())

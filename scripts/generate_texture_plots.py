#!/usr/bin/env python3
"""Generate texture plots from batch scans."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from scipy.special import voigt_profile

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
from xrd_analysis.fitting.pseudo_voigt import TrueVoigt
from xrd_analysis.methods.texture import analyze_texture
from xrd_analysis.visualization.texture_plots import (
    plot_tc_evolution,
    plot_texture_fraction_single,
)

TARGET_PEAKS = {
    (1, 1, 1): 43.316,
    (2, 0, 0): 50.448,
    (2, 2, 0): 74.124,
}


def calculate_area(amplitude: float, fwhm: float, eta: float) -> float:
    """Calculate integrated area from fitted peak parameters."""
    sigma, gamma = TrueVoigt.params_from_fwhm(fwhm, eta)
    v_max = voigt_profile(0, sigma, gamma)
    if v_max <= 0:
        return 0.0
    return float(amplitude / v_max)


def _hkl_label(hkl: tuple[int, int, int]) -> str:
    return f"({hkl[0]}{hkl[1]}{hkl[2]})"


def generate_texture_plots(
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

    files = sorted(data_dir.glob("*.txt"))
    print(f"Found {len(files)} samples in {data_dir}")
    if not files:
        return 1

    results_summary: list[dict] = []
    for idx, filepath in enumerate(files, 1):
        two_theta, intensity = load_bruker_txt(str(filepath))
        file_info = parse_filename(str(filepath))
        sample_name = filepath.stem

        intensities_area: dict[tuple[int, int, int], float] = {}
        for hkl, center in TARGET_PEAKS.items():
            try:
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

                amp = fit_res["amplitude"]
                fwhm = fit_res["fwhm"]
                eta = fit_res["eta"]
                area_ka1 = calculate_area(amp, fwhm, eta)
                intensities_area[hkl] = area_ka1 * 1.5  # Kalpha total (Ka1 + Ka2)
            except Exception as exc:
                print(f"  Skip peak {sample_name} {hkl}: {exc}")

        if len(intensities_area) < 2:
            print(f"  [SKIP] {sample_name}: not enough valid peaks ({len(intensities_area)}/3)")
            continue

        tc_result = analyze_texture(intensities_area, use_area=True)
        result_entry = {
            "name": sample_name,
            "concentration": file_info.get("concentration_ml", 0.0),
            "time": file_info.get("time_hours", 0.0),
            "tc_values": {_hkl_label(k): v for k, v in tc_result.tc_values.items()},
            "dominant": _hkl_label(tc_result.dominant_hkl) if tc_result.dominant_hkl else "None",
            "is_random": tc_result.is_random,
            "sigma": tc_result.degree_of_texture,
            "high_quality": True,
        }
        results_summary.append(result_entry)

        if idx % 10 == 0:
            print(f"  Processed {idx}/{len(files)}")

    print(f"Successfully processed {len(results_summary)} samples.")
    if not results_summary:
        return 2

    plot_tc_evolution(
        results_summary,
        x_param="time",
        output_path=str(output_dir / "texture_evolution_TC_by_direction.png"),
        normalize=False,
        dpi=300,
        show=False,
    )
    plot_tc_evolution(
        results_summary,
        x_param="time",
        output_path=str(output_dir / "texture_evolution_Fraction_by_direction.png"),
        normalize=True,
        dpi=300,
        show=False,
    )
    plot_texture_fraction_single(
        results_summary,
        x_param="time",
        output_path=str(output_dir / "texture_fraction_by_concentration.png"),
        metric="fraction",
        dpi=300,
        show=False,
    )
    plot_texture_fraction_single(
        results_summary,
        x_param="time",
        output_path=str(output_dir / "texture_TC_by_concentration.png"),
        metric="tc",
        dpi=300,
        show=False,
    )
    print(f"Done! Plots saved to: {output_dir}")
    return 0


def main() -> int:
    root = get_project_root()
    parser = argparse.ArgumentParser(description="Generate texture analysis plots.")
    parser.add_argument("--data-dir", type=Path, default=None, help="Directory containing raw .txt scans.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/plots/microstructure/texture"),
        help="Output plot directory.",
    )
    parser.add_argument(
        "--window",
        type=positive_float,
        default=2.5,
        help="Peak fitting window in degrees (>0).",
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
    return generate_texture_plots(
        data_dir=data_dir,
        output_dir=output_dir,
        window=args.window,
        min_r2=args.min_r2,
        doublet_max_iterations=args.doublet_max_iterations,
    )


if __name__ == "__main__":
    raise SystemExit(main())

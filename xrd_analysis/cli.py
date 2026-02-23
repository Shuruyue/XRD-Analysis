#!/usr/bin/env python3
"""
XRD-Analysis Command Line Interface
=============================

Unified entry point for all analysis operations.
所有分析操作的統一入口點。
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List

from xrd_analysis.__version__ import __version__


# =============================================================================
# LaB6 Standard Reference Material Peak Positions
# NIST SRM 660c - Lanthanum Hexaboride
# =============================================================================
# Reference: NIST Standard Reference Material 660c
#            Certificate of Analysis, National Institute of Standards and Technology
#            a = 4.156826 Å (certified lattice parameter at 22.5 °C)
#            Cu Kα₁ wavelength: λ = 1.540562 Å (Bearden 1967)
#
# Peak positions calculated using Bragg's Law for cubic LaB6
# These values are used for instrument calibration (Caglioti parameters)
# ═══════════════════════════════════════════════════════════════════════════
LAB6_STANDARD_PEAKS = {
    (1, 0, 0): 21.358,   # First reflection
    (1, 1, 0): 30.385,   # Second reflection
    (1, 1, 1): 37.442,   # Third reflection
    (2, 0, 0): 43.507,   # Fourth reflection
    (2, 1, 0): 48.957,   # Fifth reflection
    (2, 1, 1): 53.989,   # Sixth reflection
    (2, 2, 0): 63.218,   # Seventh reflection
    (3, 0, 0): 67.548,   # Eighth reflection
    (3, 1, 0): 71.745,   # Ninth reflection
    (3, 1, 1): 75.844,   # Tenth reflection
}

# Minimum peaks required for stable Caglioti regression
MIN_CALIBRATION_PEAKS = 5


def main(argv: Optional[List[str]] = None) -> int:
    """
    Main CLI entry point.
    主 CLI 入口點。
    """
    parser = argparse.ArgumentParser(
        prog="xrd-analysis",
        description="Advanced XRD Crystallite Size Analysis System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  xrd-analysis analyze data/sample.txt -o outputs/
  xrd-analysis calibrate data/LaB6_standard.txt
  xrd-analysis report outputs/results.csv
        """
    )
    
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands 可用指令")
    
    # analyze command
    analyze_parser = subparsers.add_parser(
        "analyze",
        help="Analyze XRD data 分析 XRD 資料"
    )
    analyze_parser.add_argument(
        "input",
        type=Path,
        help="Input file or directory 輸入檔案或目錄"
    )
    analyze_parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("outputs"),
        help="Output directory 輸出目錄"
    )
    analyze_parser.add_argument(
        "-c", "--config",
        type=Path,
        help="Config file path 配置檔案路徑"
    )
    analyze_parser.add_argument(
        "--batch",
        action="store_true",
        help="Batch processing mode 批次處理模式"
    )
    
    # calibrate command
    cal_parser = subparsers.add_parser(
        "calibrate",
        help="Calibrate instrument 校準儀器"
    )
    cal_parser.add_argument(
        "standard",
        type=Path,
        help="Standard material data file 標準樣品資料檔案"
    )
    cal_parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("calibration.yaml"),
        help="Output calibration file 輸出校正檔案"
    )
    
    # report command
    report_parser = subparsers.add_parser(
        "report",
        help="Generate analysis report 產生分析報告"
    )
    report_parser.add_argument(
        "results",
        type=Path,
        help="Results CSV file 結果 CSV 檔案"
    )
    report_parser.add_argument(
        "-f", "--format",
        choices=["markdown", "html", "pdf"],
        default="markdown",
        help="Report format 報告格式"
    )
    
    args = parser.parse_args(argv)
    
    if args.command is None:
        parser.print_help()
        return 0
    
    try:
        if args.command == "analyze":
            return _run_analyze(args)
        elif args.command == "calibrate":
            return _run_calibrate(args)
        elif args.command == "report":
            return _run_report(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    return 0


def _build_analysis_config(config_path: Optional[Path]):
    """
    Build AnalysisConfig from YAML config (best-effort mapping).
    """
    from xrd_analysis.core.config_loader import load_config
    from xrd_analysis.analysis.pipeline import AnalysisConfig

    loaded = load_config(config_path)
    config = AnalysisConfig()

    if not isinstance(loaded, dict):
        return config

    physical = loaded.get("physical_constants", {})
    if isinstance(physical, dict):
        wavelength = physical.get("wavelength")
        if isinstance(wavelength, (int, float)):
            config.wavelength = float(wavelength)

    instrument = loaded.get("instrument", {})
    if isinstance(instrument, dict):
        caglioti = instrument.get("caglioti", {})
        if isinstance(caglioti, dict):
            u = caglioti.get("U")
            v = caglioti.get("V")
            w = caglioti.get("W")
            if isinstance(u, (int, float)):
                config.caglioti_u = float(u)
            if isinstance(v, (int, float)):
                config.caglioti_v = float(v)
            if isinstance(w, (int, float)):
                config.caglioti_w = float(w)

    fitting = loaded.get("fitting", {})
    if isinstance(fitting, dict):
        peak_detection = fitting.get("peak_detection", {})
        if isinstance(peak_detection, dict):
            min_height = peak_detection.get("min_height")
            if isinstance(min_height, (int, float)):
                config.min_intensity = float(min_height)

    peak_fitting = loaded.get("peak_fitting", {})
    if isinstance(peak_fitting, dict):
        peak_window = peak_fitting.get("peak_window")
        min_intensity = peak_fitting.get("min_intensity")
        if isinstance(peak_window, (int, float)):
            config.peak_window = float(peak_window)
        if isinstance(min_intensity, (int, float)):
            config.min_intensity = float(min_intensity)

    preprocessing = loaded.get("preprocessing", {})
    if isinstance(preprocessing, dict):
        smoothing = preprocessing.get("smoothing", {})
        if isinstance(smoothing, dict):
            enable = smoothing.get("enable")
            window_size = smoothing.get("window_size")
            poly_order = smoothing.get("poly_order")
            if isinstance(enable, bool):
                config.enable_smoothing = enable
            if isinstance(window_size, int):
                config.smoothing_window = window_size
            if isinstance(poly_order, int):
                config.smoothing_poly_order = poly_order

        background = preprocessing.get("background", {})
        if isinstance(background, dict):
            enable = background.get("enable")
            method = background.get("method")
            poly_degree = background.get("poly_degree")
            if isinstance(enable, bool):
                config.enable_background = enable
            if isinstance(method, str):
                config.background_method = method
            if isinstance(poly_degree, int):
                config.background_degree = poly_degree

        kalpha = preprocessing.get("kalpha_strip", {})
        if isinstance(kalpha, dict):
            enable = kalpha.get("enable")
            if isinstance(enable, bool):
                config.enable_kalpha_strip = enable

    return config


def _collect_input_files(input_path: Path, batch: bool) -> List[Path]:
    """
    Resolve analyze input into concrete files.
    """
    input_str = str(input_path)

    if input_path.is_dir():
        files: List[Path] = []
        for pattern in ("*.txt", "*.xy", "*.csv"):
            files.extend(sorted(input_path.glob(pattern)))
        return files

    if any(ch in input_str for ch in "*?[]"):
        base_dir = input_path.parent if str(input_path.parent) not in ("", ".") else Path(".")
        return sorted(base_dir.glob(input_path.name))

    if batch and input_path.exists():
        return [input_path]

    return [input_path]


def _write_summary_csv(results, output_dir: Path) -> Optional[Path]:
    """
    Write one-row-per-sample summary CSV from pipeline results.
    """
    from xrd_analysis.analysis.report_generator import generate_csv_summary

    header = None
    rows = []

    for result in results:
        if result.comprehensive is None:
            continue
        csv_text = generate_csv_summary(result.comprehensive).strip().splitlines()
        if len(csv_text) < 2:
            continue
        if header is None:
            header = csv_text[0]
        rows.append(csv_text[1])

    if not header or not rows:
        return None

    summary_path = output_dir / "summary.csv"
    summary_path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return summary_path


def _run_analyze(args) -> int:
    """
    Run analysis command.
    執行分析指令。
    """
    from xrd_analysis.analysis.pipeline import XRDAnalysisPipeline
    
    print(f"Loading configuration... 載入配置...")
    config = _build_analysis_config(args.config)
    
    print(f"Analyzing: {args.input}")
    print(f"Output to: {args.output}")
    
    # Ensure output directory exists
    args.output.mkdir(parents=True, exist_ok=True)
    
    pipeline = XRDAnalysisPipeline(config)

    files = _collect_input_files(args.input, args.batch)
    if not files:
        print("Error: No input files found.")
        return 1

    results = []
    for f in files:
        if not f.exists():
            print(f"  Skip missing file: {f}")
            continue
        print(f"  Processing: {f.name}")
        result = pipeline.process_file(str(f), str(args.output))
        results.append(result)

    if not results:
        print("Error: No valid input files were processed.")
        return 1

    summary_path = _write_summary_csv(results, args.output)
    if summary_path:
        print(f"Summary CSV written: {summary_path}")
    
    print("Analysis complete. 分析完成。")
    return 0


def _run_calibrate(args) -> int:
    """
    Run calibration command.
    執行校正指令。
    
    Calibrates Caglioti parameters (U, V, W) using a standard material
    like LaB6 (NIST SRM 660c) or Si.
    使用標準樣品（如 LaB6 或 Si）校準 Caglioti 參數。
    """
    from datetime import datetime

    import numpy as np
    import yaml
    from xrd_analysis.analysis.pipeline import load_bruker_txt, find_peak_in_range
    from xrd_analysis.methods.caglioti import CagliotiCorrection
    
    print(f"Calibrating with standard: {args.standard}")
    
    # Load standard data
    try:
        two_theta, intensity = load_bruker_txt(str(args.standard))
    except Exception as e:
        print(f"Error loading file: {e}")
        return 1
    
    if len(two_theta) == 0:
        print("Error: No data found in file")
        return 1
    
    print(f"  Loaded {len(two_theta)} data points")
    print(f"  2θ range: {two_theta.min():.1f}° - {two_theta.max():.1f}°")
    
    # Use module-level LAB6_STANDARD_PEAKS constant (defined at top of file)
    
    # Find peaks and measure FWHM
    peaks_data = []
    matched_peaks = []
    max_position_error_deg = 0.40
    print("\nFitting standard peaks...")
    
    for hkl, expected_pos in LAB6_STANDARD_PEAKS.items():
        if two_theta.min() <= expected_pos <= two_theta.max():
            peak = find_peak_in_range(
                two_theta,
                intensity,
                expected_pos,
                window=2.0,
                use_doublet_fitting=True,
                doublet_max_iterations=8000,
            )
            if peak is not None and peak.fwhm > 0:
                delta_deg = float(peak.two_theta - expected_pos)
                peak_info = {
                    'hkl': hkl,
                    'expected_two_theta': float(expected_pos),
                    'two_theta': peak.two_theta,
                    'fwhm': peak.fwhm,
                    'intensity': peak.intensity,
                    'delta_two_theta': delta_deg,
                }
                peaks_data.append(peak_info)
                if abs(delta_deg) <= max_position_error_deg:
                    matched_peaks.append(peak_info)
                print(
                    f"  - {hkl}: expected={expected_pos:.3f} deg, "
                    f"fitted={peak.two_theta:.3f} deg, "
                    f"delta={delta_deg:+.3f} deg, FWHM={peak.fwhm:.4f} deg"
                )
    
    if len(matched_peaks) < MIN_CALIBRATION_PEAKS:
        print(
            f"\nError: Only {len(matched_peaks)} peaks match the standard position tolerance "
            f"(±{max_position_error_deg:.2f} deg). Need at least {MIN_CALIBRATION_PEAKS}."
        )
        print("Calibration aborted. Verify that input is a true standard scan (LaB6/Si).")
        return 1

    # Fit Caglioti using calibrated class
    two_theta_used = np.array([p["two_theta"] for p in matched_peaks], dtype=float)
    fwhm_used = np.array([p["fwhm"] for p in matched_peaks], dtype=float)

    caglioti = CagliotiCorrection()
    params = caglioti.calibrate(two_theta_used, fwhm_used)
    U, V, W = params.U, params.V, params.W

    # Fit diagnostics on FWHM² regression
    theta_rad = np.radians(two_theta_used / 2.0)
    tan_theta = np.tan(theta_rad)
    fwhm_sq_measured = fwhm_used**2
    fwhm_sq_pred = U * tan_theta**2 + V * tan_theta + W
    residual_sq = fwhm_sq_measured - fwhm_sq_pred

    ss_res = float(np.sum(residual_sq**2))
    ss_tot = float(np.sum((fwhm_sq_measured - np.mean(fwhm_sq_measured))**2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse_fwhm_sq = float(np.sqrt(np.mean(residual_sq**2)))
    max_abs_residual_sq = float(np.max(np.abs(residual_sq)))

    # Sanity check over measured angular window
    tt_grid = np.linspace(float(two_theta.min()), float(two_theta.max()), 200)
    tt_grid = tt_grid[(tt_grid >= 15.0) & (tt_grid <= 120.0)]
    grid_fwhm_sq = U * np.tan(np.radians(tt_grid / 2.0))**2 + V * np.tan(np.radians(tt_grid / 2.0)) + W
    has_negative_region = bool(np.any(grid_fwhm_sq <= 0))

    print(f"\nCaglioti parameters fitted from {len(matched_peaks)} matched peaks:")
    print(f"  U = {U:.6f}")
    print(f"  V = {V:.6f}")
    print(f"  W = {W:.6f}")
    fwhm_43 = np.sqrt(max(U * np.tan(np.radians(21.5))**2 + V * np.tan(np.radians(21.5)) + W, 0.0))
    print(f"  FWHM at 43 deg ~ {fwhm_43:.4f} deg")
    print(f"\nCalibration diagnostics:")
    print(f"  R² (FWHM² fit) = {r_squared:.6f}")
    print(f"  RMSE(FWHM²) = {rmse_fwhm_sq:.6e} deg²")
    print(f"  Max |residual| (FWHM²) = {max_abs_residual_sq:.6e} deg²")
    if has_negative_region:
        print("  WARNING: FWHM² becomes non-positive in part of scan range; verify standard data quality.")
    
    # Save calibration
    calibration = {
        'caglioti': {
            'U': float(U),
            'V': float(V),
            'W': float(W),
        },
        'standard': str(args.standard),
        'n_peaks_used': len(matched_peaks),
        'n_peaks_found': len(peaks_data),
        'position_tolerance_deg': max_position_error_deg,
        'calibrated_at': datetime.now().isoformat(timespec="seconds"),
        'fit_quality': {
            'r_squared_fwhm_sq': r_squared,
            'rmse_fwhm_sq_deg2': rmse_fwhm_sq,
            'max_abs_residual_fwhm_sq_deg2': max_abs_residual_sq,
            'has_non_positive_fwhm_sq_region': has_negative_region,
        },
        'peaks': [
            {
                'hkl': f"({p['hkl'][0]}{p['hkl'][1]}{p['hkl'][2]})",
                'expected_two_theta_deg': float(p['expected_two_theta']),
                'two_theta_deg': float(p['two_theta']),
                'delta_two_theta_deg': float(p['delta_two_theta']),
                'fwhm_deg': float(p['fwhm']),
            }
            for p in matched_peaks
        ],
    }
    
    with open(args.output, 'w', encoding='utf-8') as f:
        yaml.safe_dump(calibration, f, default_flow_style=False, sort_keys=False)
    
    print(f"\nCalibration saved to: {args.output}")
    print("Calibration complete. 校正完成。")
    return 0


def _run_report(args) -> int:
    """
    Run report generation command.
    執行報告生成指令。
    
    Generates analysis report from results CSV file.
    從結果 CSV 檔案生成分析報告。
    """
    import pandas as pd
    from datetime import datetime
    from xrd_analysis.analysis.report_generator import generate_comprehensive_report
    
    print(f"Generating {args.format} report from: {args.results}")
    
    # Load results
    if not args.results.exists():
        print(f"Error: Results file not found: {args.results}")
        return 1
    
    try:
        df = pd.read_csv(args.results)
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return 1
    
    print(f"  Loaded {len(df)} analysis records")
    
    # Generate report content
    report_lines = [
        "# XRD-Analysis Analysis Report",
        f"# XRD-Analysis 分析報告",
        "",
        f"**Generated / 生成時間**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        "---",
        "",
        "## Summary / 摘要",
        "",
    ]
    
    # Add summary statistics
    if 'size_nm' in df.columns:
        valid_sizes = df['size_nm'].dropna()
        if len(valid_sizes) > 0:
            report_lines.extend([
                f"- **Total samples / 總樣品數**: {len(df)}",
                f"- **Average crystallite size / 平均晶粒尺寸**: {valid_sizes.mean():.1f} nm",
                f"- **Size range / 尺寸範圍**: {valid_sizes.min():.1f} - {valid_sizes.max():.1f} nm",
                "",
            ])
    
    # Add data table
    report_lines.extend([
        "## Data Table / 資料表",
        "",
        df.to_markdown(index=False) if hasattr(df, 'to_markdown') else df.to_string(),
        "",
    ])
    
    report_content = "\n".join(report_lines)
    
    # Determine output path
    output_path = args.results.parent / f"report_{args.results.stem}.{args.format}"
    if args.format == "markdown":
        output_path = output_path.with_suffix(".md")
    
    # Save report
    if args.format == "markdown":
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
        print(f"\nMarkdown report saved to: {output_path}")
    elif args.format == "html":
        try:
            import markdown
            html_content = markdown.markdown(report_content, extensions=['tables'])
            html_output = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>XRD-Analysis Report</title>
    <style>
        body {{ font-family: 'Segoe UI', sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #4a90d9; color: white; }}
    </style>
</head>
<body>
{html_content}
</body>
</html>"""
            output_path = output_path.with_suffix(".html")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_output)
            print(f"\nHTML report saved to: {output_path}")
        except ImportError:
            print("Warning: 'markdown' package not installed. Saving as markdown instead.")
            output_path = output_path.with_suffix(".md")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(report_content)
    else:
        print(f"Warning: {args.format} format not yet implemented. Saving as markdown.")
        output_path = output_path.with_suffix(".md")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report_content)
    
    print("Report generated. 報告已生成。")
    return 0


if __name__ == "__main__":
    sys.exit(main())

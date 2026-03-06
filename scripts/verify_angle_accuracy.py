#!/usr/bin/env python3
"""Verify peak-angle correctness against expected Cu reference peaks.

This script runs the project pipeline and reports:
- mean peak offset
- RMSE of peak offsets
- max absolute peak offset
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from _script_utils import positive_float

from xrd_analysis.analysis.pipeline import XRDAnalysisPipeline


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify peak-angle correctness (measured vs expected Cu peak positions)."
    )
    parser.add_argument("input", type=Path, help="XRD data file (.txt)")
    parser.add_argument(
        "--max-offset-threshold",
        type=positive_float,
        default=0.10,
        help="Maximum acceptable absolute peak offset in degrees (default: 0.10)",
    )
    args = parser.parse_args()

    pipeline = XRDAnalysisPipeline()
    result = pipeline.analyze(str(args.input))

    if result.angle_validation_peaks == 0:
        print("No valid peaks found for angle verification.")
        return 1

    print(f"Sample: {result.sample_name}")
    print(f"Angle validation peaks: {result.angle_validation_peaks}")
    print(f"Mean offset (deg): {result.angle_offset_mean_deg:+.4f}")
    print(f"RMSE (deg): {result.angle_offset_rmse_deg:.4f}")
    print(f"Max |offset| (deg): {result.angle_offset_max_abs_deg:.4f}")

    if (
        result.angle_offset_max_abs_deg is not None
        and result.angle_offset_max_abs_deg > args.max_offset_threshold
    ):
        print(
            f"FAIL: max |offset| exceeds threshold ({args.max_offset_threshold:.3f} deg). "
            "Check zero-shift/specimen displacement."
        )
        return 2

    print("PASS: angle offsets are within threshold.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

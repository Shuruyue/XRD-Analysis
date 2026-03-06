#!/usr/bin/env python3
"""Verify background-correction behavior on one XRD scan.

Reports basic quantitative diagnostics:
- mean raw intensity
- mean estimated background
- mean corrected intensity
- background/raw mean ratio
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from xrd_analysis.analysis.pipeline import load_bruker_txt
from xrd_analysis.preprocessing.pipeline import PreprocessingPipeline


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify background-correction quality on one scan."
    )
    parser.add_argument("input", type=Path, help="XRD data file (.txt)")
    args = parser.parse_args()

    two_theta, intensity = load_bruker_txt(str(args.input))
    if len(two_theta) == 0:
        print("No data points loaded.")
        return 1

    pipeline = PreprocessingPipeline(
        enable_smoothing=True,
        enable_background=True,
        enable_kalpha_strip=True,
    )
    pre = pipeline.run(two_theta, intensity)

    raw_mean = float(np.mean(pre.raw_intensity))
    corrected_mean = float(np.mean(pre.intensity))
    print(f"Sample points: {len(two_theta)}")
    print(f"Mean raw intensity: {raw_mean:.3f}")
    print(f"Mean corrected intensity: {corrected_mean:.3f}")

    if pre.background is None:
        print("Background estimation is not available.")
        return 2

    bg_mean = float(np.mean(pre.background))
    ratio = bg_mean / max(raw_mean, 1e-12)
    print(f"Mean background intensity: {bg_mean:.3f}")
    print(f"Background/raw mean ratio: {ratio:.4f}")

    negative_after = int(np.sum(pre.intensity < 0))
    print(f"Negative points after correction: {negative_after}")
    print("PASS: background correction executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

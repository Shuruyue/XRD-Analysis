#!/usr/bin/env python3
"""Verify core XRD constants and Cu reference peaks.

This script checks:
1) Cu Kalpha constants consistency
2) Cu FCC d-spacing table consistency
3) Cu FCC 2theta table consistency via Bragg's law
4) Kalpha2 position shift direction (Kalpha2 must be at higher angle)
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from xrd_analysis.core.constants import CU_KA1, CU_KA2, CU_KA_AVG, KA2_KA1_RATIO
from xrd_analysis.core.copper_crystal import CU_CRYSTAL, CU_JCPDS_EXTENDED
from xrd_analysis.fitting.ka_doublet import calculate_ka2_position


def calc_d_spacing(hkl: tuple[int, int, int], a: float) -> float:
    h, k, l = hkl
    return a / math.sqrt(h * h + k * k + l * l)


def calc_two_theta(hkl: tuple[int, int, int], a: float, wavelength: float) -> float:
    d = calc_d_spacing(hkl, a)
    theta = math.asin(wavelength / (2.0 * d))
    return 2.0 * math.degrees(theta)


def check_constants() -> None:
    expected_avg = (2.0 * CU_KA1 + CU_KA2) / 3.0
    if abs(CU_KA_AVG - expected_avg) > 1e-6:
        raise AssertionError(
            f"CU_KA_AVG mismatch: got {CU_KA_AVG:.9f}, expected {expected_avg:.9f}"
        )

    if abs(KA2_KA1_RATIO - 0.5) > 1e-9:
        raise AssertionError(f"KA2_KA1_RATIO should be 0.5, got {KA2_KA1_RATIO}")


def check_cu_jcpds_table() -> dict[str, float]:
    max_d_err = 0.0
    max_tt_err = 0.0

    for hkl, data in CU_JCPDS_EXTENDED.items():
        d_calc = calc_d_spacing(hkl, CU_CRYSTAL.lattice_constant)
        tt_calc = calc_two_theta(hkl, CU_CRYSTAL.lattice_constant, CU_KA1)

        d_err = abs(d_calc - data["d_spacing"])
        tt_err = abs(tt_calc - data["two_theta"])

        max_d_err = max(max_d_err, d_err)
        max_tt_err = max(max_tt_err, tt_err)

    # Code table is rounded to 3 decimals; 0.002 A / 0.01 deg are strict enough.
    if max_d_err > 0.002:
        raise AssertionError(f"Max d-spacing error too large: {max_d_err:.6f} A")
    if max_tt_err > 0.01:
        raise AssertionError(f"Max 2theta error too large: {max_tt_err:.6f} deg")

    return {"max_d_err": max_d_err, "max_tt_err": max_tt_err}


def check_ka2_shift_direction() -> float:
    min_shift = float("inf")
    for hkl, data in CU_JCPDS_EXTENDED.items():
        tt1 = data["two_theta"]
        tt2 = calculate_ka2_position(tt1)
        shift = tt2 - tt1
        min_shift = min(min_shift, shift)
    return min_shift


def main() -> int:
    check_constants()
    table_stats = check_cu_jcpds_table()
    min_shift = check_ka2_shift_direction()
    if min_shift <= 0:
        raise AssertionError(f"Found invalid Kalpha2 shift: {min_shift:.6f} deg")

    print("Physics verification passed")
    print(f"  Cu Kalpha1 = {CU_KA1:.6f} A")
    print(f"  Cu Kalpha2 = {CU_KA2:.6f} A")
    print(f"  Cu Kalpha(avg) = {CU_KA_AVG:.6f} A")
    print(f"  Kalpha2/Kalpha1 ratio = {KA2_KA1_RATIO:.3f}")
    print(f"  Cu lattice constant a0 = {CU_CRYSTAL.lattice_constant:.4f} A")
    print(f"  Max d-spacing table error = {table_stats['max_d_err']:.6f} A")
    print(f"  Max 2theta table error = {table_stats['max_tt_err']:.6f} deg")
    print(f"  Min Kalpha2 shift = {min_shift:.6f} deg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

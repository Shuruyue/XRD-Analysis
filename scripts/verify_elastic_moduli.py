#!/usr/bin/env python3
"""
Verify directional elastic constants used in stress analysis.

This script checks:
1) Directional Young's modulus from stiffness tensor vs mapped values
2) Directional Poisson ratio consistency nu_hkl ~= -S12 * E_hkl
3) Zener anisotropy ratio from C11/C12/C44
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, Tuple

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from xrd_analysis.core.copper_crystal import (
    CU_ELASTIC,
    CU_POISSON,
    calculate_youngs_modulus_from_stiffness,
    get_poisson_ratio,
    get_youngs_modulus,
)


def compliance_constants(c11: float, c12: float, c44: float) -> Tuple[float, float, float]:
    denom = (c11 - c12) * (c11 + 2.0 * c12)
    s11 = (c11 + c12) / denom
    s12 = -c12 / denom
    s44 = 1.0 / c44
    return s11, s12, s44


def verify_directional_modulus() -> Dict[str, float]:
    expected_map = {
        (1, 1, 1): CU_ELASTIC.E_111,
        (1, 0, 0): CU_ELASTIC.E_100,
        (1, 1, 0): CU_ELASTIC.E_110,
    }

    max_err = 0.0
    for hkl, expected in expected_map.items():
        calc = calculate_youngs_modulus_from_stiffness(*hkl)
        max_err = max(max_err, abs(calc - expected))

    # 1 GPa tolerance is strict enough for rounded constants in docs/code.
    if max_err > 1.0:
        raise AssertionError(f"Directional Young's modulus mismatch too large: {max_err:.3f} GPa")

    return {"max_modulus_err_gpa": max_err}


def verify_poisson_consistency() -> Dict[str, float]:
    c11, c12, c44 = 168.4, 121.4, 75.4
    # Closed-form cubic expressions for principal directions.
    nu_100_calc = c12 / (c11 + c12)
    nu_111_calc = (c11 + 2.0 * c12 - 2.0 * c44) / (2.0 * (c11 + 2.0 * c12 + c44))

    checks = [
        abs(nu_100_calc - CU_POISSON.nu_200),
        abs(nu_111_calc - CU_POISSON.nu_111),
        abs(get_poisson_ratio(2, 0, 0, use_directional=True) - CU_POISSON.nu_200),
        abs(get_poisson_ratio(1, 1, 1, use_directional=True) - CU_POISSON.nu_111),
        abs(get_poisson_ratio(2, 2, 0, use_directional=True) - CU_POISSON.nu_220),
    ]
    max_err = max(checks)

    if max_err > 0.02:
        raise AssertionError(f"Poisson ratio mismatch too large: {max_err:.4f}")

    # Basic physical sanity check.
    nu_values = [
        CU_POISSON.nu_111,
        CU_POISSON.nu_200,
        CU_POISSON.nu_220,
        CU_POISSON.nu_311,
        CU_POISSON.nu_poly,
    ]
    if any((nu <= 0.0 or nu >= 0.5) for nu in nu_values):
        raise AssertionError(f"Invalid Poisson ratio range: {nu_values}")

    return {"max_poisson_err": max_err}


def verify_zener_ratio() -> float:
    c11, c12, c44 = 168.4, 121.4, 75.4
    a_zener = 2.0 * c44 / (c11 - c12)
    if abs(a_zener - 3.21) > 0.02:
        raise AssertionError(f"Zener ratio mismatch: {a_zener:.4f}")
    return a_zener


def main() -> int:
    modulus_stats = verify_directional_modulus()
    poisson_stats = verify_poisson_consistency()
    zener = verify_zener_ratio()

    print("Elastic-moduli verification passed")
    print(f"  E<111> = {CU_ELASTIC.E_111:.1f} GPa")
    print(f"  E<100> = {CU_ELASTIC.E_100:.1f} GPa")
    print(f"  E<110> = {CU_ELASTIC.E_110:.1f} GPa")
    print(f"  nu_111 = {CU_POISSON.nu_111:.3f}")
    print(f"  nu_200 = {CU_POISSON.nu_200:.3f}")
    print(f"  nu_220 = {CU_POISSON.nu_220:.3f}")
    print(f"  Max modulus error = {modulus_stats['max_modulus_err_gpa']:.3f} GPa")
    print(f"  Max poisson error = {poisson_stats['max_poisson_err']:.4f}")
    print(f"  Zener anisotropy ratio = {zener:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

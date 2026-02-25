"""Peak Finding with Pseudo-Voigt Fitting 偽Voigt擬合峰值查找
=============================================================

Find and fit peaks in XRD data using Kα doublet fitting with
Pseudo-Voigt fallback.
使用 Kα 雙峰擬合（含偽Voigt回退）在 XRD 數據中查找和擬合峰值。
"""

import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from xrd_analysis.fitting.lm_optimizer import LMOptimizer
from xrd_analysis.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams

logger = logging.getLogger(__name__)


@dataclass
class PeakData:
    """Single peak data."""

    hkl: Tuple[int, int, int]
    two_theta: float
    intensity: float
    fwhm: float
    area: float = 0.0
    eta: float = 0.5  # Pseudo-Voigt mixing parameter (0=Gaussian, 1=Lorentzian)


def _estimate_fwhm_simple(
    theta_range: np.ndarray, int_range: np.ndarray, idx_max: int, peak_int: float
) -> float:
    """Estimate FWHM using half-maximum method."""
    half_max = peak_int / 2

    left_idx = idx_max
    while left_idx > 0 and int_range[left_idx] > half_max:
        left_idx -= 1

    right_idx = idx_max
    while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
        right_idx += 1

    return max(theta_range[right_idx] - theta_range[left_idx], 0.1)


def _fit_peak_pseudo_voigt(
    theta_range: np.ndarray,
    int_range: np.ndarray,
    peak_theta: float,
    peak_int: float,
    initial_fwhm: float,
) -> Tuple[bool, float, float, float, float, float]:
    """Try Pseudo-Voigt fitting. Returns (success, theta, intensity, fwhm, eta, area)."""
    try:
        optimizer = LMOptimizer(max_iterations=500, tolerance=1e-6)

        initial_guess = PseudoVoigtParams(
            center=peak_theta, amplitude=peak_int, fwhm=initial_fwhm, eta=0.5
        )

        fit_result = optimizer.fit_single_peak(
            theta_range, int_range, initial_guess=initial_guess
        )

        if fit_result.success and fit_result.r_squared > 0.8:
            fitted_curve = PseudoVoigt.profile(
                theta_range,
                fit_result.params.center,
                fit_result.params.amplitude,
                fit_result.params.fwhm,
                fit_result.params.eta,
            )
            try:
                area = np.trapezoid(fitted_curve, theta_range)
            except AttributeError:
                area = np.trapz(fitted_curve, theta_range)

            return (
                True,
                fit_result.params.center,
                fit_result.params.amplitude,
                fit_result.params.fwhm,
                fit_result.params.eta,
                area,
            )
    except (ValueError, RuntimeError, np.linalg.LinAlgError):
        pass

    return False, peak_theta, peak_int, initial_fwhm, 0.5, 0.0


def _calculate_peak_area_simple(
    theta_range: np.ndarray, int_range: np.ndarray
) -> float:
    """Calculate peak area using trapezoidal integration."""
    try:
        return np.trapezoid(int_range, theta_range)
    except AttributeError:
        return np.trapz(int_range, theta_range)


def find_peak_in_range(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    center: float,
    window: float = 2.5,
    use_doublet_fitting: bool = True,
    doublet_max_iterations: int = 20000,
    min_fit_r_squared: float = 0.80,
) -> Optional[PeakData]:
    """Find peak near expected position using Kα doublet fitting.

    Args:
        two_theta: 2θ array
        intensity: Intensity array
        center: Expected peak center position
        window: Search window (degrees)
        use_doublet_fitting: If True, use DoubletFitter (recommended)
        doublet_max_iterations: Max iterations for doublet optimizer
        min_fit_r_squared: Minimum accepted fit R²

    Returns:
        PeakData with fitted parameters (Kα₁ only), or None if no peak found

    """
    # Select range
    mask = (two_theta >= center - window) & (two_theta <= center + window)
    if not np.any(mask):
        return None

    theta_range = two_theta[mask]
    int_range = intensity[mask]

    # Find maximum
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]

    if peak_int < 50:
        return None

    # Estimate initial FWHM
    initial_fwhm = _estimate_fwhm_simple(theta_range, int_range, idx_max, peak_int)

    # Try Kα doublet fitting using unified function
    if use_doublet_fitting:
        try:
            from xrd_analysis.fitting.peak_fitter import fit_peak_with_diagnosis

            fit_result = fit_peak_with_diagnosis(
                two_theta,
                intensity,
                center,
                window=window,
                use_doublet=True,
                doublet_max_iterations=doublet_max_iterations,
            )

            if (
                fit_result["success"]
                and fit_result.get("r_squared", 0) > min_fit_r_squared
            ):
                from xrd_analysis.fitting.pv_area import calculate_pv_area

                area = calculate_pv_area(
                    fit_result["amplitude"],
                    fit_result["fwhm"],
                    fit_result.get("eta", 0.5),
                )
                return PeakData(
                    hkl=(0, 0, 0),
                    two_theta=fit_result["center"],
                    intensity=fit_result["amplitude"],
                    fwhm=fit_result["fwhm"],
                    area=area,
                    eta=fit_result.get("eta", 0.5),
                )
        except (ValueError, RuntimeError, ImportError):
            logger.debug("Doublet fitting failed for center=%.2f, using fallback", center)

    # Fallback: simple Pseudo-Voigt fitting
    success, peak_theta, peak_int, fwhm, eta, area = _fit_peak_pseudo_voigt(
        theta_range, int_range, peak_theta, peak_int, initial_fwhm
    )
    if success:
        return PeakData(
            hkl=(0, 0, 0),
            two_theta=peak_theta,
            intensity=peak_int,
            fwhm=fwhm,
            area=area,
            eta=eta,
        )

    # Final fallback: simple method
    fwhm = max(initial_fwhm, 0.05)
    area = _calculate_peak_area_simple(theta_range, int_range)

    return PeakData(
        hkl=(0, 0, 0),
        two_theta=peak_theta,
        intensity=peak_int,
        fwhm=fwhm,
        area=area,
        eta=0.5,  # Default to intermediate
    )

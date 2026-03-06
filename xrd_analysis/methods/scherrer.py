"""Scherrer Crystallite Size Analysis.
=======================================================

Implements the Scherrer equation with dynamic K values, validity flags,
and complete unit conversion safety.

Reference:
    Langford, J. I., & Wilson, A. J. C. (1978).
    Scherrer after sixty years. J. Appl. Cryst., 11, 102-113.
"""

from dataclasses import dataclass
from enum import Enum

import numpy as np

from xrd_analysis.core.constants import (
    CU_KA1,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
    MIN_RELIABLE_SIZE,
    SCHERRER_K,
)
from xrd_analysis.core.copper_crystal import get_k_for_hkl
from xrd_analysis.fitting.hkl_assignment import assign_hkl
from xrd_analysis.fitting.pseudo_voigt import tch_components_from_eta

# =============================================================================
# Constants
# =============================================================================

# Wavelength: use CU_KA1 directly
FWHM_RATIO_THRESHOLD = MIN_BROADENING_RATIO

# =============================================================================
# Validity Flag System
# =============================================================================


class ValidityFlag(Enum):
    """Scherrer calculation validity flags."""

    VALID = "VALID"  # Normal calculation
    UNRELIABLE = "UNRELIABLE"  # FWHM ratio below threshold
    WARNING = "WARNING"  # Size exceeds limits
    ERROR = "ERROR"  # Calculation failed


class GrainShape(Enum):
    """Grain shape for Scherrer constant selection."""

    SPHERICAL = "spherical"
    CUBIC = "cubic"
    CUSTOM = "custom"


# =============================================================================
# Result Dataclass
# =============================================================================


@dataclass
class ScherrerResult:
    """Result from Scherrer crystallite size calculation.

    Includes validity flags and complete metadata.

    Attributes:
        size_nm: Crystallite size in nanometers
        size_angstrom: Crystallite size in Ångströms
        two_theta: Peak position (2θ) in degrees
        hkl: Miller indices if assigned
        k_factor: Scherrer constant used
        fwhm_observed: Observed FWHM in degrees
        fwhm_instrumental: Instrumental FWHM in degrees
        fwhm_sample: Sample broadening in degrees
        fwhm_sample_rad: Sample broadening in radians
        validity_flag: Calculation validity status
        warning_message: Warning or error message
        is_reliable: True if result is reliable

    """

    size_nm: float
    size_angstrom: float
    two_theta: float
    hkl: tuple[int, int, int] | None = None
    k_factor: float = 0.829  # L&W 1978 spherical standard
    fwhm_observed: float = 0.0
    fwhm_instrumental: float = 0.0
    fwhm_sample: float = 0.0
    fwhm_sample_rad: float = 0.0
    validity_flag: ValidityFlag = ValidityFlag.VALID
    warning_message: str = ""
    is_reliable: bool = True

    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})" if self.hkl else "N/A"
        return (
            f"ScherrerResult(hkl={hkl_str}, size={self.size_nm:.1f} nm, "
            f"K={self.k_factor:.3f}, flag={self.validity_flag.value})"
        )


# =============================================================================
# Scherrer Calculator
# =============================================================================


class ScherrerCalculator:
    """Scherrer equation for crystallite size calculation.

    D = K × λ / (β × cos θ)

    Features:
        - Dynamic K values based on hkl (cubic habit)
        - Validity flag system (VALID/UNRELIABLE/WARNING)
        - Unit conversion safety (degrees ↔ radians)
        - Caglioti instrumental broadening integration

    Args:
        wavelength: X-ray wavelength in Å (default: Cu Kα1 = 1.540562, Bearden 1967)
        use_cubic_habit: Use ED-Cu specific K values (default: True)
        caglioti_params: (U, V, W) tuple for instrumental broadening

    Example:
        >>> calc = ScherrerCalculator()
        >>> result = calc.calculate(43.32, 0.25, fwhm_instrumental=0.08)
        >>> print(f"D = {result.size_nm:.1f} nm")
        D = 49.0 nm

    """

    def __init__(
        self,
        wavelength: float = CU_KA1,
        use_cubic_habit: bool = True,
        caglioti_params: tuple[float, float, float] | None = None,
        deconvolution_method: str = "auto",
    ) -> None:
        self.wavelength = wavelength
        self.use_cubic_habit = use_cubic_habit
        self.caglioti_params = caglioti_params
        self.deconvolution_method = deconvolution_method

    def calculate(
        self,
        two_theta: float,
        fwhm_observed: float,
        fwhm_instrumental: float | None = None,
        hkl: tuple[int, int, int] | None = None,
        correction_method: str | None = None,
        **kwargs,
    ) -> ScherrerResult:
        """Calculate crystallite size using Scherrer equation.

        Args:
            two_theta: Peak position in degrees (2θ).
            fwhm_observed: Observed FWHM in degrees.
            fwhm_instrumental: Instrumental FWHM in degrees (optional).
                If None, Caglioti parameters are used.
            hkl: Miller indices for K selection (auto-assigned if None).
            correction_method:
                Broadening subtraction method:
                - "quadratic": β = sqrt(β_obs² - β_inst²)
                - "voigt": component-wise G/L subtraction
                - "auto": use "voigt" only when eta_observed is provided,
                          otherwise use "quadratic"

        Returns:
            ScherrerResult with calculated size and metadata.

        Raises:
            ValueError: If fwhm_observed <= 0.

        """
        warnings = []

        # Input validation
        if fwhm_observed <= 0:
            raise ValueError(f"fwhm_observed must be positive, got {fwhm_observed}")
        if not (0 < two_theta < 180):
            raise ValueError(f"two_theta must be between 0° and 180°, got {two_theta}")

        # Auto-assign hkl if not provided
        if hkl is None:
            hkl = assign_hkl(two_theta)

        # Get K factor
        if self.use_cubic_habit and hkl:
            k_factor = get_k_for_hkl(hkl[0], hkl[1], hkl[2])
        else:
            k_factor = SCHERRER_K.default

        # Calculate instrumental FWHM if Caglioti params available
        if fwhm_instrumental is None and self.caglioti_params:
            fwhm_instrumental = self._calculate_caglioti(two_theta)
        elif fwhm_instrumental is None:
            fwhm_instrumental = 0.0

        # Check validity threshold
        validity_flag = ValidityFlag.VALID
        is_reliable = True

        # =========================================================================
        # Advanced Voigt Component Deconvolution (Accuracy: Highest)
        # =========================================================================
        # Algorithm References:
        #
        # 1. Deconvolution Principle (G & L separation):
        #    Keijser, Th.H., Mittemeijer, E.J. & Rozendaal, H.C.F. (1983).
        #    "The determination of crystallite-size and lattice-strain parameters
        #     in conjunction with the profile-refinement method".
        #    J. Appl. Cryst. 16, 309-316.
        #    - Eq 8: β_C_structure = β_C_profile - β_C_standard
        #    - Eq 9: β_G_structure = sqrt(β_G_profile^2 - β_G_standard^2)
        #
        # 2. Recombination (Voigt FWHM Approximation):
        #    Olivero, J.J. & Longbothum, R.L. (1977).
        #    "Empirical fits to the Voigt line width: A brief review".
        #    J. Quant. Spectrosc. Radiat. Transfer, 17, 233.
        #    - Formula: f_V ≈ 0.5346 f_L + √(0.2166 f_L² + f_G²)
        #    - Error < 0.02%
        #
        # =========================================================================
        # Deconvolve Gaussian and Lorentzian components separately
        #
        # Assumptions:
        # 1. Instrument profile is largely Gaussian (η_inst ≈ 0) or user provided
        # 2. Sample profile is Voigt (convolution of G and L)
        #
        # Algorithm:
        # 1. Convert Obs and Inst (FWHM, η) -> (fG, fL) components
        # 2. Subtract: fL_samp = fL_obs - fL_inst
        #              fG_samp² = fG_obs² - fG_inst²
        # 3. Recombine: FWHM_samp ≈ 0.5346 fL + √(0.2166 fL² + fG²)
        # =========================================================================

        # Respect pipeline argument names first, then backward-compatible aliases.
        eta_obs = kwargs.get("eta_observed", kwargs.get("eta"))
        eta_inst = kwargs.get(
            "eta_instrumental", 0.0
        )  # Instrument assumed Gaussian by default
        method = (correction_method or self.deconvolution_method or "auto").lower()

        if fwhm_instrumental > 0:
            ratio = fwhm_observed / fwhm_instrumental

            if ratio < FWHM_RATIO_THRESHOLD:
                validity_flag = ValidityFlag.UNRELIABLE
                warnings.append(f"FWHM ratio {ratio:.2f} < {FWHM_RATIO_THRESHOLD}")
                is_reliable = False

            if method == "auto":
                # If eta is unavailable, use textbook quadratic subtraction
                # to remain consistent with documented examples.
                method = "voigt" if eta_obs is not None else "quadratic"

            if method == "quadratic":
                fwhm_sq_diff = fwhm_observed**2 - fwhm_instrumental**2
                if fwhm_sq_diff <= 0:
                    fwhm_sample = 0.001  # Prevent divide-by-zero / NaN
                    if fwhm_observed < fwhm_instrumental:
                        warnings.append(
                            f"FWHM_obs ({fwhm_observed:.4f}°) < FWHM_inst ({fwhm_instrumental:.4f}°)"
                        )
                else:
                    fwhm_sample = float(np.sqrt(fwhm_sq_diff))
            elif method == "voigt":
                if eta_obs is None:
                    eta_obs = 0.5
                    warnings.append(
                        "eta_observed missing; fallback eta=0.5 for Voigt subtraction"
                    )

                # 1. Decompose Observed
                fG_obs, fL_obs = self._get_voigt_components(fwhm_observed, eta_obs)

                # 2. Decompose Instrumental
                fG_inst, fL_inst = self._get_voigt_components(
                    fwhm_instrumental, eta_inst
                )

                # 3. Component Subtraction
                # Lorentzian: Linear subtraction
                fL_sample = fL_obs - fL_inst

                # Gaussian: Quadratic subtraction
                fG_sq_diff = fG_obs**2 - fG_inst**2

                # Validity Check (Physically impossible cases)
                if fL_sample < 0 and fG_sq_diff < 0:
                    # Both components smaller than instrument -> Completely unreliable
                    validity_flag = ValidityFlag.UNRELIABLE
                    warnings.append(
                        f"FWHM_obs ({fwhm_observed:.4f}°) < FWHM_inst ({fwhm_instrumental:.4f}°)"
                    )
                    is_reliable = False
                    fwhm_sample = 0.001  # prevent nan
                else:
                    # Handle cases where one component might be slightly negative due to noise/fitting error
                    # We clamp negative correlations to 0 but warn if significant
                    if fG_sq_diff < 0:
                        fG_sample = 0
                        if abs(fG_sq_diff) > 0.0001:  # Tolerance
                            warnings.append(
                                "Gaussian component smaller than instrument"
                            )
                    else:
                        fG_sample = np.sqrt(fG_sq_diff)

                    if fL_sample < 0:
                        fL_sample_raw = fL_sample
                        fL_sample = 0
                        if abs(fL_sample_raw) > 0.0001:
                            warnings.append(
                                "Lorentzian component smaller than instrument"
                            )

                    # 4. Recombine (Olivero-Longbothum approximation)
                    fwhm_sample = 0.5346 * fL_sample + np.sqrt(
                        0.2166 * fL_sample**2 + fG_sample**2
                    )
            else:
                raise ValueError(f"Unknown correction_method: {method}")
        else:
            # No instrumental correction
            fwhm_sample = fwhm_observed

        # CRITICAL: Convert to radians
        theta_rad = np.radians(two_theta / 2)
        fwhm_sample_rad = np.radians(fwhm_sample)

        # Scherrer equation: D = K × λ / (β × cos θ)
        cos_theta = np.cos(theta_rad)
        size_angstrom = (
            (k_factor * self.wavelength) / (fwhm_sample_rad * cos_theta)
            if fwhm_sample_rad > 0
            else 0
        )
        size_nm = size_angstrom / 10

        # Check size limits (only if valid calculation)
        if not np.isnan(size_nm) and size_nm > 0:
            if size_nm > MAX_RELIABLE_SIZE:
                if validity_flag == ValidityFlag.VALID:
                    validity_flag = ValidityFlag.WARNING
                warnings.append(f"Size {size_nm:.1f} nm exceeds detection limit")
            elif size_nm < MIN_RELIABLE_SIZE:
                if validity_flag == ValidityFlag.VALID:
                    validity_flag = ValidityFlag.WARNING
                warnings.append(f"Size {size_nm:.1f} nm below precision limit")
        else:
            if is_reliable:  # If marked reliable but size is nan/0
                size_nm = 0
                validity_flag = ValidityFlag.ERROR
                warnings.append("Calculation resulted in invalid size")

        return ScherrerResult(
            size_nm=size_nm,
            size_angstrom=size_angstrom,
            two_theta=two_theta,
            hkl=hkl,
            k_factor=k_factor,
            fwhm_observed=fwhm_observed,
            fwhm_instrumental=fwhm_instrumental,
            fwhm_sample=fwhm_sample,
            fwhm_sample_rad=fwhm_sample_rad,
            validity_flag=validity_flag,
            warning_message="; ".join(warnings),
            is_reliable=is_reliable,
        )

    def _get_voigt_components(self, fwhm: float, eta: float) -> tuple[float, float]:
        """Convert Pseudo-Voigt (FWHM, η) to Gaussian and Lorentzian FWHM components.

        Uses TCH inverse mapping (Thompson-Cox-Hastings 1987) for accurate
        decomposition. The previous linear approximation (fL = η*FWHM,
        fG = (1-η)*FWHM) introduced systematic errors up to ~15% at
        intermediate η values.

        Reference:
            Thompson, Cox & Hastings (1987), J. Appl. Cryst. 20, 79-83.

        """
        return tch_components_from_eta(fwhm, eta)

    def _calculate_caglioti(self, two_theta: float) -> float:
        """Calculate instrumental FWHM using Caglioti equation.

        FWHM²_inst = U·tan²θ + V·tanθ + W
        """
        if not self.caglioti_params:
            return 0.0

        U, V, W = self.caglioti_params
        theta_rad = np.radians(two_theta / 2)
        tan_theta = np.tan(theta_rad)

        fwhm_sq = U * tan_theta**2 + V * tan_theta + W
        return np.sqrt(max(fwhm_sq, 0.0))

    def batch_calculate(
        self,
        peaks: list[tuple[float, float]],
        fwhm_instrumental: float | None = None,
    ) -> list[ScherrerResult]:
        """Calculate crystallite sizes for multiple peaks.

        Args:
            peaks: List of (two_theta, fwhm_observed) tuples.
            fwhm_instrumental: Common instrumental FWHM (optional).

        Returns:
            List of ScherrerResult objects.

        """
        return [
            self.calculate(two_theta, fwhm, fwhm_instrumental)
            for two_theta, fwhm in peaks
        ]

    def average_size(
        self, results: list[ScherrerResult], include_unreliable: bool = False
    ) -> tuple[float, float]:
        """Calculate average crystallite size from multiple peaks.

        Args:
            results: List of ScherrerResult objects.
            include_unreliable: Include UNRELIABLE results in average.

        Returns:
            Tuple of (average_size_nm, std_dev_nm).

        """
        sizes = [
            r.size_nm
            for r in results
            if (r.is_reliable or include_unreliable)
            and r.validity_flag != ValidityFlag.ERROR
        ]

        if not sizes:
            return 0.0, 0.0

        return float(np.mean(sizes)), float(np.std(sizes))


# =============================================================================
# Convenience Functions
# =============================================================================


def calculate_crystallite_size(
    two_theta: float,
    fwhm: float,
    wavelength: float = CU_KA1,
    k_factor: float = 0.829,  # L&W 1978 spherical standard
    fwhm_instrumental: float = 0.0,
) -> float:
    """Quick crystallite size calculation.

    Args:
        two_theta: Peak position (degrees).
        fwhm: FWHM in degrees.
        wavelength: X-ray wavelength (Å), default Cu Kα1.
        k_factor: Scherrer constant, default 0.829 (L&W 1978).
        fwhm_instrumental: Instrumental FWHM (optional).

    Returns:
        Crystallite size in nanometers.

    """
    calc = ScherrerCalculator(wavelength=wavelength, use_cubic_habit=False)
    result = calc.calculate(two_theta, fwhm, fwhm_instrumental)
    return result.size_nm


def calculate_scherrer(
    two_theta: float,
    fwhm_observed: float,
    fwhm_instrumental: float = 0.0,
    use_cubic_habit: bool = True,
    correction_method: str = "auto",
) -> ScherrerResult:
    """Convenience function for Scherrer calculation with full metadata.

    Example:
        >>> result = calculate_scherrer(43.32, 0.25, 0.08)
        >>> print(f"D = {result.size_nm:.1f} nm")
        D = 49.0 nm

    """
    calc = ScherrerCalculator(use_cubic_habit=use_cubic_habit)
    return calc.calculate(
        two_theta,
        fwhm_observed,
        fwhm_instrumental,
        correction_method=correction_method,
    )

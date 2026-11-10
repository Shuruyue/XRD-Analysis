"""Caglioti Equation Module.
==========================================
Implements instrumental broadening correction using Caglioti equation.

Reference:
    Caglioti, G., Paoletti, A., & Ricci, F. P. (1958).
    Nucl. Instr., 3, 223-228.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class CagliotiParams:
    """Caglioti equation parameters."""

    U: float  # deg²
    V: float  # deg²
    W: float  # deg²

    def __post_init__(self):
        if any(p is None for p in [self.U, self.V, self.W]):
            raise ValueError("All Caglioti parameters must be calibrated")


class CagliotiCorrection:
    """Instrumental broadening correction using Caglioti equation.

    FWHM²_inst = U·tan²θ + V·tanθ + W

    The Caglioti equation describes how instrumental broadening
    varies with diffraction angle 2θ.

    Reference: Caglioti, Paoletti, and Ricci (1958)
    """

    def __init__(
        self,
        U: float | None = None,
        V: float | None = None,
        W: float | None = None,
    ):
        """Initialize Caglioti correction.

        Args:
            U, V, W: Caglioti parameters (calibrated from standard)

        """
        self.U = U
        self.V = V
        self.W = W
        self._calibrated = all(p is not None for p in [U, V, W])

    @property
    def is_calibrated(self) -> bool:
        """Check if instrument has been calibrated."""
        return self._calibrated

    def calculate_fwhm_inst(self, two_theta: float) -> float:
        """Calculate instrumental FWHM at given 2θ.

        Args:
            two_theta: Diffraction angle in degrees

        Returns:
            Instrumental FWHM in degrees

        """
        if not self._calibrated:
            raise ValueError(
                "Instrument not calibrated. "
                "Run calibrate() with standard data first."
            )

        theta_rad = np.radians(two_theta / 2)
        tan_theta = np.tan(theta_rad)

        fwhm_sq = self.U * tan_theta**2 + self.V * tan_theta + self.W

        if fwhm_sq < 0:
            raise ValueError(f"Negative FWHM² at 2θ={two_theta}°. Check calibration.")

        return np.sqrt(fwhm_sq)

    def correct_broadening(
        self, fwhm_observed: float, two_theta: float, method: str = "geometric"
    ) -> tuple[float, bool, str]:
        """Correct observed FWHM for instrumental broadening.

        Two textbook limits are available:

        - "geometric": β_sample = β_obs − β_inst (linear; exact when both profiles are Lorentzian)
        - "quadratic": β_sample = √(β_obs² − β_inst²) (exact when both profiles are Gaussian)

        Real profiles lie between the two; ScherrerCalculator's "voigt" method separates the
        Gaussian and Lorentzian parts instead. Results with β_obs < 1.2 β_inst are flagged
        unreliable, because the sample width is then a small difference of two similar numbers.

        Note:
            Ignoring instrumental broadening underestimates the crystallite size, more so for
            large crystallites (D > 50 nm).

        Args:
            fwhm_observed: Observed FWHM (degrees)
            two_theta: Diffraction angle (degrees)
            method: Correction method ("geometric" or "quadratic")

        Returns:
            Tuple of (corrected_fwhm, is_reliable, warning_message). corrected_fwhm is NaN
            when the observed width is not larger than the instrument width, since no sample
            broadening can be resolved.

        """
        fwhm_inst = self.calculate_fwhm_inst(two_theta)

        # Check reliability threshold
        ratio = fwhm_observed / fwhm_inst
        is_reliable = ratio >= 1.2

        # Generate appropriate warning message
        if ratio < 1.0:
            warning = (
                f"CRITICAL: β_obs ({fwhm_observed:.4f}°) < β_inst ({fwhm_inst:.4f}°). "
                "Observed width smaller than instrument - check data or calibration."
            )
        elif ratio < 1.2:
            warning = (
                f"WARNING: β_obs/β_inst = {ratio:.2f} (< 1.2). "
                "Result may significantly underestimate crystallite size (D > 50 nm range). "
                "Consider this result unreliable."
            )
        elif ratio < 1.5:
            warning = (
                f"CAUTION: β_obs/β_inst = {ratio:.2f} (borderline). "
                "Result acceptable but approaching detection limit."
            )
        else:
            warning = ""

        if method == "quadratic":
            # Quadratic subtraction (exact for pure Gaussian convolution)
            # β_sample = √(β_obs² − β_inst²)
            # Reference: Keijser, Mittemeijer & Rozendaal (1983), Eq. 9
            fwhm_sq_diff = fwhm_observed**2 - fwhm_inst**2
            fwhm_sample = float(np.sqrt(fwhm_sq_diff)) if fwhm_sq_diff > 0 else np.nan
        elif method == "geometric":
            # Linear subtraction (exact for pure Lorentzian convolution)
            # β_sample = β_obs − β_inst
            # Reference: Scherrer (1918), original assumption
            fwhm_sample = fwhm_observed - fwhm_inst
            if fwhm_sample <= 0:
                fwhm_sample = np.nan
        else:
            raise ValueError(f"Unknown method: {method}")

        # No resolvable sample broadening: NaN, never a small placeholder width that would
        # turn into a large, finite crystallite size downstream.
        if np.isnan(fwhm_sample):
            is_reliable = False

        return fwhm_sample, is_reliable, warning

    def calibrate(
        self, two_theta_std: np.ndarray, fwhm_std: np.ndarray
    ) -> CagliotiParams:
        """Calibrate Caglioti parameters from standard sample data.

        Args:
            two_theta_std: 2θ positions of standard peaks (degrees)
            fwhm_std: FWHM of standard peaks (degrees)

        Returns:
            CagliotiParams with fitted U, V, W values

        """
        # Convert to radians and calculate tan(θ)
        theta_rad = np.radians(two_theta_std / 2)
        tan_theta = np.tan(theta_rad)

        # Fit: FWHM² = U·tan²θ + V·tanθ + W
        # Design matrix: [tan²θ, tanθ, 1]
        A = np.column_stack([tan_theta**2, tan_theta, np.ones_like(tan_theta)])

        fwhm_sq = fwhm_std**2

        # Least squares fit
        coeffs, residuals, rank, s = np.linalg.lstsq(A, fwhm_sq, rcond=None)

        self.U, self.V, self.W = coeffs
        self._calibrated = True

        return CagliotiParams(U=self.U, V=self.V, W=self.W)

    def get_params(self) -> CagliotiParams:
        """Return current Caglioti parameters."""
        return CagliotiParams(U=self.U, V=self.V, W=self.W)


def calculate_instrumental_broadening(
    two_theta: float, U: float, V: float, W: float
) -> float:
    """Convenience function to calculate instrumental FWHM.

    Args:
        two_theta: Diffraction angle (degrees)
        U, V, W: Caglioti parameters

    Returns:
        Instrumental FWHM (degrees)

    """
    correction = CagliotiCorrection(U, V, W)
    return correction.calculate_fwhm_inst(two_theta)

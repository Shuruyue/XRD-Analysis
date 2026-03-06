"""Pseudo-Voigt Function Module.
=============================================
Implements the Pseudo-Voigt profile function for XRD peak fitting.
XRDVoigt

Also includes True Voigt profile using scipy.special.voigt_profile.
scipy.special.voigt_profileVoigt
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import voigt_profile


@dataclass
class VoigtParams:
    """Parameters for a true Voigt peak."""

    center: float  # 2θ center position (degrees)
    amplitude: float  # Peak amplitude (intensity)
    sigma: float  # Gaussian width parameter
    gamma: float  # Lorentzian width parameter (half-width)

    @property
    def fwhm_gaussian(self) -> float:
        """FWHM of Gaussian component."""
        return 2.0 * self.sigma * np.sqrt(2.0 * np.log(2.0))

    @property
    def fwhm_lorentzian(self) -> float:
        """FWHM of Lorentzian component."""
        return 2.0 * self.gamma

    @property
    def fwhm_total(self) -> float:
        """Total FWHM using Thompson-Cox-Hastings (1987) 5th-order polynomial.

        Reference:
            Thompson, Cox & Hastings (1987),
            J. Appl. Cryst. 20, 79-83. DOI: 10.1107/S0021889887087090
        """
        fG = self.fwhm_gaussian
        fL = self.fwhm_lorentzian
        return tch_fwhm_from_components(fG, fL)

    def to_array(self) -> np.ndarray:
        """Convert to numpy array [center, amplitude, sigma, gamma]."""
        return np.array([self.center, self.amplitude, self.sigma, self.gamma])

    @classmethod
    def from_array(cls, arr: np.ndarray) -> "VoigtParams":
        """Create from numpy array."""
        return cls(center=arr[0], amplitude=arr[1], sigma=arr[2], gamma=arr[3])


class TrueVoigt:
    """True Voigt profile function (convolution of Gaussian and Lorentzian).

    Mathematical Definition (Faddeeva Function):
    --------------------------------------------
    V(x; σ, γ) = Re[w(z)] / (σ√(2π))

    where w(z) is the Faddeeva function (scaled complex error function):
        w(z) = exp(-z²) * erfc(-iz)

    and z is the complex argument:
        z = (x + iγ) / (σ√2)

    Physical Interpretation in XRD:
    -------------------------------
    - x (Real part): Deviation from Bragg angle (2θ - 2θ₀).
      Analogous to "Frequency Drift" in spectroscopy.
    - γ (Imaginary part): Lorentzian HWHM (Size Broadening + Lifetime).
      Analogous to "Damping Coefficient" in spectroscopy.
    - σ (Normalization): Gaussian Width (Strain + Instrument).

    Note on Instrument Parameters:
    ------------------------------
    This function fits the *Total* profile (Sample ⊗ Instrument).
    You do NOT need instrument parameters to perform this fit. The optimizer
    will find the effective σ_total and γ_total from the raw data.
    Instrumental correction is performed *after* fitting during the
    analysis stage (e.g., in Williamson-Hall or Scherrer methods).

    References
    ----------
    1. Armstrong, B. H. (1967).
       "Spectrum Line Profiles: The Voigt Function".
       J. Quant. Spectrosc. Radiat. Transfer, 7, 61-88.

    2. Poppe, G. P. M., & Wijers, C. M. J. (1990).
       "More efficient computation of the complex error function".
       ACM Trans. Math. Softw., 16, 38-46.
       (Standard algorithm used by scipy.special.wofz)

    """

    @staticmethod
    def profile(
        x: np.ndarray, center: float, amplitude: float, sigma: float, gamma: float
    ) -> np.ndarray:
        """Calculate true Voigt profile.

        Args:
            x: 2θ array
            center: Peak center position
            amplitude: Peak amplitude
            sigma: Gaussian width parameter (σ)
            gamma: Lorentzian half-width parameter (γ)

        Returns:
            Intensity array

        """
        # Shift x to be centered at 0
        x_shifted = x - center

        # Use scipy's voigt_profile (normalized to unit area)
        # voigt_profile(x, sigma, gamma) computes the Voigt function
        v = voigt_profile(x_shifted, sigma, gamma)

        # Normalize and scale by amplitude
        # The voigt_profile is normalized, so we scale to match peak amplitude
        v_max = voigt_profile(0, sigma, gamma)
        if v_max > 0:
            return amplitude * (v / v_max)
        else:
            return np.zeros_like(x)

    @staticmethod
    def fwhm_from_params(sigma: float, gamma: float) -> float:
        """Calculate total FWHM from sigma and gamma.

        Uses Thompson-Cox-Hastings (1987) 5th-order polynomial.

        Reference:
            Thompson, Cox & Hastings (1987),
            J. Appl. Cryst. 20, 79-83.
        """
        fG = 2.0 * sigma * np.sqrt(2.0 * np.log(2.0))
        fL = 2.0 * gamma
        return tch_fwhm_from_components(fG, fL)

    @staticmethod
    def params_from_fwhm(fwhm: float, eta: float = 0.5) -> tuple[float, float]:
        """Estimate sigma and gamma from total FWHM and mixing ratio.

        Args:
            fwhm: Total FWHM
            eta: Lorentzian fraction (0 = Gaussian, 1 = Lorentzian)

        Returns:
            (sigma, gamma) tuple

        """
        # Approximate split between Gaussian and Lorentzian FWHM
        fL = eta * fwhm
        fG = (1 - eta) * fwhm

        sigma = fG / (2.0 * np.sqrt(2.0 * np.log(2.0))) if fG > 0 else 0.01
        gamma = fL / 2.0 if fL > 0 else 0.01

        return sigma, gamma


@dataclass
class PseudoVoigtParams:
    """Parameters for a Pseudo-Voigt peak."""

    center: float  # 2θ center position (degrees)
    amplitude: float  # Peak amplitude (intensity)
    fwhm: float  # Full Width at Half Maximum (degrees)
    eta: float  # Mixing parameter (0 = Gaussian, 1 = Lorentzian)

    def to_array(self) -> np.ndarray:
        """Convert to numpy array [center, amplitude, fwhm, eta]."""
        return np.array([self.center, self.amplitude, self.fwhm, self.eta])

    @classmethod
    def from_array(cls, arr: np.ndarray) -> "PseudoVoigtParams":
        """Create from numpy array."""
        return cls(center=arr[0], amplitude=arr[1], fwhm=arr[2], eta=arr[3])


class PseudoVoigt:
    """Pseudo-Voigt profile function.

    The Pseudo-Voigt is a linear combination of Gaussian and Lorentzian:

    I(2θ) = I₀ × [η × L(2θ) + (1-η) × G(2θ)]

    where:
    - L(2θ): Lorentzian component (size broadening)
    - G(2θ): Gaussian component (strain + instrumental)
    - η: Mixing parameter (0 < η < 1)
    """

    @staticmethod
    def gaussian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
        """Normalized Gaussian function.

        G(x) = exp(-4 ln(2) × ((x - center) / fwhm)²)
        """
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        return np.exp(-0.5 * ((x - center) / sigma) ** 2)

    @staticmethod
    def lorentzian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
        """Normalized Lorentzian function.

        L(x) = 1 / (1 + 4 × ((x - center) / fwhm)²)
        """
        gamma = fwhm / 2
        return 1 / (1 + ((x - center) / gamma) ** 2)

    @classmethod
    def profile(
        cls, x: np.ndarray, center: float, amplitude: float, fwhm: float, eta: float
    ) -> np.ndarray:
        """Calculate Pseudo-Voigt profile.

        Args:
            x: 2θ array
            center: Peak center position
            amplitude: Peak amplitude
            fwhm: Full Width at Half Maximum
            eta: Mixing parameter (0 = pure Gaussian, 1 = pure Lorentzian)

        Returns:
            Intensity array

        """
        # Ensure eta is bounded
        eta = np.clip(eta, 0, 1)

        gaussian = cls.gaussian(x, center, fwhm)
        lorentzian = cls.lorentzian(x, center, fwhm)

        return amplitude * (eta * lorentzian + (1 - eta) * gaussian)

    @classmethod
    def multi_peak(
        cls, x: np.ndarray, params_list: list, background: float = 0
    ) -> np.ndarray:
        """Calculate sum of multiple Pseudo-Voigt peaks.

        Args:
            x: 2θ array
            params_list: List of PseudoVoigtParams or parameter arrays
            background: Constant background level

        Returns:
            Total intensity array

        """
        result = np.full_like(x, background, dtype=float)

        for params in params_list:
            if isinstance(params, PseudoVoigtParams):
                center, amplitude, fwhm, eta = (
                    params.center,
                    params.amplitude,
                    params.fwhm,
                    params.eta,
                )
            else:
                center, amplitude, fwhm, eta = params[:4]

            result += cls.profile(x, center, amplitude, fwhm, eta)

        return result


def pseudo_voigt_function(
    x: np.ndarray, center: float, amplitude: float, fwhm: float, eta: float
) -> np.ndarray:
    """Convenience function for Pseudo-Voigt calculation.

    Args:
        x: 2θ array
        center: Peak center (degrees)
        amplitude: Peak amplitude
        fwhm: Full Width at Half Maximum (degrees)
        eta: Mixing parameter (0-1)

    Returns:
        Intensity array

    """
    return PseudoVoigt.profile(x, center, amplitude, fwhm, eta)


# =============================================================================
# Thompson-Cox-Hastings (1987) Parameterization
# =============================================================================
#
# Reference:
#     Thompson, P., Cox, D. E., & Hastings, J. B. (1987).
#     "Rietveld refinement of Debye-Scherrer synchrotron X-ray data from Al₂O₃."
#     J. Appl. Cryst. 20, 79-83.
#     DOI: 10.1107/S0021889887087090
#
# The TCH parameterization gives a more accurate total Voigt FWHM and
# pseudo-Voigt mixing parameter η from Gaussian and Lorentzian components
# than the Olivero-Longbothum (1977) approximation.
#
# Used by FullProf, GSAS-II, and other modern Rietveld codes.
# =============================================================================


def tch_fwhm_from_components(fG: float, fL: float) -> float:
    """Compute total Voigt FWHM from Gaussian and Lorentzian FWHM components.

    Thompson-Cox-Hastings (1987) 5th-order polynomial:
        f⁵ = fG⁵ + 2.69269·fG⁴·fL + 2.42843·fG³·fL²
             + 4.47163·fG²·fL³ + 0.07842·fG·fL⁴ + fL⁵

    Args:
        fG: Gaussian FWHM component
        fL: Lorentzian FWHM component

    Returns:
        Total Voigt FWHM

    """
    f5 = (
        fG**5
        + 2.69269 * fG**4 * fL
        + 2.42843 * fG**3 * fL**2
        + 4.47163 * fG**2 * fL**3
        + 0.07842 * fG * fL**4
        + fL**5
    )
    return f5**0.2 if f5 > 0 else 0.0


def tch_eta_from_components(fG: float, fL: float) -> tuple[float, float]:
    """Compute pseudo-Voigt η and total FWHM from G/L components via TCH.

    Thompson-Cox-Hastings (1987) parameterization.

    The η formula:
        r = fL / f_total
        η = 1.36603·r − 0.47719·r² + 0.11116·r³

    This gives the Lorentzian fraction of the pseudo-Voigt that best
    approximates the true Voigt convolution.

    Args:
        fG: Gaussian FWHM component
        fL: Lorentzian FWHM component

    Returns:
        Tuple of (total_fwhm, eta)

    """
    f_total = tch_fwhm_from_components(fG, fL)
    if f_total <= 0:
        return 0.0, 0.5

    r = fL / f_total
    eta = 1.36603 * r - 0.47719 * r**2 + 0.11116 * r**3
    eta = float(np.clip(eta, 0.0, 1.0))

    return f_total, eta


def tch_components_from_eta(fwhm_total: float, eta: float) -> tuple[float, float]:
    """Inverse TCH mapping: recover (fG, fL) from total FWHM and η.

    Given a pseudo-Voigt described by (FWHM, η), compute the Gaussian and
    Lorentzian FWHM components that reproduce these values through the
    Thompson-Cox-Hastings parameterization.

    Algorithm:
        1. Invert η(r) = 1.36603r − 0.47719r² + 0.11116r³ via Newton's method
           to obtain r = fL / fwhm_total.
        2. fL = r × fwhm_total.
        3. Solve tch_fwhm_from_components(fG, fL) = fwhm_total for fG using
           bisection (the TCH polynomial is monotonically increasing in fG).

    Args:
        fwhm_total: Total pseudo-Voigt FWHM
        eta: Pseudo-Voigt mixing parameter (0 = Gaussian, 1 = Lorentzian)

    Returns:
        (fG, fL) — Gaussian and Lorentzian FWHM components

    Reference:
        Thompson, Cox & Hastings (1987), J. Appl. Cryst. 20, 79-83.

    """
    if fwhm_total <= 0:
        return 0.0, 0.0

    eta = float(np.clip(eta, 0.0, 1.0))

    # Step 1: Invert η(r) using Newton's method
    # η = 1.36603*r - 0.47719*r² + 0.11116*r³
    # η'(r) = 1.36603 - 0.95438*r + 0.33348*r²
    r = eta  # Initial guess (good for small η)
    for _ in range(20):
        f_r = 1.36603 * r - 0.47719 * r**2 + 0.11116 * r**3 - eta
        df_r = 1.36603 - 0.95438 * r + 0.33348 * r**2
        if abs(df_r) < 1e-15:
            break
        r_new = r - f_r / df_r
        r_new = max(0.0, min(1.0, r_new))
        if abs(r_new - r) < 1e-12:
            break
        r = r_new

    r = float(np.clip(r, 0.0, 1.0))
    fL = r * fwhm_total

    # Step 2: Solve for fG via bisection on tch_fwhm_from_components(fG, fL) = fwhm_total
    fG_lo, fG_hi = 0.0, fwhm_total
    for _ in range(60):
        fG_mid = (fG_lo + fG_hi) / 2.0
        f_mid = tch_fwhm_from_components(fG_mid, fL)
        if f_mid < fwhm_total:
            fG_lo = fG_mid
        else:
            fG_hi = fG_mid
        if (fG_hi - fG_lo) < 1e-14 * fwhm_total:
            break

    fG = (fG_lo + fG_hi) / 2.0
    return fG, fL

"""Williamson-Hall Analysis Module.
===============================

Separates crystallite size and microstrain contributions to peak broadening.

Reference:
    Williamson, G. K., & Hall, W. H. (1953).
    X-ray line broadening from filed aluminium and wolfram.
    Acta Metallurgica, 1(1), 22-31.
"""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from scipy.stats import linregress

from xrd_analysis.core.constants import CU_KA1

# =============================================================================
# Constants
# =============================================================================

# Wavelength: use CU_KA1 directly

# Williamson-Hall K factor
# ═══════════════════════════════════════════════════════════════════════════
# Reference:
#     Williamson, G. K., & Hall, W. H. (1953).
#     "X-ray line broadening from filed aluminium and wolfram."
#     Acta Metallurgica, 1(1), 22-31.
#     DOI: 10.1016/0001-6160(53)90006-6
#
# Note:
#     K ≈ 0.9 is commonly used average for Scherrer/W-H analysis,
#     depending on grain shape and definition (FWHM or integral breadth)
# ═══════════════════════════════════════════════════════════════════════════
WH_K_FACTOR = 0.9

# R² quality thresholds
# ═══════════════════════════════════════════════════════════════════════════
# USER-DEFINED PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════
# These thresholds are user-defined for linear regression quality assessment.
# Can be adjusted based on requirements (e.g., 0.99 for stricter criteria).
# ═══════════════════════════════════════════════════════════════════════════
R2_EXCELLENT = 0.95
R2_ACCEPTABLE = 0.85
MIN_PEAKS = 3  # Minimum 3 peaks for linear regression

# Copper elastic anisotropy
MODULUS_MAP: dict[tuple[int, int, int], float] = {
    (1, 1, 1): 191.0,  # GPa, hardest direction
    (2, 0, 0): 67.0,  # GPa, softest direction
    (2, 2, 0): 130.0,  # GPa, intermediate
    (3, 1, 1): 128.0,  # GPa
}

ZENER_ANISOTROPY = 3.21

# =============================================================================
# Quality Level Enum
# =============================================================================


class WHQualityLevel(Enum):
    """Williamson-Hall analysis quality classification."""

    EXCELLENT = "excellent"  # R² > 0.95
    ACCEPTABLE = "acceptable"  # 0.85 ≤ R² ≤ 0.95
    POOR = "poor"  # R² < 0.85


# =============================================================================
# Result Dataclass
# =============================================================================


@dataclass
class WHResult:
    """Result from Williamson-Hall analysis.

    Attributes:
        crystallite_size_nm: Crystallite size in nm.
        microstrain: Dimensionless microstrain ε.
        r_squared: Goodness of fit R².
        intercept: y-intercept (Kλ/D).
        slope: Slope (4ε).
        intercept_stderr: Standard error of intercept.
        slope_stderr: Standard error of slope.
        size_error_nm: Error in size estimate.
        strain_error: Error in strain estimate.
        quality_level: Quality classification.
        is_reliable: Whether fit is reliable.
        n_peaks: Number of peaks used.
        peak_hkls: List of hkl indices.
        warning_message: Warning message.
        anisotropy_note: Anisotropy diagnostic.

    """

    crystallite_size_nm: float
    microstrain: float
    r_squared: float
    intercept: float
    slope: float
    intercept_stderr: float = 0.0
    slope_stderr: float = 0.0
    size_error_nm: float = 0.0
    strain_error: float = 0.0
    quality_level: WHQualityLevel = WHQualityLevel.ACCEPTABLE
    is_reliable: bool = True
    n_peaks: int = 0
    peak_hkls: list[tuple[int, int, int]] = field(default_factory=list)
    warning_message: str = ""
    anisotropy_note: str = ""

    def __repr__(self) -> str:
        return (
            f"WHResult(D={self.crystallite_size_nm:.1f} nm, "
            f"ε={self.microstrain:.2e}, R²={self.r_squared:.3f} "
            f"[{self.quality_level.value}])"
        )


# =============================================================================
# Williamson-Hall Analyzer
# =============================================================================


class WilliamsonHallAnalyzer:
    """Williamson-Hall analysis for separating size and strain broadening.

    W-H Equation:
        β cos θ = (K λ / D) + 4 ε sin θ

    where:
        β: Total peak broadening (FWHM in radians)
        θ: Bragg angle
        K: Scherrer constant
        λ: X-ray wavelength
        D: Crystallite size
        ε: Microstrain

    Features:
        - R² quality assessment with thresholds
        - Anisotropy diagnostics using MODULUS_MAP
        - Error propagation for D and ε

    Args:
        wavelength: X-ray wavelength in Å (default: Cu Kα1)
        k_factor: Scherrer constant for W-H (default: 0.9)

    Example:
        >>> analyzer = WilliamsonHallAnalyzer()
        >>> result = analyzer.analyze(two_theta, fwhm_sample)
        >>> print(f"D = {result.crystallite_size_nm:.1f} nm")

    """

    def __init__(
        self, wavelength: float = CU_KA1, k_factor: float = WH_K_FACTOR
    ) -> None:
        self.wavelength = wavelength
        self.k_factor = k_factor

    def analyze(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        hkl_list: list[tuple[int, int, int] | None] = None,
        fwhm_in_radians: bool = False,
    ) -> WHResult:
        """Perform Williamson-Hall analysis.

        Args:
            two_theta: Array of 2θ peak positions (degrees).
            fwhm_sample: Array of sample FWHM values.
            hkl_list: Optional list of (h,k,l) for each peak.
            fwhm_in_radians: If True, FWHM is in radians.

        Returns:
            WHResult with size, strain, and quality metrics.

        """
        two_theta = np.asarray(two_theta)
        fwhm_sample = np.asarray(fwhm_sample)
        n_peaks = len(two_theta)

        # Validate input
        if n_peaks < MIN_PEAKS:
            return self._create_failed_result(
                f"Need at least {MIN_PEAKS} peaks for W-H analysis", n_peaks
            )

        if len(fwhm_sample) != n_peaks:
            return self._create_failed_result(
                "2θ and FWHM arrays length mismatch", n_peaks
            )

        # Convert units
        theta_rad = np.radians(two_theta / 2.0)
        beta_rad = fwhm_sample if fwhm_in_radians else np.radians(fwhm_sample)

        # W-H coordinates: X = sin(θ), Y = β × cos(θ)
        x_data = np.sin(theta_rad)
        y_data = beta_rad * np.cos(theta_rad)

        # Linear regression
        reg_result = linregress(x_data, y_data)
        slope = reg_result.slope
        intercept = reg_result.intercept
        r_squared = reg_result.rvalue**2
        slope_stderr = reg_result.stderr
        intercept_stderr = getattr(reg_result, "intercept_stderr", 0.0)

        # Calculate physical quantities: D = Kλ / intercept, ε = slope / 4
        if intercept > 0:
            size_angstrom = self.k_factor * self.wavelength / intercept
            size_nm = size_angstrom / 10.0
            size_error_nm = (
                (size_angstrom / 10.0) * (intercept_stderr / intercept)
                if intercept_stderr > 0
                else 0.0
            )
        else:
            size_nm = float("inf")
            size_error_nm = 0.0

        microstrain = slope / 4.0
        strain_error = slope_stderr / 4.0 if slope_stderr else 0.0

        # Quality assessment
        quality_level, warning = self._assess_quality(r_squared)

        # Anisotropy diagnostics
        anisotropy_note = self._generate_anisotropy_note(r_squared, hkl_list)

        return WHResult(
            crystallite_size_nm=size_nm,
            microstrain=microstrain,
            r_squared=r_squared,
            intercept=intercept,
            slope=slope,
            intercept_stderr=intercept_stderr,
            slope_stderr=slope_stderr,
            size_error_nm=size_error_nm,
            strain_error=strain_error,
            quality_level=quality_level,
            is_reliable=(quality_level != WHQualityLevel.POOR),
            n_peaks=n_peaks,
            peak_hkls=hkl_list or [],
            warning_message=warning,
            anisotropy_note=anisotropy_note,
        )

    def analyze_with_correction(
        self,
        two_theta: np.ndarray,
        fwhm_observed: np.ndarray,
        fwhm_instrumental: np.ndarray,
    ) -> WHResult:
        """Perform W-H analysis with instrumental broadening correction.

        Uses quadratic (Gaussian) subtraction: β_sample = √(β_obs² − β_inst²)
        This is the physically correct formula for Gaussian convolution:
            G_obs = G_sample ⊗ G_inst  ⟹  FWHM_obs² = FWHM_sample² + FWHM_inst²

        Reference:
            Keijser, Mittemeijer & Rozendaal (1983), J. Appl. Cryst. 16, 309-316.
            Eq. 9: β_G_structure = sqrt(β_G_profile² - β_G_standard²)
        """
        fwhm_obs = np.asarray(fwhm_observed, dtype=float)
        fwhm_inst = np.asarray(fwhm_instrumental, dtype=float)
        # Quadratic Gaussian subtraction (exact for Gaussian profiles)
        fwhm_sq_diff = fwhm_obs**2 - fwhm_inst**2
        fwhm_corrected = np.sqrt(np.maximum(fwhm_sq_diff, 1e-8))
        return self.analyze(two_theta, fwhm_corrected)

    def get_plot_data(
        self,
        two_theta: np.ndarray,
        fwhm_sample: np.ndarray,
        fwhm_in_radians: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get data for W-H plot.

        Returns:
            Tuple of (x_data, y_data, x_fit_line, y_fit_line).

        """
        theta_rad = np.radians(np.asarray(two_theta) / 2.0)
        beta_rad = (
            fwhm_sample if fwhm_in_radians else np.radians(np.asarray(fwhm_sample))
        )

        x_data = np.sin(theta_rad)
        y_data = beta_rad * np.cos(theta_rad)

        reg = linregress(x_data, y_data)
        x_fit = np.linspace(x_data.min() * 0.9, x_data.max() * 1.1, 100)
        y_fit = reg.intercept + reg.slope * x_fit

        return x_data, y_data, x_fit, y_fit

    def _assess_quality(self, r_squared: float) -> tuple[WHQualityLevel, str]:
        """Assess analysis quality based on R²."""
        if r_squared > R2_EXCELLENT:
            return WHQualityLevel.EXCELLENT, ""
        elif r_squared > R2_ACCEPTABLE:
            return WHQualityLevel.ACCEPTABLE, f"R² = {r_squared:.3f} (acceptable)"
        else:
            return WHQualityLevel.POOR, (
                f"R² = {r_squared:.3f} < {R2_ACCEPTABLE} - "
                "significant anisotropy detected"
            )

    def _generate_anisotropy_note(
        self, r_squared: float, hkl_list: list[tuple[int, int, int] | None]
    ) -> str:
        """Generate anisotropy diagnostic note."""
        if r_squared > R2_ACCEPTABLE:
            return ""

        lines = [
            "Elastic Anisotropy Diagnostic",
            f"Cu Zener ratio A = {ZENER_ANISOTROPY} (extreme anisotropy)",
            "",
            "Direction-dependent Young's modulus:",
        ]

        for hkl, E in sorted(MODULUS_MAP.items(), key=lambda x: -x[1]):
            hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
            lines.append(f"  E{hkl_str} = {E:.0f} GPa")

        lines.extend(
            [
                "",
                f"E(111)/E(200) = {MODULUS_MAP[(1,1,1)]/MODULUS_MAP[(2,0,0)]:.1f}",
                "  → (200) strain broadening may be overestimated ~3x",
            ]
        )

        # Report which detected peaks have strongly divergent moduli
        if hkl_list:
            detected_moduli = [
                (hkl, MODULUS_MAP[hkl]) for hkl in hkl_list if hkl in MODULUS_MAP
            ]
            if len(detected_moduli) >= 2:
                lines.append("")
                lines.append("Detected peaks with divergent moduli:")
                E_vals = [E for _, E in detected_moduli]
                E_range = max(E_vals) - min(E_vals)
                for hkl, E in sorted(detected_moduli, key=lambda x: -x[1]):
                    hkl_str = f"({hkl[0]}{hkl[1]}{hkl[2]})"
                    lines.append(f"  {hkl_str}: E = {E:.0f} GPa")
                lines.append(f"  E-range = {E_range:.0f} GPa → consider Modified W-H")

        return "\n".join(lines)

    def _create_failed_result(self, message: str, n_peaks: int) -> WHResult:
        """Create a failed result object."""
        return WHResult(
            crystallite_size_nm=0.0,
            microstrain=0.0,
            r_squared=0.0,
            intercept=0.0,
            slope=0.0,
            quality_level=WHQualityLevel.POOR,
            is_reliable=False,
            n_peaks=n_peaks,
            warning_message=message,
        )


# =============================================================================
# Convenience Functions
# =============================================================================


def analyze_williamson_hall(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    hkl_list: list[tuple[int, int, int] | None] = None,
) -> WHResult:
    """Convenience function for Williamson-Hall analysis.

    Example:
        >>> two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        >>> fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        >>> result = analyze_williamson_hall(two_theta, fwhm)
        >>> print(f"D = {result.crystallite_size_nm:.1f} nm")

    """
    analyzer = WilliamsonHallAnalyzer()
    return analyzer.analyze(two_theta, fwhm_sample, hkl_list)


def get_modulus_for_hkl(hkl: tuple[int, int, int]) -> float:
    """Get Young's modulus for a given hkl direction.

    Returns:
        Young's modulus in GPa, or 120 GPa (average) if not found.

    """
    return MODULUS_MAP.get(hkl, 120.0)


# =============================================================================
# Modified Williamson-Hall (MWH) Analysis
# =============================================================================
#
# Reference:
#     Ungár, T., & Borbély, A. (1996).
#     "The effect of dislocation contrast on x-ray line broadening:
#     A new approach to line profile analysis."
#     Appl. Phys. Lett. 69(21), 3173-3175.
#     DOI: 10.1063/1.117951
#
# For elastically anisotropic materials (e.g. Cu, Zener ratio A=3.21),
# classical W-H gives poor R² because strain broadening depends on the
# dislocation contrast factor C̄_hkl for each reflection. MWH accounts
# for this by plotting ΔK vs K√C̄ instead of β cos θ vs sin θ.
# =============================================================================

# Default contrast factor parameters for FCC Cu
# C̄_hkl = C̄_h00 × (1 − q × H²)
# where H² = (h²k² + k²l² + l²h²) / (h² + k² + l²)²
#
# The q parameter depends on the dislocation character (edge vs screw)
# and the elastic constants.
# For Cu FCC edge dislocations: q ≈ 1.72
# For Cu FCC screw dislocations: q ≈ 2.37
# Average: q ≈ 1.9 (mixed character)
#
# C̄_h00 depends on elastic constants:
# Edge: C̄_h00 ≈ 0.307, Screw: C̄_h00 ≈ 0.279
# Average: C̄_h00 ≈ 0.293
MWH_C_H00_DEFAULT = 0.293  # Average contrast factor for <h00>
MWH_Q_DEFAULT = 1.9  # Average q for mixed dislocations


@dataclass
class MWHResult:
    """Modified Williamson-Hall analysis result.

    Attributes:
        crystallite_size_nm: Crystallite size from intercept
        dislocation_density_m2: Dislocation density (m⁻²)
        r_squared: R² of the MWH linear fit
        slope: Slope of ΔK vs K√C̄ plot
        intercept: Intercept (= 0.9/D)
        quality_level: Fit quality classification
        is_reliable: Whether the result is trustworthy
        n_peaks: Number of peaks used
        contrast_factors: C̄_hkl for each peak
        q_parameter: The q parameter used
        warning_message: Any diagnostic messages

    """

    crystallite_size_nm: float
    dislocation_density_m2: float
    r_squared: float
    slope: float
    intercept: float
    quality_level: WHQualityLevel
    is_reliable: bool
    n_peaks: int
    contrast_factors: dict[tuple[int, int, int], float] = field(default_factory=dict)
    q_parameter: float = MWH_Q_DEFAULT
    warning_message: str = ""


def compute_H_squared(hkl: tuple[int, int, int]) -> float:
    """Compute anisotropy factor H² for a given (hkl).

    H² = (h²k² + k²l² + l²h²) / (h² + k² + l²)²

    This factor equals 0 for <100>, 1/4 for <110>, and 1/3 for <111>.

    """
    h, k, l = hkl
    h2, k2, l2 = h**2, k**2, l**2
    numerator = h2 * k2 + k2 * l2 + l2 * h2
    denominator = (h2 + k2 + l2) ** 2
    return numerator / denominator if denominator > 0 else 0.0


def compute_contrast_factor(
    hkl: tuple[int, int, int],
    C_h00: float = MWH_C_H00_DEFAULT,
    q: float = MWH_Q_DEFAULT,
) -> float:
    """Compute average dislocation contrast factor C̄_hkl.

    C̄_hkl = C̄_h00 × (1 − q × H²)

    Reference:
        Ungár & Borbély (1996), Appl. Phys. Lett. 69, 3173.

    Args:
        hkl: Miller indices
        C_h00: Average contrast factor for <h00> reflections
        q: Anisotropy parameter (depends on dislocation character)

    Returns:
        Average dislocation contrast factor (always >= 0)

    """
    H2 = compute_H_squared(hkl)
    C_bar = C_h00 * (1.0 - q * H2)
    return max(C_bar, 1e-6)  # Ensure positive


def analyze_modified_wh(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    hkl_list: list[tuple[int, int, int]],
    wavelength: float = CU_KA1,
    C_h00: float = MWH_C_H00_DEFAULT,
    q: float = MWH_Q_DEFAULT,
) -> MWHResult:
    """Perform Modified Williamson-Hall analysis.

    Plots ΔK vs K√C̄_hkl where:
        K = 2 sin θ / λ  (reciprocal space)
        ΔK = β cos θ / λ  (broadening in reciprocal space)
        C̄_hkl = dislocation contrast factor

    Points should collapse onto a straight line even for anisotropic metals.

    Intercept → 0.9/D (crystallite size)
    Slope → proportional to √ρ (dislocation density)

    Args:
        two_theta: Array of 2θ positions (degrees)
        fwhm_sample: Array of sample FWHM values (degrees, instrument-corrected)
        hkl_list: List of (h,k,l) Miller indices for each peak
        wavelength: X-ray wavelength (Å)
        C_h00: Average contrast factor for <h00>
        q: Anisotropy parameter

    Returns:
        MWHResult with crystallite size and dislocation density

    Reference:
        Ungár, T. & Borbély, A. (1996). Appl. Phys. Lett. 69, 3173-3175.

    """
    two_theta_arr = np.asarray(two_theta, dtype=float)
    fwhm_arr = np.asarray(fwhm_sample, dtype=float)

    n_peaks = len(two_theta_arr)
    if n_peaks < MIN_PEAKS:
        return MWHResult(
            crystallite_size_nm=0.0,
            dislocation_density_m2=0.0,
            r_squared=0.0,
            slope=0.0,
            intercept=0.0,
            quality_level=WHQualityLevel.POOR,
            is_reliable=False,
            n_peaks=n_peaks,
            warning_message=f"Need ≥ {MIN_PEAKS} peaks, got {n_peaks}",
        )

    # Convert to reciprocal space coordinates
    theta_rad = np.radians(two_theta_arr / 2.0)
    fwhm_rad = np.radians(fwhm_arr)

    K = 2.0 * np.sin(theta_rad) / wavelength  # 1/Å
    delta_K = fwhm_rad * np.cos(theta_rad) / wavelength  # ΔK (1/Å)

    # Compute contrast factors and MWH x-axis
    contrast_factors = {}
    x_mwh = np.zeros(n_peaks)
    for i, hkl in enumerate(hkl_list):
        C_bar = compute_contrast_factor(hkl, C_h00, q)
        contrast_factors[hkl] = C_bar
        x_mwh[i] = K[i] * np.sqrt(C_bar)  # K√C̄

    y_mwh = delta_K  # ΔK

    # Linear regression: ΔK = intercept + slope × K√C̄
    reg = linregress(x_mwh, y_mwh)
    slope = float(reg.slope)
    intercept = float(reg.intercept)
    r_squared = float(reg.rvalue**2)

    # Extract crystallite size from intercept: intercept = 0.9/D
    # D is in Å (same units as wavelength)
    if intercept > 0:
        D_angstrom = 0.9 / intercept
        D_nm = D_angstrom / 10.0
    else:
        D_nm = 0.0

    # Extract dislocation density from slope
    # slope ≈ (π M² b² / 2)^(1/2) × √ρ
    # For FCC Cu: b = a₀/√2 ≈ 2.556 Å, M ≈ 2 (Wilkens parameter)
    b_angstrom = 3.615 / np.sqrt(2.0)  # Burgers vector for Cu FCC (Å)
    M = 2.0  # Wilkens cut-off parameter (typical)
    prefactor = np.sqrt(np.pi * M**2 * b_angstrom**2 / 2.0)
    if prefactor > 0 and slope > 0:
        rho = (slope / prefactor) ** 2  # Å⁻²
        rho_m2 = rho * 1e20  # Convert from Å⁻² to m⁻²
    else:
        rho_m2 = 0.0

    # Quality assessment
    if r_squared >= R2_EXCELLENT:
        quality = WHQualityLevel.EXCELLENT
    elif r_squared >= R2_ACCEPTABLE:
        quality = WHQualityLevel.ACCEPTABLE
    else:
        quality = WHQualityLevel.POOR

    is_reliable = r_squared >= R2_ACCEPTABLE and D_nm > 0

    return MWHResult(
        crystallite_size_nm=D_nm,
        dislocation_density_m2=rho_m2,
        r_squared=r_squared,
        slope=slope,
        intercept=intercept,
        quality_level=quality,
        is_reliable=is_reliable,
        n_peaks=n_peaks,
        contrast_factors=contrast_factors,
        q_parameter=q,
    )


# =============================================================================
# Size-Strain Plot (SSP) / Halder-Wagner Method
# =============================================================================
#
# Reference:
#     Halder, N. C., & Wagner, C. N. J. (1966).
#     "Separation of particle size and lattice strain in integral breadth
#      measurements."
#     Acta Crystallographica, 20(2), 312-313.
#     DOI: 10.1107/S0365110X66000628
#
# The SSP method assumes a Voigt profile (Gaussian strain + Lorentzian size)
# and plots (d_hkl × β × cosθ)² vs (d_hkl² × β × cosθ) / sinθ.
# This yields a straight line whose slope gives size and intercept gives strain.
#
# Advantage over classical W-H:
#     Less weight on high-angle reflections which often have larger errors.
# =============================================================================


@dataclass
class SSPResult:
    """Size-Strain Plot (Halder-Wagner) result.

    Attributes:
        crystallite_size_nm: Crystallite size from slope
        microstrain: Microstrain from intercept
        r_squared: Goodness of linear fit
        slope: Slope of the SSP plot
        intercept: Intercept of the SSP plot
        quality_level: Fit quality classification
        is_reliable: Whether result is trustworthy
        n_peaks: Number of peaks used
        warning_message: Diagnostic messages

    """

    crystallite_size_nm: float
    microstrain: float
    r_squared: float
    slope: float
    intercept: float
    quality_level: WHQualityLevel = WHQualityLevel.ACCEPTABLE
    is_reliable: bool = True
    n_peaks: int = 0
    warning_message: str = ""


def analyze_ssp(
    two_theta: np.ndarray,
    fwhm_sample: np.ndarray,
    hkl_list: list[tuple[int, int, int]],
    wavelength: float = CU_KA1,
    k_factor: float = WH_K_FACTOR,
    lattice_constant: float = 3.615,
) -> SSPResult:
    """Perform Size-Strain Plot (Halder-Wagner) analysis.

    Plots (d²βcosθ)² vs (d²βcosθ)/(sinθ):
        Y = (d_hkl × β_rad × cosθ)²
        X = (d_hkl² × β_rad × cosθ) / sinθ

    The linear fit yields:
        slope = Kλ / D  → D = Kλ / slope
        intercept = (ε/2)²  → ε = 2√intercept

    Args:
        two_theta: 2θ peak positions (degrees)
        fwhm_sample: Sample FWHM values (degrees, instrument-corrected)
        hkl_list: (h,k,l) Miller indices for each peak
        wavelength: X-ray wavelength in Å
        k_factor: Scherrer constant (default: 0.9)
        lattice_constant: Lattice parameter in Å (default: Cu 3.615 Å)

    Returns:
        SSPResult with crystallite size and microstrain

    Reference:
        Halder & Wagner (1966), Acta Cryst. 20, 312-313.

    """
    two_theta_arr = np.asarray(two_theta, dtype=float)
    fwhm_arr = np.asarray(fwhm_sample, dtype=float)
    n_peaks = len(two_theta_arr)

    if n_peaks < MIN_PEAKS:
        return SSPResult(
            crystallite_size_nm=0.0,
            microstrain=0.0,
            r_squared=0.0,
            slope=0.0,
            intercept=0.0,
            quality_level=WHQualityLevel.POOR,
            is_reliable=False,
            n_peaks=n_peaks,
            warning_message=f"Need ≥ {MIN_PEAKS} peaks, got {n_peaks}",
        )

    theta_rad = np.radians(two_theta_arr / 2.0)
    beta_rad = np.radians(fwhm_arr)

    # Compute d-spacings from hkl and lattice constant (cubic)
    d_hkl = np.array(
        [lattice_constant / np.sqrt(h**2 + k**2 + l**2) for h, k, l in hkl_list]
    )

    # SSP coordinates
    # Y = (d_hkl * β * cosθ)²
    # X = (d_hkl² * β * cosθ) / sinθ
    beta_cos = beta_rad * np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)

    y_ssp = (d_hkl * beta_cos) ** 2
    x_ssp = (d_hkl**2 * beta_cos) / sin_theta

    # Linear regression
    reg = linregress(x_ssp, y_ssp)
    slope = float(reg.slope)
    intercept = float(reg.intercept)
    r_squared = float(reg.rvalue**2)

    # Extract physical quantities
    # slope = Kλ/D → D = Kλ/slope
    if slope > 0:
        D_angstrom = k_factor * wavelength / slope
        D_nm = D_angstrom / 10.0
    else:
        D_nm = 0.0

    # intercept = (ε/2)² → ε = 2√|intercept|
    if intercept >= 0:
        microstrain = 2.0 * np.sqrt(intercept)
    else:
        # Negative intercept can occur with noisy data; report zero strain
        microstrain = 0.0

    # Quality
    if r_squared >= R2_EXCELLENT:
        quality = WHQualityLevel.EXCELLENT
    elif r_squared >= R2_ACCEPTABLE:
        quality = WHQualityLevel.ACCEPTABLE
    else:
        quality = WHQualityLevel.POOR

    warning = ""
    if r_squared < R2_ACCEPTABLE:
        warning = f"R² = {r_squared:.3f} < {R2_ACCEPTABLE} — SSP fit unreliable"

    return SSPResult(
        crystallite_size_nm=D_nm,
        microstrain=microstrain,
        r_squared=r_squared,
        slope=slope,
        intercept=intercept,
        quality_level=quality,
        is_reliable=(r_squared >= R2_ACCEPTABLE and D_nm > 0),
        n_peaks=n_peaks,
        warning_message=warning,
    )


# =============================================================================
# Backward Compatibility Aliases
# =============================================================================

WHResultEnhanced = WHResult
WilliamsonHallEnhanced = WilliamsonHallAnalyzer
analyze_williamson_hall_enhanced = analyze_williamson_hall
williamson_hall_analysis = analyze_williamson_hall

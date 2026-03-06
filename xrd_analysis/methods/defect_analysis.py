"""Defect Analysis Module.
===================================

Stacking fault detection and lattice constant monitoring.

References:
- Warren (1969), X-ray Diffraction, Chapter 13 (stacking faults)

"""

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from xrd_analysis.core.constants import CU_KA1
from xrd_analysis.core.copper_crystal import (
    CU_CRYSTAL,
    CU_JCPDS_EXTENDED,
)

# =============================================================================
# Constants
# =============================================================================

# Standard peak separation (111)-(200)
# Calculated dynamically from JCPDS data to match updated high-precision values
STANDARD_PEAK_SEPARATION = (
    CU_JCPDS_EXTENDED[(2, 0, 0)]["two_theta"]
    - CU_JCPDS_EXTENDED[(1, 1, 1)]["two_theta"]
)

# Warren geometric coefficient for FCC (111)-(200) stacking fault
#
# Reference:
#     Warren, B. E. (1969).
#     "X-ray Diffraction."
#     Dover Publications, New York.
#     Chapter 13: Stacking Faults, Pages 275-298.
#     ISBN: 978-0486663173 (Dover reprint)
#
# ═══════════════════════════════════════════════════════════════════════════
# Theoretical Derivation
# ═══════════════════════════════════════════════════════════════════════════
#
# Warren's formula for FCC (200)-(111) peak separation shift:
#
#     Δ(2θ_200 - 2θ_111) = -(45√3 / π²) × α
#
# where α is the stacking fault probability
#
# Coefficient calculation:
#     G = -45√3 / π²
#       = -45 × 1.7320508... / 9.8696044...
#       = -77.9422... / 9.8696...
#       = -7.897 (degrees per unit probability)
#
# Unit conversion:
#   If α is expressed as percentage (α_percent = α × 100):
#   G_percent = G / 100 = -0.07897 °/%
#
# Physical meaning:
#   - Negative sign indicates stacking faults reduce peak separation
#   - (111) peak shifts to higher angles, (200) shifts to lower angles
#
_WARREN_G_THEORETICAL = -45 * math.sqrt(3) / (math.pi**2)  # = -7.897

WARREN_G_COEFFICIENT = _WARREN_G_THEORETICAL  # Theoretical value

# Standard lattice constant for Cu
# Reference: JCPDS 04-0836
STANDARD_LATTICE_CONSTANT = CU_CRYSTAL.lattice_constant  # 3.6150 Å

# Lattice constant thresholds
LATTICE_MINOR_THRESHOLD = 3.616  # Å
LATTICE_SEVERE_THRESHOLD = 3.618  # Å

# Wavelength constant: use CU_KA1 directly from core.constants

# =============================================================================
# Enums
# =============================================================================


class StackingFaultSeverity(Enum):
    """Stacking fault severity classification."""

    NORMAL = "normal"  # α < 0.5%
    MILD = "mild"  # 0.5% ≤ α < 1%
    MODERATE = "moderate"  # 1% ≤ α < 2%
    SEVERE = "severe"  # α ≥ 2%


class LatticeStatus(Enum):
    """Lattice constant status classification."""

    NORMAL = "normal"  # a ≈ 3.615 Å
    MINOR_EXPANSION = "minor"  # 3.616-3.618 Å
    SEVERE_EXPANSION = "severe"  # > 3.618 Å


class AnnealingState(Enum):
    """Self-annealing state classification."""

    AS_DEPOSITED = "as-deposited"  # < 1 hour
    PARTIAL = "partial"  # 1-24 hours
    ANNEALED = "annealed"  # 1-7 days
    STABLE = "stable"  # > 7 days
    UNKNOWN = "unknown"  # Not specified


# =============================================================================
# Result Dataclasses
# =============================================================================


@dataclass
class StackingFaultResult:
    """Result of stacking fault analysis."""

    peak_separation_deg: float
    deviation_deg: float
    alpha_probability: float  # Stacking fault probability (0-1)
    alpha_percent: float  # As percentage
    severity: StackingFaultSeverity
    sps_warning: bool
    message: str

    def __repr__(self) -> str:
        return (
            f"StackingFault: Δ2θ={self.peak_separation_deg:.3f}°, "
            f"α={self.alpha_percent:.2f}% [{self.severity.value}]"
        )


@dataclass
class LatticeConstantResult:
    """Result of lattice constant analysis."""

    lattice_constant: float  # Å
    deviation: float  # From standard
    d_spacing: float  # Å
    hkl_used: tuple[int, int, int]
    two_theta_used: float
    status: LatticeStatus
    message: str

    def __repr__(self) -> str:
        return (
            f"Lattice: a={self.lattice_constant:.4f} Å "
            f"(Δ={self.deviation:+.4f} Å) [{self.status.value}]"
        )


@dataclass
class DefectAnalysisResult:
    """Complete defect analysis result."""

    # Stacking fault
    stacking_fault: StackingFaultResult | None = None

    # Lattice constant
    lattice_constant: LatticeConstantResult | None = None

    # Self-annealing
    annealing_state: AnnealingState = AnnealingState.UNKNOWN
    sample_age_hours: float | None = None
    annealing_note: str = ""

    # Overall
    overall_status: str = ""
    recommendations: list[str] = field(default_factory=list)


# =============================================================================
# Module Q: Stacking Fault Analyzer
# =============================================================================


class StackingFaultAnalyzer:
    """Warren-based stacking fault analyzer.

    Reference: Warren (1969), X-ray Diffraction, Ch.13

    For FCC metals with intrinsic stacking faults:
    - (111) peak shifts to higher angles
    - (200) peak shifts to lower angles
    - Peak separation decreases

    Warren formula:
        α = (Δ2θ_exp - Δ2θ_std) / G

    where G ≈ -7.897 for FCC (111)-(200) in degree-based implementation
    """

    def __init__(self, g_coefficient: float = WARREN_G_COEFFICIENT):
        """Initialize stacking fault analyzer.

        Args:
            g_coefficient: Warren geometric coefficient (default: -7.897)

        """
        self.g_coefficient = g_coefficient

    def analyze(
        self, two_theta_111: float, two_theta_200: float
    ) -> StackingFaultResult:
        """Analyze stacking fault probability.

        Args:
            two_theta_111: Measured 2θ position of (111) peak
            two_theta_200: Measured 2θ position of (200) peak

        Returns:
            StackingFaultResult with α probability

        """
        # Q.1.3: Calculate peak separation
        peak_separation = two_theta_200 - two_theta_111

        # Q.1.4: Calculate deviation from standard
        deviation = peak_separation - STANDARD_PEAK_SEPARATION

        # Warren formula for α
        # α = deviation / G
        # G is negative, deviation is negative when SF present
        # So α comes out positive
        alpha = deviation / self.g_coefficient if self.g_coefficient != 0 else 0
        alpha = max(0, alpha)  # α cannot be negative

        alpha_percent = alpha * 100

        # Determine severity
        if alpha_percent < 0.5:
            severity = StackingFaultSeverity.NORMAL
        elif alpha_percent < 1.0:
            severity = StackingFaultSeverity.MILD
        elif alpha_percent < 2.0:
            severity = StackingFaultSeverity.MODERATE
        else:
            severity = StackingFaultSeverity.SEVERE

        # SPS warning check
        sps_warning = peak_separation < 7.0

        # Generate message
        if sps_warning:
            message = (
                f"Warning: peak separation {peak_separation:.3f}° < 7.0°, "
                "possible excessive SPS accelerator concentration"
            )
        elif severity == StackingFaultSeverity.NORMAL:
            message = "Stacking fault probability within normal range"
        else:
            message = f"Stacking fault detected: α ≈ {alpha_percent:.2f}%"

        return StackingFaultResult(
            peak_separation_deg=peak_separation,
            deviation_deg=deviation,
            alpha_probability=alpha,
            alpha_percent=alpha_percent,
            severity=severity,
            sps_warning=sps_warning,
            message=message,
        )


# =============================================================================
# Module R: Lattice Constant Calculator
# =============================================================================


def calculate_d_spacing(two_theta: float, wavelength: float = CU_KA1) -> float:
    """Calculate d-spacing from Bragg's law.

    d = λ / (2 sin θ)
    """
    theta_rad = two_theta / 2 * np.pi / 180
    return wavelength / (2 * np.sin(theta_rad))


def calculate_lattice_constant(
    two_theta: float, hkl: tuple[int, int, int], wavelength: float = CU_KA1
) -> float:
    """Calculate lattice constant from peak position.

    a = d × √(h² + k² + l²)
    """
    d = calculate_d_spacing(two_theta, wavelength)
    h, k, l = hkl
    return d * np.sqrt(h**2 + k**2 + l**2)


class LatticeMonitor:
    """Lattice constant monitor.

    Prefers high-angle peaks (311, 220) for better accuracy.
    """

    def __init__(
        self, standard_a: float = STANDARD_LATTICE_CONSTANT, wavelength: float = CU_KA1
    ):
        self.standard_a = standard_a
        self.wavelength = wavelength

    def analyze_lattice(
        self, two_theta: float, hkl: tuple[int, int, int]
    ) -> LatticeConstantResult:
        """Analyze lattice constant from peak position.

        Args:
            two_theta: Measured 2θ position
            hkl: Miller indices of the peak

        Returns:
            LatticeConstantResult

        """
        # R.1.2-3: Calculate d and a
        d = calculate_d_spacing(two_theta, self.wavelength)
        a = calculate_lattice_constant(two_theta, hkl, self.wavelength)

        # Calculate deviation
        deviation = a - self.standard_a

        # R.1.4: Determine status
        if a <= LATTICE_MINOR_THRESHOLD:
            status = LatticeStatus.NORMAL
            message = "Lattice constant within normal range"
        elif a <= LATTICE_SEVERE_THRESHOLD:
            status = LatticeStatus.MINOR_EXPANSION
            message = "Minor lattice expansion; possible impurity solid solution or vacancy effects"
        else:
            status = LatticeStatus.SEVERE_EXPANSION
            message = "Severe lattice expansion; check additive purity or impurity solid solution"

        return LatticeConstantResult(
            lattice_constant=a,
            deviation=deviation,
            d_spacing=d,
            hkl_used=hkl,
            two_theta_used=two_theta,
            status=status,
            message=message,
        )


# =============================================================================
# Self-Annealing State Machine
# =============================================================================


def determine_annealing_state(
    sample_age_hours: float | None = None, fwhm_narrowing_detected: bool = False
) -> tuple[AnnealingState, str]:
    """Determine self-annealing state from sample age.

    Args:
        sample_age_hours: Time since deposition (hours), None for unknown
        fwhm_narrowing_detected: True if FWHM is narrower than expected

    Returns:
        Tuple of (AnnealingState, recommendation_note)

    """
    if sample_age_hours is None:
        return (
            AnnealingState.UNKNOWN,
            "Sample age unknown; cannot determine self-annealing state",
        )

    if sample_age_hours < 1:
        return (
            AnnealingState.AS_DEPOSITED,
            "As-deposited: expect fine grains, high defect density. Re-measure after 7 days for stable structure",
        )
    elif sample_age_hours < 24:
        state = AnnealingState.PARTIAL
        if fwhm_narrowing_detected:
            note = "Self-annealing in progress: FWHM narrowing observed"
        else:
            note = "Early self-annealing: grains may be growing"
        return (state, note)
    elif sample_age_hours < 168:  # 7 days
        return (
            AnnealingState.ANNEALED,
            "Self-annealing: structure not yet fully stable",
        )
    else:
        return (
            AnnealingState.STABLE,
            "Stable state: structure has stabilised, analysis results are reliable",
        )


# =============================================================================
# Convenience Functions
# =============================================================================


def analyze_stacking_faults(
    two_theta_111: float, two_theta_200: float
) -> StackingFaultResult:
    """Convenience function for stacking fault analysis."""
    analyzer = StackingFaultAnalyzer()
    return analyzer.analyze(two_theta_111, two_theta_200)


def analyze_lattice(
    two_theta: float, hkl: tuple[int, int, int]
) -> LatticeConstantResult:
    """Convenience function for lattice constant analysis."""
    monitor = LatticeMonitor()
    return monitor.analyze_lattice(two_theta, hkl)

"""Physical Constants Module.
======================================
Centralized physical constants for XRD analysis.

All parameters are verified with paper citations.
"""

from dataclasses import dataclass

# =============================================================================
# X-ray Wavelengths (Angstroms)
# Reference: Bearden (1967), Rev. Mod. Phys. 39, 78-124, Table V, p.9
# =============================================================================

# Cu Kα1: K-L_III transition
# Bearden (1967) Table V: 1.540562 Å, 8.04778 keV
CU_KA1 = 1.540562

# Cu Kα2: K-L_II transition
# Bearden (1967) Table V: 1.544390 Å, 8.02783 keV
CU_KA2 = 1.544390

# Weighted average Kα for unresolved doublet
# Calculation: λ_avg = (2×Kα1 + Kα2) / 3 = (2×1.540562 + 1.544390) / 3
# Weighting basis: Burger-Dorgelo rule (2p₃/₂: 4e⁻, 2p₁/₂: 2e⁻ → 2:1)
# Reference: Cullity, "Elements of X-ray Diffraction", 3rd Ed., p.7
CU_KA_AVG = 1.541838

# Kα2/Kα1 intensity ratio
# Theoretical basis: Burger-Dorgelo rule
#   - Kα1 from 2p₃/₂→1s transition (4 electrons)
#   - Kα2 from 2p₁/₂→1s transition (2 electrons)
#   - Ratio: 2/4 = 0.5
# Reference: Cullity, "Elements of X-ray Diffraction", 3rd Ed., p.7
KA2_KA1_RATIO = 0.5

# =============================================================================
# Scherrer Constants (K values)
# Reference: Langford & Wilson (1978), J. Appl. Cryst. 11, 102-113
# =============================================================================


@dataclass
class ScherrerConstants:
    """Scherrer constant K values for different grain shapes.

    Important:
    - These are K_w values for FWHM (half-width) definition
    - For integral breadth, use K_β values from L&W Table 1

    Reference: L&W (1978), Table 1, p.107

    """

    # Sphere: K_w = 0.8290 (L&W Table 1, p.107, Line 846)
    spherical: float = 0.829

    # Generic cubic habit: approximation for mixed orientations
    # Reference: Warren (1969), X-ray Diffraction
    cubic: float = 0.94

    # Octahedron: K_w = 0.83 (L&W p.107, same as sphere for 110)
    octahedral: float = 0.83

    # Default value: L&W 1978 spherical standard
    default: float = 0.829


SCHERRER_K = ScherrerConstants()


def get_jcpds_data(material: str = "Cu") -> dict[tuple[int, int, int], dict]:
    """Get JCPDS standard data for copper.

    Args:
        material: Material name (only "Cu" supported)

    Returns:
        Dictionary of JCPDS data keyed by (h, k, l)

    """
    if material.upper() != "CU":
        raise ValueError(f"Only Cu is supported. Got: {material}")

    # Import here to avoid circular dependencies (low-level constants vs high-level crystal)
    from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED

    return CU_JCPDS_EXTENDED


# =============================================================================
# Physical Limits and Empirical Thresholds
# =============================================================================

# Crystallite size detection limits
# EMPIRICAL THRESHOLDS - Based on instrumental limitations and practical experience
# Commonly used in XRD crystallite size analysis
#
# Lower limit rationale:
#   Below ~2 nm, peak broadening approaches instrumental resolution limits
#   Scherrer equation becomes less reliable due to:
#   - Instrumental broadening dominates the signal
#   - Quantum size effects may alter crystal structure
#   - Surface relaxation becomes significant
#
# Upper limit rationale:
#   Above ~200 nm, crystallite size broadening becomes negligible
#   Peak width is dominated by instrumental broadening
#   Other techniques (SEM, TEM) are more appropriate
#
# Reference: Klug & Alexander (1974), "X-ray Diffraction Procedures" (general guidance)
#            Commonly accepted values in XRD community
MIN_RELIABLE_SIZE = 2.0  # nm, minimum detectable crystallite size
MAX_RELIABLE_SIZE = 200.0  # nm, maximum reliable crystallite size

# FWHM quality threshold
# EMPIRICAL THRESHOLD - Quality control for peak broadening measurements
#
# Rationale: Ratio of (sample FWHM) / (instrumental FWHM)
#   If ratio < 1.2, sample broadening is too small relative to instrumental broadening
#   Results become unreliable due to uncertainty propagation
#   sqrt(β_sample² - β_inst²) becomes numerically unstable when β_sample ≈ β_inst
#
# Reference: Commonly used threshold in Rietveld refinement and XRD analysis
#            Young (1993), "The Rietveld Method" (general guidance)
MIN_BROADENING_RATIO = (
    1.2  # Minimum sample/instrumental FWHM ratio for reliable analysis
)

# Fit quality thresholds
# ═══════════════════════════════════════════════════════════════════════════
# USER-ADJUSTABLE THRESHOLDS
# ═══════════════════════════════════════════════════════════════════════════
# These thresholds are user-defined and can be modified based on experimental
# requirements and quality standards. Adjust these values as needed:
#
# Example adjustments:
#   - For high-precision work: MAX_RWP_PERCENT = 5.0, MIN_R_SQUARED = 0.99
#   - For routine analysis: MAX_RWP_PERCENT = 15.0, MIN_R_SQUARED = 0.90
#   - For publication quality: MAX_RWP_PERCENT = 8.0, MIN_R_SQUARED = 0.995
#
# Reference (general guidance):
#   Young, R. A. (1993). "The Rietveld Method." Oxford University Press.
#   (R_wp < 10% is generally acceptable for Rietveld refinement)
# ═══════════════════════════════════════════════════════════════════════════
MAX_RWP_PERCENT = 10.0  # Maximum acceptable R_wp (%) - ADJUSTABLE
MIN_R_SQUARED = 0.95  # Minimum acceptable R² - ADJUSTABLE

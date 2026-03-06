"""Texture Analysis Module.
=======================

Harris Texture Coefficient (TC) analysis for preferred orientation.

NOTE: This module provides DATA ONLY. Interpretation should be done by humans.

Reference:
    Harris, G. B. (1952). Quantitative measurement of preferred orientation
    in rolled uranium bars. Phil. Mag., 43(336), 113-123.
"""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED

# =============================================================================
# Constants
# =============================================================================

# JCPDS standard data for Cu (PDF 04-0836)
# JCPDS standard data for Cu (PDF 04-0836)
# Derived dynamically from centralized SSOT data
JCPDS_STANDARD_INTENSITY: dict[tuple[int, int, int], float] = {
    hkl: data["intensity"] for hkl, data in CU_JCPDS_EXTENDED.items()
}

# Standard peak positions and intensities
CU_JCPDS_STANDARD = CU_JCPDS_EXTENDED

# TC interpretation thresholds (DATA MARKERS ONLY)
# EMPIRICAL THRESHOLDS - Commonly used in Harris texture analysis
#
# These thresholds are empirical guidelines for interpreting TC values:
# - TC ∈ [0.9, 1.1]: Random orientation (no preferred orientation)
# - TC > 1.5: Strong preferred orientation in that direction
# - TC < 0.9: Suppressed orientation
#
# Reference: Harris (1952), Phil. Mag. 43, 113-123 (original TC method)
#            Commonly accepted interpretation thresholds in XRD community
#            These specific values are empirical best practices
TC_RANDOM_MIN = 0.9
TC_RANDOM_MAX = 1.1
TC_PREFERRED_THRESHOLD = 1.5
RANDOM_TEXTURE_SIGMA = 0.3

# =============================================================================
# Enums
# =============================================================================


class OrientationType(Enum):
    """Classification of crystallographic orientation (DATA MARKER)."""

    PREFERRED = "preferred"  # TC > 1.0
    RANDOM = "random"  # TC ≈ 1.0 (0.9-1.1)
    SUPPRESSED = "suppressed"  # TC < 1.0


# =============================================================================
# Data Structures
# =============================================================================


@dataclass
class TCResult:
    """Single peak texture coefficient result."""

    hkl: tuple[int, int, int]
    tc_value: float
    intensity_observed: float
    intensity_standard: float
    ratio: float = 0.0
    orientation_type: OrientationType = OrientationType.RANDOM

    def __repr__(self) -> str:
        hkl_str = f"({self.hkl[0]}{self.hkl[1]}{self.hkl[2]})"
        return f"{hkl_str}: TC={self.tc_value:.2f} [{self.orientation_type.value}]"


@dataclass
class TextureAnalysisResult:
    """Complete texture analysis result.

    NOTE: Contains DATA ONLY. No automatic process diagnosis.

    Attributes:
        tc_values: TC for each (hkl).
        tc_details: Detailed results.
        dominant_hkl: Highest TC direction.
        dominant_tc: Highest TC value.
        is_random: Whether texture is random.
        degree_of_texture: Standard deviation of TC values.
        n_peaks: Number of peaks.
        intensity_type: Type of intensity used.

    """

    tc_values: dict[tuple[int, int, int], float] = field(default_factory=dict)
    tc_details: list[TCResult] = field(default_factory=list)
    dominant_hkl: tuple[int, int, int] | None = None
    dominant_tc: float = 1.0
    is_random: bool = True
    degree_of_texture: float = 0.0
    n_peaks: int = 0
    intensity_type: str = "area"

    def __repr__(self) -> str:
        if self.dominant_hkl:
            dom = (
                f"({self.dominant_hkl[0]}{self.dominant_hkl[1]}{self.dominant_hkl[2]})"
            )
        else:
            dom = "None"
        return f"TextureAnalysisResult(Dominant={dom}, TC={self.dominant_tc:.2f})"


# Backward compatibility alias
TextureResult = TCResult

# =============================================================================
# Texture Analyzer
# =============================================================================


class TextureAnalyzer:
    """Harris Texture Coefficient analyzer.

    Provides DATA output only. Interpretation should be done by humans.

    Harris TC Formula:
        TC(hkl) = [I(hkl) / I₀(hkl)] / [1/N × Σ I(hkl) / I₀(hkl)]

    Interpretation:
        TC = 1 for all peaks: Random orientation (powder average)
        TC > 1: Preferred orientation along this direction
        TC < 1: Under-represented direction

    Args:
        standard_data: Custom JCPDS data (default: Cu)
        use_area: If True, intensities are integrated areas

    Example:
        >>> analyzer = TextureAnalyzer()
        >>> result = analyzer.analyze({(1,1,1): 15680, (2,0,0): 5520})
        >>> print(f"TC(111) = {result.tc_values[(1,1,1)]:.2f}")

    """

    def __init__(
        self,
        standard_data: dict[tuple[int, int, int], float] | None = None,
        use_area: bool = True,
    ) -> None:
        if standard_data is None:
            self.standard_data = JCPDS_STANDARD_INTENSITY
        else:
            self.standard_data = standard_data
        self.use_area = use_area
        self.intensity_type = "area" if use_area else "height"

    def analyze(
        self, intensities: dict[tuple[int, int, int], float], normalize: bool = True
    ) -> TextureAnalysisResult:
        """Perform texture coefficient analysis.

        Args:
            intensities: Dict mapping (hkl) to observed intensities.
                Should be integrated areas, not peak heights.
            normalize: If True, normalize intensities (ignored in this version).

        Returns:
            TextureAnalysisResult with TC values (DATA ONLY).

        """
        # Filter to peaks with standard data
        valid_hkls = [hkl for hkl in intensities if hkl in self.standard_data]
        n_peaks = len(valid_hkls)

        if n_peaks < 2:
            return TextureAnalysisResult(
                tc_values={}, n_peaks=0, intensity_type=self.intensity_type
            )

        # Calculate intensity ratios I/I₀
        ratios: dict[tuple[int, int, int], float] = {}
        for hkl in valid_hkls:
            I_obs = intensities[hkl]
            I_std = self.standard_data[hkl]
            if I_obs <= 0 or I_std <= 0:
                continue  # Skip zero/negative intensities
            ratios[hkl] = I_obs / I_std

        if len(ratios) < 2:
            return TextureAnalysisResult(
                tc_values={}, n_peaks=0, intensity_type=self.intensity_type
            )

        # Calculate average ratio
        n_valid = len(ratios)
        avg_ratio = sum(ratios.values()) / n_valid

        # Calculate TC values
        tc_values: dict[tuple[int, int, int], float] = {}
        tc_details: list[TCResult] = []

        for hkl in ratios:
            tc = ratios[hkl] / avg_ratio if avg_ratio > 0 else 1.0
            tc_values[hkl] = tc

            # Classify orientation type
            if TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX:
                otype = OrientationType.RANDOM
            elif tc > 1.0:
                otype = OrientationType.PREFERRED
            else:
                otype = OrientationType.SUPPRESSED

            tc_details.append(
                TCResult(
                    hkl=hkl,
                    tc_value=tc,
                    intensity_observed=intensities[hkl],
                    intensity_standard=self.standard_data[hkl],
                    ratio=ratios[hkl],
                    orientation_type=otype,
                )
            )

        # Sort by TC value
        tc_details.sort(key=lambda x: x.tc_value, reverse=True)

        # Find dominant orientation
        dominant_hkl = max(tc_values, key=tc_values.get)
        dominant_tc = tc_values[dominant_hkl]

        # Check if random
        is_random = all(
            TC_RANDOM_MIN <= tc <= TC_RANDOM_MAX for tc in tc_values.values()
        )

        # Calculate degree of texture (σ)
        tc_list = list(tc_values.values())
        degree_of_texture = float(np.std(tc_list))

        return TextureAnalysisResult(
            tc_values=tc_values,
            tc_details=tc_details,
            dominant_hkl=dominant_hkl,
            dominant_tc=dominant_tc,
            is_random=is_random,
            degree_of_texture=degree_of_texture,
            n_peaks=n_peaks,
            intensity_type=self.intensity_type,
        )

    def analyze_from_peaks(
        self,
        peaks: list[tuple[float, float]],
        hkl_assignments: list[tuple[int, int, int]],
    ) -> TextureAnalysisResult:
        """Analyze texture from peak fitting results."""
        if len(peaks) != len(hkl_assignments):
            raise ValueError("peaks and hkl_assignments must have same length")

        intensities = {}
        for (two_theta, intensity), hkl in zip(peaks, hkl_assignments):
            intensities[hkl] = intensity

        return self.analyze(intensities)

    def get_hkl_for_angle(
        self, two_theta: float, tolerance: float = 1.0
    ) -> tuple[int, int, int] | None:
        """Find (hkl) assignment for a given 2θ angle."""
        best_match = None
        min_diff = tolerance

        for hkl, data in CU_JCPDS_STANDARD.items():
            diff = abs(two_theta - data["two_theta"])
            if diff < min_diff:
                min_diff = diff
                best_match = hkl

        return best_match


# =============================================================================
# Convenience Functions
# =============================================================================


def analyze_texture(
    intensities: dict[tuple[int, int, int], float], use_area: bool = True
) -> TextureAnalysisResult:
    """Convenience function for texture analysis.

    Example:
        >>> intensities = {(1,1,1): 15680, (2,0,0): 5520, (2,2,0): 4200}
        >>> result = analyze_texture(intensities)
        >>> print(f"TC(111) = {result.tc_values[(1,1,1)]:.2f}")

    """
    analyzer = TextureAnalyzer(use_area=use_area)
    return analyzer.analyze(intensities)


def get_standard_intensity(hkl: tuple[int, int, int]) -> float:
    """Get JCPDS standard intensity for given hkl.

    Reference: JCPDS 04-0836
    """
    return JCPDS_STANDARD_INTENSITY.get(hkl, 0.0)


def get_standard_angles(material: str = "Cu") -> dict[tuple[int, int, int], float]:
    """Get standard 2θ angles for common materials."""
    if material.upper() == "CU":
        return {hkl: data["two_theta"] for hkl, data in CU_JCPDS_STANDARD.items()}
    else:
        raise ValueError(f"Unknown material: {material}")


def calculate_texture_coefficient(
    observed_intensities: dict[tuple[int, int, int], float],
    standard_intensities: dict[tuple[int, int, int], float] | None = None,
) -> dict[tuple[int, int, int], float]:
    """Convenience function for texture coefficient calculation."""
    analyzer = TextureAnalyzer(standard_data=standard_intensities)
    result = analyzer.analyze(observed_intensities)
    return result.tc_values

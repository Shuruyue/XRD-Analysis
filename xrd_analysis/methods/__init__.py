"""
xrd_analysis Methods Module
=====================

Analysis methods for XRD data.
XRD 資料分析方法。
"""

from xrd_analysis.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    GrainShape,
    calculate_scherrer,
    calculate_crystallite_size,
)

from xrd_analysis.methods.williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
    WHQualityLevel,
    analyze_williamson_hall,
    get_modulus_for_hkl,
)

from xrd_analysis.methods.texture import (
    TextureAnalyzer,
    TextureAnalysisResult,
    TCResult,
    OrientationType,
    analyze_texture,
    get_standard_intensity,
    calculate_texture_coefficient,
)

from xrd_analysis.methods.caglioti import (
    CagliotiCorrection,
    CagliotiParams,
    calculate_instrumental_broadening,
)

__all__ = [
    # Scherrer
    "ScherrerCalculator",
    "ScherrerResult",
    "ValidityFlag",
    "GrainShape",
    "calculate_scherrer",
    "calculate_crystallite_size",
    # W-H
    "WilliamsonHallAnalyzer",
    "WHResult",
    "WHQualityLevel",
    "analyze_williamson_hall",
    "get_modulus_for_hkl",
    # Texture
    "TextureAnalyzer",
    "TextureAnalysisResult",
    "TCResult",
    "OrientationType",
    "analyze_texture",
    "get_standard_intensity",
    "calculate_texture_coefficient",
    # Caglioti
    "CagliotiCorrection",
    "CagliotiParams",
    "calculate_instrumental_broadening",
]


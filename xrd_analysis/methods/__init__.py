"""xrd_analysis Methods Module.
=====================

Analysis methods for XRD data.
"""

from xrd_analysis.methods.caglioti import (
    CagliotiCorrection,
    CagliotiParams,
    calculate_instrumental_broadening,
)
from xrd_analysis.methods.scherrer import (
    GrainShape,
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    calculate_crystallite_size,
    calculate_scherrer,
)
from xrd_analysis.methods.texture import (
    OrientationType,
    TCResult,
    TextureAnalysisResult,
    TextureAnalyzer,
    analyze_texture,
    calculate_texture_coefficient,
    get_standard_intensity,
)
from xrd_analysis.methods.williamson_hall import (
    WHQualityLevel,
    WHResult,
    WilliamsonHallAnalyzer,
    analyze_williamson_hall,
    get_modulus_for_hkl,
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

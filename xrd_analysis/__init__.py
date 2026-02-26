"""XRD-Analysis: Advanced XRD Crystallite Size Analysis System
=====================================================

A comprehensive toolkit for XRD data analysis including:
- Scherrer crystallite size calculation
- Harris Texture Coefficient analysis
- Pseudo-Voigt peak fitting

Example:
    >>> import xrd_analysis
    >>> print(xrd_analysis.__version__)
    0.1.0

"""

import logging

# Library best practice: add NullHandler to prevent "No handlers" warnings.
# Consuming applications configure their own handlers.
logging.getLogger(__name__).addHandler(logging.NullHandler())

from xrd_analysis.__version__ import __version__, __version_info__
from xrd_analysis.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
    calculate_crystallite_size,
    calculate_scherrer,
)
from xrd_analysis.methods.texture import (
    TextureAnalysisResult,
    TextureAnalyzer,
    analyze_texture,
)

__all__ = [
    # Version
    "__version__",
    "__version_info__",
    # Scherrer
    "ScherrerCalculator",
    "ScherrerResult",
    "ValidityFlag",
    "calculate_scherrer",
    "calculate_crystallite_size",
    # Texture
    "TextureAnalyzer",
    "TextureAnalysisResult",
    "analyze_texture",
]

"""Validation Module
Contains validation functions for XRD analysis results.
"""

from .error_analysis import ErrorAnalyzer, check_broadening_ratio, validate_size_range
from .goodness_of_fit import (
    FitQualityResult,
    calculate_chi_squared,
    calculate_r_squared,
    calculate_rwp,
)

__all__ = [
    "calculate_rwp",
    "calculate_r_squared",
    "calculate_chi_squared",
    "FitQualityResult",
    "ErrorAnalyzer",
    "validate_size_range",
    "check_broadening_ratio",
]

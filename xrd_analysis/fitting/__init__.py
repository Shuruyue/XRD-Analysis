# Fitting Module
"""Peak fitting module for XRD analysis.
Includes peak detection, Pseudo-Voigt fitting, quality metrics, and hkl assignment.
"""

from .hkl_assignment import (
    JCPDS_COPPER_PEAKS,
    PeakAssignment,
    assign_all_peaks,
    assign_hkl,
    assign_hkl_detailed,
    format_hkl,
)
from .lm_optimizer import FitResult, LMOptimizer, fit_peaks
from .peak_detection import PeakDetector, PeakInfo, find_peaks
from .pseudo_voigt import PseudoVoigt, PseudoVoigtParams, pseudo_voigt_function
from .quality_metrics import (
    FitQualityReport,
    QualityLevel,
    calculate_r_squared,
    calculate_r_wp,
    calculate_rss,
    generate_quality_report,
    validate_fit_parameters,
)

__all__ = [
    # Peak detection
    "PeakDetector",
    "PeakInfo",
    "find_peaks",
    # Pseudo-Voigt
    "PseudoVoigt",
    "PseudoVoigtParams",
    "pseudo_voigt_function",
    # LM Optimizer
    "LMOptimizer",
    "FitResult",
    "fit_peaks",
    # Quality metrics
    "QualityLevel",
    "FitQualityReport",
    "calculate_r_wp",
    "calculate_rss",
    "calculate_r_squared",
    "validate_fit_parameters",
    "generate_quality_report",
    # HKL assignment
    "PeakAssignment",
    "assign_hkl",
    "assign_hkl_detailed",
    "assign_all_peaks",
    "format_hkl",
    "JCPDS_COPPER_PEAKS",
]

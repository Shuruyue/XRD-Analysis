"""Custom Exception Hierarchy.
==========================================

Centralised exception classes for XRD-Analysis.

Usage::

    from xrd_analysis.core.exceptions import FittingError

    raise FittingError("Doublet fitting failed to converge")
"""


class XRDAnalysisError(Exception):
    """Base exception for all XRD-Analysis errors."""


class DataLoadError(XRDAnalysisError):
    """Failed to load or parse XRD data file."""


class CalibrationError(XRDAnalysisError):
    """Instrument calibration error (e.g. invalid Caglioti parameters)."""


class FittingError(XRDAnalysisError):
    """Peak fitting convergence or numerical error."""


class FittingConvergenceError(FittingError):
    """Fitting failed to converge within maximum iterations."""


class ValidationError(XRDAnalysisError):
    """Data or result validation failed."""


class ConfigurationError(XRDAnalysisError):
    """Invalid configuration parameter."""


class PreprocessingError(XRDAnalysisError):
    """Error during data preprocessing."""

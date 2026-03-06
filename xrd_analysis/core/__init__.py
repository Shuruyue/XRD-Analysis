"""xrd_analysis Core Module.
==================

Core domain logic including constants, configuration, and crystal data.
"""

from xrd_analysis.core.config_loader import load_config
from xrd_analysis.core.constants import (
    CU_KA1,
    CU_KA2,
    CU_KA_AVG,
    KA2_KA1_RATIO,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
    MIN_RELIABLE_SIZE,
    SCHERRER_K,
)
from xrd_analysis.core.copper_crystal import (
    CU_ELASTIC,
    SCHERRER_CUBIC_K,
    get_k_for_hkl,
    get_standard_peaks,
)
from xrd_analysis.core.exceptions import (
    CalibrationError,
    ConfigurationError,
    DataLoadError,
    FittingConvergenceError,
    FittingError,
    PreprocessingError,
    ValidationError,
    XRDAnalysisError,
)
from xrd_analysis.core.units import (
    angstrom_to_nm,
    deg_to_rad,
    nm_to_angstrom,
    rad_to_deg,
)

__all__ = [
    # Constants - X-ray wavelengths
    "CU_KA1",
    "CU_KA2",
    "CU_KA_AVG",
    "KA2_KA1_RATIO",
    # Constants - Scherrer
    "SCHERRER_K",
    # Constants - Limits
    "MIN_RELIABLE_SIZE",
    "MAX_RELIABLE_SIZE",
    "MIN_BROADENING_RATIO",
    # Constants - JCPDS
    # CU_JCPDS removed (Use copper_crystal.CU_JCPDS_EXTENDED)
    # Copper crystal
    "get_k_for_hkl",
    "get_standard_peaks",
    "SCHERRER_CUBIC_K",
    "CU_ELASTIC",
    # Config
    "load_config",
    # Units
    "deg_to_rad",
    "rad_to_deg",
    "angstrom_to_nm",
    "nm_to_angstrom",
    # Exceptions
    "XRDAnalysisError",
    "DataLoadError",
    "CalibrationError",
    "FittingError",
    "FittingConvergenceError",
    "ValidationError",
    "ConfigurationError",
    "PreprocessingError",
]

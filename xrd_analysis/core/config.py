"""Unified Parameter Configuration System.
======================================

Centralized configuration management for all xrd_analysis parameters.

This module provides a hierarchical configuration system:
1. Physical constants (constants.py) - Read-only
2. Default configuration (ParameterConfig) - Programmatic
3. User configuration (config.yaml) - File-based override
"""

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# Import physical constants
from xrd_analysis.core.constants import (
    CU_KA1,
    MAX_RELIABLE_SIZE,
    MIN_BROADENING_RATIO,
    MIN_RELIABLE_SIZE,
)

logger = logging.getLogger(__name__)

# =============================================================================
# Configuration Dataclasses
# =============================================================================


@dataclass
class InstrumentConfig:
    """Instrument parameters configuration.

    Attributes:
        caglioti_u: Caglioti U parameter (FWHM²_inst = U·tan²θ + V·tanθ + W)
        caglioti_v: Caglioti V parameter
        caglioti_w: Caglioti W parameter (default: 0.003 → FWHM ≈ 0.055°)
        wavelength: X-ray wavelength in Angstroms (default: Cu Kα₁)

    """

    caglioti_u: float = 0.0
    caglioti_v: float = 0.0
    caglioti_w: float = 0.003  # → FWHM_inst ≈ 0.055°
    wavelength: float = field(default_factory=lambda: CU_KA1)

    def validate(self) -> None:
        """Validate instrument parameters."""
        if self.caglioti_w < 0:
            raise ValueError("Caglioti W must be non-negative")
        if self.caglioti_u < 0:
            raise ValueError("Caglioti U must be non-negative")
        if self.wavelength <= 0:
            raise ValueError("Wavelength must be positive")


@dataclass
class PeakDetectionConfig:
    """Peak detection parameters configuration.

    Attributes:
        peak_window: Search window around expected peak position (degrees)
        fitting_window: Window for detailed fitting diagnosis (degrees)
        min_intensity: Minimum peak intensity threshold (counts)

    """

    peak_window: float = 2.0  # degrees
    fitting_window: float = 2.5  # degrees (for diagnosis)
    min_intensity: float = 100  # counts

    def validate(self) -> None:
        """Validate peak detection parameters."""
        if self.peak_window <= 0:
            raise ValueError("peak_window must be positive")
        if self.fitting_window <= 0:
            raise ValueError("fitting_window must be positive")
        if self.min_intensity < 0:
            raise ValueError("min_intensity must be non-negative")


@dataclass
class ValidationConfig:
    """Validation thresholds configuration.

    Attributes:
        max_rwp: Maximum R_wp percentage for acceptable fit
        min_r_squared: Minimum R² for acceptable fit
        min_broadening_ratio: Minimum β_obs/β_inst ratio for reliable size
        min_reliable_size: Minimum detectable crystallite size (nm)
        max_reliable_size: Maximum detectable crystallite size (nm)

    """

    max_rwp: float = 10.0  # %
    min_r_squared: float = 0.95
    min_broadening_ratio: float = field(default_factory=lambda: MIN_BROADENING_RATIO)
    min_reliable_size: float = field(default_factory=lambda: MIN_RELIABLE_SIZE)
    max_reliable_size: float = field(default_factory=lambda: MAX_RELIABLE_SIZE)

    def validate(self) -> None:
        """Validate validation thresholds."""
        if self.max_rwp <= 0:
            raise ValueError("max_rwp must be positive")
        if not (0 <= self.min_r_squared <= 1):
            raise ValueError("min_r_squared must be between 0 and 1")
        if self.min_broadening_ratio <= 0:
            raise ValueError("min_broadening_ratio must be positive")


@dataclass
class VisualizationConfig:
    """Visualization parameters configuration.

    Attributes:
        dpi: Output resolution (dots per inch)
        figure_format: Output file format (png, svg, pdf)
        colormap: Color scheme name

    """

    dpi: int = 600
    figure_format: str = "png"
    colormap: str = "colorblind_safe"

    def validate(self) -> None:
        """Validate visualization parameters."""
        if self.dpi <= 0:
            raise ValueError("dpi must be positive")
        if self.figure_format not in ["png", "svg", "pdf", "jpg"]:
            warnings.warn(f"Uncommon figure format: {self.figure_format}")


@dataclass
class ParameterConfig:
    """Unified xrd_analysis parameter configuration.

    This is the main configuration class that aggregates all parameter categories.

    Attributes:
        instrument: Instrument-related parameters
        peak_detection: Peak detection parameters
        validation: Validation thresholds
        visualization: Visualization parameters

    Example:
        >>> # Use default configuration
        >>> config = ParameterConfig()
        >>>
        >>> # Customize specific parameters
        >>> config = ParameterConfig(
        ...     instrument=InstrumentConfig(caglioti_w=0.005),
        ...     validation=ValidationConfig(min_r_squared=0.99)
        ... )
        >>>
        >>> # Load from YAML file
        >>> config = ParameterConfig.from_yaml("custom_config.yaml")

    """

    instrument: InstrumentConfig = field(default_factory=InstrumentConfig)
    peak_detection: PeakDetectionConfig = field(default_factory=PeakDetectionConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)

    @classmethod
    def from_yaml(cls, filepath: Path | None = None) -> "ParameterConfig":
        """Load configuration from YAML file.

        Args:
            filepath: Path to YAML config file (default: "config.yaml")

        Returns:
            ParameterConfig instance with values from YAML

        Example:
            >>> config = ParameterConfig.from_yaml("my_config.yaml")

        """
        from xrd_analysis.core.config_loader import load_config

        if filepath is None:
            filepath = Path("config.yaml")

        try:
            # Fix: load_config expects Path object, not string
            yaml_config = load_config(filepath)
        except (FileNotFoundError, AttributeError):
            warnings.warn(f"Config file {filepath} not found, using defaults")
            return cls()

        # Extract instrument parameters
        instrument_dict = yaml_config.get("instrument", {})
        caglioti = instrument_dict.get("caglioti", {})

        instrument = InstrumentConfig(
            caglioti_u=caglioti.get("U", 0.0),
            caglioti_v=caglioti.get("V", 0.0),
            caglioti_w=caglioti.get("W", 0.003),
            wavelength=yaml_config.get("physical", {}).get("wavelength_ka1", CU_KA1),
        )

        # Extract peak detection parameters
        peak_fitting = yaml_config.get("peak_fitting", {})
        peak_detection = PeakDetectionConfig(
            peak_window=peak_fitting.get("peak_window", 2.0),
            min_intensity=peak_fitting.get("min_intensity", 100),
        )

        # Extract validation parameters
        validation_dict = yaml_config.get("validation", {})
        validation = ValidationConfig(
            max_rwp=validation_dict.get("max_rwp", 10.0),
            min_r_squared=validation_dict.get("min_r_squared", 0.95),
        )

        # Extract visualization parameters
        vis_dict = yaml_config.get("visualization", {})
        visualization = VisualizationConfig(
            dpi=vis_dict.get("dpi", 600),
            figure_format=vis_dict.get("format", "png"),
        )

        return cls(
            instrument=instrument,
            peak_detection=peak_detection,
            validation=validation,
            visualization=visualization,
        )

    def validate_all(self) -> None:
        """Validate all configuration parameters.

        Raises:
            ValueError: If any parameter is invalid

        """
        self.instrument.validate()
        self.peak_detection.validate()
        self.validation.validate()
        self.visualization.validate()

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to dictionary.

        Returns:
            Nested dictionary representation

        """
        return {
            "instrument": {
                "caglioti_u": self.instrument.caglioti_u,
                "caglioti_v": self.instrument.caglioti_v,
                "caglioti_w": self.instrument.caglioti_w,
                "wavelength": self.instrument.wavelength,
            },
            "peak_detection": {
                "peak_window": self.peak_detection.peak_window,
                "fitting_window": self.peak_detection.fitting_window,
                "min_intensity": self.peak_detection.min_intensity,
            },
            "validation": {
                "max_rwp": self.validation.max_rwp,
                "min_r_squared": self.validation.min_r_squared,
                "min_broadening_ratio": self.validation.min_broadening_ratio,
                "min_reliable_size": self.validation.min_reliable_size,
                "max_reliable_size": self.validation.max_reliable_size,
            },
            "visualization": {
                "dpi": self.visualization.dpi,
                "figure_format": self.visualization.figure_format,
                "colormap": self.visualization.colormap,
            },
        }


# =============================================================================
# Convenience Functions
# =============================================================================


def get_default_config() -> ParameterConfig:
    """Get default xrd_analysis configuration.

    Returns:
        ParameterConfig with default values

    """
    return ParameterConfig()


def load_config_from_file(filepath: str = "config.yaml") -> ParameterConfig:
    """Load configuration from file with fallback to defaults.

    Args:
        filepath: Path to config file

    Returns:
        ParameterConfig instance

    """
    try:
        return ParameterConfig.from_yaml(Path(filepath))
    except (FileNotFoundError, OSError, ValueError) as e:
        logger.warning(
            "Failed to load config from %s: %s. Using defaults.", filepath, e
        )
        return get_default_config()

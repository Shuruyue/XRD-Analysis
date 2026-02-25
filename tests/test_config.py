"""Tests for core/config.py configuration system."""

import warnings
import pytest
from pathlib import Path

from xrd_analysis.core.config import (
    InstrumentConfig,
    PeakDetectionConfig,
    ValidationConfig,
    VisualizationConfig,
    ParameterConfig,
    get_default_config,
    load_config_from_file,
)


class TestInstrumentConfig:
    """Tests for InstrumentConfig validation."""

    def test_default_values(self):
        cfg = InstrumentConfig()
        assert cfg.caglioti_w == 0.003
        assert cfg.wavelength > 0

    def test_validate_negative_w(self):
        cfg = InstrumentConfig(caglioti_w=-1)
        with pytest.raises(ValueError, match="non-negative"):
            cfg.validate()

    def test_validate_negative_wavelength(self):
        cfg = InstrumentConfig(wavelength=-1)
        with pytest.raises(ValueError, match="positive"):
            cfg.validate()

    def test_validate_negative_u(self):
        cfg = InstrumentConfig(caglioti_u=-1)
        with pytest.raises(ValueError, match="non-negative"):
            cfg.validate()


class TestPeakDetectionConfig:
    """Tests for PeakDetectionConfig validation."""

    def test_validate_negative_window(self):
        cfg = PeakDetectionConfig(peak_window=-1)
        with pytest.raises(ValueError):
            cfg.validate()

    def test_validate_negative_intensity(self):
        cfg = PeakDetectionConfig(min_intensity=-1)
        with pytest.raises(ValueError):
            cfg.validate()


class TestValidationConfig:
    """Tests for ValidationConfig validation."""

    def test_validate_r_squared_out_of_range(self):
        cfg = ValidationConfig(min_r_squared=2.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            cfg.validate()

    def test_validate_negative_rwp(self):
        cfg = ValidationConfig(max_rwp=-1)
        with pytest.raises(ValueError, match="positive"):
            cfg.validate()


class TestVisualizationConfig:
    """Tests for VisualizationConfig validation."""

    def test_validate_negative_dpi(self):
        cfg = VisualizationConfig(dpi=-1)
        with pytest.raises(ValueError, match="positive"):
            cfg.validate()

    def test_validate_uncommon_format(self):
        cfg = VisualizationConfig(figure_format="bmp")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cfg.validate()
            assert len(w) == 1
            assert "Uncommon" in str(w[0].message)


class TestParameterConfig:
    """Tests for ParameterConfig."""

    def test_validate_all_defaults(self):
        """Default config should pass validation."""
        cfg = ParameterConfig()
        cfg.validate_all()

    def test_to_dict(self):
        """Should produce nested dict with expected keys."""
        d = ParameterConfig().to_dict()
        assert "instrument" in d
        assert "peak_detection" in d
        assert "validation" in d
        assert "visualization" in d
        assert "caglioti_w" in d["instrument"]


class TestConvenienceFunctions:
    """Tests for module-level convenience functions."""

    def test_get_default_config(self):
        cfg = get_default_config()
        assert isinstance(cfg, ParameterConfig)

    def test_load_config_nonexistent_file(self):
        """Should return defaults when file doesn't exist."""
        cfg = load_config_from_file("nonexistent_config.yaml")
        assert isinstance(cfg, ParameterConfig)

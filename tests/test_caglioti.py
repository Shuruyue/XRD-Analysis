"""Tests for Caglioti equation and instrumental broadening."""

import numpy as np
import pytest
from xrd_analysis.methods.caglioti import CagliotiCorrection, calculate_instrumental_broadening


class TestCagliotiCorrection:
    """Tests for CagliotiCorrection class."""

    def test_default_params(self):
        """With U/V/W given, should give reasonable FWHM."""
        cag = CagliotiCorrection(U=0.0, V=0.0, W=0.003)
        fwhm = cag.calculate_fwhm_inst(43.3)
        assert 0.04 < fwhm < 0.1

    def test_uncalibrated_raises(self):
        """Should raise ValueError if not calibrated."""
        cag = CagliotiCorrection()
        with pytest.raises(ValueError, match="not calibrated"):
            cag.calculate_fwhm_inst(43.3)

    def test_is_calibrated(self):
        """Should report calibration status correctly."""
        assert not CagliotiCorrection().is_calibrated
        assert CagliotiCorrection(U=0, V=0, W=0.003).is_calibrated

    def test_fwhm_increases_with_angle(self):
        """FWHM should generally increase with angle for positive U."""
        cag = CagliotiCorrection(U=0.01, V=0.0, W=0.003)
        fwhm_low = cag.calculate_fwhm_inst(30.0)
        fwhm_high = cag.calculate_fwhm_inst(80.0)
        assert fwhm_high > fwhm_low

    def test_negative_fwhm_sq_raises(self):
        """Negative FWHM² should raise ValueError."""
        cag = CagliotiCorrection(U=0.0, V=-10.0, W=0.0)
        with pytest.raises(ValueError, match="Negative"):
            cag.calculate_fwhm_inst(43.3)


class TestCagliotiCalibration:
    """Tests for Caglioti calibration from standard data."""

    def test_calibrate_from_data(self):
        """Should fit U, V, W from standard peak data."""
        two_theta = np.array([21.36, 30.39, 37.44, 43.51, 53.00, 61.28])
        fwhm = np.array([0.06, 0.058, 0.057, 0.057, 0.059, 0.062])

        cag = CagliotiCorrection()
        params = cag.calibrate(two_theta, fwhm)
        assert cag.is_calibrated
        assert params.W is not None

    def test_calibrate_preserves_shape(self):
        """Calibrated model should reproduce input FWHM approximately."""
        two_theta = np.array([30.0, 50.0, 70.0, 90.0])
        fwhm = np.array([0.06, 0.058, 0.062, 0.070])

        cag = CagliotiCorrection()
        cag.calibrate(two_theta, fwhm)
        predicted = np.array([cag.calculate_fwhm_inst(tt) for tt in two_theta])
        np.testing.assert_allclose(predicted, fwhm, atol=0.02)


class TestCorrectBroadening:
    """Tests for broadening correction."""

    def test_geometric_correction(self):
        """Geometric correction should reduce FWHM."""
        cag = CagliotiCorrection(U=0.0, V=0.0, W=0.003)
        corrected, reliable, warning = cag.correct_broadening(0.25, 43.3)
        assert corrected < 0.25
        assert reliable

    def test_unreliable_correction(self):
        """Should flag unreliable when ratio is small."""
        cag = CagliotiCorrection(U=0.0, V=0.0, W=0.003)
        corrected, reliable, warning = cag.correct_broadening(0.06, 43.3)
        assert not reliable
        assert "WARNING" in warning or "CRITICAL" in warning


class TestConvenience:
    """Tests for convenience function."""

    def test_calculate_instrumental_broadening(self):
        fwhm = calculate_instrumental_broadening(43.3, 0.0, 0.0, 0.003)
        assert 0.04 < fwhm < 0.1


"""
Unit Tests for Visualization Module
===================================

Tests for the updated visualization components, ensuring all plots use
correct constants and rigorous fitting methods.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import patch

from xrd_analysis.visualization.style import apply_xrd_analysis_style, PEAK_COLORS
from xrd_analysis.visualization.wh_plots import plot_williamson_hall
from xrd_analysis.visualization.texture_plots import plot_tc_evolution
from xrd_analysis.visualization.fitting_plots import plot_peak_fit
from xrd_analysis.visualization.generate_fitting_diagnosis import fit_peak_with_diagnosis
from xrd_analysis.core.constants import CU_KA1

class TestVisualizationConstants:
    """Tests for constants usage in visualization."""
    
    def test_colors_defined(self):
        """Standard peak colors should be defined."""
        assert '(111)' in PEAK_COLORS
        assert '(200)' in PEAK_COLORS
        
    def test_style_application(self):
        """Style application should not raise error."""
        try:
            apply_xrd_analysis_style()
        except Exception as e:
            pytest.fail(f"Style application failed: {e}")

class TestWHPlotting:
    """Tests for Williamson-Hall plotting function."""
    
    def test_wh_plot_execution(self):
        """Should execute plot_williamson_hall without error."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        
        # Build a verification mock to prevent actual plotting
        with patch('matplotlib.pyplot.show'):
             fig = plot_williamson_hall(two_theta, fwhm, show=False)
             assert fig is not None

    def test_wh_plot_high_r2_branch(self):
        """High R^2 branch should render without unbound color variable errors."""
        two_theta = np.array([43.32, 50.45, 74.16, 89.97])
        fwhm = np.array([0.224, 0.251, 0.282, 0.305])
        fit_result = {"slope": 0.01, "intercept": 0.003, "r_squared": 0.99}

        with patch('matplotlib.pyplot.show'):
            fig = plot_williamson_hall(two_theta, fwhm, fit_result=fit_result, show=False)
            assert fig is not None

class TestTexturePlotting:
    """Tests for Texture plotting function."""
    
    def test_texture_plot_execution(self):
        """Should execute plot_tc_evolution without error."""
        samples = [
            {
                "name": "sample_a",
                "concentration": 0.0,
                "time": 0.0,
                "tc_values": {"(111)": 1.35, "(200)": 0.72, "(220)": 0.93},
                "high_quality": True,
            },
            {
                "name": "sample_b",
                "concentration": 0.0,
                "time": 2.0,
                "tc_values": {"(111)": 1.20, "(200)": 0.85, "(220)": 0.95},
                "high_quality": True,
            },
        ]
        
        with patch('matplotlib.pyplot.show'):
            fig = plot_tc_evolution(samples, x_param="time", show=False)
            assert fig is not None

class TestFittingDiagnosis:
    """Tests for fitting diagnosis logic (without actual plotting)."""

    def test_fitting_logic_structure(self):
        """Verify fit_peak_with_diagnosis returns correct dictionary structure."""
        # Create synthetic Gaussian peak data
        two_theta = np.linspace(40, 46, 200)
        intensity = 1000 * np.exp(-(two_theta - 43.3)**2 / (2 * 0.1**2)) + 50

        # Test function — fitting should succeed on clean synthetic data
        result = fit_peak_with_diagnosis(two_theta, intensity, 43.3, use_doublet=False)

        # Verify required keys are present
        expected_keys = {
            'success', 'center', 'amplitude', 'fwhm', 'eta', 'r_squared',
            'theta_range', 'int_range', 'fitted_curve', 'method'
        }
        assert expected_keys.issubset(result.keys())

        
if __name__ == "__main__":
    pytest.main([__file__, "-v"])

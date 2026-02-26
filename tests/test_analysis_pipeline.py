"""Tests for analysis pipeline coupling and implementation constraints."""

import numpy as np

from xrd_analysis.analysis import pipeline as pipeline_module
from xrd_analysis.analysis.pipeline import (
    AnalysisConfig,
    XRDAnalysisPipeline,
)


def test_find_peak_in_range_respects_window(monkeypatch):
    """Peak fitter should receive caller-provided search window."""
    captured = {}

    def fake_fit_peak_with_diagnosis(
        two_theta,
        intensity,
        expected_center,
        window=2.5,
        use_doublet=False,
        doublet_max_iterations=20000,
    ):
        captured["window"] = window
        captured["doublet_max_iterations"] = doublet_max_iterations
        return {
            "success": True,
            "r_squared": 0.99,
            "amplitude": 1000.0,
            "fwhm": 0.2,
            "center": expected_center,
            "eta": 0.5,
        }

    monkeypatch.setattr(
        "xrd_analysis.fitting.peak_fitter.fit_peak_with_diagnosis",
        fake_fit_peak_with_diagnosis,
    )

    two_theta = np.linspace(40.0, 50.0, 1000)
    intensity = np.full_like(two_theta, 100.0)

    peak = pipeline_module.find_peak_in_range(
        two_theta=two_theta,
        intensity=intensity,
        center=45.0,
        window=0.9,
        use_doublet_fitting=True,
        doublet_max_iterations=1234,
        min_fit_r_squared=0.8,
    )
    assert peak is not None
    assert captured["window"] == 0.9
    assert captured["doublet_max_iterations"] == 1234


def test_find_peaks_uses_configured_fitting_controls(monkeypatch):
    """Pipeline should forward fitting controls from AnalysisConfig."""
    calls = []

    def fake_find_peak_in_range(
        two_theta,
        intensity,
        center,
        window=2.5,
        use_doublet_fitting=True,
        doublet_max_iterations=20000,
        min_fit_r_squared=0.8,
    ):
        calls.append(
            {
                "window": window,
                "use_doublet_fitting": use_doublet_fitting,
                "doublet_max_iterations": doublet_max_iterations,
                "min_fit_r_squared": min_fit_r_squared,
            }
        )
        return None

    monkeypatch.setattr(pipeline_module, "find_peak_in_range", fake_find_peak_in_range)

    cfg = AnalysisConfig(
        peak_window=1.25,
        use_doublet_fitting=False,
        doublet_max_iterations=6789,
        min_fit_r_squared=0.93,
        EXPECTED_PEAKS={(1, 1, 1): 43.316},
    )
    pipeline = XRDAnalysisPipeline(cfg)
    _ = pipeline._find_peaks_from_data(
        two_theta=np.linspace(30.0, 90.0, 100),
        intensity=np.full(100, 100.0),
    )

    assert len(calls) == 1
    assert calls[0]["window"] == 1.25
    assert calls[0]["use_doublet_fitting"] is False
    assert calls[0]["doublet_max_iterations"] == 6789
    assert calls[0]["min_fit_r_squared"] == 0.93

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


def test_double_correction_guard_disables_kalpha_strip_for_doublet_fits():
    """Doublet fitting models Kα2, so Kα2 must not also be stripped beforehand."""
    cfg = AnalysisConfig(use_doublet_fitting=True, enable_kalpha_strip=True)
    pipeline = XRDAnalysisPipeline(cfg)
    assert pipeline.preprocessing.enable_kalpha_strip is False
    assert cfg.enable_kalpha_strip is True  # caller's config is not modified

    cfg = AnalysisConfig(use_doublet_fitting=False, enable_kalpha_strip=True)
    assert XRDAnalysisPipeline(cfg).preprocessing.enable_kalpha_strip is True


def test_average_size_skips_unreliable_error_and_nan():
    from xrd_analysis.analysis.pipeline import PeakData

    pipeline = XRDAnalysisPipeline(
        AnalysisConfig(caglioti_u=0.0, caglioti_v=0.0, caglioti_w=0.0064)  # β_inst = 0.08°
    )
    peaks = [
        PeakData(hkl=(1, 1, 1), two_theta=43.316, intensity=1000.0, fwhm=0.25, eta=0.5),
        # narrower than the instrument and Gaussian: nothing left after Voigt deconvolution
        PeakData(hkl=(2, 0, 0), two_theta=50.448, intensity=500.0, fwhm=0.05, eta=0.0),
    ]
    results, avg = pipeline._run_scherrer_analysis(peaks)
    assert np.isnan(results[1].size_nm)
    assert avg == results[0].size_nm

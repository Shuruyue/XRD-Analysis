"""Tests for CLI config-to-AnalysisConfig mapping."""

from pathlib import Path

from xrd_analysis.cli import _build_analysis_config


def test_build_analysis_config_maps_fitting_controls(tmp_path: Path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "\n".join(
            [
                "physical_constants:",
                "  wavelength: 1.540562",
                "fitting:",
                "  max_iterations: 32100",
                "  peak_detection:",
                "    min_height: 123",
                "    peak_window: 1.7",
                "validation:",
                "  min_r_squared: 0.91",
            ]
        ),
        encoding="utf-8",
    )

    result = _build_analysis_config(cfg)
    assert result.doublet_max_iterations == 32100
    assert result.min_intensity == 123.0
    assert result.peak_window == 1.7
    assert result.min_fit_r_squared == 0.91


def test_build_analysis_config_sanitizes_savgol_settings(tmp_path: Path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "\n".join(
            [
                "preprocessing:",
                "  smoothing:",
                "    enable: true",
                "    window_size: 10",
                "    poly_order: 20",
            ]
        ),
        encoding="utf-8",
    )

    result = _build_analysis_config(cfg)
    # Even window size should be shifted to odd.
    assert result.smoothing_window == 11
    # poly_order should be clipped to window_size - 1.
    assert result.smoothing_poly_order == 10

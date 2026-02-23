"""Tests for analysis pipeline coupling and methodology constraints."""

import numpy as np

from xrd_analysis.analysis.pipeline import (
    AnalysisConfig,
    PeakData,
    XRDAnalysisPipeline,
)
from xrd_analysis.methods.scherrer import ScherrerResult


def test_prepare_wh_input_uses_scherrer_corrected_fwhm():
    """W-H input must use Scherrer-corrected sample broadening."""
    pipeline = XRDAnalysisPipeline(AnalysisConfig())

    peaks = [
        PeakData(hkl=(1, 1, 1), two_theta=43.3, intensity=1000.0, fwhm=0.30),
        PeakData(hkl=(2, 0, 0), two_theta=50.4, intensity=700.0, fwhm=0.28),
        PeakData(hkl=(2, 2, 0), two_theta=74.1, intensity=500.0, fwhm=0.26),
    ]
    scherrer_results = [
        ScherrerResult(size_nm=20.0, size_angstrom=200.0, two_theta=43.3, fwhm_sample=0.21),
        ScherrerResult(size_nm=21.0, size_angstrom=210.0, two_theta=50.4, fwhm_sample=0.19),
        ScherrerResult(size_nm=22.0, size_angstrom=220.0, two_theta=74.1, fwhm_sample=0.17),
    ]

    two_theta_arr, fwhm_arr, hkl_list = pipeline._prepare_wh_input(peaks, scherrer_results)

    assert np.allclose(two_theta_arr, [43.3, 50.4, 74.1])
    assert np.allclose(fwhm_arr, [0.21, 0.19, 0.17])
    assert hkl_list == [(1, 1, 1), (2, 0, 0), (2, 2, 0)]


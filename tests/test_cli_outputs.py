"""Tests for CLI analysis output serialization/writing."""

from pathlib import Path

from xrd_analysis.analysis.pipeline import PeakData, PipelineResult
from xrd_analysis.cli import _serialize_pipeline_result, _write_analysis_outputs
from xrd_analysis.methods.scherrer import ScherrerResult, ValidityFlag


def _build_min_result() -> PipelineResult:
    result = PipelineResult(
        filepath="data/sample_a.txt",
        sample_name="sample_a",
        leveler_concentration=1.5,
    )
    result.peaks = [
        PeakData(
            hkl=(1, 1, 1),
            two_theta=43.316,
            intensity=1200.0,
            fwhm=0.24,
            area=2500.0,
            eta=0.42,
        )
    ]
    result.scherrer_results = [
        ScherrerResult(
            size_nm=34.2,
            size_angstrom=342.0,
            two_theta=43.316,
            hkl=(1, 1, 1),
            k_factor=0.855,
            fwhm_observed=0.24,
            fwhm_instrumental=0.08,
            fwhm_sample=0.226,
            validity_flag=ValidityFlag.VALID,
            is_reliable=True,
        )
    ]
    result.average_size_nm = 34.2
    return result


def test_serialize_pipeline_result_has_expected_structure():
    payload = _serialize_pipeline_result(_build_min_result())

    assert payload["sample_name"] == "sample_a"
    assert payload["peaks"][0]["hkl"] == "(111)"
    assert payload["scherrer"]["results"][0]["validity_flag"] == "VALID"
    assert payload["scherrer"]["average_size_nm"] == 34.2


def test_write_analysis_outputs_creates_summary_files(tmp_path: Path):
    outputs = _write_analysis_outputs([_build_min_result()], tmp_path)

    assert outputs["summary_json"].exists()
    assert outputs["summary_csv"].exists()
    assert (tmp_path / "sample_a_analysis.json").exists()

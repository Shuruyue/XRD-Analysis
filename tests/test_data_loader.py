"""Tests for data_loader module."""

import numpy as np
import pytest

from xrd_analysis.analysis.data_loader import load_bruker_txt, parse_filename


class TestLoadBrukerTxt:
    """Tests for Bruker TXT file loading."""

    def test_load_valid_file(self, tmp_path):
        """Should parse angle and intensity from [Data] section."""
        content = (
            "[Header]\ntype=XRD\n"
            "[Data]\nAngle,PSD\n"
            "10.000,100\n20.000,200\n30.000,300\n"
        )
        f = tmp_path / "sample.txt"
        f.write_text(content, encoding="utf-8")

        two_theta, intensity = load_bruker_txt(str(f))
        assert len(two_theta) == 3
        np.testing.assert_allclose(two_theta, [10, 20, 30])
        np.testing.assert_allclose(intensity, [100, 200, 300])

    def test_load_empty_file(self, tmp_path):
        """Should return empty arrays for empty file."""
        f = tmp_path / "empty.txt"
        f.write_text("", encoding="utf-8")

        two_theta, intensity = load_bruker_txt(str(f))
        assert len(two_theta) == 0
        assert len(intensity) == 0

    def test_load_no_data_section(self, tmp_path):
        """Should return empty arrays if no [Data] section."""
        f = tmp_path / "nodata.txt"
        f.write_text("[Header]\nfoo=bar\n", encoding="utf-8")

        two_theta, intensity = load_bruker_txt(str(f))
        assert len(two_theta) == 0

    def test_load_with_spaces(self, tmp_path):
        """Should handle space-separated data."""
        content = "[Data]\n10.0 100\n20.0 200\n"
        f = tmp_path / "spaces.txt"
        f.write_text(content, encoding="utf-8")

        two_theta, intensity = load_bruker_txt(str(f))
        assert len(two_theta) == 2

    def test_load_skips_invalid_lines(self, tmp_path):
        """Should skip lines that can't be parsed as floats."""
        content = "[Data]\n10.0,100\nabc,def\n20.0,200\n"
        f = tmp_path / "bad_lines.txt"
        f.write_text(content, encoding="utf-8")

        two_theta, intensity = load_bruker_txt(str(f))
        assert len(two_theta) == 2

    def test_file_not_found(self):
        """Should raise OSError for missing file."""
        with pytest.raises(OSError):
            load_bruker_txt("nonexistent_file.txt")


class TestParseFilename:
    """Tests for filename parsing."""

    def test_parse_standard_format(self):
        """Should extract concentration and time."""
        result = parse_filename("20251125_0ml_2h.txt")
        assert result["name"] == "20251125_0ml_2h"
        assert result["concentration_ml"] == 0.0
        assert result["time_hours"] == 2

    def test_parse_with_minutes(self):
        """Should handle hour + minute format."""
        result = parse_filename("20251125_4.5ml_0h_30min.txt")
        assert result["concentration_ml"] == 4.5
        assert result["time_hours"] == 0.5

    def test_parse_no_concentration(self):
        """Should return None for missing concentration."""
        result = parse_filename("sample_2h.txt")
        assert result["concentration_ml"] is None

    def test_parse_decimal_concentration(self):
        """Should parse decimal concentration."""
        result = parse_filename("20251125_4.5ml_2h.txt")
        assert result["concentration_ml"] == 4.5

    def test_parse_no_time(self):
        """Should return 0 for missing time."""
        result = parse_filename("random_file.txt")
        assert result["time_hours"] == 0

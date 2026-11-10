"""End-to-end run of the Quick Start on the synthetic patterns in examples/.

The patterns are simulated by scripts/make_example_data.py, so the instrument widths and the
crystallite size are known. Both commands run in subprocesses with an ASCII-only console
encoding, the strictest case (Windows consoles with cp950/cp1252 pipes failed on θ and ²).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
sys.path.insert(0, str(ROOT / "scripts"))

import make_example_data as truth  # noqa: E402

from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED  # noqa: E402


def _run(args: list[str], cwd: Path) -> str:
    env = dict(os.environ, PYTHONIOENCODING="ascii")
    out = subprocess.run(
        [sys.executable, *args], cwd=cwd, env=env, capture_output=True, text=True, check=False
    )
    assert out.returncode == 0, out.stdout + out.stderr
    return out.stdout


@pytest.fixture(scope="module")
def quick_start(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("quick_start")
    config = tmp / "config.yaml"
    config.write_text((ROOT / "config.yaml").read_text(encoding="utf-8"), encoding="utf-8")
    _run(
        [
            str(ROOT / "scripts" / "calibrate_instrument.py"),
            str(EXAMPLES / "lab6_synthetic.txt"),
            "-o",
            str(tmp / "calibration.yaml"),
            "--update-config",
            str(config),
        ],
        cwd=ROOT,
    )
    _run(
        [
            "-m",
            "xrd_analysis.cli",
            "analyze",
            str(EXAMPLES / "cu_synthetic.txt"),
            "--calibration",
            str(tmp / "calibration.yaml"),
            "-o",
            str(tmp / "out"),
        ],
        cwd=ROOT,
    )
    calibration = yaml.safe_load((tmp / "calibration.yaml").read_text(encoding="utf-8"))
    analysis = json.loads((tmp / "out" / "cu_synthetic_analysis.json").read_text(encoding="utf-8"))
    return calibration, analysis, yaml.safe_load(config.read_text(encoding="utf-8"))


def test_calibration_recovers_instrument_width(quick_start):
    calibration, _, config = quick_start
    c = calibration["caglioti"]
    assert calibration["n_peaks_used"] == 10
    for two_theta in (25.0, 45.0, 75.0):
        t = math.tan(math.radians(two_theta / 2))
        fitted = math.sqrt(c["U"] * t * t + c["V"] * t + c["W"])
        assert fitted == pytest.approx(truth.instrument_fwhm(two_theta), rel=0.02)
    assert config["instrument"]["caglioti"] == c  # --update-config wrote the same U, V, W


def test_analysis_finds_all_peaks_at_their_positions(quick_start):
    """Stripping Kα2 before a doublet fit shifted (111) by -0.04°; positions must stay put."""
    _, analysis, _ = quick_start
    results = analysis["scherrer"]["results"]
    assert len(results) == 5
    for r in results:
        hkl = tuple(int(ch) for ch in r["hkl"].strip("()"))
        assert r["two_theta"] == pytest.approx(CU_JCPDS_EXTENDED[hkl]["two_theta"], abs=0.01)
        assert r["validity_flag"] == "VALID"


def test_analysis_recovers_crystallite_size(quick_start):
    """Default preprocessing (smoothing, global background) biases the size a few % low."""
    _, analysis, _ = quick_start
    sizes = [r["size_nm"] for r in analysis["scherrer"]["results"]]
    assert all(s == pytest.approx(truth.D_TRUE_NM, rel=0.15) for s in sizes)
    assert analysis["scherrer"]["average_size_nm"] == pytest.approx(truth.D_TRUE_NM, rel=0.10)

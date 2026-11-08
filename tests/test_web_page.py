"""The interactive page web/crystallite-size.html computes the same numbers as the package."""

import importlib.util
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.skipif(shutil.which("node") is None, reason="Node.js is needed to run the page's JavaScript")
def test_page_physics_matches_package():
    spec = importlib.util.spec_from_file_location("verify_web_page", ROOT / "scripts" / "verify_web_page.py")
    verify = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verify)
    worst, flags_differ, table = verify.run(n_cases=500)
    assert max(worst.values()) <= verify.TOLERANCE, worst
    assert flags_differ == 0
    assert table < 0.002

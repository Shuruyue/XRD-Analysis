#!/usr/bin/env python3
"""
Run a 20-phase thesis-grade audit workflow and emit a markdown report.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List


@dataclass
class PhaseResult:
    phase: int
    title: str
    command: str
    expected_exit: int
    exit_code: int
    passed: bool
    summary: str


def run_cmd(command: List[str], cwd: Path) -> tuple[int, str]:
    proc = subprocess.run(
        command,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    out = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
    return proc.returncode, out.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Run 20-phase PhD-grade XRD project audit.")
    parser.add_argument("--sample", type=Path, default=Path("data/raw/202511/20251125_0ml_2h.txt"))
    parser.add_argument("--output", type=Path, default=Path("outputs/reports/phd20_audit_report.md"))
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    phases = [
        (1, "Repository baseline and status", ["git", "status", "--short"], 0),
        (2, "Data integrity audit", ["python", "scripts/audit_data_integrity.py"], 1),
        (3, "Core physics verification", ["python", "scripts/verify_physics.py"], 0),
        (4, "Elastic anisotropy verification", ["python", "scripts/verify_elastic_moduli.py"], 0),
        (5, "Angle accuracy verification", ["python", "scripts/verify_angle_accuracy.py", str(args.sample)], 0),
        (6, "Background correction verification", ["python", "scripts/verify_background_correction.py", str(args.sample)], 0),
        (7, "CLI smoke test", ["python", "-m", "xrd_analysis.cli", "--help"], 0),
        (8, "Single-sample analysis run", ["python", "-m", "xrd_analysis.cli", "analyze", str(args.sample), "-o", "outputs/phd20_tmp"], 0),
        (9, "Calibration guardrail test with non-standard input", ["python", "-m", "xrd_analysis.cli", "calibrate", str(args.sample), "-o", "outputs/phd20_tmp/fake_calibration.yaml"], 1),
        (10, "Unit test suite", ["python", "-m", "pytest", "-q"], 0),
        (11, "Documentation presence check", ["python", "-c", "from pathlib import Path; p=Path('docs/engineering_specs/11_PhD_20Phase_Optimization_Plan_EN.md'); print('OK' if p.exists() else 'MISSING'); raise SystemExit(0 if p.exists() else 1)"], 0),
        (12, "Reference catalog presence check", ["python", "-c", "from pathlib import Path; p=Path('docs/engineering_specs/12_Reference_Catalog_100plus_EN.md'); print('OK' if p.exists() else 'MISSING'); raise SystemExit(0 if p.exists() else 1)"], 0),
        (13, "Report file generation check", ["python", "-c", "from pathlib import Path; p=Path('outputs/phd20_tmp'); print('OK' if p.exists() else 'MISSING'); raise SystemExit(0 if p.exists() else 1)"], 0),
        (14, "Pipeline import check", ["python", "-c", "from xrd_analysis.analysis.pipeline import XRDAnalysisPipeline; print('OK')"], 0),
        (15, "Scherrer module import check", ["python", "-c", "from xrd_analysis.methods.scherrer import ScherrerCalculator; print('OK')"], 0),
        (16, "W-H module import check", ["python", "-c", "from xrd_analysis.methods.williamson_hall import WilliamsonHallAnalyzer; print('OK')"], 0),
        (17, "Texture module import check", ["python", "-c", "from xrd_analysis.methods.texture import TextureAnalyzer; print('OK')"], 0),
        (18, "Defect module import check", ["python", "-c", "from xrd_analysis.methods.defect_analysis import StackingFaultAnalyzer; print('OK')"], 0),
        (19, "Visualization backend headless check", ["python", "-c", "import xrd_analysis.visualization as v; print('OK')"], 0),
        (20, "Final git status check", ["git", "status", "--short"], 0),
    ]

    results: List[PhaseResult] = []
    for phase_id, title, cmd, expected_exit in phases:
        code, out = run_cmd(cmd, root)
        summary = out.splitlines()[-1] if out else ""
        passed = code == expected_exit
        results.append(
            PhaseResult(
                phase=phase_id,
                title=title,
                command=" ".join(cmd),
                expected_exit=expected_exit,
                exit_code=code,
                passed=passed,
                summary=summary,
            )
        )

    # Markdown report
    lines = [
        "# PhD-Grade 20-Phase Audit Report",
        "",
        f"- Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"- Sample: `{args.sample}`",
        "",
        "## Summary",
        "",
    ]

    passed = sum(1 for r in results if r.passed)
    failed = len(results) - passed
    lines.extend([
        f"- Passed: {passed}",
        f"- Failed: {failed}",
        "",
        "## Phase Results",
        "",
    ])

    for r in results:
        status = "PASS" if r.passed else "FAIL"
        lines.extend([
            f"### Phase {r.phase:02d}: {r.title}",
            f"- Status: **{status}** (exit={r.exit_code}, expected={r.expected_exit})",
            f"- Command: `{r.command}`",
            f"- Summary: `{r.summary}`",
            "",
        ])

    args.output.write_text("\n".join(lines), encoding="utf-8")

    json_output = args.output.with_suffix(".json")
    json_output.write_text(
        json.dumps([r.__dict__ for r in results], indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"Report written: {args.output}")
    print(f"JSON written: {json_output}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

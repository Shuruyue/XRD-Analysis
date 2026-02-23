#!/usr/bin/env python3
"""
Instrument calibration helper for XRD-Analysis.

Workflow:
1) Run CLI calibration on a standard scan (LaB6/Si)
2) Optionally write fitted U/V/W back to config.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from types import SimpleNamespace

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from xrd_analysis.cli import _run_calibrate


def update_config(config_path: Path, calibration_path: Path) -> None:
    """Write calibrated Caglioti parameters into config YAML."""
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    if not calibration_path.exists():
        raise FileNotFoundError(f"Calibration file not found: {calibration_path}")

    with config_path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    with calibration_path.open("r", encoding="utf-8") as f:
        calibration = yaml.safe_load(f) or {}

    caglioti = calibration.get("caglioti")
    if not isinstance(caglioti, dict):
        raise ValueError("Invalid calibration file: missing 'caglioti' section")

    instrument = config.setdefault("instrument", {})
    instrument_caglioti = instrument.setdefault("caglioti", {})
    instrument_caglioti["U"] = float(caglioti["U"])
    instrument_caglioti["V"] = float(caglioti["V"])
    instrument_caglioti["W"] = float(caglioti["W"])

    with config_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, default_flow_style=False, sort_keys=False)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calibrate instrument using standard XRD scan and optionally update config.yaml."
    )
    parser.add_argument("standard", type=Path, help="Path to standard scan file (e.g. LaB6)")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("outputs/calibration.yaml"),
        help="Calibration output YAML path",
    )
    parser.add_argument(
        "--update-config",
        type=Path,
        default=None,
        help="Config YAML to update in-place with calibrated U/V/W",
    )
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    run_args = SimpleNamespace(standard=args.standard, output=args.output)
    code = _run_calibrate(run_args)
    if code != 0:
        return code

    if args.update_config is not None:
        update_config(args.update_config, args.output)
        print(f"Updated config with calibrated U/V/W: {args.update_config}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

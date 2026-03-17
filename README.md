# XRD-Analysis

[![CI](https://github.com/Shuruyue/XRD-Analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/Shuruyue/XRD-Analysis/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## XRD Crystallite Size Analysis System

Automated crystallite size analysis system using True Voigt + Kalpha doublet fitting and instrumental correction.

---

## Project Overview

xrd_analysis is an automated XRD analysis system that addresses common issues in manual Scherrer size calculations:

- Subjective baseline selection errors
- Inaccurate instrumental broadening correction
- Incorrect peak profile function selection

Default constants and examples in this repository target Cu scans, but the workflow itself is not limited to electroplated Cu.

### Core Features

- **True Voigt + Kalpha Doublet Fitting**: Primary path uses physically rigorous line shape with Cu Kalpha1/Kalpha2 coupling
- **Pseudo-Voigt Fallback**: Used only when strict doublet fitting cannot converge
- **Double-Correction Guard**: When doublet fitting is enabled, Kalpha2 pre-stripping is automatically disabled to avoid over-correction
- **Caglioti Instrumental Correction**: Full-angle instrumental width correction using NIST SRM 660c (LaB6)
- **High Precision**: Applicable to crystallite sizes in the 2-100 nm range

---

## Quick Start

### Installation

```bash
pip install -e .
```

Or for development:

```bash
pip install -e ".[dev]"
```

### CLI Usage

```bash
# Run analysis on a single file
xrd-analysis analyze data/raw/202511/20251125_0ml_2h.txt -o outputs/

# Instrument calibration (LaB6)
xrd-analysis calibrate data/standards/lab6_standard.txt -o calibration.yaml

# Practical calibration + update config
python scripts/calibrate_instrument.py data/standards/lab6_standard.txt \
  -o outputs/calibration.yaml \
  --update-config config.yaml
```

`xrd-analysis analyze` now writes:
- Per-sample machine-readable result: `outputs/<sample>_analysis.json`
- Batch summaries: `outputs/analysis_summary.json` and `outputs/analysis_summary.csv`

### Validation Scripts

```bash
# Verify 2-theta peak positions (Bragg's Law)
python scripts/verify_physics.py

# Verify peak-angle correctness (zero-shift sanity check)
python scripts/verify_angle_accuracy.py data/raw/202511/20251125_0ml_2h.txt

# Verify background correction statistics
python scripts/verify_background_correction.py data/raw/202511/20251125_0ml_2h.txt

# Instrument calibration helper (with fit-quality diagnostics)
python scripts/calibrate_instrument.py data/standards/lab6_standard.txt -o outputs/calibration.yaml

```

---

## Project Structure

xrd_analysis/
├── config.yaml              # Global configuration
├── pyproject.toml           # Python project configuration
├── data/                    # Data directory
│   ├── raw/                 # Raw XRD data
│   ├── standards/           # NIST standard data
│   └── processed/           # Preprocessed data
├── xrd_analysis/                  # Core source code
│   ├── core/                # Physical constants and crystallography
│   ├── preprocessing/       # Data preprocessing modules
│   ├── fitting/             # Peak fitting core
│   ├── methods/             # Analysis methods (Scherrer, W-H, Texture)
│   ├── validation/          # Error analysis and validation
│   └── visualization/       # Plotting
├── scripts/                 # Verification & Utility scripts
├── tests/                   # Unit tests
├── references/              # Literature references (PDFs)
└── outputs/                 # Output directory (Ignored by Git)

---

## Theoretical Background

### True Voigt Peak Profile

I(2theta) = A * V(2theta; sigma, gamma) + Background

where V is the Voigt profile (Gaussian-Lorentzian convolution via Faddeeva function).

For Cu radiation, Kalpha2 is constrained from Kalpha1 by Bragg shift and fixed intensity ratio 0.5.

### Caglioti Equation

FWHM_inst^2 = U * tan^2(theta) + V * tan(theta) + W

### Scherrer Equation

D = K * lambda / (beta * cos(theta))

---

## References

1. Bearden, J.A. (1967). X-Ray Wavelengths. Rev. Mod. Phys. 39, 78-124.
2. Langford, J.I. & Wilson, A.J.C. (1978). Scherrer after Sixty Years. J. Appl. Cryst. 11, 102-113.
3. NIST Standard Reference Material 660c (LaB6)

---

## License

MIT License

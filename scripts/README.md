# XRD-Analysis Verification & Utility Scripts
## XRD-Analysis 驗證與公用腳本

This directory contains standalone scripts for validating physical calculations, testing instrument calibration, and generating debug plots.
此目錄包含用於驗證物理計算、測試儀器校準和生成除錯圖表的獨立腳本。

### Verification Scripts (驗證腳本)

#### `calibrate_instrument.py`
**Purpose**: Practical instrument calibration using standard scans.
**Uses**:
- Runs Caglioti fit (`U`, `V`, `W`) from LaB6/Si standard peaks.
- Exports calibration diagnostics to YAML.
- Optionally updates `config.yaml` in-place.

**Example**:
```bash
python scripts/calibrate_instrument.py data/standards/lab6_standard.txt \
  -o outputs/calibration.yaml \
  --update-config config.yaml
```

#### `verify_physics.py`
**Purpose**: Validates basic crystallographic formulas.
**Uses**:
- Checks Bragg's Law calculations.
- Verifies interplanar spacing ($d_{hkl}$) for Cubic systems.
- Confirms wavelength constants (Cu K$\alpha$).

#### `verify_elastic_moduli.py`
**Purpose**: Verification of Young's Modulus calculations.
**Uses**:
- Calculates directional Young's Modulus ($E_{hkl}$) for Copper.
- Verifies Reuss/Voigt/Hill averaging methods.

#### `verify_angle_accuracy.py`
**Purpose**: Verifies 2θ correctness from detected peak positions.
**Uses**:
- Compares measured Cu peak positions against expected references.
- Reports mean/RMSE/max peak-angle offsets.
- Flags zero-shift/specimen displacement risk when offset is too large.

#### `verify_background_correction.py`
**Purpose**: Verifies preprocessing background subtraction behavior.
**Uses**:
- Reports raw/background/corrected intensity statistics.
- Quantifies background-to-signal ratio.
- Confirms no negative-intensity artifacts remain after correction.

### Plot Generation (繪圖工具)

#### `generate_fitting_diagnosis.py`
**Purpose**: Generates visual diagnostics for peak fitting.
**Output**:
- `outputs/fitting_diagnosis/`: Individual peak fits and residual plots.

#### `generate_fwhm_plots.py`
**Purpose**: Plots Full Width at Half Maximum (FWHM) evolution.
**Output**:
- `outputs/fwhm_plots/`: Williamson-Hall style FWHM vs $sin(\theta)$ plots.

### Usage Example

Run scripts from the project root directory:

```bash
# Run physics verification
python scripts/verify_physics.py

# Generate diagnosis plots
python scripts/generate_fitting_diagnosis.py
```

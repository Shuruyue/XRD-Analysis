# Thesis 5-Phase Audit and Optimization

**Document ID**: DOC-ENG-09  
**Scope**: `xrd_analysis` constants, fitting path, methodology formulas, and cross-module consistency  
**Date**: 2026-02-23

---

## Phase 1: Core Constants and Material Data Verification

### Checked items
- Cu Kalpha constants: `CU_KA1`, `CU_KA2`, `CU_KA_AVG`, `KA2_KA1_RATIO`
- Cu FCC reference data: `CU_CRYSTAL.lattice_constant`, `CU_JCPDS_EXTENDED`
- LaB6 calibration peak list: `LAB6_STANDARD_PEAKS`
- Elastic constants and directional properties: `E_hkl`, `nu_hkl`

### Implemented optimizations
- Added executable verification scripts:
  - `scripts/verify_physics.py`
  - `scripts/verify_elastic_moduli.py`
- Updated LaB6 calibration table to Bragg-calculated values using SRM 660c lattice parameter and Cu Kalpha1:
  - file: `xrd_analysis/cli.py`

### Outcome
- Constant linkage is now testable by script, not only by static comments.
- Calibration references are internally consistent with Bragg-law calculation.

---

## Phase 2: True Voigt / Kalpha Doublet Fitting Path Hardening

### Target requirement
Primary fitting should be True Voigt + Kalpha doublet, with lower-accuracy methods only as fallback.

### Implemented optimizations
- Clarified workflow in docs/config:
  - `README.md`
  - `config.yaml` (`fitting.function: true_voigt_doublet`)
- Confirmed code path:
  - primary: `xrd_analysis/fitting/ka_doublet.py` + `TrueVoigt`
  - fallback: pseudo-Voigt in `xrd_analysis/fitting/peak_fitter.py`

### Outcome
- Method hierarchy is explicit and aligned with thesis intent.

---

## Phase 3: Methodology Formula and Constraint Corrections

### Scherrer
- Added explicit correction mode control (`auto`, `quadratic`, `voigt`) in `ScherrerCalculator`.
- Fixed `eta` parameter plumbing (`eta_observed`, `eta_instrumental`) so pipeline-provided values are truly used.
- Kept textbook quadratic subtraction available for documented hand-calculation scenarios.
- File: `xrd_analysis/methods/scherrer.py`

### Pipeline parameter coupling
- Scherrer analyzer now receives:
  - `wavelength`
  - `use_cubic_habit`
  - full Caglioti tuple `(U, V, W)`
- Instrumental broadening now computed angle-dependently in Scherrer (instead of using only `sqrt(W)`).
- File: `xrd_analysis/analysis/pipeline.py`

### Defect / stress method consistency
- Corrected outdated Warren G wording (`-20` -> `-7.897`) in stack-fault docs/comments.
- Corrected residual stress formula wording to match implemented plane-stress model:
  - `sigma = -E/(2nu) * (d-d0)/d0`
- Files:
  - `xrd_analysis/methods/defect_analysis.py`
  - `xrd_analysis/core/copper_crystal.py`

---

## Phase 4: Code-Doc Cross-Link Consistency

### Implemented optimizations
- Updated API constant note for Scherrer K meaning:
  - `docs/API_REFERENCE.md`
- Fixed CLI command examples to match actual parser behavior:
  - `docs/API_REFERENCE.md`
  - `README.md`
- Corrected grain-size engineering spec from "Integral Breadth implementation" to actual current code path (corrected FWHM + Voigt-aware subtraction):
  - `docs/engineering_specs/02_Grain_Size_Algorithm.md`
  - `docs/engineering_specs/02_Grain_Size_Algorithm_EN.md`

### Outcome
- Formula descriptions, CLI examples, and implementation are now mutually consistent.

---

## Phase 5: Regression Validation and Reproducibility

### Recommended validation commands
```bash
python scripts/verify_physics.py
python scripts/verify_elastic_moduli.py
python scripts/calibrate_instrument.py data/standards/lab6_standard.txt -o outputs/calibration.yaml --update-config config.yaml
python -m pytest -q
python -m xrd_analysis.cli --help
```

### Known residual risk
- Existing Scherrer-document-example tests historically failed in this repo before this audit due to formula-path mismatch.
- This audit introduces explicit correction-mode control to remove hidden ambiguity.
- Practical instrument calibration still requires an actual standard scan file in `data/standards/`.
- CLI calibration no longer falls back to placeholder Caglioti values when data are insufficient.

---

## Practical Method Constraints (for thesis chapter)

- **Scherrer**: reliable only when `beta_obs / beta_inst >= 1.2`
- **W-H**: minimum 3 peaks required for linear regression
- **Texture TC**: at least 2 valid peaks required
- **Stacking fault alpha**: requires both (111) and (200) peak positions
- **Residual stress**: current formula assumes plane stress and equi-biaxial in-plane stress
- **Background correction**: verify mean background ratio and negative-point handling (`scripts/verify_background_correction.py`)
- **Angle correctness**: verify peak-position offsets vs reference (`scripts/verify_angle_accuracy.py`, recommended max |Δ2θ| <= 0.10°)

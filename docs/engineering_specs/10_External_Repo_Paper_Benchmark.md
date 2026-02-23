# External Repo/Paper Benchmark and Practical Calibration

**Document ID**: DOC-ENG-10  
**Date**: 2026-02-23  
**Scope**: cross-check core algorithms against major XRD toolchains and primary literature; define practical calibration workflow.

---

## 1. Benchmark Targets

This project was cross-checked against:

- **GSAS-II** docs and repo  
  - https://gsas-ii.readthedocs.io/en/latest/objvarorg.html  
  - https://github.com/AdvancedPhotonSource/GSAS-II
- **FullProf** docs  
  - https://www.ill.eu/sites/fullprof/php/programs.html  
  - https://www.ill.eu/sites/fullprof/php/fullprof_News_2003.htm
- **MAUD** repo  
  - https://github.com/luttero/maud
- **HEXRD** docs  
  - https://hexrd.readthedocs.io/en/stable/
- **xrayutilities** docs  
  - https://xrayutilities.sourceforge.io/simulations.html#x-ray-line-profile-analysis

Primary papers/data:

- Scherrer review: Langford & Wilson (1978), J. Appl. Cryst., DOI: 10.1107/S0021889878012856  
  - https://journals.iucr.org/paper?a17087=
- Size/strain profile-refinement basis: de Keijser et al. (1983), J. Appl. Cryst., DOI: 10.1107/S0021889883010493  
  - https://journals.iucr.org/paper?a22481=
- TCH pseudo-Voigt profile model: Thompson, Cox, Hastings (1987), DOI: 10.1107/S0021889887087090  
  - https://www.osti.gov/etdeweb/biblio/6227400
- GSAS-II profile-term mapping from FPA (NIST/JAC): DOI: 10.1107/S1600576722001169  
  - https://www.nist.gov/publications/determination-physically-based-pseudo-voigt-powder-diffraction-profile-terms
- Cu elastic constants source: Ledbetter & Naimon (1974), DOI: 10.1063/1.3253150  
  - https://pubs.aip.org/jpcrd/article/3/4/897/1046636/Elastic-Properties-of-Metals-and-Alloys-II-Copper
- NIST SRM 660c reference landing/archive  
  - https://www-s.nist.gov/srmors/view_cert.cfm?srm=660c  
  - https://www-s.nist.gov/srmors/archive/view_detail.cfm?srm=660c

---

## 2. Cross-Check Summary

### 2.1 True Voigt + Kα doublet
- Project primary path (`DoubletFitter` + `TrueVoigt`) is aligned with modern practice that avoids pseudo-Voigt-only approximation when possible.
- Kα2/Kα1 fixed intensity ratio `0.5` is consistent with common refinement parameterization (GSAS-II docs include `I(L2)/I(L1)` as explicit ratio term).
- Kα2 position is constrained from Bragg shift using wavelength ratio; this is physically consistent.

### 2.2 Instrumental broadening model
- Project uses Caglioti form: `FWHM_inst^2 = U tan^2(theta) + V tan(theta) + W`.
- This form is consistent with mainstream profile parameterization used in GSAS-II/FullProf-style workflows.

### 2.3 Scherrer / size-strain path
- Project Scherrer path supports explicit correction mode and reliability thresholding.
- Deconvolution logic is conceptually aligned with de Keijser component-separation philosophy (G/L treatment and recombination assumptions explicitly documented).

### 2.4 Pseudo-Voigt fallback
- TCH-style pseudo-Voigt remains fallback-only in this project, not the primary path.
- This hierarchy is suitable for thesis claims emphasizing precision.

### 2.5 Stress/elastic anisotropy
- Directional `E_hkl` and `nu_hkl` use Ledbetter/Naimon-class stiffness basis and are internally validated by scripts.

---

## 3. Practical Calibration Workflow (Required for Thesis-Grade Results)

### 3.1 Data requirement
- Use a **standard material scan** (LaB6 SRM 660c or equivalent Si standard) measured under the same optics/slits/scan conditions as samples.
- Place file under `data/standards/` (currently repository has no real standard scan file committed).

### 3.2 Command
```bash
python scripts/calibrate_instrument.py data/standards/lab6_standard.txt \
  -o outputs/calibration.yaml \
  --update-config config.yaml
```

### 3.3 Acceptance criteria
- Peaks used: `>= 5`
- Fit quality: `R^2(FWHM^2)` should be high (recommended `>= 0.98`)
- No non-positive `FWHM^2` region in practical 2θ range
- Calibration file stores traceable diagnostics and peak list

### 3.4 Current status (2026-02-23)
- `data/standards/` does not include a real standard scan yet.
- Therefore true instrument calibration cannot be finalized until standard data is provided.

---

## 4. Core Analysis Flow (Thesis Layout)

Implemented report flow:

1. Load data
2. Preprocessing (smoothing, background correction, Kα2 strip)
3. Peak finding around expected positions
4. Scherrer
5. Williamson-Hall
6. Texture
7. Defect/stress diagnostics
8. Text report output

Methodology quality checkpoints:

- Background correction metrics in report + `scripts/verify_background_correction.py`
- Angle correctness metrics in report + `scripts/verify_angle_accuracy.py`
- If angle max offset exceeds recommended bound (0.10°), report adds warning

---

## 5. Implemented Optimizations in This Repo

- Calibration command now **fails explicitly** when peak count is insufficient (no silent fallback to placeholder U/V/W).
- Calibration YAML now includes:
  - timestamp
  - fit-quality metrics (`R^2`, RMSE, max residual)
  - peak list used in regression
- Added script: `scripts/calibrate_instrument.py` (optional direct config update).

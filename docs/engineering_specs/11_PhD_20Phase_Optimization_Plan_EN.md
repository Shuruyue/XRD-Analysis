# PhD-Grade XRD Optimization Plan (20 Phases)

**Document ID**: DOC-ENG-11  
**Date**: 2026-02-23  
**Objective**: Upgrade this bachelor-thesis XRD codebase toward PhD-grade methodology rigor, software robustness, and reproducible reporting.

---

## Scope

This plan treats the project as a full thesis workflow audit, not only a coding task.  
Core flow under review:

`Load data -> peak search around expected positions -> Scherrer -> Williamson-Hall -> Texture -> Defect -> text report`

Additional mandatory review tracks:

- background correction methodology
- angle correctness verification
- instrument calibration reliability
- cross-module consistency
- English-only technical documentation artifacts

---

## 20 Phases

1. **Baseline Freeze and Risk Register**  
   Freeze repository baseline, collect known issues, define high-risk scientific assumptions.

2. **Thesis QA Rubric Construction**  
   Define acceptance criteria for data quality, model validity, code quality, and reproducibility.

3. **External Repo Reconnaissance**  
   Map feature coverage against open-source XRD ecosystems, including `XRD-AutoAnalyzer`.

4. **Primary Literature Expansion**  
   Expand reference corpus with method-defining papers and standards.

5. **Reference Quality Screening**  
   Keep primary sources; classify by method area (line profile, strain, texture, stress, calibration).

6. **Data Integrity Audit**  
   Check raw/standard data availability, naming, metadata consistency, and provenance readiness.

7. **Preprocessing Validation**  
   Validate smoothing/background/Kalpha workflows and edge-case behavior.

8. **Peak Model Audit**  
   Validate True Voigt + Kalpha doublet primary path and fallback boundaries.

9. **Calibration Audit**  
   Enforce calibration guardrails, diagnostics, and rejection of non-standard input.

10. **Scherrer Pathway Audit**  
    Verify deconvolution assumptions, K constants, unit conversions, and validity flags.

11. **Williamson-Hall Pathway Audit**  
    Ensure W-H uses sample broadening after instrumental correction, not raw FWHM.

12. **Texture Pathway Audit**  
    Validate TC formula assumptions, normalization, and minimum-peak requirements.

13. **Defect/Stress Pathway Audit**  
    Verify stacking-fault and stress formulas, limitations, and directional elasticity usage.

14. **Cross-Module Coupling Audit**  
    Verify parameter flow across preprocessing, fitting, methods, reporting.

15. **Numerical Reliability Hardening**  
    Tighten bounds/tolerances, expose diagnostics, reduce silent failure paths.

16. **Testing Architecture Upgrade**  
    Add/extend unit and integration checks for methodology-coupling behavior.

17. **Maintainability Refactor**  
    Improve architecture cohesion, API clarity, and failure messages.

18. **English Documentation Rewrite**  
    Produce English-only technical artifacts for thesis audit and reproducibility.

19. **Reproducible Benchmark Packaging**  
    Build scripted multi-phase audit runner and machine-readable outputs.

20. **Final Readiness Gate**  
    Produce pass/fail matrix, unresolved risks, and publication-readiness checklist.

---

## Acceptance Milestones

- M1: Pipeline scientifically consistent with declared methodology.
- M2: Calibration cannot silently succeed on non-standard data.
- M3: Report includes preprocessing, angle-verification, and method-limit diagnostics.
- M4: Reference catalog expanded by at least 100 additional entries.
- M5: One-command audit run produces reproducible artifacts.


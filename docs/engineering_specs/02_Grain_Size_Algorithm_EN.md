# Chapter 2: Crystallite Size Analysis

**Document ID**: DOC-ENG-02-EN
**Module Path**: `xrd_analysis.methods.scherrer.ScherrerCalculator`
**Physics Core**: FWHM-based Scherrer Equation with Voigt-aware broadening correction

---

## 2.1 Theoretical Gold Standard: Fourier Analysis (Warren-Averbach)
The most rigorous definition of "Crystallite Size" is derived from the Fourier transform coefficients of the diffraction peak.

**Definition**:
The derivative of the Fourier coefficient $A_n$ at $n \to 0$ gives the **Area-Weighted Column Length** $\langle L \rangle_{area}$:
$$ -\frac{dA_n^{size}}{dn} \bigg|_{n=0} = \frac{1}{\langle L \rangle_{area}} $$

*   **Pros**: Perfectly separates size and strain, independent of peak shape assumptions.
*   **Cons**: Requires multiple diffraction orders (e.g., (111) **AND** (222)); extremely sensitive to tail noise.

---

## 2.2 Implementation: FWHM Scherrer Method (Code Path)
**Project Algorithm**: `xrd_analysis.methods.scherrer.ScherrerCalculator`

Current implementation uses Scherrer equation based on **corrected FWHM**.
When `eta_observed` is available from peak fitting, the code uses Voigt component subtraction.
If `eta_observed` is unavailable, it falls back to textbook quadratic subtraction.

**Formula**:
$$ D_{v} = \frac{K \lambda}{\beta \cos\theta} $$
*   $D_v$: Volume-weighted average crystallite size.
*   $\beta$: Corrected peak width (FWHM) in radians.

### Correction Strategy Used in Code
1.  **Voigt component subtraction (preferred)**: Decompose observed and instrumental widths into Gaussian/Lorentzian components, subtract separately, then recombine with Olivero-Longbothum approximation.
2.  **Quadratic subtraction (fallback)**: $\beta = \sqrt{\beta_{obs}^2-\beta_{inst}^2}$ when no reliable `eta` is available.

**Method Comparison from Image**:
*   **Method 1 (Single PV)**: Ignores $K\alpha_2$ splitting. FWHM is artificially widened by $K\alpha_2$ $\to$ Grain size severely underestimated (Artifact).
*   **Method 2 (Stripping)**: Mathematically subtracts $K\alpha_2$. High noise, prone to peak distortion.
*   **Method 3 (Doublet Fit - Adopted)**: Uses True Voigt Doublet Model for precise separation and provides `eta_observed` for Scherrer correction.

---

## 2.3 Geometric Shape Factor (K)

The Scherrer constant $K$ depends on **Grain Shape** and **Diffraction Direction (hkl)**.
This project rejects the "Spherical Assumption ($K=0.9$)" because electroplated copper grows in columns.

**Implementation: Cubic Habit Model**
Based on **Langford & Wilson (1978)**, assuming cubic grains:

*   **(111) Diffraction**: Line of sight along the body diagonal. Projection is hexagonal.
    *   **$K = 0.855$**
*   **(200) Diffraction**: Line of sight along the edge. Projection is square.
    *   **$K = 0.886$**

### Numerical Trace

To demonstrate the impact of $K$ value selection, we simulate a real scenario.
**Input**:
*   $\beta = 0.3^\circ$ (corrected FWHM) $\approx 0.005236$ rad
*   $\theta = 21.65^\circ$ (for 111)
*   $\lambda = 1.5406 \mathring{A}$

**Scenario A: Using Generic K=0.9 (Spherical Assumption)**
$$ D = \frac{0.9 \times 1.5406}{0.005236 \times \cos(21.65^\circ)} = \frac{1.3865}{0.00486} = \mathbf{285.2 \mathring{A}} $$

**Scenario B: Using Correct K=0.855 (Cubic (111) Correction)**
$$ D = \frac{0.855 \times 1.5406}{0.005236 \times \cos(21.65^\circ)} = \frac{1.3172}{0.00486} = \mathbf{271.0 \mathring{A}} $$

**Engineering Conclusion**:
Using the generic K value introduces a **5.2% Systematic Error** (285 vs 271).
For precision process monitoring, this 5% error could lead to misjudgment of the Leveler's grain refining efficiency. Therefore, we insist on using $K=0.855$.

---

## 2.4 Constants Data Table

| Parameter | Value | Unit | Description | Source / DOI |
| :--- | :--- | :--- | :--- | :--- |
| $\lambda$ | **1.540562** | Å | Cu Kα1 Wavelength | Bearden (1967) [10.1103/RevModPhys.39.78](https://doi.org/10.1103/RevModPhys.39.78) |
| $K_{111}$ | **0.855** | - | Cubic (111) Shape Factor | Langford (1978) [10.1107/S0021889878012844](https://doi.org/10.1107/S0021889878012844) |
| $K_{200}$ | **0.886** | - | Cubic (200) Shape Factor | Langford (1978) |
| Method | **FWHM (corrected)** | - | Width definition used in code path | Scherrer + Voigt correction |

**References**:
1.  **Langford, J. I., & Wilson, A. J. C. (1978)**. "Scherrer after sixty years". *J. Appl. Cryst.*, 11, 102-113.
2.  **Olivero, J. J., & Longbothum, R. L. (1977)**. "Empirical fits to the Voigt line width". *JQSRT*, 17, 233.

# Chapter 8: Algorithm Principles & Original References

**Document ID**: DOC-ENG-08-EN
**Abstract**: Per engineering requirements, this document strictly traces the "Original Source" of all mathematical models used, providing DOIs for academic verification. We cite only the original papers or seminal textbooks that proposed the mathematical forms, avoiding secondary sources or presentations.

---

## 8.1 Voigt Function & Faddeeva Algorithm

### Mathematical Form
The Voigt profile $V(x; \sigma, \gamma)$ is calculated via the real part of the Complex Error Function $w(z)$:

$$ V(x) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}}, \quad z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$

### Original Sources

#### 1. Physics Definition
*   **Paper**: "Spectrum Line Profiles: The Voigt Function"
*   **Author**: B. H. Armstrong
*   **Journal**: *Journal of Quantitative Spectroscopy and Radiative Transfer*, Vol 7, Issue 1, pp 61-88.
*   **Year**: 1967
*   **DOI**: [10.1016/0022-4073(67)90057-X](https://doi.org/10.1016/0022-4073(67)90057-X)
*   **Core Contribution**: Provided the complete mathematical properties and expansion of the Voigt function, establishing it as the standard for spectral line shapes.

#### 2. Numerical Algorithm
*   **Paper**: "Algorithm 680: Evaluation of the Complex Error Function"
*   **Author**: G. P. M. Poppe & C. M. J. Wijers
*   **Journal**: *ACM Transactions on Mathematical Software*, Vol 16, Issue 1, pp 47.
*   **Year**: 1990
*   **DOI**: [10.1145/77626.77630](https://doi.org/10.1145/77626.77630)
*   **Core Contribution**: Proposed the fast and high-precision algorithm (Algorithm 680) for evaluating $w(z)$, which forms the basis of Python's `scipy.special.wofz` implementation.

---

## 8.2 X-ray Penetration Depth

### Mathematical Form
Relation between Cumulative Intensity Fraction $G_x$ and depth $x$:
$$ G_x = 1 - e^{-(2\mu x / \sin\theta)} $$

### Original Sources

#### 1. Derivation
*   **Book**: *Elements of X-Ray Diffraction* (2nd Edition)
*   **Author**: B. D. Cullity
*   **Publisher**: Addison-Wesley
*   **Year**: 1978
*   **Page**: p. 292, Equation 9-7.
*   **Proof**: The authoritative source for deriving absorption in asymmetric Bragg geometry. Cited as a standard textbook definition.

#### 2. Mass Attenuation Coefficients
*   **Paper**: "Photon Mass Attenuation and Energy-absorption Coefficients"
*   **Author**: J. H. Hubbell
*   **Journal**: *International Journal of Applied Radiation and Isotopes*, Vol 33, Issue 11, pp 1269-1290.
*   **Year**: 1982
*   **DOI**: [10.1016/0020-708X(82)90248-4](https://doi.org/10.1016/0020-708X(82)90248-4)
*   **Core Contribution**: Provided the standard standard values for $\mu/\rho$ at 8.04 keV (Basis of NIST database).

---

## 8.3 Levenberg-Marquardt Optimization Form

### Mathematical Form
Parameter update iteration formula introducing the damping factor $\lambda$:
$$ \mathbf{P}_{new} = \mathbf{P}_{old} - (\mathbf{J}^T \mathbf{J} + \lambda \mathbf{I})^{-1} \mathbf{J}^T \mathbf{r} $$

### Original Source

*   **Paper**: "An Algorithm for Least-Squares Estimation of Nonlinear Parameters"
*   **Author**: Donald W. Marquardt
*   **Journal**: *Journal of the Society for Industrial and Applied Mathematics*, Vol 11, No 2, pp 431-441.
*   **Year**: 1963
*   **DOI**: [10.1137/0111030](https://doi.org/10.1137/0111030)
*   **Core Contribution**: Proposed the method to dynamically switch between Gauss-Newton (Fast convergence) and Gradient Descent (Guaranteed convergence), solving instability in nonlinear fitting. The `least_squares` core in this project is based on this paper.

---

## 8.4 Doublet Separation

### Mathematical Form
Total intensity is a superposition of $K\alpha_1$ and $K\alpha_2$, following the Rachinger correction:
$$ I_{total} = I_{\alpha1}(2\theta) + \frac{1}{2} I_{\alpha1}(2\theta - \Delta 2\theta) $$

### Original Source

*   **Paper**: "A Correction for the $\alpha_1 \alpha_2$ Doublet in the Measurement of Widths of X-ray Diffraction Lines"
*   **Author**: W. A. Rachinger
*   **Journal**: *Journal of Scientific Instruments*, Vol 25, Issue 7, pp 254-255.
*   **Year**: 1948
*   **DOI**: [10.1088/0950-7671/25/7/305](https://doi.org/10.1088/0950-7671/25/7/305)
*   **Core Contribution**: First proposed the method to mathematically separate the pure $K\alpha_1$ signal from a mixed pattern using geometric relations and the 0.5 intensity ratio (Rachinger Correction).

---

## 8.5 Summary Reference Table

All core algorithms in this project correspond directly to the following classic literature:

| Algorithm Module | Key Formula | Original Source (Year) | DOI |
| :--- | :--- | :--- | :--- |
| **Voigt Profile** | $V(x) = \text{Re}[w(z)]$ | Armstrong (1967) | 10.1016/0022-4073(67)90057-X |
| **Faddeeva** | Algorithm 680 | Poppe (1990) | 10.1145/77626.77630 |
| **Optimization** | $(J^TJ + \lambda I)^{-1}$ | Marquardt (1963) | 10.1137/0111030 |
| **Doublet** | Rachinger Correction | Rachinger (1948) | 10.1088/0950-7671/25/7/305 |
| **Absorption** | $G_x = 1 - e^{-2\mu x}$ | Cullity (1978) | (Textbook) |

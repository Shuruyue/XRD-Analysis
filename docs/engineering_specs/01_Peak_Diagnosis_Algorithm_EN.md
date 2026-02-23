# Chapter 1: Peak Diagnosis & Doublet Deconvolution Algorithm

**Document ID**: DOC-ENG-01-EN
**Module Path**: `xrd_analysis.fitting.peak_fitter.fit_peak_with_diagnosis`
**Physics Core**: Quantum Mechanical Spin-Orbit Coupling & True Voigt Convolution

---

## 1.1 Quantum Origin: Spin-Orbit Coupling
In X-ray diffraction patterns, every Bragg Reflection is split into a doublet consisting of $K\alpha_1$ and $K\alpha_2$ lines. This phenomenon originates from the **Spin-Orbit Coupling (L-S Coupling)** of electrons within the copper atom.

### A. Selection Rules
Characteristic X-rays originate from electron transitions from the $L$ shell ($n=2$) back to the $K$ shell ($n=1$).
According to quantum mechanical selection rules:
1.  Principal quantum number change: arbitrary ($\Delta n \neq 0$).
2.  Orbital angular momentum change: $\Delta l = \pm 1$.
3.  Total angular momentum change: $\Delta j = 0, \pm 1$.

Therefore, the $2s$ orbital ($l=0$) in the $L$ shell cannot transition to $1s$ ($l=0$) because $\Delta l=0$ constitutes a **Forbidden Transition**. Only electrons from the $2p$ orbital ($l=1$) can transition.

### B. Fine Structure Splitting
The spin ($s=1/2$) of the $2p$ electron couples with its orbital angular momentum ($l=1$), causing the energy level to split into two states:
1.  **$K\alpha_1$ Source ($2P_{3/2}$)**: $j = l + s = 1 + 1/2 = 3/2$. Higher energy, shorter wavelength.
2.  **$K\alpha_2$ Source ($2P_{1/2}$)**: $j = l - s = 1 - 1/2 = 1/2$. Lower energy, longer wavelength.

### C. Derivation of 0.5 Intensity Ratio
The spectral line intensity is proportional to the **Degeneracy** of the energy level, given by $N = 2j+1$.
*   **$2P_{3/2}$ ($K\alpha_1$)**: $N = 2(3/2) + 1 = 4$.
*   **$2P_{1/2}$ ($K\alpha_2$)**: $N = 2(1/2) + 1 = 2$.

$$ \text{Intensity Ratio} = \frac{I(K\alpha_2)}{I(K\alpha_1)} = \frac{2}{4} = 0.5 $$

This is the physical law origin for `intensity_ratio = 0.5` in the code, not an empirical constant.

---

## 1.2 Mathematical Model: True Voigt Convolution
The fitting function $V(x)$ adopted by this system is the **True Voigt**, defined as the mathematical convolution of a Gaussian function ($G$) and a Lorentzian function ($L$):

$$ V(x; \sigma, \gamma) = G(x; \sigma) \otimes L(x; \gamma) = \int_{-\infty}^{\infty} G(x') L(x-x') dx' $$

This project uses the **Faddeeva function $w(z)$** for calculation, which provides the analytical solution for the Voigt function in the complex plane:
$$ V(x) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}}, \quad z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$
*   $x$: $2\theta - 2\theta_0$ (Deviation angle)
*   $\sigma$: Gaussian standard deviation (Corresponds to Microstrain)
*   $\gamma$: Lorentzian HWHM (Corresponds to Crystallite Size inverse)

The total fitting function is a doublet superposition:
$$ I_{total}(2\theta) = A \cdot V(2\theta - \theta_1) + 0.5A \cdot V(2\theta - \theta_2) + Bkg $$

## 1.2.1 Step-by-Step Substitution Trace

You asked **"How to substitute values into this?"**. Here we demonstrate how to plug real numbers into the complex formula to calculate a single point at the peak apex.

**Hypothetical Case**:
*   Target: Calculate intensity at **Peak Center** ($x=0$).
*   Fitting Parameters: $\sigma = 0.05$ (Gaussian Width), $\gamma = 0.05$ (Lorentzian Width).

**Step 1: Calculate Complex Variable z**
$$ z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$
Substitution:
$$ z = \frac{0 + 0.05i}{0.05 \times 1.414} = \frac{0.05i}{0.0707} = \mathbf{0.707i} $$
(This is a pure imaginary number)

**Step 2: Calculate Faddeeva Function w(z)**
This requires a lookup table or software (e.g., Python `scipy.special.wofz`):
$$ w(0.707i) = \text{erfcx}(0.707) \approx \mathbf{0.5231} $$
(For pure imaginary input, the output of $w(z)$ is real)

**Step 3: Calculate Voigt Function Value V(x)**
$$ V(0) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}} $$
Substitution:
$$ V(0) = \frac{0.5231}{0.05 \times 2.5066} = \frac{0.5231}{0.1253} = \mathbf{4.174} $$

**Step 4: Substitute into Total Intensity Formula**
Assume Amplitude $A=1000$ (Counts), Main Peak Center $\theta_1=43^\circ$.
We want to calculate intensity at $\theta=43^\circ$:
*   **$K\alpha_1$ Term**: $1000 \times V(0) = 1000 \times 4.174 = 4174$
*   **$K\alpha_2$ Term**: Located $0.1^\circ$ away. Assume calculated $V(-0.1) \approx 0.5$. Intensity $= 500 \times 0.5 = 250$.
*   **Total Intensity**: $4174 + 250 + \text{Bkg}$.

This details how the code converts **$z$ (Complex Coordinate)** into **$I$ (Real Intensity)**.

---

## 1.3 Numerical Trace

To verify the algorithm's correctness, we perform an "End-to-End" numerical trace using the data from the provided image.

**Input Data from Image**:
*   Plane (HKL): (111)
*   Fitting Center ($2\theta_{fit}$): **43.353°** ($K\alpha_1$ position)
*   Wavelength Constants: $\lambda_1 = 1.540562 \mathring{A}, \lambda_2 = 1.544390 \mathring{A}$

### Step 1: Calculate Theoretical $K\alpha_2$ Position
Using Bragg's Law $n\lambda = 2d\sin\theta$, for the same interplanar spacing $d$, the diffraction angle is proportional to the wavelength:
$$ \sin(\frac{2\theta_2}{2}) = \frac{\lambda_2}{\lambda_1} \sin(\frac{2\theta_1}{2}) $$

1.  **Calculate $\theta_1$ for $K\alpha_1$**:
    $$ \theta_1 = \frac{43.353}{2} = 21.6765^\circ $$
    $$ \sin(21.6765^\circ) = 0.369369 $$

2.  **Calculate $\sin\theta_2$ for $K\alpha_2$**:
    $$ \text{Ratio} = \frac{1.544390}{1.540562} = 1.0024848 $$
    $$ \sin\theta_2 = 0.369369 \times 1.0024848 = 0.370287 $$

3.  **Derive $2\theta_2$**:
    $$ \theta_2 = \arcsin(0.370287) = 21.7331^\circ $$
    $$ 2\theta_2 = 21.7331 \times 2 = \mathbf{43.466^\circ} $$

**Verification Result**: Compare with the red text box `Ka-43.466` in your image. This is an **Exact Match**.

### Step 2: Intensity Summation
At $2\theta = 43.353^\circ$ (Primary Peak Apex), the total intensity $I_{total}$ consists of two parts:

1.  **$K\alpha_1$ Contribution**:
    At the center ($x=0$). Voigt function $V(0)$ is at its maximum. Let intensity be $I_{max}$.
    Contribution $\approx I_{max}$.

2.  **$K\alpha_2$ Contribution**:
    Located $0.113^\circ$ away.
    Since (111) FWHM is $0.248^\circ$, HWHM $\approx 0.124^\circ$.
    The $K\alpha_2$ center is approximately 1 HWHM distance from the observation point.
    According to Voigt function properties, intensity decays to ~50% at 1 HWHM.
    Combined with the quantum mechanical intensity scaling of 0.5.
    Contribution $\approx I_{max} \times 0.5 \times 0.5 = 0.25 I_{max}$.

**Engineering Conclusion**:
At the main peak apex, approximately **20-25%** of the signal actually comes from the tail extension of the "invisible $K\alpha_2$ peak" to the right.
If the code does not perform this **Doublet Deconvolution** and measures peak height directly, it results in a 25% error, and the peak width is incorrectly broadened, leading to underestimated grain size.

---

## 1.4 Engineering Diagnosis

### A. Coefficient of Determination ($R^2$)
$$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$
*   **Value**: Image shows $R^2 = 0.9997$.
*   **Interpretation**: The model explains 99.97% of signal variance. Physically, this implies the **internal strain field is extremely uniform**. Residuals would show wave-like patterns if lattice gradients or secondary phases existed.

### B. Mixing Parameter ($\eta$)
$$ \eta = \frac{\Gamma_L}{\Gamma_L + \Gamma_G} $$
*   **Value**: Image shows $\eta = 0.916$ (Close to 1.0).
*   **Interpretation**: The peak shape is highly **Lorentzian**.
    *   **Lorentzian Dominant**: Broadening mechanism is mainly **Size Broadening** (Finite crystallite size).
    *   **Gaussian Dominant**: Broadening mechanism is mainly **Microstrain**.
    *   **Conclusion**: This electroplated copper layer has very fine (nanocrystalline) (111) grains, but dislocation distribution is relatively simple, causing minimal Gaussian lattice distortion.

### C. FWHM
*   **Value**: (111) $0.248^\circ$ vs (200) $0.534^\circ$.
*   **Interpretation**: The (200) peak width is more than double that of (111). This confirms the plating has **(111) Preferred Orientation**, allowing (111) grains to coalesce and grow along the growth direction (Large Size $\to$ Narrow Peak), while (200) grain growth is suppressed (Small Size $\to$ Broad Peak).

---

## 1.5 Code Module Mapping
*   **TrueVoigt**: `xrd_analysis.fitting.pseudo_voigt.py` (Performs `wofz` complex error function calculation)

---

## 1.6 References

1.  **Cullity, B. D. (1978)**. *Elements of X-Ray Diffraction*, 2nd Ed., Addison-Wesley.
    *   Chapter 2: Physics of X-rays (Spin-Orbit Coupling & Selection Rules).
2.  **Armstrong, B. H. (1967)**. "Spectrum Line Profiles: The Voigt Function". *Journal of Quantitative Spectroscopy and Radiative Transfer*, 7(1), 61-88.
    *   Mathematical definition and properties of the Voigt function.
3.  **Poppe, G. P. M., & Wijers, C. M. J. (1990)**. "Algorithm 680: Evaluation of the Complex Error Function". *ACM Transactions on Mathematical Software*, 16(1), 47.
    *   Source algorithm for `scipy.special.wofz` (Faddeeva function).

# Chapter 4: Stacking Fault Analysis Algorithm

**Document ID**: DOC-ENG-04-EN
**Module Path**: `xrd_analysis.methods.defect_analysis.StackingFaultAnalyzer`
**Physics Core**: Warren-Averbach Planar Defect Theory

---

## 4.1 Crystal Physics Origin: Symmetry Breaking
Face-Centered Cubic (FCC) metals like copper have an atomic stacking sequence of **ABCABC...**.
When errors occur during electroplating growth or plastic deformation, the sequence becomes **ABC|BC|ABC...** (Intrinsic) or **ABC|B|ABC...** (Extrinsic). These 2D planar defects are called **Stacking Faults**.

### Why Peaks Shift?
In Reciprocal Space, stacking faults break the translational symmetry of the lattice, causing Reciprocal Lattice Points to elongate and shift along specific directions.
According to derivations by **Paterson (1952)** and **Warren (1969)**, this causes regular shifts in X-ray diffraction peaks:
1.  **(111) Plane**: Shifts to **Higher Angle** ($+$).
2.  **(200) Plane**: Shifts to **Lower Angle** ($-$).

These two peaks act like they are "pushed together" by an invisible hand, causing their **Peak Separation** to shrink.

---

## 4.2 Mathematical Model: Warren Geometric Coefficient G
Instead of relying on absolute peak shifts (which are highly susceptible to sample height errors), we measure the **contraction of the (200)-(111) Peak Separation**.

**Definition**:
$$ \alpha = \frac{ \Delta SEP_{obs} }{ G } = \frac{ (2\theta_{200}-2\theta_{111})_{obs} - (2\theta_{200}-2\theta_{111})_{std} }{ G } $$
Where $\alpha$ is the Stacking Fault Probability (ratio of faulted layers to total layers).

### Derivation of G
Warren proved the relationship between peak shift $\Delta 2\theta$ and $\alpha$:
$$ \Delta 2\theta_{hkl} = \pm \frac{90\sqrt{3} \tan\theta}{\pi^2 h_0^2 (u+b)} \sum (\pm L_0) \alpha $$

Substituting geometric parameters for Copper (111) and (200):
1.  **(111) Shift**: $\Delta 2\theta_{111} = + \frac{45\sqrt{3}}{4\pi^2} \alpha \tan\theta_{111}$ (Corrected for polycrystalline averaging)
2.  **(200) Shift**: $\Delta 2\theta_{200} = - \frac{90\sqrt{3}}{2\pi^2} \alpha \tan\theta_{200}$

**Peak Separation Change $\Delta SEP$**:
$$ \Delta SEP = \Delta 2\theta_{200} - \Delta 2\theta_{111} = - \frac{45\sqrt{3}}{\pi^2} (\tan\theta_{111} + 2\tan\theta_{200}) \alpha $$

Plugging in standard angles for Copper ($\theta_{111} \approx 21.66^\circ, \theta_{200} \approx 25.22^\circ$):
$$ G \approx -7.897 \text{ deg/unit } \alpha $$

This is the physical origin of the code constant `WARREN_G_COEFFICIENT = -7.897`.

---

## 4.3 Numerical Trace

To verify algorithm sensitivity, we simulate a realistic scenario with **"0.5% Fault Probability"**. This is a "mild but measurable" defect concentration in electroplated copper.

**Standard State ($\alpha=0$)**:
*   $\lambda = 1.540562 \mathring{A}$
*   $d_{111} = 2.0871 \mathring{A} \rightarrow 2\theta_{111} = 43.316^\circ$
*   $d_{200} = 1.8075 \mathring{A} \rightarrow 2\theta_{200} = 50.448^\circ$
*   **Standard Separation $SEP_{std} = 50.448 - 43.316 = 7.132^\circ$**

**Defected State ($\alpha=0.005$)**:
Assume the sample has 0.5% stacking faults (1 fault every 200 layers).

1.  **Calculate Peak Shift**:
    $$ \text{Shift} = G \times \alpha = -7.897 \times 0.005 = -0.0395^\circ $$
    Ideally, peak separation shrinks by about 0.04 degrees.

2.  **Observed Separation**:
    $$ SEP_{obs} = SEP_{std} + \text{Shift} = 7.132 - 0.0395 = 7.0925^\circ $$

3.  **Algorithm Reverse Calculation**:
    When our code reads $SEP_{obs} = 7.0925^\circ$:
    $$ \text{Deviation} = 7.0925 - 7.132 = -0.0395^\circ $$
    $$ \alpha = \frac{-0.0395}{-7.897} = 0.005001 $$
    $$ \alpha_{\%} = 0.50 \% $$

**Engineering Conclusion**:
Even if the peak separation shrinks by a mere **0.04 degrees** (requiring high-resolution instruments to distinguish), our algorithm can precisely back-calculate a **0.5%** fault probability. This proves the high sensitivity of this algorithm to microstructural changes.

---

## 4.4 Engineering Interpretation

What does the output $\alpha$ value mean physically?

| Fault Probability $\alpha$ | State Description | Impact on Properties |
| :--- | :--- | :--- |
| **< 0.1%** | Near Perfect Crystal | Excellent conductivity, but mechanical strength might be low (too soft). |
| **0.2% - 0.5%** | Mild Defects | **Optimal Balance**. Moderate faults impede dislocation slip (hardening) without significantly degrading conductivity. |
| **> 1.0%** | High Defect Density | Extremely hard (Work Hardened), but resistivity rises significantly and **Elongation** drops drastically (Brittleness). |

### Relation to Residual Stress (Stress Decoupling)
Stacking faults cause (111) to shift to a higher angle, which is the same direction as a shift caused by "Compressive Stress".
**If fault correction is ignored**, and stress is calculated directly from the (111) position, it yields a **False Giant Compressive Stress**.
This module prioritizes calculating $\alpha$, warning users that the (111) position contains non-stress factors. This is why in `01_Peak_Diagnosis`, we recommend using the (220) peak (less affected by faults) for stress monitoring.

---

## 4.5 References
1.  **Warren, B. E. (1969)**. *X-ray Diffraction*. Dover Publications. (Chapter 13).
2.  **Paterson, M. S. (1952)**. "X-ray Diffraction by Face-Centered Cubic Crystals with Deformation Faults." *Journal of Applied Physics*, 23(8), 805-811.

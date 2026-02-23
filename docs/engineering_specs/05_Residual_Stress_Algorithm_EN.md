# Chapter 5: Residual Stress Analysis Algorithm

**Document ID**: DOC-ENG-05-EN
**Module Path**: `xrd_analysis.core.copper_crystal.get_youngs_modulus`
**Physics Core**: Anisotropic Elasticity Theory (Hooke's Law in Textured Films)

---

## 5.1 Physics Origin: Elastic Anisotropy
General mechanics textbooks often assume materials are "Isotropic", meaning Young's Modulus $E$ and Poisson's Ratio $\nu$ are the same in all directions.
**This is completely wrong for electroplated copper.** Copper is one of the most elastically anisotropic metals:

*   **<111> Direction (Hardest)**: Densest atomic packing, $E_{111} \approx 191$ GPa.
*   **<100> Direction (Softest)**: Least dense packing, $E_{100} \approx 67$ GPa.

The difference is nearly **3-fold**. If you calculate stress using a generic $E_{bulk} \approx 120$ GPa, your stress data could have **huge errors of 50%~100%**.

### Algorithm Innovation
We implemented direction-dependent stress calculation based on single-crystal data from **Ledbetter & Naimon (1974)**:
$$ E_{hkl} = \left( S_{11} - 2(S_{11}-S_{12}-\frac{1}{2}S_{44})\Gamma \right)^{-1} $$
$$ \nu_{hkl} = - \frac{S_{12}}{S_{11} - 2(S_{11}-S_{12}-\frac{1}{2}S_{44})\Gamma} $$
Where $\Gamma$ is the direction factor.

---

## 5.2 Mathematical Model: Equi-biaxial Stress
Electroplated copper films are typically in a **Plane Stress State**, meaning normal stress $\sigma_z = 0$.
Assuming in-plane stress is uniform (Equi-biaxial): $\sigma_x = \sigma_y = \sigma$.

**Hooke's Law Derivation**:
The relationship between normal strain $\epsilon_z$ (which we measure) and in-plane stress $\sigma$ is:
$$ \epsilon_z = \frac{1}{E} [\sigma_z - \nu(\sigma_x + \sigma_y)] $$
Substitute $\sigma_z=0, \sigma_x=\sigma_y=\sigma$:
$$ \epsilon_z = -\frac{2\nu}{E} \sigma $$
Rearranging gives the stress formula:
$$ \sigma = -\frac{E_{hkl}}{2\nu_{hkl}} \cdot \epsilon_z = -\frac{E_{hkl}}{2\nu_{hkl}} \cdot \frac{d_{obs} - d_{0}}{d_{0}} $$

*   **Sign Convention**: If $d$ expands (Lattice stretch $\epsilon_z > 0$), then $\sigma$ is negative (Compressive Stress).

---

## 5.3 Numerical Trace

To demonstrate the power of **Anisotropy Correction**, we calculate using the **(200) plane**, which is most sensitive to stress.

**Step 0: Prepare Constants**
The code automatically loads constants specific to the (200) direction (Non-polycrystalline average):
*   $E_{200} = 66.7$ GPa (Softest direction, highest sensitivity)
*   $\nu_{200} = 0.419$ (High Poisson's ratio)
*   $d_{0, 200} = 1.8080 \mathring{A}$ (Corresponding to $2\theta = 50.448^\circ$)

**Step 1: Input Observed Data**
Assume we observe the (200) peak shifted slightly left to **$50.400^\circ$**.
1.  **Calculate Observed d-spacing $d_{obs}$**:
    $$ \theta_{obs} = 25.200^\circ $$
    $$ d_{obs} = \frac{1.540562}{2 \sin(25.200^\circ)} = 1.8091 \mathring{A} $$

**Step 2: Calculate Strain**
$$ \epsilon_z = \frac{1.8091 - 1.8080}{1.8080} = +6.13 \times 10^{-4} $$
The lattice is stretched vertically by about 0.06%.

**Step 3: Calculate Stress**
$$ \sigma = - \frac{66.7 \times 10^3 \text{ MPa}}{2 \times 0.419} \times (6.13 \times 10^{-4}) $$
$$ \sigma = - 79594 \times 0.000613 = \mathbf{-48.8 \text{ MPa}} $$

**(Comparison) Consequence of using wrong Polycrystalline Constants ($E=130, \nu=0.34$)**
$$ \sigma_{wrong} = - \frac{130000}{0.68} \times 0.000613 = \mathbf{-117.2 \text{ MPa}} $$
**Error up to 140%!** This explains why stress calculated by traditional generic XRD often mismatches Wafer Bow measurements.

---

## 5.4 Engineering Interpretation

| Stress Type | Physics Phenomenon | Engineering Risk |
| :--- | :--- | :--- |
| **Compressive (< 0)** | Lattice squeezed laterally (Vertical expansion) | **Stress Migration (SM)**. Atoms are squeezed and migrate, forming **Hillocks** after long storage, causing shorts. |
| **Tensile (> 0)** | Lattice pulled laterally (Vertical contraction) | **Voiding**. Grain boundaries pulled apart, leading to resistance increase or open circuits. |

### Why Recommend (200) for Stress?
Although (111) has the strongest signal, it is insensitive to stress ($E=191$ GPa is too hard; shift is only 1/3 of (200) for the same stress). This makes (111) shifts easily overwhelmed by instrument error.
**(200) is the softest direction**, acting like a spring with the largest response to stress, making it the best indicator for thin film stress monitoring.

---

## 5.5 References
1.  **Ledbetter, H. M., & Naimon, E. R. (1974)**. "Elastic Properties of Metals and Alloys. II. Copper". *J. Phys. Chem. Ref. Data*, 3(4).
    *   Provides $C_{11}, C_{12}, C_{44}$ for calculating $E_{hkl}$.
2.  **Nye, J. F. (1957)**. *Physical Properties of Crystals*. Oxford University Press.
    *   The bible of anisotropic elasticity.

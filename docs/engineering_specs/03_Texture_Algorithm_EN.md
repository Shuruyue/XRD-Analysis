# Chapter 3: Texture Analysis Algorithm

**Document ID**: DOC-ENG-03-EN
**Module Path**: `xrd_analysis.methods.texture.TextureAnalyzer`
**Physics Core**: Harris Texture Coefficient (Inverse Pole Figure Method)

---

## 3.1 Physics Origin: Preferred Orientation
In an ideal polycrystalline powder, crystallite orientations are randomly distributed. However, in electroplated copper processes, grains tend to grow along specific crystallographic directions (e.g., <111> or <200>) to minimize surface energy.
This non-random alignment is called **Preferred Orientation** or **Texture**.

### Why is this important?
*   **<111> Texture**: Closest atomic packing, highest resistance to Electromigration. The gold standard for semiconductor interconnects.
*   **<200> Texture**: Lowest surface energy, common in recrystallized grains after annealing.
*   **<220> Texture**: Often seen in extremely fine grain structures formed under excessive Suppressor action.

---

## 3.2 Mathematical Model: Harris Texture Coefficient
We use the **Texture Coefficient (TC)** proposed by Harris (1952) to quantify this effect. This is a simplified "Inverse Pole Figure" algorithm suitable for $\theta-2\theta$ scans.

**Definition**:
$$ TC(hkl) = \frac{ \frac{I_{obs}(hkl)}{I_{std}(hkl)} }{ \frac{1}{N} \sum_{i=1}^{N} \frac{I_{obs}(h_i k_i l_i)}{I_{std}(h_i k_i l_i)} } $$

*   $I_{obs}$: Experimentally measured Integrated Intensity (Note: Must be area, not peak height, to account for broadening).
*   $I_{std}$: JCPDS Standard Powder Diffraction intensity (Random orientation).
*   $N$: Total number of diffraction peaks included in the calculation.

**Interpretation Criteria**:
*   $TC = 1$: Random orientation (No texture).
*   $TC > 1$: Preferred orientation.
*   $TC < 1$: Suppressed orientation.

---

## 3.3 Numerical Trace

To thoroughly understand how the Harris algorithm works, we manually calculate a set of simulated data with **High <111> Texture**.

**Standard Intensity (JCPDS 04-0836)**:
Read from `xrd_analysis.core.copper_crystal` module:
*   $I_{std}(111) = 100$
*   $I_{std}(200) = 46$
*   $I_{std}(220) = 20$

**Observed Data**:
Assume we measured the following background-subtracted intensities:
*   $I_{obs}(111) = 15000$
*   $I_{obs}(200) = 2000$
*   $I_{obs}(220) = 500$

### Step 1: Calculate Normalized Ratios
$$ R_{111} = \frac{15000}{100} = 150.0 $$
$$ R_{200} = \frac{2000}{46} = 43.48 $$
$$ R_{220} = \frac{500}{20} = 25.0 $$

### Step 2: Calculate Average Ratio
$$ R_{avg} = \frac{150.0 + 43.48 + 25.0}{3} = \frac{218.48}{3} = 72.827 $$

### Step 3: Calculate Texture Coefficients (TC)
$$ TC_{111} = \frac{R_{111}}{R_{avg}} = \frac{150.0}{72.827} = \mathbf{2.06} $$
$$ TC_{200} = \frac{R_{200}}{R_{avg}} = \frac{43.48}{72.827} = \mathbf{0.60} $$
$$ TC_{220} = \frac{R_{220}}{R_{avg}} = \frac{25.0}{72.827} = \mathbf{0.34} $$

**Verify Sum**:
In the Harris method, the mean value of TC must be normalized to 1 by definition.
$$ \text{Mean TC} = \frac{2.06 + 0.60 + 0.34}{3} = \frac{3.0}{3} = 1.0 \quad \text{(Correct)} $$

**Engineering Conclusion**:
The $TC(111) = 2.06$ for this sample is significantly greater than 1.
Physically, this means that among the grains participating in diffraction at the surface, the number of grains with (111) planes parallel to the surface is **2 times** that of a random distribution. This is a very strong (111) Preferred Orientation signal.

---

## 3.4 Engineering Diagnosis & Application

You often see trends corresponding to FWHM in `compare_integral_methods` and `fwhm_evolution` charts.

| Parameter | Physical Interpretation | Engineering Diagnosis |
| :--- | :--- | :--- |
| **TC(111) > 1.5** | Strong (111) Texture | **Preferred**. High electromigration resistance, suitable for fine line processes. Usually accompanied by narrower (111) FWHM. |
| **TC(200) > 1.2** | (200) Recrystallization | **Warning**. Requires attention. Often occurs after long annealing or with incorrect plating parameters (e.g., excessively high current density). Usually accompanied by a sharp drop in (200) FWHM. |
| **Degree of Texture (σ)** | $\sqrt{\frac{\sum (TC-1)^2}{N}}$ | Higher value indicates stronger texture. $\sigma < 0.2$ is considered random structure. |

### Why is this related to FWHM charts?
The `fwhm_evolution` is essentially a mirror of texture.
Under strong (111) texture, (111) grains preferentially coalesce to reduce surface energy. This causes the **(111) FWHM to start decreasing earlier than other planes**. We can confirm this by comparing the TC curve with the FWHM curve.

---

## 3.5 References
1.  **Harris, G. B. (1952)**. "Quantitative Measurement of Preferred Orientation in Rolled Uranium Bars." *Philosophical Magazine*, 43(336), 113-123.
2.  **JCPDS Card 04-0836**. Joint Committee on Powder Diffraction Standards.

# 第八章：算法核心原理與原始文獻 (Algorithm Principles & Original References)

**文件編號**: DOC-ENG-08
**摘要**: 應工程規範要求，本文件嚴格溯源所有數學模型的「原始出處 (Original Source)」，並提供 DOI 以供學術查證。我們不引用間接文獻或簡報，僅引用提出該數學形式的原始論文或聖經級教科書。

---

## 8.1 Voigt 函數與 Faddeeva 計算 (Voigt Profile & Faddeeva)

### 數學形式 (Mathematical Form)
我們使用的 $V(x; \sigma, \gamma)$ 是透過複數誤差函數 (Complex Error Function) $w(z)$ 的實部來計算：

$$ V(x) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}}, \quad z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$

### 原始出處 (Original Sources)

#### 1. 物理定義 (Physics Definition)
*   **論文**: "Spectrum Line Profiles: The Voigt Function"
*   **作者**: B. H. Armstrong
*   **期刊**: *Journal of Quantitative Spectroscopy and Radiative Transfer*, Vol 7, Issue 1, pp 61-88.
*   **年份**: 1967
*   **DOI**: [10.1016/0022-4073(67)90057-X](https://doi.org/10.1016/0022-4073(67)90057-X)
*   **核心貢獻**: 給出了 Voigt 函數的完整數學性質與展開式，確立了其作為光譜線形標準的地位。

#### 2. 數值演算法 (Numerical Algorithm)
*   **論文**: "Algorithm 680: Evaluation of the Complex Error Function"
*   **作者**: G. P. M. Poppe & C. M. J. Wijers
*   **期刊**: *ACM Transactions on Mathematical Software*, Vol 16, Issue 1, pp 47.
*   **年份**: 1990
*   **DOI**: [10.1145/77626.77630](https://doi.org/10.1145/77626.77630)
*   **核心貢獻**: 提出了計算 $w(z)$ 的快速且高精度演算法 (Algorithm 680)，這是 Python `scipy.special.wofz` 函數的底層實作依據。

---

## 8.2 X-ray 穿透深度 (X-ray Penetration Depth)

### 數學形式 (Mathematical Form)
累積強度分率 $G_x$ 與深度 $x$ 的關係：
$$ G_x = 1 - e^{-(2\mu x / \sin\theta)} $$

### 原始出處 (Original Sources)

#### 1. 公式推導 (Derivation)
*   **書籍**: *Elements of X-Ray Diffraction* (2nd Edition)
*   **作者**: B. D. Cullity
*   **出版社**: Addison-Wesley
*   **年份**: 1978
*   **頁數**: p. 292, Equation 9-7.
*   **證明**: 這是推導不對稱 Bragg 幾何吸收的權威來源。教科書級定義，無須引用單篇論文。

#### 2. 質量吸收係數 (Mass Attenuation Coefficients)
*   **論文**: "Photon Mass Attenuation and Energy-absorption Coefficients"
*   **作者**: J. H. Hubbell
*   **期刊**: *International Journal of Applied Radiation and Isotopes*, Vol 33, Issue 11, pp 1269-1290.
*   **年份**: 1982
*   **DOI**: [10.1016/0020-708X(82)90248-4](https://doi.org/10.1016/0020-708X(82)90248-4)
*   **核心貢獻**: 提供了 8.04 keV 下 $\mu/\rho$ 的標準值 (NIST 數據庫基礎)。

---

## 8.3 Levenberg-Marquardt 優化算法 (L-M Algorithm)

### 數學形式 (Mathematical Form)
參數更新迭代公式，引入了阻尼因子 $\lambda$：
$$ \mathbf{P}_{new} = \mathbf{P}_{old} - (\mathbf{J}^T \mathbf{J} + \lambda \mathbf{I})^{-1} \mathbf{J}^T \mathbf{r} $$

### 原始出處 (Original Source)

*   **論文**: "An Algorithm for Least-Squares Estimation of Nonlinear Parameters"
*   **作者**: Donald W. Marquardt
*   **期刊**: *Journal of the Society for Industrial and Applied Mathematics*, Vol 11, No 2, pp 431-441.
*   **年份**: 1963
*   **DOI**: [10.1137/0111030](https://doi.org/10.1137/0111030)
*   **核心貢獻**: 提出了在 Gauss-Newton (快速收斂) 與 Gradient Descent (保證收斂) 之間動態切換的方法，解決了非線性擬合的不穩定性。本專案的 `least_squares` 核心即基於此論文。

---

## 8.4 雙峰分離與疊加 (Doublet Separation)

### 數學形式 (Mathematical Form)
總強度由 $K\alpha_1$ 與 $K\alpha_2$ 疊加而成，且遵循 Rachinger 修正：
$$ I_{total} = I_{\alpha1}(2\theta) + \frac{1}{2} I_{\alpha1}(2\theta - \Delta 2\theta) $$

### 原始出處 (Original Source)

*   **論文**: "A Correction for the $\alpha_1 \alpha_2$ Doublet in the Measurement of Widths of X-ray Diffraction Lines"
*   **作者**: W. A. Rachinger
*   **期刊**: *Journal of Scientific Instruments*, Vol 25, Issue 7, pp 254-255.
*   **年份**: 1948
*   **DOI**: [10.1088/0950-7671/25/7/305](https://doi.org/10.1088/0950-7671/25/7/305)
*   **核心貢獻**: 最早提出利用幾何關係與 0.5 強度比，從混合圖譜中數學分離出純 $K\alpha_1$ 訊號的方法 (即 Rachinger Correction)。

---

## 8.5 總結 (Summary)

本專案的所有算法核心，均直接對應至以下經典文獻：

| 算法模組 | 關鍵公式 | 原始文獻 (Year) | DOI |
| :--- | :--- | :--- | :--- |
| **Voigt Profile** | $V(x) = \text{Re}[w(z)]$ | Armstrong (1967) | 10.1016/0022-4073(67)90057-X |
| **Faddeeva** | Algorithm 680 | Poppe (1990) | 10.1145/77626.77630 |
| **Optimization** | $(J^TJ + \lambda I)^{-1}$ | Marquardt (1963) | 10.1137/0111030 |
| **Doublet** | Rachinger Correction | Rachinger (1948) | 10.1088/0950-7671/25/7/305 |
| **Absorption** | $G_x = 1 - e^{-2\mu x}$ | Cullity (1978) | (Textbook) |

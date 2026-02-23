# 第二章：晶粒尺寸分析 (Crystallite Size Analysis)

**文件編號**: DOC-ENG-02
**模組位置**: `xrd_analysis.methods.scherrer.ScherrerCalculator`
**物理核心**: 基於 FWHM 的 Scherrer 方程與 Voigt 展寬校正

---

## 2.1 理論金標準：傅立葉分析 (Warren-Averbach)
定義「晶粒尺寸」最嚴謹的方法是透過繞射峰形的傅立葉變換係數。

**定義式**:
傅立葉係數 $A_n$ 在 $n \to 0$ 處的導數，給出了 **面積加權平均柱長 (Area-Weighted Column Length)** $\langle L \rangle_{area}$：
$$ -\frac{dA_n^{size}}{dn} \bigg|_{n=0} = \frac{1}{\langle L \rangle_{area}} $$

*   **優點**: 能完美區分尺寸與應變，不依賴峰形假設。
*   **缺點**: 需要多階數衍射峰 (如 (111) **與** (222))；對峰尾 (Tails) 雜訊極度敏感。

---

## 2.2 程式實作：FWHM Scherrer 法 (Implementation)
**本專案演算法**: `xrd_analysis.methods.scherrer.ScherrerCalculator`

目前程式主路徑使用 **校正後 FWHM** 的 Scherrer 方程式。
當峰擬合可提供 `eta_observed` 時，使用 Voigt 分量扣除；
若無 `eta_observed`，則自動回退到教科書的二次扣除。

**實作公式**:
$$ D_{v} = \frac{K \lambda}{\beta \cos\theta} $$
*   $D_v$: 體積加權平均晶粒尺寸。
*   $\beta$: 校正後峰寬 (FWHM)，單位為弧度。

### 程式中的寬化校正策略
1.  **Voigt 分量扣除（優先）**：將觀測峰與儀器峰分解為高斯/洛倫茲分量，分別扣除後再重組。
2.  **二次扣除（回退）**：$\beta = \sqrt{\beta_{obs}^2-\beta_{inst}^2}$，用於 `eta` 不可用時。

**算法差異比對 (Method Comparison from Image)**:
*   **Method 1 (Single PV)**: 忽略 $K\alpha_2$ 分裂。FWHM 被 $K\alpha_2$ 撐大 $\to$ 晶粒尺寸被嚴重低估 (假象)。
*   **Method 2 (Stripping)**: 數學減去 $K\alpha_2$。雜訊大，峰形易扭曲。
*   **Method 3 (Doublet Fit - 本專案採用)**: 使用 True Voigt 雙峰模型，並提供 `eta_observed` 供 Scherrer 校正使用。

---

## 2.3 幾何形狀因子 (Shape Factor K)

Scherrer 常數 $K$ 取決於 **晶粒形狀** 與 **衍射方向 (hkl)**。
本專案拒絕使用「球形假設 ($K=0.9$)」，因為電鍍銅是柱狀生長。

**本程式實作：立方晶習模型 (Cubic Habit Model)**
依據 **Langford & Wilson (1978)**，假設晶粒為立方體：

*   **(111) 衍射**: 視線沿著體對角線。投影為六邊形。
    *   **$K = 0.855$**
*   **(200) 衍射**: 視線沿著邊長。投影為正方形。
    *   **$K = 0.886$**

### 數值演算示範 (Numerical Trace)

為證明 $K$ 值選擇的影響，我們代入一個實際情境。
**輸入**:
*   $\beta = 0.3^\circ$ (校正後 FWHM) $\approx 0.005236$ rad
*   $\theta = 21.65^\circ$ (for 111)
*   $\lambda = 1.5406 \mathring{A}$

**情境 A：使用通用 K=0.9 (球形假設)**
$$ D = \frac{0.9 \times 1.5406}{0.005236 \times \cos(21.65^\circ)} = \frac{1.3865}{0.00486} = \mathbf{285.2 \mathring{A}} $$

**情境 B：使用正確 K=0.855 (立方體 (111) 修正)**
$$ D = \frac{0.855 \times 1.5406}{0.005236 \times \cos(21.65^\circ)} = \frac{1.3172}{0.00486} = \mathbf{271.0 \mathring{A}} $$

**工程結論**:
使用通用 K 值會造成 **5.2% 的系統性誤差** (285 vs 271)。
對於精密製程監控，這 5% 的誤差可能導致誤判添加劑 (Leveler) 的晶粒細化效能。因此我們堅持使用 $K=0.855$。

---

## 2.4 物理常數與參數表 (Constants Data Table)

| 參數 | 數值 | 單位 | 描述 | 來源 / DOI |
| :--- | :--- | :--- | :--- | :--- |
| $\lambda$ | **1.540562** | Å | Cu Kα1 波長 | Bearden (1967) [10.1103/RevModPhys.39.78](https://doi.org/10.1103/RevModPhys.39.78) |
| $K_{111}$ | **0.855** | - | 立方體 (111) 形狀因子 | Langford (1978) [10.1107/S0021889878012844](https://doi.org/10.1107/S0021889878012844) |
| $K_{200}$ | **0.886** | - | 立方體 (200) 形狀因子 | Langford (1978) |
| Method | **FWHM (校正後)** | - | 程式目前使用的寬度定義 | Scherrer + Voigt 校正 |

**引用文獻**:
1.  **Langford, J. I., & Wilson, A. J. C. (1978)**. "Scherrer after sixty years". *J. Appl. Cryst.*, 11, 102-113.
2.  **Olivero, J. J., & Longbothum, R. L. (1977)**. "Empirical fits to the Voigt line width". *JQSRT*, 17, 233.

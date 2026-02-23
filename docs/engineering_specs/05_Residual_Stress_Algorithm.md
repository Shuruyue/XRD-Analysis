# 第五章：殘留應力分析算法 (Residual Stress Analysis Algorithm)

**文件編號**: DOC-ENG-05
**模組位置**: `xrd_analysis.core.copper_crystal.get_youngs_modulus`
**物理核心**: Anisotropic Elasticity Theory (Hooke's Law in Textured Films)

---

## 5.1 物理根源：彈性各向異性 (Elastic Anisotropy)
在一般材料力學教科書中，常假設材料是「各向同性 (Isotropic)」的，即楊氏模數 $E$ 和泊松比 $\nu$ 在所有方向都相同。
**這對電鍍銅完全錯誤。** 銅是彈性各向異性最強的金屬之一：

*   **<111> 方向 (最硬)**: 原子排列最密，$E_{111} \approx 191$ GPa。
*   **<100> 方向 (最軟)**: 原子排列最疏，$E_{100} \approx 67$ GPa。

兩者相差近 **3倍**。如果您在計算應力時使用通用的 $E_{bulk} \approx 120$ GPa，您的應力數據可能會有 **50%~100% 的巨大誤差**。

### 本專案算法創新
我們根據 **Ledbetter & Naimon (1974)** 的單晶數據，實作了方向相依的應力計算：
$$ E_{hkl} = \left( S_{11} - 2(S_{11}-S_{12}-\frac{1}{2}S_{44})\Gamma \right)^{-1} $$
$$ \nu_{hkl} = - \frac{S_{12}}{S_{11} - 2(S_{11}-S_{12}-\frac{1}{2}S_{44})\Gamma} $$
其中 $\Gamma$ 為方向因子。

---

## 5.2 數學模型：等雙軸應力 (Mathematical Model)
電鍍銅薄膜通常處於 **平面應力狀態 (Plane Stress)**，即表面法向應力 $\sigma_z = 0$。
假設面內應力是均勻的（等雙軸應力）：$\sigma_x = \sigma_y = \sigma$。

**胡克定律推導 (Hooke's Law Derivation)**:
法向應變 $\epsilon_z$ (我們可以測量的也是這個) 與面內應力 $\sigma$ 的關係為：
$$ \epsilon_z = \frac{1}{E} [\sigma_z - \nu(\sigma_x + \sigma_y)] $$
代入 $\sigma_z=0, \sigma_x=\sigma_y=\sigma$:
$$ \epsilon_z = -\frac{2\nu}{E} \sigma $$
移項得應力計算公式：
$$ \sigma = -\frac{E_{hkl}}{2\nu_{hkl}} \cdot \epsilon_z = -\frac{E_{hkl}}{2\nu_{hkl}} \cdot \frac{d_{obs} - d_{0}}{d_{0}} $$

*   **符號意義**: 若 $d$ 變大 (晶格拉長 $\epsilon_z > 0$)，則 $\sigma$ 為負 (壓縮應力)。

---

## 5.3 數值演算示範 (Numerical Trace)

為展示**各向異性修正**的威力，我們使用對應力最敏感的 **(200) 晶面** 進行演算。

**步驟 0：準備物理常數**
程式碼自動載入 (200) 方向專用常數 (非多晶平均值):
*   $E_{200} = 66.7$ GPa (銅最軟方向，靈敏度最高)
*   $\nu_{200} = 0.419$ (泊松比極大)
*   $d_{0, 200} = 1.8080 \mathring{A}$ (對應 $2\theta = 50.448^\circ$)

**步驟 1：輸入觀測數據**
假設我們觀測到 (200) 峰位稍微向左偏移到 **$50.400^\circ$**。
1.  **計算觀測面間距 $d_{obs}$**:
    $$ \theta_{obs} = 25.200^\circ $$
    $$ d_{obs} = \frac{1.540562}{2 \sin(25.200^\circ)} = 1.8091 \mathring{A} $$

**步驟 2：計算應變 (Strain)**
$$ \epsilon_z = \frac{1.8091 - 1.8080}{1.8080} = +6.13 \times 10^{-4} $$
晶格垂直方向被拉長了約 0.06%。

**步驟 3：計算應力 (Stress)**
$$ \sigma = - \frac{66.7 \times 10^3 \text{ MPa}}{2 \times 0.419} \times (6.13 \times 10^{-4}) $$
$$ \sigma = - 79594 \times 0.000613 = \mathbf{-48.8 \text{ MPa}} $$

**(比較) 若錯誤使用多晶常數 ($E=130, \nu=0.34$)會有什麼後果？**
$$ \sigma_{wrong} = - \frac{130000}{0.68} \times 0.000613 = \mathbf{-117.2 \text{ MPa}} $$
**誤差高達 140%！** 這解釋了為什麼傳統 XRD 算出的應力常常與晶圓翹曲 (Wafer Bow) 量測對不上。

---

## 5.4 工程診斷意義 (Engineering Interpretation)

| 應力類型 | 物理現象 | 工程風險 |
| :--- | :--- | :--- |
| **壓縮應力 (< 0)** | 晶格側向被擠壓 (垂直膨脹) | **Stress Migration (SM)**。原子被擠得想跑，容易在長時間存放後長成 **小山丘 (Hillocks)**，造成線路短路。 |
| **拉伸應力 (> 0)** | 晶格側向被拉扯 (垂直收縮) | **Voiding (孔洞)**。晶界容易被拉開，導致電阻升高或開路。 |

### 為什麼推薦用 (200) 算應力？
雖然 (111) 訊號最強，但 (111) 對應力很不敏感 ($E=191$ GPa 太硬了，同樣應力下位移量只有 (200) 的 1/3)。這使得 (111) 的位移容易被儀器誤差淹沒。
**(200) 是最軟的方向**，像彈簧一樣，對應力反應最大最準，是監控薄膜應力的最佳指標。

---

## 5.5 引用文獻
1.  **Ledbetter, H. M., & Naimon, E. R. (1974)**. "Elastic Properties of Metals and Alloys. II. Copper". *J. Phys. Chem. Ref. Data*, 3(4).
    *   提供 $C_{11}, C_{12}, C_{44}$ 用於計算 $E_{hkl}$。
2.  **Nye, J. F. (1957)**. *Physical Properties of Crystals*. Oxford University Press.
    *   各向異性彈性力學的聖經。

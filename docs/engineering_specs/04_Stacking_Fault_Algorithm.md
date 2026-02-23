# 第四章：堆疊層錯分析算法 (Stacking Fault Analysis Algorithm)

**文件編號**: DOC-ENG-04
**模組位置**: `xrd_analysis.methods.defect_analysis.StackingFaultAnalyzer`
**物理核心**: Warren-Averbach Planar Defect Theory

---

## 4.1 晶體物理根源：對稱性破缺 (Crystal Physics Origin)
面心立方 (FCC) 金屬如銅，其原子堆疊順序為 **ABCABC...**。
當晶體在電鍍生長或塑性變形過程中發生錯誤，順序變成 **ABC|BC|ABC...** (抽取型層錯) 或 **ABC|B|ABC...** (插入型層錯)，這種二維的面缺陷稱為 **堆疊層錯 (Stacking Faults)**。

### 為什麼這會導致峰位移？(Why Peaks Shift?)
在倒置空間 (Reciprocal Space) 中，層錯破壞了晶格的平移對稱性，導致倒格點 (Reciprocal Lattice Points) 沿特定方向拉長並移動。
根據 **Paterson (1952)** 與 **Warren (1969)** 的推導，這會導致 X 光繞射峰發生規律性位移：
1.  **(111) 晶面**: 向 **高角度** 位移 ($+$)。
2.  **(200) 晶面**: 向 **低角度** 位移 ($-$)。

這兩個峰就像被一隻隱形的手「互相推近」，導致它們的 **峰間距 (Peak Separation)** 縮小。

---

## 4.2 數學模型：Warren 幾何係數 G (Mathematical Model)
我們不依賴單一峰的絕對位移 (受樣品高度誤差干擾大)，而是測量 **(200)-(111) 峰間距的收縮量**。

**定義式**:
$$ \alpha = \frac{ \Delta SEP_{obs} }{ G } = \frac{ (2\theta_{200}-2\theta_{111})_{obs} - (2\theta_{200}-2\theta_{111})_{std} }{ G } $$
其中 $\alpha$ 為層錯機率 (Stacking Fault Probability)，即每多少層原子中出現一層錯誤的比率。

### 係數 G 的嚴格推導 (Derivation of G)
Warren 證明了峰位移量 $\Delta 2\theta$ 與 $\alpha$ 的關係：
$$ \Delta 2\theta_{hkl} = \pm \frac{90\sqrt{3} \tan\theta}{\pi^2 h_0^2 (u+b)} \sum (\pm L_0) \alpha $$

針對銅的 (111) 與 (200) 晶面，代入幾何參數：
1.  **(111) 位移**: $\Delta 2\theta_{111} = + \frac{45\sqrt{3}}{4\pi^2} \alpha \tan\theta_{111}$ (係數修正：考慮多晶平均)
2.  **(200) 位移**: $\Delta 2\theta_{200} = - \frac{90\sqrt{3}}{2\pi^2} \alpha \tan\theta_{200}$
3.  **(220) 位移**: $\Delta 2\theta_{220} = + \frac{45\sqrt{3}}{4\pi^2} \alpha \tan\theta_{220}$

**峰間距縮小量 $\Delta SEP$**:
$$ \Delta SEP = \Delta 2\theta_{200} - \Delta 2\theta_{111} = - \frac{45\sqrt{3}}{\pi^2} (\tan\theta_{111} + 2\tan\theta_{200}) \alpha $$

代入銅的標準角度 ($\theta_{111} \approx 21.66^\circ, \theta_{200} \approx 25.22^\circ, \theta_{220} \approx 37.06^\circ$)：
$$ G \approx -7.897 \text{ deg/unit } \alpha $$

此即程式碼常數 `WARREN_G_COEFFICIENT = -7.897` 的物理來源。

---

## 4.3 數值演算示範 (Numerical Trace)

為驗證算法靈敏度，我們模擬一個 **"0.5% 層錯率"** 的真實情境。這在電鍍銅中屬於「輕微但可測」的缺陷濃度。

**標準狀態 (Standard State, $\alpha=0$)**:
*   $\lambda = 1.540562 \mathring{A}$
*   $d_{111} = 2.0871 \mathring{A} \rightarrow 2\theta_{111} = 43.316^\circ$
*   $d_{111} = 2.0871 \mathring{A} \rightarrow 2\theta_{111} = 43.316^\circ$
*   $d_{200} = 1.8075 \mathring{A} \rightarrow 2\theta_{200} = 50.448^\circ$
*   $d_{220} = 1.2781 \mathring{A} \rightarrow 2\theta_{220} = 74.126^\circ$
*   **標準間距 $SEP_{std} = 50.448 - 43.316 = 7.132^\circ$**

**含缺陷狀態 (Defected State, $\alpha=0.005$)**:
假設樣品含有 0.5% 的層錯 (每 200 層錯 1 層)。

1.  **計算位移量 (Peak Shift)**:
    $$ \text{Shift} = G \times \alpha = -7.897 \times 0.005 = -0.0395^\circ $$
    也就是說，峰間距會縮小約 0.04 度。

2.  **觀測間距 (Observed Separation)**:
    $$ SEP_{obs} = SEP_{std} + \text{Shift} = 7.132 - 0.0395 = 7.0925^\circ $$

3.  **程式逆運算 (Algorithm Execution)**:
    當我們的程式碼讀到 $SEP_{obs} = 7.0925^\circ$ 時：
    $$ \text{Deviation} = 7.0925 - 7.132 = -0.0395^\circ $$
    $$ \alpha = \frac{-0.0395}{-7.897} = 0.005001 $$
    $$ \alpha_{\%} = 0.50 \% $$

**(220) 峰位檢驗 (220 Check)**:
雖然我們不使用 (220) 計算層錯 (因為高角峰強度較低且可能與 (311) 重疊)，但可以預測其行為：
$$ \text{Shift}_{220} = + \frac{45\sqrt{3}}{4\pi^2} \alpha \tan\theta_{220} $$
$$ \text{Shift}_{220} \approx + 1.974 \times 0.005 \times \tan(37.06^\circ) \approx + 0.007^\circ $$
預測 (220) 峰位將微幅移至 $74.133^\circ$。由於位移量極小 (僅約 (200) 位移量的 1/5)，因此不適合作為主要計算依據。

**工程結論**:
即使峰間距只縮小了微不足道的 **0.04 度** (這需要高解析度儀器才能分辨)，我們的算法也能精確反推出樣品含有 **0.5%** 的層錯。這證明了此算法對微觀結構變化具有極高的靈敏度。

---

## 4.4 工程診斷意義 (Engineering Interpretation)

當程式輸出 $\alpha$ 值時，代表什麼物理意義？

| 層錯機率 $\alpha$ | 狀態描述 | 工程影響 (Impact on Properties) |
| :--- | :--- | :--- |
| **< 0.1%** | 近似完美晶體 | 導電率極佳，但機械強度可能較低 (太軟)。 |
| **0.2% - 0.5%** | 輕微缺陷 | **最佳平衡區**。適量的層錯能阻礙差排滑移，提升硬度，同時不顯著降低導電率。 |
| **> 1.0%** | 高缺陷密度 | 硬度極高 (加工硬化)，但電阻率顯著上升，且延展性 (Elongation) 大幅下降，容易脆斷。 |

### 與殘留應力的關係 (Stress Decoupling)
層錯會導致 (111) 向高角度位移，這與「壓縮應力」造成的位移方向相同。
**如果忽略層錯修正**，直接用 (111) 位置算應力，會算出 **虛假的巨大壓縮應力**。
本模組優先計算 $\alpha$，提醒使用者 (111) 峰位包含非應力因素，這就是為什麼我們在 `01_Peak_Diagnosis` 中建議使用受層錯影響較小的 (220) 峰來監控應力。

---

## 4.5 引用文獻
1.  **Warren, B. E. (1969)**. *X-ray Diffraction*. Dover Publications. (Chapter 13).
2.  **Paterson, M. S. (1952)**. "X-ray Diffraction by Face-Centered Cubic Crystals with Deformation Faults." *Journal of Applied Physics*, 23(8), 805-811.

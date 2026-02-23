# 第三章：織構分析算法 (Texture Analysis Algorithm)

**文件編號**: DOC-ENG-03
**模組位置**: `xrd_analysis.methods.texture.TextureAnalyzer`
**物理核心**: Harris Texture Coefficient (Inverse Pole Figure Method)

---

## 3.1 物理根源：擇優取向 (Preferred Orientation)
理想的多晶粉末，其晶粒取向是完全隨機分佈的。但在電鍍銅製程中，晶粒傾向於沿著特定的結晶學方向（如 <111> 或 <200>）生長以降低表面能。
這種非隨機的排列稱為 **擇優取向 (Preferred Orientation)** 或 **織構 (Texture)**。

### 為什麼這很重要？
*   **<111> 織構**: 原子排列最緊密，抗電遷移 (Electromigration) 能力最強，是半導體導線的黃金標準。
*   **<200> 織構**: 表面能最低，常見於退火後的再結晶晶粒。
*   **<220> 織構**: 常見於使用了過量抑制劑 (Suppressor) 的極細晶粒結構。

---

## 3.2 數學模型：Harris 紋理係數 (Mathematical Model)
我們使用 Harris (1952) 提出的 **紋理係數 (Texture Coefficient, TC)** 來量化這種效應。這是一種「反極圖 (Inverse Pole Figure)」的簡化算法，適用於 $ \theta-2\theta $ 掃瞄。

**定義式**:
$$ TC(hkl) = \frac{ \frac{I_{obs}(hkl)}{I_{std}(hkl)} }{ \frac{1}{N} \sum_{i=1}^{N} \frac{I_{obs}(h_i k_i l_i)}{I_{std}(h_i k_i l_i)} } $$

*   $I_{obs}$: 實驗測量的積分強度 (注意：必須是積分面積，非峰高，以涵蓋寬化效應)。
*   $I_{std}$: JCPDS 標準粉末衍射卡片 (無織構狀態) 的強度。
*   $N$: 參與計算的衍射峰總數。

**判定標準**:
*   $TC = 1$: 隨機取向 (無織構)。
*   $TC > 1$: 擇優取向 (Preferred)。
*   $TC < 1$:受抑取向 (Suppressed)。

---

## 3.3 數值演算示範 (Numerical Trace)

為了透徹理解 Harris 算法如何運作，我們代入一組 **高 <111> 織構** 的模擬數據進行手算。

**標準強度 (Standard Intensity, JCPDS 04-0836)**:
從 `xrd_analysis.core.copper_crystal` 模組讀取：
*   $I_{std}(111) = 100$
*   $I_{std}(200) = 46$
*   $I_{std}(220) = 20$

**觀測數據 (Observed Data)**:
假設我們測量到一組數據 (強度值已扣除背景)：
*   $I_{obs}(111) = 15000$
*   $I_{obs}(200) = 2000$
*   $I_{obs}(220) = 500$

### 步驟 1：計算歸一化比率 (Ratio Calculation)
$$ R_{111} = \frac{15000}{100} = 150.0 $$
$$ R_{200} = \frac{2000}{46} = 43.48 $$
$$ R_{220} = \frac{500}{20} = 25.0 $$

### 步驟 2：計算平均比率 (Average Ratio)
$$ R_{avg} = \frac{150.0 + 43.48 + 25.0}{3} = \frac{218.48}{3} = 72.827 $$

### 步驟 3：計算紋理係數 (TC Calculation)
$$ TC_{111} = \frac{R_{111}}{R_{avg}} = \frac{150.0}{72.827} = \mathbf{2.06} $$
$$ TC_{200} = \frac{R_{200}}{R_{avg}} = \frac{43.48}{72.827} = \mathbf{0.60} $$
$$ TC_{220} = \frac{R_{220}}{R_{avg}} = \frac{25.0}{72.827} = \mathbf{0.34} $$

**驗證總和**:
在 Harris 方法中，TC 的平均值必須由定義歸一化為 1。
$$ \text{Mean TC} = \frac{2.06 + 0.60 + 0.34}{3} = \frac{3.0}{3} = 1.0 \quad \text{(Correct)} $$

**工程結論**:
此樣品的 $TC(111) = 2.06$，遠大於 1。
物理上，這代表在樣品表面參與衍射的晶粒中，(111) 面平行的晶粒數量是隨機分佈的 **2 倍**。這是一個非常強的 (111) 擇優取向訊號。

---

## 3.4 工程診斷與應用 (Diagnosis & Application)

在 `compare_integral_methods` 與 `fwhm_evolution` 圖表中，您經常看到與 FWHM 相對應的趨勢。

| 參數 | 物理意義 | 工程判讀 |
| :--- | :--- | :--- |
| **TC(111) > 1.5** | 強 (111) 織構 | **優選**。抗電遷移能力強，適合細線路製程。通常伴隨著較窄的 (111) FWHM。 |
| **TC(200) > 1.2** | (200) 再結晶 | **警示**。這通常發生在長時間退火後，或者是使用了錯誤的電鍍參數 (如電流密度過高)。通常伴隨著 (200) FWHM 的急劇下降。 |
| **織構度 (σ)** | $\sqrt{\frac{\sum (TC-1)^2}{N}}$ | 數值越高，代表織構越強烈。$\sigma < 0.2$ 視為隨機結構。 |

### 為什麼這與 FWHM 圖表相關？
您提到的 `fwhm_evolution` 其實是織構的一面鏡子。
在強 (111) 織構下，(111) 晶粒會為了降低表面能而優先長大 (Coalescence)，這導致 **(111) 的 FWHM 會比其他晶面更早開始下降**。我們可以通過比對 TC 曲線與 FWHM 曲線來確認這一點。

---

## 3.5 引用文獻
1.  **Harris, G. B. (1952)**. "Quantitative Measurement of Preferred Orientation in Rolled Uranium Bars." *Philosophical Magazine*, 43(336), 113-123.
2.  **JCPDS Card 04-0836**. Joint Committee on Powder Diffraction Standards.

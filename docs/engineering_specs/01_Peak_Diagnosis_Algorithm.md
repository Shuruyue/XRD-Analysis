# 第一章：峰形診斷與雙峰去卷積算法 (Peak Diagnosis & Doublet Deconvolution Algorithm)

**文件編號**: DOC-ENG-01
**模組位置**: `xrd_analysis.fitting.peak_fitter.fit_peak_with_diagnosis`
**物理核心**: Quantum Mechanical Spin-Orbit Coupling & True Voigt Convolution

---

## 1.1 量子力學根源：自旋-軌道耦合 (Quantum Origin: Spin-Orbit Coupling)
在 X 光繞射圖譜中，每一個布拉格回波 (Bragg Reflection) 皆分裂為 $K\alpha_1$ 與 $K\alpha_2$ 雙重譜線。此現象源於銅原子內部電子的 **自旋-軌道耦合 (Spin-Orbit Coupling, L-S Coupling)**。

### A. 能級躍遷選擇定則 (Selection Rules)
X 光特徵輻射源自電子從 $L$ 殼層 ($n=2$) 躍遷回 $K$ 殼層 ($n=1$)。
根據量子力學選擇定則：
1.  主量子數變化任意 ($\Delta n \neq 0$)。
2.  角動量量子數變化 $\Delta l = \pm 1$。
3.  總角動量量子數變化 $\Delta j = 0, \pm 1$。

因此，$L$ 殼層中的 $2s$ 軌域 ($l=0$) 無法躍遷至 $1s$ ($l=0$)，因為 $\Delta l=0$ 是 **禁戒躍遷 (Forbidden Transition)**。僅有 $2p$ 軌域 ($l=1$) 的電子可躍遷。

### B. 精細結構分裂 (Fine Structure Splitting)
$2p$ 電子的自旋 ($s=1/2$) 與軌道角動量 ($l=1$) 發生磁矩耦合，導致能級分裂為兩個狀態：
1.  **$K\alpha_1$ 來源 ($2P_{3/2}$)**: $j = l + s = 1 + 1/2 = 3/2$。能量較高，波長較短。
2.  **$K\alpha_2$ 來源 ($2P_{1/2}$)**: $j = l - s = 1 - 1/2 = 1/2$。能量較低，波長較長。

### C. 強度比 0.5 的嚴格推導 (Derivation of Intensity Ratio)
譜線強度正比於該能級的 **簡併度 (Degeneracy)**，公式為 $N = 2j+1$。
*   **$2P_{3/2}$ ($K\alpha_1$)**: $N = 2(3/2) + 1 = 4$。
*   **$2P_{1/2}$ ($K\alpha_2$)**: $N = 2(1/2) + 1 = 2$。

$$ \text{Intensity Ratio} = \frac{I(K\alpha_2)}{I(K\alpha_1)} = \frac{2}{4} = 0.5 $$

此即為程式碼中 `intensity_ratio = 0.5` 的物理定律來源，非經驗常數。

---

## 1.2 數學模型：True Voigt 卷積 (Mathematical Model)
本系統採用的擬合函數 $V(x)$ 為 **True Voigt**，定義為高斯函數 ($G$) 與勞倫茲函數 ($L$) 的數學卷積：

$$ V(x; \sigma, \gamma) = G(x; \sigma) \otimes L(x; \gamma) = \int_{-\infty}^{\infty} G(x') L(x-x') dx' $$

本專案使用 **Faddeeva 函數 $w(z)$** 進行計算，這是 Voigt 函數的複數解析解：
$$ V(x) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}}, \quad z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$
*   $x$: $2\theta - 2\theta_0$ (偏離角度)
*   $\sigma$: 高斯標準差 (對應微觀應變)
*   $\gamma$: 勞倫茲半寬 (對應晶粒尺寸倒數)

擬合總函數為雙峰疊加：
$$ I_{total}(2\theta) = A \cdot V(2\theta - \theta_1) + 0.5A \cdot V(2\theta - \theta_2) + Bkg $$

---

## 1.2.1 代入運算詳解 (Step-by-Step Substitution Trace)

您詢問 **「這東西要怎麼代？」**。以下我們示範如何將具體數值代入上述複數公式，計算出峰頂的一個點。

**情境假設 (Hypothetical Case)**:
*   我們要計算 **峰頂中心點** 的強度 ($x=0$)。
*   假設擬合參數：$\sigma = 0.05$ (高斯寬), $\gamma = 0.05$ (勞倫茲寬)。

**步驟 1：計算複數變數 z**
$$ z = \frac{x + i\gamma}{\sigma\sqrt{2}} $$
代入數值：
$$ z = \frac{0 + 0.05i}{0.05 \times 1.414} = \frac{0.05i}{0.0707} = \mathbf{0.707i} $$
(這是一個純虛數)

**步驟 2：計算 Faddeeva 函數 w(z)**
這一步需要查表或程式計算 (如 Python `scipy.special.wofz`)：
$$ w(0.707i) = \text{erfcx}(0.707) \approx \mathbf{0.5231} $$
(對於純虛數輸入，$w(z)$ 的輸出是實數)

**步驟 3：計算 Voigt 函數值 V(x)**
$$ V(0) = \frac{\text{Re}[w(z)]}{\sigma\sqrt{2\pi}} $$
代入數值：
$$ V(0) = \frac{0.5231}{0.05 \times 2.5066} = \frac{0.5231}{0.1253} = \mathbf{4.174} $$

**步驟 4：代入雙峰總強度公式**
假設振幅 $A=1000$ (計數值)，主峰中心 $\theta_1=43^\circ$。
我們想算 $\theta=43^\circ$ 處的強度：
*   **$K\alpha_1$ 項**: $1000 \times V(0) = 1000 \times 4.174 = 4174$
*   **$K\alpha_2$ 項**: 距離 $0.1^\circ$ 處。假設算出 $V(-0.1) \approx 0.5$。強度 $= 500 \times 0.5 = 250$。
*   **總強度**: $4174 + 250 + \text{Bkg}$。

這就是程式碼內部如何將 **$z$ (複數座標)** 轉換為 **$I$ (實數強度)** 的過程。

---

## 1.3 數值演算示範 (Numerical Trace)

為驗證算法的正確性，我們代入您提供的圖片數據進行「從頭到尾」的實數演算。

**輸入數據 (Input Data from Image)**:
*   晶面 (HKL): (111)
*   擬合中心 ($2\theta_{fit}$): **43.353°** ($K\alpha_1$ 位置)
*   波長常數: $\lambda_1 = 1.540562 \mathring{A}, \lambda_2 = 1.544390 \mathring{A}$

### 第一步：計算 $K\alpha_2$ 的理論位置
利用布拉格定律 $n\lambda = 2d\sin\theta$，在相同晶面间距 $d$ 下，衍射角與波長成正比關係：
$$ \sin(\frac{2\theta_2}{2}) = \frac{\lambda_2}{\lambda_1} \sin(\frac{2\theta_1}{2}) $$

1.  **計算 $K\alpha_1$ 的 $\theta_1$**:
    $$ \theta_1 = \frac{43.353}{2} = 21.6765^\circ $$
    $$ \sin(21.6765^\circ) = 0.369369 $$

2.  **計算 $K\alpha_2$ 的 $\sin\theta_2$**:
    $$ \text{Ratio} = \frac{1.544390}{1.540562} = 1.0024848 $$
    $$ \sin\theta_2 = 0.369369 \times 1.0024848 = 0.370287 $$

3.  **反推 $2\theta_2$**:
    $$ \theta_2 = \arcsin(0.370287) = 21.7331^\circ $$
    $$ 2\theta_2 = 21.7331 \times 2 = \mathbf{43.466^\circ} $$

**驗證結果**: 請比對您圖片右側紅色文字框中的標註：`Ka-43.466`。與本手算結果 **完全吻合 (Exact Match)**。

### 第二步：總強度合成 (Intensity Summation)
在 $2\theta = 43.353^\circ$ (主峰頂點) 處，總強度 $I_{total}$ 由兩部分組成：

1.  **$K\alpha_1$ 貢獻**:
    位於中心 ($x=0$)。Voigt 函數 $V(0)$ 處於最大值。設強度為 $I_{max}$。
    貢獻量 $\approx I_{max}$。

2.  **$K\alpha_2$ 貢獻**:
    位於 $0.113^\circ$ 之外。
    由於 (111) 峰寬 FWHM 為 $0.248^\circ$，半寬 HWHM $\approx 0.124^\circ$。
    此時 $K\alpha_2$ 的中心距離觀測點約 1 倍 HWHM。
    根據 Voigt 函數特性，在 1 倍 HWHM 處強度衰減至約 50%。
    加上原本強度已被量子力學限制為 0.5 倍。
    貢獻量 $\approx I_{max} \times 0.5 \times 0.5 = 0.25 I_{max}$。

**工程結論**:
在主峰頂點處，有約 **20-25%** 的信號其實是來自右邊那個「看不見的 $K\alpha_2$ 峰」的尾部延伸。
若程式碼沒有執行此雙峰去卷積 (Doublet Deconvolution)，直接測量峰高會導致 25% 的誤差，且峰寬會被錯誤地拉寬，導致計算出的晶粒尺寸偏小。

---

## 1.4 工程診斷指標解釋 (Engineering Diagnosis)

### A. 擬合優度 (Coefficient of Determination, $R^2$)
$$ R^2 = 1 - \frac{SS_{res}}{SS_{tot}} $$
*   **數值**: 圖片顯示 $R^2 = 0.9997$。
*   **工程解釋**: 模型解釋了 99.97% 的訊號變異。這並非僅僅代表「擬合得好」，更深層的物理意義是：**樣品內部應變場極度均勻**。若存在晶格梯度或第二相，殘差 (Residuals) 會呈現波浪狀，導致 $R^2$ 下降。

### B. 混合參數 (Mixing Parameter, $\eta$)
$$ \eta = \frac{\Gamma_L}{\Gamma_L + \Gamma_G} $$
*   **數值**: 圖片顯示 $\eta = 0.916$ (接近 1.0)。
*   **工程解釋**: 此峰形呈現高度 **Lorentzian** 特徵。
    *   **Lorentzian 主導**: 代表寬化機制主要來自 **晶粒尺寸有限 (Size Broadening)**。
    *   **Gaussian 主導**: 代表寬化機制主要來自 **微觀應變 (Microstrain)**。
    *   **結論**: 此電鍍銅層的 (111) 晶粒非常細小 (奈米晶)，但內部的差排 (Dislocation) 分佈相對單純，沒有造成劇烈的高斯型晶格扭曲。

### C. 半高寬 (FWHM)
*   **數值**: (111) $0.248^\circ$ vs (200) $0.534^\circ$。
*   **工程解釋**: (200) 峰寬是 (111) 的兩倍以上。這證實了該鍍層具有 **(111) 擇優取向 (Preferred Orientation)**，導致 (111) 晶粒沿生長方向融合長大 (尺寸大 $\to$ 峰窄)，而 (200) 晶粒生長受抑 (尺寸小 $\to$ 峰寬)。

---

## 1.5 程式模組對照
*   **TrueVoigt**: `xrd_analysis.fitting.pseudo_voigt.py` (執行 `wofz` 複數誤差函數運算)

---

## 1.6 引用文獻 (References)

1.  **Cullity, B. D. (1978)**. *Elements of X-Ray Diffraction*, 2nd Ed., Addison-Wesley.
    *   Chapter 2: Physics of X-rays (Spin-Orbit Coupling & Selection Rules).
2.  **Armstrong, B. H. (1967)**. "Spectrum Line Profiles: The Voigt Function". *Journal of Quantitative Spectroscopy and Radiative Transfer*, 7(1), 61-88.
    *   Voigt 函數的數學定義與性質。
3.  **Poppe, G. P. M., & Wijers, C. M. J. (1990)**. "Algorithm 680: Evaluation of the Complex Error Function". *ACM Transactions on Mathematical Software*, 16(1), 47.
    *   `scipy.special.wofz` 背後的 Fadeeva 函數快速演算法來源。

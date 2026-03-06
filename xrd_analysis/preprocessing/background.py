"""Background Subtraction Module.
Implements Sonneveld-Visser and Chebyshev polynomial methods.
"""

import numpy as np
from scipy import sparse
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve


class BackgroundSubtractor:
    """Background subtraction for XRD data.

    Three methods available:
    1. Sonneveld-Visser: Second derivative method
    2. Chebyshev: Polynomial fitting to peak-free regions
    3. SNIP: Statistics-sensitive Non-linear Iterative Peak-clipping

    Reference (SNIP):
        Ryan, C. G., et al. (1988). NIM-B 34, 396-402.
        Morháč, M. & Matoušek, V. (2008). Appl. Spectrosc. 62, 91-106.
    """

    def __init__(self, method: str = "chebyshev", **kwargs):
        """Initialize background subtractor.

        Args:
            method: "chebyshev", "sonneveld_visser", or "snip"
            **kwargs: Method-specific parameters
                - poly_degree: For Chebyshev method (default: 5)
                - iterations: For Sonneveld-Visser (default: 100)
                - snip_iterations: For SNIP method (default: 24)

        """
        self.method = method.lower()
        self.poly_degree = kwargs.get("poly_degree", 5)
        self.iterations = kwargs.get("iterations", 100)
        self.snip_iterations = kwargs.get("snip_iterations", 24)

        self.als_lambda = kwargs.get("als_lambda", 1e6)
        self.als_p = kwargs.get("als_p", 0.01)
        self.als_iterations = kwargs.get("als_iterations", 15)

        if self.method not in ["chebyshev", "sonneveld_visser", "snip", "als"]:
            raise ValueError(f"Unknown method: {method}")

    def subtract(
        self, two_theta: np.ndarray, intensity: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Subtract background from XRD data.

        Args:
            two_theta: 2θ angle array
            intensity: Intensity array

        Returns:
            Tuple of (corrected_intensity, background)

        """
        if self.method == "chebyshev":
            return self._chebyshev_background(two_theta, intensity)
        elif self.method == "snip":
            return self._snip_background(two_theta, intensity)
        elif self.method == "als":
            return self._als_background(two_theta, intensity)
        else:
            return self._sonneveld_visser_background(two_theta, intensity)

    def _chebyshev_background(
        self, two_theta: np.ndarray, intensity: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Chebyshev polynomial background fitting.

        Iteratively fits polynomial, excluding points above the fit.
        """
        # Normalize x for numerical stability
        x_norm = (
            2 * (two_theta - two_theta.min()) / (two_theta.max() - two_theta.min()) - 1
        )

        # Initial fit
        mask = np.ones(len(intensity), dtype=bool)

        for _ in range(self.iterations):
            # Fit Chebyshev polynomial to masked data
            coeffs = np.polynomial.chebyshev.chebfit(
                x_norm[mask], intensity[mask], self.poly_degree
            )
            background = np.polynomial.chebyshev.chebval(x_norm, coeffs)

            # Update mask: keep points below or near background
            new_mask = intensity <= background * 1.05

            if np.array_equal(mask, new_mask):
                break
            mask = new_mask

        corrected = intensity - background
        corrected[corrected < 0] = 0  # Ensure non-negative

        return corrected, background

    def _sonneveld_visser_background(
        self, two_theta: np.ndarray, intensity: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Sonneveld-Visser second derivative method.

        Uses second derivative to identify background regions.
        """
        # Compute smoothed second derivative
        window = min(21, len(intensity) // 10)
        if window % 2 == 0:
            window += 1

        deriv2 = savgol_filter(intensity, window, 3, deriv=2)

        # Background points have near-zero second derivative
        threshold = np.std(deriv2) * 0.1
        bg_mask = np.abs(deriv2) < threshold

        if np.sum(bg_mask) < 10:
            # Fallback: use minimum envelope
            bg_mask = intensity < np.percentile(intensity, 20)

        # Interpolate background from identified points
        bg_indices = np.where(bg_mask)[0]
        bg_values = intensity[bg_mask]

        background = np.interp(np.arange(len(intensity)), bg_indices, bg_values)

        corrected = intensity - background
        corrected[corrected < 0] = 0

        return corrected, background

    def _als_background(
        self, two_theta: np.ndarray, intensity: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Asymmetric Least Squares (ALS) background estimation.

        Minimizes: Σ w_i (y_i - z_i)² + λ Σ (Δ²z_i)²
        where w_i is asymmetric: p for y_i > z_i, (1-p) for y_i ≤ z_i.

        Reference:
            Eilers, P. H. C., & Boelens, H. F. M. (2005).
            "Baseline estimation with asymmetric least squares smoothing."
            Leiden University Medical Centre Report.

        """
        background = als_baseline(
            intensity,
            lam=self.als_lambda,
            p=self.als_p,
            n_iter=self.als_iterations,
        )
        corrected = intensity - background
        corrected[corrected < 0] = 0
        return corrected, background

    def _snip_background(
        self, two_theta: np.ndarray, intensity: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """SNIP: Statistics-sensitive Non-linear Iterative Peak-clipping.

        Operates in log-log-sqrt space to linearize background, then
        applies an iterative clipping operator to remove peaks.

        Reference:
            Ryan, C. G., Clayton, E., Griffin, W. L., Sie, S. H., & Cousens, D. R.
            (1988). NIM-B 34, 396-402. DOI: 10.1016/0168-583X(88)90063-8

            Morháč, M. & Matoušek, V. (2008). Appl. Spectrosc. 62, 91-106.
        """
        background = snip_background(intensity, n_iterations=self.snip_iterations)
        corrected = intensity - background
        corrected[corrected < 0] = 0
        return corrected, background


def snip_background(
    y: np.ndarray,
    n_iterations: int = 24,
    decreasing: bool = True,
) -> np.ndarray:
    """SNIP background estimation (vectorized).

    Statistics-sensitive Non-linear Iterative Peak-clipping algorithm.
    Operates in LLS (log-log-sqrt) space to handle Poisson counting statistics.

    Args:
        y: 1-D intensity array (counts)
        n_iterations: Number of clipping iterations (typically 15-40 for XRD)
        decreasing: If True, iterate from wide to narrow window (more robust)

    Returns:
        Estimated background array, same shape as y

    Reference:
        Ryan et al. (1988), NIM-B 34, 396-402.
        Morháč & Matoušek (2008), Appl. Spectrosc. 62, 91-106.

    """
    # Transform to LLS space: log(log(sqrt(y+1)+1)+1)
    # Guard against zeros/negatives
    y_safe = np.maximum(y, 0.0)
    v = np.log(np.log(np.sqrt(y_safe + 1.0) + 1.0) + 1.0)

    # Iterative peak clipping with decreasing (or increasing) window
    window_range: range
    if decreasing:
        window_range = range(n_iterations, 0, -1)
    else:
        window_range = range(1, n_iterations + 1)

    for w in window_range:
        # Vectorized clipping: v[i] = min(v[i], (v[i-w] + v[i+w]) / 2)
        left = np.empty_like(v)
        right = np.empty_like(v)
        left[w:] = v[:-w]
        left[:w] = v[:w]
        right[:-w] = v[w:]
        right[-w:] = v[-w:]
        avg = (left + right) / 2.0
        v = np.minimum(v, avg)

    # Inverse LLS transform: y = (exp(exp(v) - 1) - 1)² - 1
    background = (np.exp(np.exp(v) - 1.0) - 1.0) ** 2 - 1.0
    background = np.maximum(background, 0.0)

    return background


def als_baseline(
    y: np.ndarray,
    lam: float = 1e6,
    p: float = 0.01,
    n_iter: int = 15,
) -> np.ndarray:
    """Asymmetric Least Squares (ALS) baseline estimation.

    Iteratively solves a penalized least squares problem with
    asymmetric weights to estimate the baseline of spectral data.

    Args:
        y: 1-D intensity array
        lam: Smoothness parameter (larger → smoother baseline, typical: 1e5–1e7)
        p: Asymmetry parameter (0 < p < 1, typically 0.001–0.05).
            Smaller p penalizes overshoot above baseline more strongly.
        n_iter: Number of reweighting iterations

    Returns:
        Estimated baseline array

    Reference:
        Eilers & Boelens (2005), Leiden University Medical Centre Report.

    """
    n = len(y)
    # Second-order difference matrix D
    D = sparse.diags([1, -2, 1], [0, 1, 2], shape=(n - 2, n))
    D_T_D = lam * D.T @ D

    w = np.ones(n)
    for _ in range(n_iter):
        W = sparse.diags(w, 0)
        Z = W + D_T_D
        z = spsolve(Z, w * y)
        # Asymmetric reweighting
        w = np.where(y > z, p, 1 - p)

    return z


def subtract_background(
    two_theta: np.ndarray, intensity: np.ndarray, method: str = "chebyshev", **kwargs
) -> tuple[np.ndarray, np.ndarray]:
    """Convenience function for background subtraction.

    Args:
        two_theta: 2θ angle array
        intensity: Intensity array
        method: "chebyshev" or "sonneveld_visser"
        **kwargs: Method-specific parameters

    Returns:
        Tuple of (corrected_intensity, background)

    """
    subtractor = BackgroundSubtractor(method, **kwargs)
    return subtractor.subtract(two_theta, intensity)

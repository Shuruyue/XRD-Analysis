"""Peak Detection Module.
==================================
Implements second derivative method for automatic peak finding.
"""

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks as scipy_find_peaks
from scipy.signal import savgol_filter


@dataclass
class PeakInfo:
    """Data class for peak information."""

    index: int
    two_theta: float
    intensity: float
    estimated_fwhm: float


class PeakDetector:
    """Automatic peak detection using second derivative method.

    The second derivative identifies inflection points where
    peaks have local minima (negative curvature).
    """

    def __init__(
        self,
        min_height: float = 100,
        min_distance: float = 0.5,
        window_size: int = 11,
        poly_order: int = 3,
    ):
        """Initialize peak detector.

        Args:
            min_height: Minimum peak intensity threshold
            min_distance: Minimum distance between peaks (degrees)
            window_size: Savitzky-Golay window size
            poly_order: Savitzky-Golay polynomial order

        """
        self.min_height = min_height
        self.min_distance = min_distance
        self.window_size = window_size
        self.poly_order = poly_order

    def detect(self, two_theta: np.ndarray, intensity: np.ndarray) -> list[PeakInfo]:
        """Detect peaks in XRD data.

        Args:
            two_theta: 2θ angle array
            intensity: Intensity array

        Returns:
            List of PeakInfo objects

        """
        # Calculate step size
        step = np.mean(np.diff(two_theta))
        min_distance_points = max(1, int(self.min_distance / step))

        # Compute second derivative for peak position refinement
        deriv2 = savgol_filter(intensity, self.window_size, self.poly_order, deriv=2)

        # Find peaks in original data
        peaks_idx, properties = scipy_find_peaks(
            intensity,
            height=self.min_height,
            distance=min_distance_points,
            prominence=self.min_height * 0.1,
        )

        # Refine peak positions using second derivative minima
        # True peak centers correspond to minima in the second derivative
        refined_idx = self._refine_peak_positions(
            two_theta, intensity, deriv2, peaks_idx
        )

        # Estimate FWHM for each peak
        peaks = []
        for idx in refined_idx:
            fwhm = self._estimate_fwhm(two_theta, intensity, idx)
            peaks.append(
                PeakInfo(
                    index=idx,
                    two_theta=two_theta[idx],
                    intensity=intensity[idx],
                    estimated_fwhm=fwhm,
                )
            )

        return peaks

    def _refine_peak_positions(
        self,
        two_theta: np.ndarray,
        intensity: np.ndarray,
        deriv2: np.ndarray,
        peaks_idx: np.ndarray,
    ) -> list[int]:
        """Refine peak positions using second derivative minima.

        True peak centers correspond to minima (most negative values)
        in the second derivative of the intensity profile.

        Args:
            two_theta: 2θ angle array
            intensity: Intensity array
            deriv2: Second derivative of intensity
            peaks_idx: Initial peak indices from scipy.find_peaks

        Returns:
            List of refined peak indices

        """
        step = np.mean(np.diff(two_theta))
        search_radius = max(1, int(0.2 / step))  # Search within ±0.2°

        refined = []
        for idx in peaks_idx:
            lo = max(0, idx - search_radius)
            hi = min(len(deriv2) - 1, idx + search_radius)

            # Find the minimum of second derivative in the search window
            local_deriv2 = deriv2[lo : hi + 1]
            min_offset = np.argmin(local_deriv2)
            refined_idx = lo + min_offset

            # Only accept if the refined position still has reasonable intensity
            if intensity[refined_idx] >= intensity[idx] * 0.5:
                refined.append(refined_idx)
            else:
                refined.append(idx)

        return refined

    def _estimate_fwhm(
        self, two_theta: np.ndarray, intensity: np.ndarray, peak_idx: int
    ) -> float:
        """Estimate FWHM for a single peak using linear interpolation.

        Accounts for local background level and uses interpolation
        between data points for sub-step accuracy.

        Args:
            two_theta: 2θ array
            intensity: Intensity array
            peak_idx: Peak index

        Returns:
            Estimated FWHM in degrees

        """
        peak_height = intensity[peak_idx]

        # Estimate local background from the edges of a window around the peak
        window = min(50, len(intensity) // 4)
        lo = max(0, peak_idx - window)
        hi = min(len(intensity) - 1, peak_idx + window)
        local_bg = min(
            np.mean(intensity[lo : lo + 5]), np.mean(intensity[hi - 4 : hi + 1])
        )

        half_max = local_bg + (peak_height - local_bg) / 2.0

        # Find left crossing with interpolation
        left_idx = peak_idx
        while left_idx > 0 and intensity[left_idx] > half_max:
            left_idx -= 1
        # Interpolate between left_idx and left_idx+1
        if left_idx < peak_idx and intensity[left_idx + 1] != intensity[left_idx]:
            frac = (half_max - intensity[left_idx]) / (
                intensity[left_idx + 1] - intensity[left_idx]
            )
            left_pos = two_theta[left_idx] + frac * (
                two_theta[left_idx + 1] - two_theta[left_idx]
            )
        else:
            left_pos = two_theta[left_idx]

        # Find right crossing with interpolation
        right_idx = peak_idx
        while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
            right_idx += 1
        # Interpolate between right_idx-1 and right_idx
        if right_idx > peak_idx and intensity[right_idx - 1] != intensity[right_idx]:
            frac = (half_max - intensity[right_idx]) / (
                intensity[right_idx - 1] - intensity[right_idx]
            )
            right_pos = two_theta[right_idx] + frac * (
                two_theta[right_idx - 1] - two_theta[right_idx]
            )
        else:
            right_pos = two_theta[right_idx]

        fwhm = right_pos - left_pos
        return max(fwhm, 0.05)  # Minimum 0.05 degrees


def find_peaks(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    min_height: float = 100,
    min_distance: float = 0.5,
) -> list[PeakInfo]:
    """Convenience function for peak detection.

    Args:
        two_theta: 2θ angle array
        intensity: Intensity array
        min_height: Minimum peak intensity
        min_distance: Minimum distance between peaks (degrees)

    Returns:
        List of PeakInfo objects

    """
    detector = PeakDetector(min_height=min_height, min_distance=min_distance)
    return detector.detect(two_theta, intensity)

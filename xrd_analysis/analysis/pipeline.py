"""xrd_analysis Complete Analysis Pipeline.
==========================================

Unified pipeline integrating all analysis phases.

- Phase 02: Preprocessing
- Phase 03: Peak Fitting
- Phase 04: Scherrer Size Calculation
- Phase 05: Williamson-Hall Strain Analysis
- Phase 06: Texture Analysis
- Phase 07: Defect Diagnosis
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from xrd_analysis.core.config import ParameterConfig
from xrd_analysis.core.constants import CU_KA1
from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED
from xrd_analysis.methods.defect_analysis import (
    AnnealingState,
    LatticeConstantResult,
    LatticeMonitor,
    StackingFaultAnalyzer,
    StackingFaultResult,
    determine_annealing_state,
)
from xrd_analysis.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
)
from xrd_analysis.methods.texture import (
    TextureAnalysisResult,
    TextureAnalyzer,
)
from xrd_analysis.methods.williamson_hall import (
    WHResult,
    WilliamsonHallAnalyzer,
)
from xrd_analysis.preprocessing.pipeline import PreprocessingPipeline

logger = logging.getLogger(__name__)

# =============================================================================
# Configuration
# =============================================================================


@dataclass
class AnalysisConfig:
    """Configuration for xrd_analysis analysis pipeline."""

    # X-ray parameters
    # Default: Cu Kα1 from constants (Bearden 1967)
    wavelength: float = CU_KA1

    # Scherrer parameters
    use_cubic_habit: bool = True

    # Peak detection
    peak_window: float = 2.0  # degrees around expected position
    min_intensity: float = 100  # minimum counts
    use_doublet_fitting: bool = True
    doublet_max_iterations: int = 20000
    min_fit_r_squared: float = 0.80

    # Preprocessing
    enable_smoothing: bool = True
    smoothing_window: int = 11
    smoothing_poly_order: int = 3
    enable_background: bool = True
    background_method: str = "chebyshev"
    background_degree: int = 5
    enable_kalpha_strip: bool = True

    # Instrumental broadening (Caglioti U, V, W)
    # FWHM_inst² = U·tan²θ + V·tanθ + W
    #
    # Default values (Empirical ranges for typical lab diffractometers):
    #   High-end lab:     W ≈ 0.001-0.002  (FWHM ≈ 0.032-0.045°)
    #   Standard lab:     W ≈ 0.002-0.005  (FWHM ≈ 0.045-0.071°)  ← Default
    #   Aging equipment:  W ≈ 0.008-0.015  (FWHM ≈ 0.089-0.122°)
    #
    # IMPORTANT: These are placeholder values for initial analysis.
    # For publication-quality results, calibrate using LaB6 (NIST SRM 660c):
    #   >>> xrd-analysis calibrate data/lab6_standard.txt
    #
    # Reference: Caglioti et al. (1958). Nucl. Instrum. 3, 223-228.
    caglioti_u: float = 0.0  # Default: 0 (typically -0.005 to +0.005)
    caglioti_v: float = 0.0  # Default: 0 (typically -0.002 to +0.002)
    caglioti_w: float = 0.003  # Default: standard lab (FWHM_inst ≈ 0.055°)

    # Cu peak positions from JCPDS 04-0836
    # Dynamically generated from CU_JCPDS constants
    EXPECTED_PEAKS: dict[tuple[int, int, int], float] = field(
        default_factory=lambda: {
            hkl: data["two_theta"] for hkl, data in CU_JCPDS_EXTENDED.items()
        }
    )

    @classmethod
    def from_parameter_config(cls, param_config: ParameterConfig) -> "AnalysisConfig":
        """Create AnalysisConfig from the unified ParameterConfig.

        Bridges the two configuration systems to avoid duplicated fields.

        Args:
            param_config: Unified parameter configuration

        Returns:
            AnalysisConfig populated from ParameterConfig values

        """
        inst = param_config.instrument
        peak = param_config.peak_detection
        return cls(
            wavelength=inst.wavelength,
            caglioti_u=inst.caglioti_u,
            caglioti_v=inst.caglioti_v,
            caglioti_w=inst.caglioti_w,
            peak_window=peak.peak_window,
            min_intensity=peak.min_intensity,
            min_fit_r_squared=param_config.validation.min_r_squared,
        )


# =============================================================================
# Data Loader (delegated to data_loader module)
# =============================================================================
# Backward-compatible re-exports
from xrd_analysis.analysis.data_loader import (  # noqa: E402,F401
    load_bruker_txt,
    parse_filename,
)

# =============================================================================
# Peak Finding (delegated to peak_finder module)
# =============================================================================
# Backward-compatible re-exports
from xrd_analysis.analysis.peak_finder import (  # noqa: E402,F401
    PeakData,
    find_peak_in_range,
)

# =============================================================================
# Result Container
# =============================================================================


@dataclass
class PipelineResult:
    """Complete pipeline analysis result."""

    # Input
    filepath: str
    sample_name: str

    # Sample metadata
    leveler_concentration: float | None = None
    plating_time_hours: float | None = None
    sample_age_hours: float | None = None

    # Peak data
    peaks: list[PeakData] = field(default_factory=list)

    # Preprocessing and angle-quality diagnostics
    preprocessing_notes: list[str] = field(default_factory=list)
    background_applied: bool = False
    background_method: str = ""
    background_fraction_mean: float | None = None
    angle_offset_mean_deg: float | None = None
    angle_offset_rmse_deg: float | None = None
    angle_offset_max_abs_deg: float | None = None
    angle_validation_peaks: int = 0

    # Phase 04: Scherrer
    scherrer_results: list[ScherrerResult] = field(default_factory=list)
    average_size_nm: float | None = None

    # Phase 05: W-H
    wh_result: WHResult | None = None

    # Phase 06: Texture
    texture_result: TextureAnalysisResult | None = None

    # Phase 07: Defects
    stacking_fault: StackingFaultResult | None = None
    lattice_result: LatticeConstantResult | None = None
    annealing_state: AnnealingState = AnnealingState.UNKNOWN


# =============================================================================
# Main Pipeline
# =============================================================================


class XRDAnalysisPipeline:
    """Complete xrd_analysis analysis pipeline.

    Integrates all Phase 04-07 analysis modules into a unified workflow.
    """

    def __init__(self, config: AnalysisConfig | None = None):
        """Initialize pipeline with configuration."""
        self.config = config or AnalysisConfig()

        # Double-Correction Guard: when doublet fitting handles Kα2 explicitly,
        # pre-stripping must be disabled to avoid subtracting Kα2 twice.
        if self.config.use_doublet_fitting and self.config.enable_kalpha_strip:
            logger.info(
                "Double-Correction Guard: disabling Kα2 pre-stripping "
                "because doublet fitting is active"
            )
            self.config.enable_kalpha_strip = False

        # Initialize analyzers
        self.scherrer = ScherrerCalculator(
            wavelength=self.config.wavelength,
            use_cubic_habit=self.config.use_cubic_habit,
            caglioti_params=(
                self.config.caglioti_u,
                self.config.caglioti_v,
                self.config.caglioti_w,
            ),
        )

        # Initialize Williamson-Hall analyzer
        self.wh = WilliamsonHallAnalyzer()

        # Initialize Texture analyzer
        self.texture = TextureAnalyzer()
        self.sf_analyzer = StackingFaultAnalyzer()
        self.lattice = LatticeMonitor()
        self.preprocessing = PreprocessingPipeline(
            window_size=self.config.smoothing_window,
            poly_order=self.config.smoothing_poly_order,
            background_method=self.config.background_method,
            background_degree=self.config.background_degree,
            enable_smoothing=self.config.enable_smoothing,
            enable_background=self.config.enable_background,
            enable_kalpha_strip=self.config.enable_kalpha_strip,
        )

    def _load_and_validate_data(self, filepath: str):
        """Load XRD data and validate. Returns (two_theta, intensity, error_message)."""
        try:
            two_theta, intensity = load_bruker_txt(filepath)
        except (OSError, ValueError, UnicodeDecodeError) as e:
            return None, None, f"Error loading file: {e}"

        if len(two_theta) == 0:
            return None, None, "No data found in file"

        return two_theta, intensity, ""

    def _find_peaks_from_data(self, two_theta, intensity):
        """Find peaks in XRD data using expected peak positions."""
        peaks = []
        for hkl, expected_pos in self.config.EXPECTED_PEAKS.items():
            peak = find_peak_in_range(
                two_theta,
                intensity,
                expected_pos,
                window=self.config.peak_window,
                use_doublet_fitting=self.config.use_doublet_fitting,
                doublet_max_iterations=self.config.doublet_max_iterations,
                min_fit_r_squared=self.config.min_fit_r_squared,
            )
            if peak:
                peak.hkl = hkl
                peaks.append(peak)
        return peaks

    def _run_preprocessing(self, two_theta: np.ndarray, intensity: np.ndarray):
        """Run configured preprocessing and return processed arrays with metadata."""
        pre_result = self.preprocessing.run(two_theta, intensity)
        notes = [
            f"{step.name}: {'ON' if step.applied else 'OFF'} ({step.notes})"
            for step in pre_result.steps
        ]

        background_fraction_mean = None
        if pre_result.background is not None:
            denom = float(np.mean(np.maximum(pre_result.raw_intensity, 1.0)))
            if denom > 0:
                background_fraction_mean = float(np.mean(pre_result.background) / denom)

        return (
            pre_result.two_theta,
            pre_result.intensity,
            notes,
            pre_result.background is not None,
            background_fraction_mean,
        )

    def _validate_peak_angles(self, peaks):
        """Validate measured peak positions against expected references."""
        offsets = []
        for peak in peaks:
            expected = self.config.EXPECTED_PEAKS.get(peak.hkl)
            if expected is None:
                continue
            offsets.append(float(peak.two_theta - expected))

        if not offsets:
            return None, None, None, 0

        offset_arr = np.array(offsets, dtype=float)
        mean_offset = float(np.mean(offset_arr))
        rmse = float(np.sqrt(np.mean(offset_arr**2)))
        max_abs = float(np.max(np.abs(offset_arr)))
        return mean_offset, rmse, max_abs, len(offsets)

    def _run_scherrer_analysis(self, peaks):
        """Run Scherrer analysis on all peaks. Returns (results, average_size)."""
        scherrer_results = []
        for peak in peaks:
            result = self.scherrer.calculate(
                two_theta=peak.two_theta,
                fwhm_observed=peak.fwhm,
                # Let ScherrerCalculator compute angle-dependent instrumental broadening
                # via Caglioti U, V, W when available.
                fwhm_instrumental=None,
                hkl=peak.hkl,
                eta_observed=peak.eta,
                eta_instrumental=0.0,  # Assume Gaussian instrument (Caglioti standard)
            )
            scherrer_results.append(result)

        # Calculate average
        valid_sizes = [
            r.size_nm
            for r in scherrer_results
            if r.validity_flag != ValidityFlag.UNRELIABLE
        ]
        avg_size = np.mean(valid_sizes) if valid_sizes else None
        return scherrer_results, avg_size

    def _prepare_wh_input(self, peaks, scherrer_results):
        """Prepare W-H input using instrument-corrected sample broadening.

        Methodology rationale:
        - W-H should use sample broadening (beta_sample), not raw observed FWHM.
        - We therefore reuse Scherrer pathway deconvolution output (fwhm_sample).
        """
        if not scherrer_results or len(peaks) != len(scherrer_results):
            return None, None, None

        two_theta_arr = np.array([p.two_theta for p in peaks], dtype=float)
        fwhm_sample_arr = np.array(
            [r.fwhm_sample for r in scherrer_results], dtype=float
        )
        hkl_list = [p.hkl for p in peaks]
        return two_theta_arr, fwhm_sample_arr, hkl_list

    def _run_defect_analysis(self, peaks):
        """Run defect and lattice analysis. Returns (stacking_fault, lattice)."""
        # Stacking fault
        peak_111 = next((p for p in peaks if p.hkl == (1, 1, 1)), None)
        peak_200 = next((p for p in peaks if p.hkl == (2, 0, 0)), None)
        sf_result = None

        if peak_111 and peak_200:
            sf_result = self.sf_analyzer.analyze(peak_111.two_theta, peak_200.two_theta)

        # Lattice constant (prefer high-angle peaks)
        high_angle_peak = next(
            (p for p in peaks if p.hkl in [(3, 1, 1), (2, 2, 0)]),
            peaks[-1] if peaks else None,
        )
        lattice_result = None
        if high_angle_peak:
            lattice_result = self.lattice.analyze_lattice(
                high_angle_peak.two_theta, high_angle_peak.hkl
            )

        return sf_result, lattice_result

    def analyze(
        self, filepath: str, sample_age_hours: float | None = None
    ) -> PipelineResult:
        """Run complete analysis on XRD data file.

        Args:
            filepath: Path to XRD data file
            sample_age_hours: Time since deposition (hours)

        Returns:
            PipelineResult with all analysis data

        """
        # Parse filename and initialize result
        file_info = parse_filename(filepath)
        effective_sample_age = (
            sample_age_hours
            if sample_age_hours is not None
            else file_info["time_hours"]
        )

        result = PipelineResult(
            filepath=filepath,
            sample_name=file_info["name"],
            leveler_concentration=file_info["concentration_ml"],
            plating_time_hours=None,
            sample_age_hours=effective_sample_age,
        )

        # Step 1: Load data
        two_theta, intensity, error = self._load_and_validate_data(filepath)
        if error:
            return result

        # Step 1.5: Preprocessing (smoothing/background/Kα2 strip)
        (
            two_theta_proc,
            intensity_proc,
            preprocessing_notes,
            background_applied,
            background_fraction_mean,
        ) = self._run_preprocessing(two_theta, intensity)
        result.preprocessing_notes = preprocessing_notes
        result.background_applied = background_applied
        result.background_method = (
            self.config.background_method if background_applied else "none"
        )
        result.background_fraction_mean = background_fraction_mean

        # Step 2: Find peaks
        result.peaks = self._find_peaks_from_data(two_theta_proc, intensity_proc)
        if len(result.peaks) < 2:
            return result

        # Step 2.5: Angle correctness check (measured vs expected positions)
        (
            result.angle_offset_mean_deg,
            result.angle_offset_rmse_deg,
            result.angle_offset_max_abs_deg,
            result.angle_validation_peaks,
        ) = self._validate_peak_angles(result.peaks)

        # Step 3: Scherrer analysis
        result.scherrer_results, result.average_size_nm = self._run_scherrer_analysis(
            result.peaks
        )

        # Step 4: W-H analysis (if enough peaks)
        if len(result.peaks) >= 3:
            two_theta_arr, fwhm_arr, hkl_list = self._prepare_wh_input(
                result.peaks, result.scherrer_results
            )
            if two_theta_arr is not None:
                result.wh_result = self.wh.analyze(two_theta_arr, fwhm_arr, hkl_list)

        # Step 5: Texture analysis
        # Use integrated area (not peak height) for proper comparison with JCPDS
        # standard intensities, which are relative integrated intensities.
        # Fall back to peak height if area is zero (e.g., simple fallback fitting).
        intensities = {
            p.hkl: (p.area if p.area > 0 else p.intensity) for p in result.peaks
        }
        result.texture_result = self.texture.analyze(intensities)

        # Step 6: Defect analysis
        result.stacking_fault, result.lattice_result = self._run_defect_analysis(
            result.peaks
        )

        # Step 7: Annealing state
        result.annealing_state, _ = determine_annealing_state(effective_sample_age)

        return result

    def process_file(
        self,
        filepath: str,
        output_dir: str,
        sample_age_hours: float | None = None,
    ) -> PipelineResult:
        """Analyze one file and ensure output directory exists."""
        result = self.analyze(filepath, sample_age_hours=sample_age_hours)

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        return result


# =============================================================================
# Convenience Functions
# =============================================================================


def run_full_analysis(
    filepath: str,
    sample_age_hours: float | None = None,
    config: AnalysisConfig | None = None,
) -> PipelineResult:
    """Run complete xrd_analysis analysis on a single file.

    Example:
        >>> result = run_full_analysis("data/raw/sample.txt")
        >>> print(result.sample_name)

    """
    pipeline = XRDAnalysisPipeline(config)
    return pipeline.analyze(filepath, sample_age_hours)


def batch_analyze(
    filepaths: list[str],
    sample_age_hours: float | None = None,
    config: AnalysisConfig | None = None,
) -> list[PipelineResult]:
    """Run analysis on multiple files."""
    pipeline = XRDAnalysisPipeline(config)
    results = []

    for fp in filepaths:
        result = pipeline.analyze(fp, sample_age_hours)
        results.append(result)

    return results

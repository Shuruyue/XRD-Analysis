"""
xrd_analysis Complete Analysis Pipeline 完整分析管道
==========================================

Unified pipeline integrating all analysis phases.
統一的分析流程，整合所有分析階段。

- Phase 02: Preprocessing 預處理
- Phase 03: Peak Fitting 峰値擬合
- Phase 04: Scherrer Size Calculation 晶粒尺寸計算
- Phase 05: Williamson-Hall Strain Analysis 應變分析
- Phase 06: Texture Analysis 織構分析
- Phase 07: Defect and Stress Diagnosis 缺陷與應力診斷
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict, Any
from pathlib import Path
import re

from xrd_analysis.core.constants import CU_KA1
from xrd_analysis.core.copper_crystal import CU_JCPDS_EXTENDED
from xrd_analysis.fitting.ka_doublet import DoubletFitter

from xrd_analysis.methods.scherrer import (
    ScherrerCalculator,
    ScherrerResult,
    ValidityFlag,
)
from xrd_analysis.methods.williamson_hall import (
    WilliamsonHallAnalyzer,
    WHResult,
)
from xrd_analysis.methods.texture import (
    TextureAnalyzer,
    TextureAnalysisResult,
)
from xrd_analysis.methods.defect_analysis import (
    StackingFaultAnalyzer,
    StackingFaultResult,
    LatticeMonitor,
    LatticeConstantResult,
    AnnealingState,
    determine_annealing_state,
)
from xrd_analysis.fitting.hkl_assignment import assign_hkl
from xrd_analysis.fitting.lm_optimizer import LMOptimizer
from xrd_analysis.fitting.pseudo_voigt import PseudoVoigt, PseudoVoigtParams
from xrd_analysis.preprocessing.pipeline import PreprocessingPipeline

# Import after to avoid circular dependency
from xrd_analysis.analysis.report_generator import (
    ComprehensiveResult,
    generate_comprehensive_report,
    generate_csv_summary,
)


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class AnalysisConfig:
    """
    Configuration for xrd_analysis analysis pipeline.
    xrd_analysis 分析管道配置。
    """
    
    # X-ray parameters / X-ray 參數
    # Default: Cu Kα1 from constants (Bearden 1967)
    wavelength: float = CU_KA1
    
    # Scherrer parameters / Scherrer 參數
    use_cubic_habit: bool = True
    
    # Peak detection / 峰值偵測
    peak_window: float = 2.0  # degrees around expected position
    min_intensity: float = 100  # minimum counts

    # Preprocessing / 預處理
    enable_smoothing: bool = True
    smoothing_window: int = 11
    smoothing_poly_order: int = 3
    enable_background: bool = True
    background_method: str = "chebyshev"
    background_degree: int = 5
    enable_kalpha_strip: bool = True
    
    # Instrumental broadening (Caglioti U, V, W) / 儀器展寬
    # FWHM_inst² = U·tan²θ + V·tanθ + W
    # 
    # Default values 默認值 (Empirical ranges for typical lab diffractometers):
    #   High-end lab:     W ≈ 0.001-0.002  (FWHM ≈ 0.032-0.045°)
    #   Standard lab:     W ≈ 0.002-0.005  (FWHM ≈ 0.045-0.071°)  ← Default
    #   Aging equipment:  W ≈ 0.008-0.015  (FWHM ≈ 0.089-0.122°)
    #
    # ⚠️ IMPORTANT 重要: These are placeholder values for initial analysis.
    # For publication-quality results, calibrate using LaB6 (NIST SRM 660c):
    #   >>> xrd-analysis calibrate data/lab6_standard.txt
    # 
    # Reference 參考: Caglioti et al. (1958). Nucl. Instrum. 3, 223-228.
    caglioti_u: float = 0.0      # Default: 0 (typically -0.005 to +0.005)
    caglioti_v: float = 0.0      # Default: 0 (typically -0.002 to +0.002)
    caglioti_w: float = 0.003    # Default: standard lab (FWHM_inst ≈ 0.055°)
    
    # Cu peak positions from JCPDS 04-0836 / 銅峰位從 JCPDS 標準
    # Dynamically generated from CU_JCPDS constants
    EXPECTED_PEAKS: Dict[Tuple[int, int, int], float] = field(
        default_factory=lambda: {
            hkl: data["two_theta"] for hkl, data in CU_JCPDS_EXTENDED.items()
        }
    )


# =============================================================================
# Result Container
# =============================================================================

@dataclass
class PeakData:
    """Single peak data."""
    hkl: Tuple[int, int, int]
    two_theta: float
    intensity: float
    fwhm: float
    area: float = 0.0
    eta: float = 0.5  # Pseudo-Voigt mixing parameter (0=Gaussian, 1=Lorentzian)


@dataclass
class PipelineResult:
    """Complete pipeline analysis result."""
    
    # Input
    filepath: str
    sample_name: str
    
    # Sample metadata
    leveler_concentration: Optional[float] = None
    plating_time_hours: Optional[float] = None
    sample_age_hours: Optional[float] = None
    
    # Peak data
    peaks: List[PeakData] = field(default_factory=list)

    # Preprocessing and angle-quality diagnostics
    preprocessing_notes: List[str] = field(default_factory=list)
    background_applied: bool = False
    background_method: str = ""
    background_fraction_mean: Optional[float] = None
    angle_offset_mean_deg: Optional[float] = None
    angle_offset_rmse_deg: Optional[float] = None
    angle_offset_max_abs_deg: Optional[float] = None
    angle_validation_peaks: int = 0
    
    # Phase 04: Scherrer
    scherrer_results: List[ScherrerResult] = field(default_factory=list)
    average_size_nm: Optional[float] = None
    
    # Phase 05: W-H
    wh_result: Optional[WHResult] = None
    
    # Phase 06: Texture
    texture_result: Optional[TextureAnalysisResult] = None
    
    # Phase 07: Defects
    stacking_fault: Optional[StackingFaultResult] = None
    lattice_result: Optional[LatticeConstantResult] = None
    annealing_state: AnnealingState = AnnealingState.UNKNOWN
    
    # Comprehensive
    comprehensive: Optional[ComprehensiveResult] = None
    report: str = ""


# =============================================================================
# Data Loader
# =============================================================================

def load_bruker_txt(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load Bruker TXT format XRD data.
    
    Returns:
        Tuple of (two_theta, intensity) arrays
    """
    two_theta = []
    intensity = []
    in_data_section = False
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            if '[Data]' in line:
                in_data_section = True
                continue
            
            if in_data_section and line:
                # Skip header row
                if 'Angle' in line or 'PSD' in line:
                    continue
                
                # Parse data
                parts = line.replace(',', ' ').split()
                if len(parts) >= 2:
                    try:
                        theta = float(parts[0])
                        counts = float(parts[1])
                        two_theta.append(theta)
                        intensity.append(counts)
                    except ValueError:
                        continue
    
    return np.array(two_theta), np.array(intensity)


def parse_filename(filepath: str) -> Dict[str, Any]:
    """
    Parse sample info from filename.
    
    Format: YYYYMMDD_Xml_Xh.txt or YYYYMMDD_Xml_Xh_Xmin.txt
    """
    name = Path(filepath).stem
    
    result = {
        'name': name,
        'concentration_ml': None,
        'time_hours': None,
    }
    
    # Extract concentration (e.g., "0ml", "4.5ml", "9ml", "18ml")
    conc_match = re.search(r'(\d+\.?\d*)ml', name)
    if conc_match:
        result['concentration_ml'] = float(conc_match.group(1))
    
    # Extract time (e.g., "2h", "0h_15min", "24h")
    time_hours = 0
    hour_match = re.search(r'(\d+)h', name)
    if hour_match:
        time_hours = int(hour_match.group(1))
    
    min_match = re.search(r'(\d+)min', name)
    if min_match:
        time_hours += int(min_match.group(1)) / 60
    
    result['time_hours'] = time_hours
    
    return result


# =============================================================================
# Peak Finding with Pseudo-Voigt Fitting
# =============================================================================

def _estimate_fwhm_simple(
    theta_range: np.ndarray,
    int_range: np.ndarray,
    idx_max: int,
    peak_int: float
) -> float:
    """Estimate FWHM using half-maximum method."""
    half_max = peak_int / 2
    
    left_idx = idx_max
    while left_idx > 0 and int_range[left_idx] > half_max:
        left_idx -= 1
    
    right_idx = idx_max
    while right_idx < len(int_range) - 1 and int_range[right_idx] > half_max:
        right_idx += 1
    
    return max(theta_range[right_idx] - theta_range[left_idx], 0.1)


def _fit_peak_pseudo_voigt(
    theta_range: np.ndarray,
    int_range: np.ndarray,
    peak_theta: float,
    peak_int: float,
    initial_fwhm: float
) -> Tuple[bool, float, float, float, float, float]:
    """Try Pseudo-Voigt fitting. Returns (success, theta, intensity, fwhm, eta, area)."""
    try:
        optimizer = LMOptimizer(max_iterations=500, tolerance=1e-6)
        
        initial_guess = PseudoVoigtParams(
            center=peak_theta,
            amplitude=peak_int,
            fwhm=initial_fwhm,
            eta=0.5
        )
        
        fit_result = optimizer.fit_single_peak(theta_range, int_range, initial_guess=initial_guess)
        
        if fit_result.success and fit_result.r_squared > 0.8:
            fitted_curve = PseudoVoigt.profile(
                theta_range, fit_result.params.center,
                fit_result.params.amplitude, fit_result.params.fwhm, fit_result.params.eta
            )
            try:
                area = np.trapezoid(fitted_curve, theta_range)
            except AttributeError:
                area = np.trapz(fitted_curve, theta_range)
            
            return (
                True,
                fit_result.params.center,
                fit_result.params.amplitude,
                fit_result.params.fwhm,
                fit_result.params.eta,
                area
            )
    except Exception:
        pass
    
    return False, peak_theta, peak_int, initial_fwhm, 0.5, 0.0


def _calculate_peak_area_simple(
    theta_range: np.ndarray,
    int_range: np.ndarray
) -> float:
    """Calculate peak area using trapezoidal integration."""
    try:
        return np.trapezoid(int_range, theta_range)
    except AttributeError:
        return np.trapz(int_range, theta_range)



def find_peak_in_range(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    center: float,
    window: float = 2.5,
    use_doublet_fitting: bool = True,
    doublet_max_iterations: int = 20000,
) -> Optional[PeakData]:
    """
    Find peak near expected position using Kα doublet fitting.
    
    Args:
        two_theta: 2θ array
        intensity: Intensity array
        center: Expected peak center position
        window: Search window (degrees)
        use_doublet_fitting: If True, use DoubletFitter (recommended)
        
    Returns:
        PeakData with fitted parameters (Kα₁ only), or None if no peak found
    """
    # Select range
    mask = (two_theta >= center - window) & (two_theta <= center + window)
    if not np.any(mask):
        return None
    
    theta_range = two_theta[mask]
    int_range = intensity[mask]
    
    # Find maximum
    idx_max = np.argmax(int_range)
    peak_theta = theta_range[idx_max]
    peak_int = int_range[idx_max]
    
    if peak_int < 50:
        return None
    
    # Estimate initial FWHM
    initial_fwhm = _estimate_fwhm_simple(theta_range, int_range, idx_max, peak_int)
    
    # Try Kα doublet fitting using unified function
    if use_doublet_fitting:
        try:
            from xrd_analysis.fitting.peak_fitter import fit_peak_with_diagnosis
            fit_result = fit_peak_with_diagnosis(
                two_theta,
                intensity,
                center,
                window=2.5,
                use_doublet=True,
                doublet_max_iterations=doublet_max_iterations,
            )
            
            if fit_result['success'] and fit_result.get('r_squared', 0) > 0.8:
                from xrd_analysis.fitting.pv_area import calculate_pv_area
                area = calculate_pv_area(
                     fit_result['amplitude'],
                     fit_result['fwhm'],
                     fit_result.get('eta', 0.5)
                )
                return PeakData(
                    hkl=(0, 0, 0),
                    two_theta=fit_result['center'],
                    intensity=fit_result['amplitude'],
                    fwhm=fit_result['fwhm'],
                    area=area,
                    eta=fit_result.get('eta', 0.5)
                )
        except Exception:
            pass
    
    # Fallback: simple Pseudo-Voigt fitting
    success, peak_theta, peak_int, fwhm, eta, area = _fit_peak_pseudo_voigt(
        theta_range, int_range, peak_theta, peak_int, initial_fwhm
    )
    if success:
        return PeakData(
            hkl=(0, 0, 0),
            two_theta=peak_theta,
            intensity=peak_int,
            fwhm=fwhm,
            area=area,
            eta=eta
        )
    
    # Final fallback: simple method
    fwhm = max(initial_fwhm, 0.05)
    area = _calculate_peak_area_simple(theta_range, int_range)
    
    return PeakData(
        hkl=(0, 0, 0),
        two_theta=peak_theta,
        intensity=peak_int,
        fwhm=fwhm,
        area=area,
        eta=0.5 # Default to intermediate
    )


# =============================================================================
# Main Pipeline
# =============================================================================

class XRDAnalysisPipeline:
    """
    Complete xrd_analysis analysis pipeline.
    
    Integrates all Phase 04-07 analysis modules into a unified workflow.
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        """Initialize pipeline with configuration."""
        self.config = config or AnalysisConfig()
        
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
        # 初始化 Williamson-Hall 分析器
        self.wh = WilliamsonHallAnalyzer()
        
        # Initialize Texture analyzer
        # 初始化紋理分析器
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
        except Exception as e:
            return None, None, f"Error loading file: {e}"
        
        if len(two_theta) == 0:
            return None, None, "No data found in file"
        
        return two_theta, intensity, ""
    
    def _find_peaks_from_data(self, two_theta, intensity):
        """Find peaks in XRD data using expected peak positions."""
        peaks = []
        for hkl, expected_pos in self.config.EXPECTED_PEAKS.items():
            peak = find_peak_in_range(
                two_theta, intensity, expected_pos, self.config.peak_window
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

        return pre_result.two_theta, pre_result.intensity, notes, pre_result.background is not None, background_fraction_mean

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
                eta_instrumental=0.0  # Assume Gaussian instrument (Caglioti standard)
            )
            scherrer_results.append(result)
        
        # Calculate average
        valid_sizes = [
            r.size_nm for r in scherrer_results
            if r.validity_flag != ValidityFlag.UNRELIABLE
        ]
        avg_size = np.mean(valid_sizes) if valid_sizes else None
        return scherrer_results, avg_size
    
    def _run_defect_analysis(self, peaks):
        """Run defect and lattice analysis. Returns (stacking_fault, lattice)."""
        # Stacking fault
        peak_111 = next((p for p in peaks if p.hkl == (1, 1, 1)), None)
        peak_200 = next((p for p in peaks if p.hkl == (2, 0, 0)), None)
        sf_result = None
        
        if peak_111 and peak_200:
            sf_result = self.sf_analyzer.analyze(
                peak_111.two_theta, peak_200.two_theta
            )
        
        # Lattice constant (prefer high-angle peaks)
        high_angle_peak = next(
            (p for p in peaks if p.hkl in [(3, 1, 1), (2, 2, 0)]),
            peaks[-1] if peaks else None
        )
        lattice_result = None
        if high_angle_peak:
            lattice_result = self.lattice.analyze_lattice(
                high_angle_peak.two_theta, high_angle_peak.hkl
            )
        
        return sf_result, lattice_result
    
    def analyze(
        self,
        filepath: str,
        sample_age_hours: Optional[float] = None
    ) -> PipelineResult:
        """
        Run complete analysis on XRD data file.
        
        Args:
            filepath: Path to XRD data file
            sample_age_hours: Time since deposition (hours)
            
        Returns:
            PipelineResult with all analysis data
        """
        # Parse filename and initialize result
        file_info = parse_filename(filepath)
        effective_sample_age = sample_age_hours if sample_age_hours is not None else file_info['time_hours']
        
        result = PipelineResult(
            filepath=filepath,
            sample_name=file_info['name'],
            leveler_concentration=file_info['concentration_ml'],
            plating_time_hours=None,
            sample_age_hours=effective_sample_age,
        )
        
        # Step 1: Load data
        two_theta, intensity, error = self._load_and_validate_data(filepath)
        if error:
            result.report = error
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
        result.background_method = self.config.background_method if background_applied else "none"
        result.background_fraction_mean = background_fraction_mean
        
        # Step 2: Find peaks
        result.peaks = self._find_peaks_from_data(two_theta_proc, intensity_proc)
        if len(result.peaks) < 2:
            result.report = f"Only {len(result.peaks)} peaks found"
            return result

        # Step 2.5: Angle correctness check (measured vs expected positions)
        (
            result.angle_offset_mean_deg,
            result.angle_offset_rmse_deg,
            result.angle_offset_max_abs_deg,
            result.angle_validation_peaks,
        ) = self._validate_peak_angles(result.peaks)
        
        # Step 3: Scherrer analysis
        result.scherrer_results, result.average_size_nm = self._run_scherrer_analysis(result.peaks)
        
        # Step 4: W-H analysis (if enough peaks)
        if len(result.peaks) >= 3:
            two_theta_arr = np.array([p.two_theta for p in result.peaks])
            fwhm_arr = np.array([p.fwhm for p in result.peaks])
            hkl_list = [p.hkl for p in result.peaks]
            result.wh_result = self.wh.analyze(two_theta_arr, fwhm_arr, hkl_list)
        
        # Step 5: Texture analysis
        intensities = {p.hkl: p.intensity for p in result.peaks}
        result.texture_result = self.texture.analyze(intensities)
        
        # Step 6: Defect analysis
        result.stacking_fault, result.lattice_result = self._run_defect_analysis(result.peaks)
        
        # Step 7: Annealing state and comprehensive report
        result.annealing_state, _ = determine_annealing_state(effective_sample_age)
        result.comprehensive = self._build_comprehensive(result)
        result.report = generate_comprehensive_report(result.comprehensive)
        
        return result
    
    def _build_comprehensive(self, result: PipelineResult) -> ComprehensiveResult:
        """Build ComprehensiveResult from pipeline result."""
        comp = ComprehensiveResult(
            sample_name=result.sample_name,
            sample_age_hours=result.sample_age_hours,
        )

        comp.preprocessing_summary = list(result.preprocessing_notes)
        comp.background_applied = result.background_applied
        comp.background_method = result.background_method
        comp.background_fraction_mean = result.background_fraction_mean
        comp.angle_offset_mean_deg = result.angle_offset_mean_deg
        comp.angle_offset_rmse_deg = result.angle_offset_rmse_deg
        comp.angle_offset_max_abs_deg = result.angle_offset_max_abs_deg
        comp.angle_validation_peaks = result.angle_validation_peaks
        if (
            result.angle_offset_max_abs_deg is not None
            and result.angle_offset_max_abs_deg > 0.10
        ):
            comp.warnings.append(
                f"Peak-angle max offset is {result.angle_offset_max_abs_deg:.3f} deg (>0.10 deg). "
                "Check specimen displacement/zero-shift."
            )
        
        # Scherrer
        if result.average_size_nm:
            comp.scherrer_size_nm = result.average_size_nm
            comp.scherrer_validity = "VALID"
        
        # W-H
        if result.wh_result:
            comp.wh_size_nm = result.wh_result.crystallite_size_nm
            comp.wh_strain = result.wh_result.microstrain
            comp.wh_r_squared = result.wh_result.r_squared
            comp.wh_quality = result.wh_result.quality_level.value
        
        # Texture (DATA ONLY, no diagnosis)
        if result.texture_result:
            comp.dominant_orientation = result.texture_result.dominant_hkl
            comp.dominant_tc = result.texture_result.dominant_tc
            comp.is_random_texture = result.texture_result.is_random
        
        # Defects
        if result.stacking_fault:
            comp.peak_separation_deg = result.stacking_fault.peak_separation_deg
            comp.stacking_fault_alpha = result.stacking_fault.alpha_percent
            comp.stacking_fault_severity = result.stacking_fault.severity.value
        
        if result.lattice_result:
            comp.lattice_constant = result.lattice_result.lattice_constant
            comp.lattice_status = result.lattice_result.status.value
        
        comp.annealing_state = result.annealing_state.value
        
        return comp

    def process_file(
        self,
        filepath: str,
        output_dir: str,
        sample_age_hours: Optional[float] = None,
    ) -> PipelineResult:
        """
        Analyze one file and save a text report to output directory.
        """
        result = self.analyze(filepath, sample_age_hours=sample_age_hours)

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        report_path = out_dir / f"{result.sample_name}_report.txt"
        report_content = result.report or "No report generated."
        report_path.write_text(report_content, encoding="utf-8")

        return result


# =============================================================================
# Convenience Functions
# =============================================================================

def run_full_analysis(
    filepath: str,
    sample_age_hours: Optional[float] = None,
    config: Optional[AnalysisConfig] = None
) -> PipelineResult:
    """
    Run complete xrd_analysis analysis on a single file.
    
    Example:
        >>> result = run_full_analysis("data/raw/sample.txt")
        >>> print(result.report)
    """
    pipeline = XRDAnalysisPipeline(config)
    return pipeline.analyze(filepath, sample_age_hours)


def batch_analyze(
    filepaths: List[str],
    sample_age_hours: Optional[float] = None,
    config: Optional[AnalysisConfig] = None
) -> List[PipelineResult]:
    """
    Run analysis on multiple files.
    """
    pipeline = XRDAnalysisPipeline(config)
    results = []
    
    for fp in filepaths:
        result = pipeline.analyze(fp, sample_age_hours)
        results.append(result)
    
    return results

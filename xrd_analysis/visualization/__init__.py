"""xrd_analysis Visualization Module
===========================

Provides structured visualization tools for XRD analysis results.
提供結構化的 XRD 分析結果視覺化工具。

Submodules:
    - style: Common styling configuration and color palettes
    - fwhm_plots: FWHM evolution and comparison plots
    - scherrer_plots: Crystallite size visualizations
    - wh_plots: Williamson-Hall analysis plots
    - texture_plots: Texture coefficient visualizations
    - fitting_plots: Peak fitting diagnosis plots
"""

from xrd_analysis.visualization.fitting_plots import (
    plot_doublet_comparison,
    plot_fit_residuals,
    plot_peak_fit,
)
from xrd_analysis.visualization.fwhm_plots import (
    plot_fwhm_by_peak,
    plot_fwhm_evolution,
    plot_fwhm_by_concentration,
)
from xrd_analysis.visualization.scherrer_plots import (
    plot_scherrer_evolution_by_peak,
    plot_scherrer_by_concentration,
)
from xrd_analysis.visualization.style import (
    XRD_ANALYSIS_STYLE,
    COLORBLIND_SAFE,
    apply_xrd_analysis_style,
    get_color_palette,
)
from xrd_analysis.visualization.texture_plots import (
    plot_tc_evolution,
    plot_texture_polar,
)
from xrd_analysis.visualization.wh_plots import (
    plot_wh_residuals,
    plot_williamson_hall,
)

__all__ = [
    # Style
    "XRD_ANALYSIS_STYLE",
    "COLORBLIND_SAFE",
    "apply_xrd_analysis_style",
    "get_color_palette",
    # FWHM
    "plot_fwhm_evolution",
    "plot_fwhm_by_peak",
    "plot_fwhm_by_concentration",
    # Scherrer
    "plot_scherrer_evolution_by_peak",
    "plot_scherrer_by_concentration",
    # Williamson-Hall
    "plot_williamson_hall",
    "plot_wh_residuals",
    # Texture
    "plot_texture_polar",
    "plot_tc_evolution",
    # Fitting
    "plot_peak_fit",
    "plot_doublet_comparison",
    "plot_fit_residuals",
]

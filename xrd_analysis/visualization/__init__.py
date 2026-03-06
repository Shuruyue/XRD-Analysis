"""xrd_analysis Visualization Module.
===========================

Provides structured visualization tools for XRD analysis results.

Submodules:
    - style: Common styling configuration and color palettes
    - fwhm_plots: FWHM evolution and comparison plots
    - scherrer_plots: Crystallite size visualizations
    - texture_plots: Texture coefficient visualizations
    - fitting_plots: Peak fitting diagnosis plots
"""

import os

import matplotlib

# Use a non-interactive backend by default in headless environments.
if not os.environ.get("MPLBACKEND"):
    matplotlib.use("Agg")

from xrd_analysis.visualization.fitting_plots import (
    plot_doublet_comparison,
    plot_fit_residuals,
    plot_peak_fit,
)
from xrd_analysis.visualization.fwhm_plots import (
    plot_fwhm_by_concentration,
    plot_fwhm_by_peak,
    plot_fwhm_evolution,
)
from xrd_analysis.visualization.scherrer_plots import (
    plot_scherrer_by_concentration,
    plot_scherrer_evolution_by_peak,
)
from xrd_analysis.visualization.style import (
    COLORBLIND_SAFE,
    XRD_ANALYSIS_STYLE,
    apply_xrd_analysis_style,
    get_color_palette,
)
from xrd_analysis.visualization.texture_plots import (
    plot_tc_evolution,
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
    # Texture
    "plot_tc_evolution",
    # Fitting
    "plot_peak_fit",
    "plot_doublet_comparison",
    "plot_fit_residuals",
]

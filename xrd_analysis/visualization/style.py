"""xrd_analysis Style Configuration.
==========================

Unified matplotlib styling and colorblind-safe color palettes.
"""

from typing import Any

import matplotlib.pyplot as plt

# =============================================================================
# XRD-Analysis Standard Style
# =============================================================================

XRD_ANALYSIS_STYLE: dict[str, Any] = {
    # Figure settings
    "figure.figsize": (10, 6),
    "figure.dpi": 100,
    "figure.facecolor": "white",
    # 'figure.edgecolor': 'black', # Some journals prefer this
    # Save settings: 300 DPI for common journal raster requirements
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "savefig.bbox": "tight",
    # Font settings - Strictly Times New Roman / Serif
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "Liberation Serif"],
    "font.size": 14,
    "text.usetex": False,  # Avoid external dependency unless requested
    "mathtext.fontset": "stix",  # LaTeX-like math font
    # Axes settings - Box style
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "axes.titleweight": "bold",
    "axes.grid": False,  # Academic standard: usually no grid or very subtle
    "axes.spines.top": True,  # Box style
    "axes.spines.right": True,
    "axes.linewidth": 1.5,
    "axes.edgecolor": "black",
    # Grid settings (if enabled manually)
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "grid.linewidth": 0.5,
    # Legend settings
    "legend.fontsize": 12,
    "legend.frameon": True,
    "legend.framealpha": 1.0,  # Opaque
    "legend.edgecolor": "black",
    "legend.fancybox": False,  # Square corners
    # Tick settings - Inward ticks are standard in physics/materials
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,  # Ticks on top
    "ytick.right": True,  # Ticks on right
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    # Line settings
    "lines.linewidth": 2.0,
    "lines.markersize": 8,
    "lines.markeredgewidth": 1.0,
}

# =============================================================================
# Colorblind-Safe Color Palette (Wong, 2011)
# =============================================================================

COLORBLIND_SAFE: list[str] = [
    "#0077BB",  # Blue - primary
    "#EE7733",  # Orange - secondary
    "#009988",  # Teal
    "#CC3311",  # Red
    "#EE3377",  # Magenta
    "#33BBEE",  # Cyan
    "#BBBBBB",  # Gray
    "#000000",  # Black
]

# Semantic color mapping for XRD peaks
PEAK_COLORS: dict[str, str] = {
    "(111)": "#0077BB",  # Blue
    "(200)": "#EE7733",  # Orange
    "(220)": "#009988",  # Teal
    "(311)": "#CC3311",  # Red
    "(222)": "#EE3377",  # Magenta
}

# =============================================================================
# Style Application Functions
# =============================================================================


def apply_xrd_analysis_style() -> None:
    """Apply xrd_analysis style to all subsequent matplotlib plots.

    Example:
        >>> from xrd_analysis.visualization.style import apply_xrd_analysis_style
        >>> apply_xrd_analysis_style()
        >>> # All subsequent plots will use xrd_analysis style

    """
    plt.rcParams.update(XRD_ANALYSIS_STYLE)


def get_color_palette(n_colors: int = 6) -> list[str]:
    """Get a colorblind-safe color palette.

    Args:
        n_colors: Number of colors needed (max 8).

    Returns:
        List of hex color strings.

    Example:
        >>> colors = get_color_palette(3)
        >>> print(colors)
        ['#0077BB', '#EE7733', '#009988']

    """
    return COLORBLIND_SAFE[: min(n_colors, len(COLORBLIND_SAFE))]


def get_peak_color(hkl: str) -> str:
    """Get color for a specific (hkl) peak.

    Args:
        hkl: Peak identifier string, e.g. "(111)".

    Returns:
        Hex color string. Returns gray if not found.

    """
    return PEAK_COLORS.get(hkl, "#BBBBBB")


def create_figure(figsize: tuple = None, apply_style: bool = True) -> tuple:
    """Create a figure with xrd_analysis styling.

    Args:
        figsize: Optional figure size (width, height) in inches.
        apply_style: Whether to apply xrd_analysis style.

    Returns:
        Tuple of (fig, ax) matplotlib objects.

    """
    if apply_style:
        apply_xrd_analysis_style()

    if figsize is None:
        figsize = XRD_ANALYSIS_STYLE["figure.figsize"]

    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def save_figure(
    fig, filepath: str, dpi: int = 300, format: str = None, transparent: bool = False
) -> None:
    """Save figure with standardized settings.

    Args:
        fig: Matplotlib figure object.
        filepath: Output file path.
        dpi: Resolution (dots per inch).
        format: Output format (png, svg, pdf). Auto-detected from filepath if None.
        transparent: Whether to use transparent background.

    """
    fig.savefig(
        filepath,
        dpi=dpi,
        format=format,
        transparent=transparent,
        bbox_inches="tight",
        facecolor="white" if not transparent else "none",
        edgecolor="none",
    )

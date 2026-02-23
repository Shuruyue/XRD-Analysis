"""
Defect Analysis Plots Module
============================

Visualization for stacking fault analysis.
"""

from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt

from xrd_analysis.visualization.style import apply_xrd_analysis_style, save_figure


def plot_stacking_fault_evolution(
    data: List[Dict[str, Any]],
    output_path: Optional[str] = None,
    dpi: int = 600,
    show: bool = True,
) -> plt.Figure:
    """Plot stacking fault probability evolution."""
    apply_xrd_analysis_style()

    high_quality_data = [s for s in data if s.get("high_quality", True)]
    concentrations = sorted(set(s.get("concentration", 0) for s in high_quality_data))

    n_conc = len(concentrations)
    n_cols = 2
    n_rows = (n_conc + 1) // 2

    figsize = (10, 10)
    if n_rows > 2:
        figsize = (10, 10 * (n_rows / 2))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=figsize, squeeze=False, constrained_layout=True
    )
    axes = axes.flatten()

    line_color = "#CC79A7"

    for idx, conc in enumerate(concentrations):
        ax = axes[idx]
        conc_samples = [s for s in high_quality_data if s.get("concentration") == conc]
        times = sorted(list(set(s["time"] for s in conc_samples)))

        x_vals = []
        y_vals = []

        for t in times:
            sample = next((s for s in conc_samples if s["time"] == t), None)
            if sample and "stacking_fault" in sample:
                sf_res = sample["stacking_fault"]
                if sf_res and sf_res.get("alpha_percent") is not None:
                    x_vals.append(t)
                    y_vals.append(sf_res.get("alpha_percent"))

        if x_vals:
            ax.plot(
                x_vals,
                y_vals,
                marker="D",
                color=line_color,
                linestyle="-",
                linewidth=2,
                markersize=7,
                label="Stacking Fault α% ((111)-(200) Shift)",
                alpha=0.9,
            )

        ax.set_title(f"{conc} mL/1.5L", fontsize=12, fontweight="bold")
        ax.set_xlabel("Annealing Time (hours)", fontsize=12)
        ax.set_ylabel("Stacking Fault Probability α (%)", fontsize=12)

        current_max = max(y_vals) if y_vals else 0.5
        y_limit = max(1.0, current_max * 1.2)
        ax.set_ylim(0, y_limit)

        ax.grid(True, alpha=0.3, linestyle=":")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_box_aspect(1)

    for idx in range(len(concentrations), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        "Stacking Fault Probability Evolution (Warren Method)",
        fontsize=16,
        fontweight="bold",
    )

    if output_path:
        save_figure(fig, output_path, dpi=dpi)

    if show:
        plt.show()

    return fig

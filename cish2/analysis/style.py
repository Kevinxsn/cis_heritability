"""Shared figure style.

Palette values are taken unchanged from the validated reference palette:
categorical slots 1-3 (blue / orange / aqua), which clear the all-pairs CVD and
normal-vision gates in light mode, and the blue sequential ramp for ordinal
encodings (starting no lighter than step 250 so the lightest step still clears
2:1 on a light surface).  Nothing is re-stepped here, so nothing needs
re-validating.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --- categorical: identity, fixed order, never cycled -----------------------
SERIES: Dict[str, str] = {
    "GCTA": "#2a78d6",   # slot 1, blue
    "SuSiE": "#eb6834",  # slot 2, orange
    "HEELS": "#1baf7a",  # slot 3, aqua
    "third": "#1baf7a",  # legacy alias for slot 3
}
# Slots 1-3 of the reference palette validate all-pairs in both modes, so three
# estimators is exactly the number this colour scheme supports without folding
# anything into "Other".
ESTIMATOR_ORDER: List[str] = ["GCTA", "SuSiE", "HEELS"]

# --- sequential / ordinal: one hue, light -> dark ---------------------------
BLUE_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
ORANGE_RAMP = ["#f5a884", "#f08b5b", "#eb6834", "#c4501f", "#9b3d15", "#7a2f10"]

# --- diverging: two poles, neutral gray midpoint ---------------------------
DIVERGING = "RdBu_r"
SEQUENTIAL = "Blues"

# --- ink & surface ---------------------------------------------------------
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8880"
GRID = "#e6e5e1"
REFERENCE = "#8a8880"
SURFACE = "#ffffff"

RC = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK_SECONDARY,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.7,
    "grid.alpha": 1.0,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "figure.dpi": 120,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,   # embed real fonts, not paths - editable in Illustrator
    "ps.fonttype": 42,
}


def apply_style() -> None:
    plt.rcParams.update(RC)


def ordinal_colors(k: int, ramp: Sequence[str] = tuple(BLUE_RAMP)) -> List[str]:
    """``k`` evenly spaced steps of an ordinal ramp, light -> dark."""
    ramp = list(ramp)
    if k <= 1:
        return [ramp[len(ramp) // 2]]
    if k <= len(ramp):
        idx = [round(i * (len(ramp) - 1) / (k - 1)) for i in range(k)]
        return [ramp[i] for i in idx]
    cmap = mpl.colors.LinearSegmentedColormap.from_list("ramp", ramp)
    return [mpl.colors.to_hex(cmap(i / (k - 1))) for i in range(k)]


def panel_label(ax, text: str) -> None:
    """Bold panel letter in the top-left corner, outside the data area."""
    ax.text(-0.16, 1.06, text, transform=ax.transAxes, fontsize=11,
            fontweight="bold", color=INK, ha="left", va="bottom")


def finish(ax, xlabel: str = "", ylabel: str = "", title: str = "") -> None:
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, color=INK, pad=8)
    ax.set_axisbelow(True)


def save(fig, out_dir: str, name: str, formats: Sequence[str] = ("png", "pdf")) -> List[str]:
    """Write one figure in every requested format; PDF is the vector master."""
    import os

    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ext in formats:
        path = os.path.join(out_dir, f"{name}.{ext}")
        fig.savefig(path)
        paths.append(path)
    plt.close(fig)
    return paths

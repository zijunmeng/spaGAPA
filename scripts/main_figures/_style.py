"""Shared style for all spaGAPA BIB main figures.

Publication conventions:
- 300 DPI, colorblind-friendly palette (Okabe-Ito derived)
- Labeled axes with units, clear legends, bold panel letters
- DejaVu Sans (universally available)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import os

# ---- Colorblind-friendly palette (Okabe-Ito) ----
# Primary triad requested by reviewer + supporting colors
BLUE   = "#0072B2"   # primary
ORANGE = "#D55E00"   # primary
GREEN  = "#009E73"   # primary
SKYBLU = "#56B4E9"
YELLOW = "#F0E442"
RED    = "#CC79A7"   # vermillion-pink
GREY   = "#999999"
BLACK  = "#000000"

# 5-method benchmark palette (consistent across figures)
METHOD_COLORS = {
    "spaGAPA-GP":  BLUE,
    "mean":        GREY,
    "spatial-KNN": GREEN,
    "stAPAminer":  ORANGE,
    "spvAPA":      RED,
}
METHOD_ORDER = ["spaGAPA-GP", "mean", "spatial-KNN", "stAPAminer", "spvAPA"]
METHOD_LABELS = {
    "spaGAPA-GP":  "spaGAPA-GP",
    "mean":        "Mean",
    "spatial-KNN": "Spatial-KNN",
    "stAPAminer":  "stAPAminer",
    "spvAPA":      "spvAPA",
}

# Noise method palette (Fig 4)
NOISE_METHOD_COLORS = {
    "A_constant":      GREY,
    "B_local_gene":    BLUE,
    "C_spatial_spot":  SKYBLU,
    "D_residual_spot": ORANGE,
}
NOISE_METHOD_LABELS = {
    "A_constant":      "A: Constant",
    "B_local_gene":    "B: Local-gene",
    "C_spatial_spot":  "C: Spatial-spot",
    "D_residual_spot": "D: Residual-spot",
}

OUT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/main_figures"

def setup_rc():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelsize": 9,
        "axes.labelweight": "bold",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "savefig.dpi": 300,
        "figure.dpi": 100,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

def panel_label(ax, letter, x=-0.14, y=1.06, fontsize=13):
    """Bold uppercase panel letter in axes-fraction coords."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved -> {path}")
    return path

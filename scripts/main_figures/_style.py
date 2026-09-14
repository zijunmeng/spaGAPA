"""Shared style for all spaGAPA BIB main figures.

Publication conventions:
- Vector PDF output (text/axes/lines stay vector; rasterize only dense
  spatial-scatter layers with ax.set_rasterized(True)) + 300-DPI PNG preview.
- Colorblind-friendly Okabe-Ito palette; DejaVu Sans throughout.
- Figures sized for BIB double-column print (~178 mm / 7.0 in wide) so all
  text remains >= 7-8 pt at final scale. Panels show data only; explanatory
  prose belongs in the figure caption, not inside the panel.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import glob as _glob

# Register system Arial fonts (msfonts not in matplotlib default search path)
for _fp in _glob.glob("/usr/share/fonts/msfonts/ARIAL*.TTF"):
    try:
        fm.fontManager.addfont(_fp)
    except Exception:
        pass
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

# NC/NAR page geometry (inches). 183 mm = 7.2 in (Nature double-column, NAR full width).
PAGE_WIDTH_IN = 7.2    # double-column width (~183 mm)
SINGLE_COL_IN = 3.5    # single column (~89 mm)

def setup_rc():
    """NC/NAR publication style: Arial, 183mm double-column width."""
    plt.rcParams.update({
        # NC/NAR mandate: Arial throughout (incl. mathtext)
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
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
        "pdf.fonttype": 42,   # TrueType embed (journal requirement)
        "ps.fonttype": 42,
    })

def panel_label(ax, letter, x=-0.14, y=1.06, fontsize=13):
    """Bold uppercase panel letter in axes-fraction coords."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    """Save figure as vector PDF (publication) + 300-DPI PNG (preview) to OUT.

    `name` may be given with or without extension; both .pdf and .png are
    written with the same basename. For dense spatial-scatter axes, call
    ax.set_rasterized(True) in the figure script before saving so the PDF
    stays small while text remains vector.
    """
    base, _ = os.path.splitext(name)
    pdf_path = os.path.join(OUT, base + ".pdf")
    png_path = os.path.join(OUT, base + ".png")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")           # vector
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")  # preview
    plt.close(fig)
    print(f"  saved -> {pdf_path} (+png)")
    return pdf_path

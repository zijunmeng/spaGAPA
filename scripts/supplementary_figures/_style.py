"""Shared style for spaGAPA BIB supplementary figures.

Mirrors main_figures/_style.py (Okabe-Ito palette, 300 DPI, DejaVu Sans).
Imported by all supp_fig*.py scripts for visual consistency.
"""
import sys, os, importlib.util
# Load the main_figures style module by file path to avoid name shadowing
# (this dir also has a _style.py which would otherwise be imported in a loop).
_MAIN_STYLE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "main_figures", "_style.py"))
_spec = importlib.util.spec_from_file_location("_main_style", _MAIN_STYLE)
_m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_m)
setup_rc = _m.setup_rc
panel_label = _m.panel_label
save = _m.save
BLUE = _m.BLUE; ORANGE = _m.ORANGE; GREEN = _m.GREEN; SKYBLU = _m.SKYBLU
YELLOW = _m.YELLOW; RED = _m.RED; GREY = _m.GREY; BLACK = _m.BLACK
METHOD_COLORS = _m.METHOD_COLORS
METHOD_ORDER = _m.METHOD_ORDER
METHOD_LABELS = _m.METHOD_LABELS
NOISE_METHOD_COLORS = _m.NOISE_METHOD_COLORS
NOISE_METHOD_LABELS = _m.NOISE_METHOD_LABELS

# Output dir for supplementary figures (PNGs write to pipeline_output, not the
# tracked scripts dir — mirrors main_figures/_style.py).
OUT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/supplementary_figures"

# Tissue color map (for dataset-level figures S1/S6/S9).
TISSUE_COLORS = {
    "kidney":       BLUE,
    "brain":        ORANGE,
    "colon":        GREEN,
    "liver":        SKYBLU,
    "glioma":       RED,
    "mouse_colon":  GREEN,
    "Visium-brain": ORANGE,
    "Visium-kidney": BLUE,
    "Visium-liver":  SKYBLU,
    "Mouse-colon":   GREEN,
}

# Method palette with Mean/StAPAminer/spvAPA already in METHOD_COLORS.
# Add additional categorical colors for non-method groupings.
EXTRA_COLORS = [YELLOW, "#CC79A7", "#0072B2", "#009E73", "#D55E00", "#56B4E9", GREY, BLACK]


def save_supp(fig, name):
    """Save a supplementary figure as vector PDF + 300-DPI PNG preview to OUT."""
    base, _ = os.path.splitext(name)
    pdf_path = os.path.join(OUT, base + ".pdf")
    png_path = os.path.join(OUT, base + ".png")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")           # vector
    fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")  # preview
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  saved -> {pdf_path} (+png)")
    return pdf_path

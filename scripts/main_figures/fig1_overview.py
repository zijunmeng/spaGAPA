"""Figure 1: spaGAPA framework overview (4-panel conceptual schematic).

Caption draft (kept here per BIB convention; not rendered inside panels):
  spaGAPA addresses three gaps in spatial APA analysis (A): the APA-usage
  matrix is severely sparse (median 4.6% of gene x spot entries observed
  across 7 calibration samples; range 3.1-26.5%, Visium 3-5% / Stereo-seq
  ~27% - source: pipeline_output/spagapa_per_sample_summary.json), prior
  tools return point estimates without calibrated uncertainty, and
  neighbor-based comparators did not complete at >=42k locations under
  benchmark limits. (B) Input pipeline:
  spatial reads -> PAS peak counts -> APA usage matrix Y in [0,1]^(G x N)
  with observed mask M -> spaGAPA imputation. (C) Sparse Gaussian Process:
  posterior mean with 68/95% credible bands over M=150 inducing points,
  O(N * M^2) complexity scales to 100k spots. (D) Split-conformal
  calibration: train / calibrate / test -> quantile q_hat -> 90% prediction
  intervals with finite-sample marginal coverage under exchangeability.
  Dataset details (8 GSE, 2 species, 2 platforms) are reported in Table 1
  and Supplementary Fig. S1, not in this framework figure.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup_rc, panel_label, save, PAGE_WIDTH_IN, SINGLE_COL_IN, \
    BLUE, ORANGE, GREEN, BLACK

setup_rc()
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

# --- Load real sparsity stats from ground-truth JSON (no hardcoded numbers) ---
_SUMMARY = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/spagapa_per_sample_summary.json"
import json
with open(_SUMMARY) as fh:
    _samples = json.load(fh)
_obs = [s["observed_fraction"] for s in _samples]
_n_samples = len(_samples)
_med_obs = float(np.median(_obs))
_min_obs = float(min(_obs))
_max_obs = float(max(_obs))
# Visium samples are everything except the binned Stereo-seq run
_visium_obs = [s["observed_fraction"] for s in _samples
               if "binned" not in s["sample"]]
_visium_med = float(np.median(_visium_obs))

# Layout: 2x2 grid within PAGE_WIDTH. Two rows of equal height.
fig = plt.figure(figsize=(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.82))
# Margin / gap in figure-fraction coords
LM, RM, TM, BM, GX, GY = 0.045, 0.985, 0.955, 0.055, 0.090, 0.105
col_w = ((RM - LM) - GX) / 2.0
row_h = ((TM - BM) - GY) / 2.0
axA = fig.add_axes([LM,            TM - row_h, col_w, row_h])  # top-left
axB = fig.add_axes([LM + col_w+GX, TM - row_h, col_w, row_h])  # top-right
axC = fig.add_axes([LM,            BM,         col_w, row_h])  # bot-left
axD = fig.add_axes([LM + col_w+GX, BM,         col_w, row_h])  # bot-right


def nice_box(ax, xy, w, h, text, fc, ec, fontsize=8, textcolor="white",
             bold=True):
    x, y = xy
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.008,rounding_size=0.02",
                         fc=fc, ec=ec, lw=1.2, mutation_aspect=ax.get_data_ratio())
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, color=textcolor,
            fontweight="bold" if bold else "normal", wrap=True)


def arrow(ax, p0, p1, color=BLACK, lw=1.4, style="->", mut=12):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=mut,
                        color=color, lw=lw, shrinkA=2, shrinkB=2)
    ax.add_patch(a)


# =====================================================================
# Panel A: Three gaps
# =====================================================================
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.set_title("Three gaps in spatial APA analysis", loc="left", fontsize=10)

gaps = [
    ("Gap 1", "Sparse APA-usage matrix",
     "Only ~%.0f%% of gene x spot entries\n"
     "observed (Visium median, %d samples)"
     % (_visium_med * 100, _n_samples),
     ORANGE),
    ("Gap 2", "No calibrated intervals",
     "Existing tools give point\n"
     "estimates - no coverage guarantee", BLUE),
    ("Gap 3", "Cannot scale",
     "Neighbor-based comparators did\n"
     "not complete at $\\geq$42k locations\n"
     "under benchmark limits", GREEN),
]
y0 = 0.82
step = 0.285
for i, (tag, title, body, col) in enumerate(gaps):
    yy = y0 - i * step
    # tag chip
    nice_box(axA, (0.03, yy - 0.135), 0.18, 0.17, tag, col, col, fontsize=8)
    # card
    card = FancyBboxPatch((0.235, yy - 0.135), 0.735, 0.17,
                          boxstyle="round,pad=0.008,rounding_size=0.015",
                          fc="#F5F5F5", ec=col, lw=1.3,
                          mutation_aspect=axA.get_data_ratio())
    axA.add_patch(card)
    axA.text(0.265, yy - 0.015, title, fontsize=8.5, fontweight="bold",
             color=col, va="top")
    axA.text(0.265, yy - 0.055, body, fontsize=7.2, color="#333333", va="top")
panel_label(axA, "A")

# =====================================================================
# Panel B: Input pipeline  (spatial reads -> PAS -> APA matrix -> spaGAPA)
# =====================================================================
axB.set_xlim(0, 1); axB.set_ylim(0, 1); axB.axis("off")
axB.set_title("Input pipeline", loc="left", fontsize=10)

# Top row: 4 stages left to right (narrower boxes -> longer visible arrows)
nice_box(axB, (0.025, 0.74), 0.175, 0.18, "Spatial\nreads\n(BAM)",
         "#444444", "#222222", 7.5)
nice_box(axB, (0.275, 0.74), 0.175, 0.18, "PAS peak\ncounts",
         "#666666", "#333333", 7.5)
nice_box(axB, (0.525, 0.74), 0.175, 0.18, "APA usage\nmatrix\n[gene x spot]",
         ORANGE, "#A04500", 7.0)
nice_box(axB, (0.775, 0.74), 0.20, 0.18, "spaGAPA", BLUE, "#004F7A", 8.5)
arrow(axB, (0.205, 0.83), (0.270, 0.83), lw=1.8, mut=14)
arrow(axB, (0.455, 0.83), (0.520, 0.83), lw=1.8, mut=14)
arrow(axB, (0.705, 0.83), (0.770, 0.83), lw=1.8, mut=14)

# Sparsity visualization: a gene x spot grid with most cells empty
axB.text(0.035, 0.63, "Observed entries (Visium ~%.0f%%):"
         % (_visium_med * 100), fontsize=7, fontweight="bold", color="#333")
np.random.seed(1)
gx, gy = 0.035, 0.40
cw, ch = 0.0305, 0.040
nrows, ncols = 6, 28
mask = np.random.rand(nrows, ncols) < _visium_med
for r in range(nrows):
    for c in range(ncols):
        cell = Rectangle((gx + c * cw, gy + (nrows - 1 - r) * ch),
                         cw * 0.90, ch * 0.90,
                         fc=(ORANGE if mask[r, c] else "#E8E8E8"),
                         ec="white", lw=0.3)
        axB.add_patch(cell)
axB.text(0.035, 0.32, "orange = observed APA usage", fontsize=6.5,
         color=ORANGE, style="italic")

# Goal statement (compact, no overlap with grid)
axB.text(0.035, 0.20, "Input:  Y $\\in$ [0,1]$^{G\\times N}$,  observed mask M",
         fontsize=7.5, color="#222", family="monospace")
axB.text(0.035, 0.105, "Goal:  impute Y$_{ij}$  +  uncertainty  for all (i, j)",
         fontsize=7.5, color=BLUE, family="monospace", fontweight="bold")
panel_label(axB, "B")

# =====================================================================
# Panel C: Sparse GP schematic (posterior mean + 68/95% CI, M inducing)
# =====================================================================
axC.set_xlim(0, 1); axC.set_ylim(0, 1); axC.axis("off")
axC.set_title("Sparse Gaussian Process", loc="left", fontsize=10)

# Inset data axis (in figure fraction coords, anchored to panel C region)
# Compute panel C bbox in figure fraction to place the inset
_cx0 = LM
_cy0 = BM
ax_gp = fig.add_axes([_cx0 + 0.055, _cy0 + 0.115, col_w - 0.115, row_h - 0.27])
x = np.linspace(0, 1, 200)
mean = 0.5 + 0.3 * np.sin(3 * np.pi * x) * np.exp(-1.2 * x)
sig = 0.05 + 0.12 * x
ax_gp.fill_between(x, mean - 1.96 * sig, mean + 1.96 * sig, color=BLUE,
                   alpha=0.15, label="95% CI")
ax_gp.fill_between(x, mean - sig, mean + sig, color=BLUE, alpha=0.30,
                   label="68% CI")
ax_gp.plot(x, mean, color=BLUE, lw=2, label="posterior mean")
xo = np.array([0.08, 0.22, 0.41, 0.63, 0.88])
yo = 0.5 + 0.3 * np.sin(3 * np.pi * xo) * np.exp(-1.2 * xo) \
    + np.random.RandomState(2).normal(0, 0.04, len(xo))
ax_gp.scatter(xo, yo, color=ORANGE, s=22, zorder=5, edgecolor="white",
              lw=0.5, label="observed")
zi = np.array([0.12, 0.33, 0.55, 0.78])
ax_gp.scatter(zi, 0.5 + 0.3 * np.sin(3 * np.pi * zi) * np.exp(-1.2 * zi),
              marker="^", s=55, color=GREEN, zorder=6, edgecolor="white",
              lw=0.5, label="inducing (M=150)")
ax_gp.set_xlim(0, 1); ax_gp.set_ylim(0, 1)
ax_gp.set_xticks([]); ax_gp.set_yticks([])
ax_gp.set_xlabel("spatial coordinate", fontsize=7)
ax_gp.set_ylabel("APA usage", fontsize=7)
ax_gp.legend(fontsize=6, loc="upper right", handlelength=1.2,
             borderpad=0.3, ncol=2, columnspacing=1.0)
for s in ax_gp.spines.values():
    s.set_linewidth(0.7)
panel_label(axC, "C")
# complexity annotation under the inset
axC.text(0.5, 0.045,
         "Complexity  O(N$\\cdot$M$^2$),  M $\\ll$ N   $\\rightarrow$  scales to 100k spots",
         ha="center", fontsize=7.3, color=GREEN, fontweight="bold")

# =====================================================================
# Panel D: Split-conformal flow (train | calibrate | test -> q_hat -> interval)
# =====================================================================
axD.set_xlim(0, 1); axD.set_ylim(0, 1); axD.axis("off")
axD.set_title("Split-conformal calibration", loc="left", fontsize=10)

# Three partitions (narrower boxes -> longer visible arrows)
nice_box(axD, (0.02, 0.80), 0.27, 0.13, "Train\n(fit GP)",
         BLUE, "#004F7A", 7.8)
nice_box(axD, (0.365, 0.80), 0.27, 0.13, "Calibrate\n(compute $\\hat{q}$)",
         ORANGE, "#A04500", 7.8)
nice_box(axD, (0.71, 0.80), 0.27, 0.13, "Test\n(evaluate)",
         GREEN, "#00734F", 7.8)
arrow(axD, (0.29, 0.865), (0.36, 0.865), lw=1.8, mut=14)
arrow(axD, (0.635, 0.865), (0.705, 0.865), lw=1.8, mut=14)

# Nonconformity score
card = FancyBboxPatch((0.025, 0.50), 0.95, 0.18,
                      boxstyle="round,pad=0.008,rounding_size=0.015",
                      fc="#F5F5F5", ec="#888", lw=1.0)
axD.add_patch(card)
axD.text(0.50, 0.625, "Nonconformity score:", fontsize=7.8, fontweight="bold",
         ha="center")
axD.text(0.50, 0.555, "s$_i$ = |y$_i$ - $\\hat{y}_i$|", fontsize=10,
         ha="center", family="monospace", color=BLACK)

# Quantile -> interval
card2 = FancyBboxPatch((0.025, 0.225), 0.95, 0.20,
                       boxstyle="round,pad=0.008,rounding_size=0.015",
                       fc="#F5F5F5", ec="#888", lw=1.0)
axD.add_patch(card2)
axD.text(0.50, 0.375, "90% prediction interval:", fontsize=7.8,
         fontweight="bold", ha="center")
axD.text(0.50, 0.295, "$\\hat{y}$ $\\pm$ $\\hat{q}_{0.90}$", fontsize=11,
         ha="center", family="monospace", color=ORANGE, fontweight="bold")

# Coverage guarantee (compact, single line, well inside panel)
axD.text(0.50, 0.135,
         "P(y $\\in$ interval) $\\geq$ 0.90  under exchangeability",
         fontsize=7.2, ha="center", color="#333", style="italic")
axD.text(0.50, 0.075, "finite-sample marginal coverage",
         fontsize=7.0, ha="center", color=GREEN)
panel_label(axD, "D")

save(fig, "fig1_framework_overview.png")
print("Figure 1 done")

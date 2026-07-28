"""Figure 1: Framework overview (6 panels, mostly schematic)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import *

setup_rc()
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import numpy as np

fig = plt.figure(figsize=(12, 7.5))
# Layout: A (top-left, wide), B (top-mid), C (top-right); D, E, F (bottom row)
# Use manual axes for schematic flexibility
axA = fig.add_axes([0.035, 0.55, 0.30, 0.40])
axB = fig.add_axes([0.36, 0.55, 0.30, 0.40])
axC = fig.add_axes([0.685, 0.55, 0.29, 0.40])
axD = fig.add_axes([0.035, 0.06, 0.30, 0.40])
axE = fig.add_axes([0.36, 0.06, 0.30, 0.40])
axF = fig.add_axes([0.685, 0.06, 0.29, 0.40])

def nice_box(ax, xy, w, h, text, fc, ec, fontsize=8, textcolor="white", bold=True):
    x, y = xy
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.02",
                         fc=fc, ec=ec, lw=1.2, mutation_aspect=ax.get_data_ratio())
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, color=textcolor, fontweight="bold" if bold else "normal",
            wrap=True)

def arrow(ax, p0, p1, color=BLACK, lw=1.4, style="->", mut=12):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=mut,
                        color=color, lw=lw, shrinkA=2, shrinkB=2)
    ax.add_patch(a)

# ---------------- Panel A: Three gaps ----------------
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.set_title("Three gaps in spatial APA analysis", loc="left", fontsize=10)
# Three stacked gap cards
gaps = [
    ("Gap 1", "Sparse PAS matrix", "APA usage is observed for\n<12% of gene x spot entries", ORANGE),
    ("Gap 2", "No intervals", "Existing tools give point\nestimates only - no calibrated\nuncertainty", BLUE),
    ("Gap 3", "Cannot scale", "O(N²) KNN methods fail at\n>15k spots (42k, 100k)", GREEN),
]
y0 = 0.70
for i, (tag, title, body, col) in enumerate(gaps):
    yy = y0 - i * 0.265
    # tag chip
    nice_box(axA, (0.02, yy-0.16), 0.17, 0.18, tag, col, col, fontsize=8)
    # card
    card = FancyBboxPatch((0.22, yy-0.16), 0.76, 0.18, boxstyle="round,pad=0.008,rounding_size=0.015",
                          fc="#F5F5F5", ec=col, lw=1.3, mutation_aspect=axA.get_data_ratio())
    axA.add_patch(card)
    axA.text(0.25, yy-0.02, title, fontsize=8.5, fontweight="bold", color=col, va="top")
    axA.text(0.25, yy-0.08, body, fontsize=7.2, color="#333333", va="top")
panel_label(axA, "A")

# ---------------- Panel B: Input pipeline ----------------
axB.set_xlim(0, 1); axB.set_ylim(0, 1); axB.axis("off")
axB.set_title("Input pipeline", loc="left", fontsize=10)
# Raw reads -> PAS counts -> APA matrix -> spaGAPA
nice_box(axB, (0.04, 0.74), 0.20, 0.16, "Spatial\nreads\n(BAM)", "#444444", "#222222", 7.5)
nice_box(axB, (0.30, 0.74), 0.22, 0.16, "PAS peak\ncounts", "#666666", "#333333", 7.5)
nice_box(axB, (0.58, 0.74), 0.20, 0.16, "APA usage\nmatrix\n[gene x spot]", ORANGE, "#A04500", 7.2)
nice_box(axB, (0.82, 0.74), 0.15, 0.16, "spa\nGAPA", BLUE, "#004F7A", 8)
arrow(axB, (0.24, 0.82), (0.30, 0.82))
arrow(axB, (0.52, 0.82), (0.58, 0.82))
arrow(axB, (0.78, 0.82), (0.82, 0.82))

# Show sparsity visually: small grid with many empty cells
axB.text(0.04, 0.62, "Observed entries (<12%):", fontsize=7, fontweight="bold", color="#333")
np.random.seed(1)
gx, gy = 0.04, 0.40
cw, ch = 0.022, 0.028
nrows, ncols = 7, 26
mask = np.random.rand(nrows, ncols) < 0.12
for r in range(nrows):
    for c in range(ncols):
        cell = Rectangle((gx + c*cw, gy + (nrows-1-r)*ch), cw*0.92, ch*0.92,
                         fc=(ORANGE if mask[r, c] else "#E8E8E8"),
                         ec="white", lw=0.3)
        axB.add_patch(cell)
axB.text(0.04, 0.34, "orange = observed APA usage", fontsize=6.5, color=ORANGE, style="italic")
axB.text(0.50, 0.20, "Input:  Y ∈ [0,1]^{G×N},  observed mask M", fontsize=7.5,
         ha="left", color="#222", family="monospace")
axB.text(0.50, 0.12, "Goal:  impute Yᵢⱼ  +  uncertainty  for all (i,j)", fontsize=7.5,
         ha="left", color=BLUE, family="monospace", fontweight="bold")
panel_label(axB, "B")

# ---------------- Panel C: Sparse GP schematic ----------------
axC.set_xlim(0, 1); axC.set_ylim(0, 1); axC.axis("off")
axC.set_title("Sparse Gaussian Process", loc="left", fontsize=10)
# Draw a GP posterior curve with inducing points
ax_gp = fig.add_axes([0.715, 0.575, 0.20, 0.32])
x = np.linspace(0, 1, 200)
# toy posterior: mean + bands
mean = 0.5 + 0.3*np.sin(3*np.pi*x) * np.exp(-1.2*x)
sig = 0.05 + 0.12*x
ax_gp.fill_between(x, mean-1.96*sig, mean+1.96*sig, color=BLUE, alpha=0.15, label="95% CI")
ax_gp.fill_between(x, mean-sig, mean+sig, color=BLUE, alpha=0.30, label="68% CI")
ax_gp.plot(x, mean, color=BLUE, lw=2, label="posterior mean")
# observed points
xo = np.array([0.08, 0.22, 0.41, 0.63, 0.88])
yo = 0.5 + 0.3*np.sin(3*np.pi*xo)*np.exp(-1.2*xo) + np.random.RandomState(2).normal(0,0.04,len(xo))
ax_gp.scatter(xo, yo, color=ORANGE, s=22, zorder=5, edgecolor="white", lw=0.5, label="observed")
# inducing points
zi = np.array([0.12, 0.33, 0.55, 0.78])
ax_gp.scatter(zi, 0.5 + 0.3*np.sin(3*np.pi*zi)*np.exp(-1.2*zi), marker="^",
              s=55, color=GREEN, zorder=6, edgecolor="white", lw=0.5, label="inducing (M=150)")
ax_gp.set_xlim(0,1); ax_gp.set_ylim(0,1)
ax_gp.set_xticks([]); ax_gp.set_yticks([])
ax_gp.set_xlabel("spatial coordinate", fontsize=7)
ax_gp.set_ylabel("APA usage", fontsize=7)
ax_gp.legend(fontsize=6, loc="upper right", handlelength=1.2, borderpad=0.3)
for s in ax_gp.spines.values(): s.set_linewidth(0.7)
panel_label(axC, "C")
# annotation under gp axis
axC.text(0.5, 0.04, "Complexity  O(N·M²),  M ≪ N   →  scales to 100k spots",
         ha="center", fontsize=7.3, color=GREEN, fontweight="bold")

# ---------------- Panel D: Split-conformal flow ----------------
axD.set_xlim(0, 1); axD.set_ylim(0, 1); axD.axis("off")
axD.set_title("Split-conformal calibration", loc="left", fontsize=10)
# Three partitions: train | calibrate | test
nice_box(axD, (0.04, 0.70), 0.26, 0.16, "Train\n(fit GP)", BLUE, "#004F7A", 7.5)
nice_box(axD, (0.34, 0.70), 0.26, 0.16, "Calibrate\n(compute qₜ)", ORANGE, "#A04500", 7.5)
nice_box(axD, (0.64, 0.70), 0.26, 0.16, "Test\n(evaluate)", GREEN, "#00734F", 7.5)
arrow(axD, (0.30, 0.78), (0.34, 0.78))
arrow(axD, (0.60, 0.78), (0.64, 0.78))
# nonconformity score formula
card = FancyBboxPatch((0.04, 0.44), 0.86, 0.18, boxstyle="round,pad=0.008,rounding_size=0.015",
                      fc="#F5F5F5", ec="#888", lw=1.0)
axD.add_patch(card)
axD.text(0.47, 0.56, "Nonconformity score:", fontsize=7.8, fontweight="bold", ha="center")
axD.text(0.47, 0.49, "sᵢ = |yᵢ − ŷᵢ|", fontsize=10, ha="center",
         family="monospace", color=BLACK)
# quantile -> interval
card2 = FancyBboxPatch((0.04, 0.18), 0.86, 0.20, boxstyle="round,pad=0.008,rounding_size=0.015",
                       fc="#F5F5F5", ec="#888", lw=1.0)
axD.add_patch(card2)
axD.text(0.47, 0.32, "90% interval:", fontsize=7.8, fontweight="bold", ha="center")
axD.text(0.47, 0.235, "ŷ ± q₀.9₀", fontsize=10, ha="center",
         family="monospace", color=ORANGE, fontweight="bold")
axD.text(0.47, 0.075, "Guarantees  P(y ∈ interval) ≥ 0.90  under exchangeability",
         fontsize=7.2, ha="center", color="#333", style="italic")
panel_label(axD, "D")

# ---------------- Panel E: Output ----------------
axE.set_xlim(0, 1); axE.set_ylim(0, 1); axE.axis("off")
axE.set_title("spaGAPA outputs", loc="left", fontsize=10)
# Four output mini-panels: mean, interval width, uncertainty map, domains
# mini-axes
positions = [(0.04, 0.40, 0.20), (0.28, 0.40, 0.20), (0.52, 0.40, 0.20), (0.76, 0.40, 0.20)]
titles = ["Posterior mean", "Interval width", "Uncertainty map", "Domains"]
cmaps = ["viridis", "magma", "plasma", "tab10"]
np.random.seed(7)
for (px, py, pw), tt, cm in zip(positions, titles, cmaps):
    axm = fig.add_axes([px, py, pw, pw])
    # synthesize a spatial blob
    g = np.linspace(-3, 3, 40)
    XX, YY = np.meshgrid(g, g)
    field = np.exp(-(XX**2 + YY**2)/4) - 0.5*np.exp(-((XX-2)**2 + (YY+1)**2)/3)
    field += np.random.RandomState(3).normal(0, 0.08, field.shape)
    if tt == "Domains":
        dom = (XX + YY > 0).astype(int) + (XX > 1).astype(int)
        axm.imshow(dom, cmap="tab10", interpolation="nearest")
    else:
        axm.imshow(field, cmap=cm, interpolation="bilinear")
    axm.set_xticks([]); axm.set_yticks([])
    for s in axm.spines.values(): s.set_linewidth(0.6)
    axE.text(px + pw/2, py - 0.045, tt, fontsize=7, ha="center", fontweight="bold")
axE.text(0.5, 0.32, "Calibrated uncertainty + downstream", ha="center", fontsize=7.5,
         style="italic", color="#444")
axE.text(0.5, 0.26, "tissue domains from imputed APA", ha="center", fontsize=7.5,
         style="italic", color="#444")
panel_label(axE, "E")

# ---------------- Panel F: Dataset atlas ----------------
axF.set_xlim(0, 1); axF.set_ylim(0, 1); axF.axis("off")
axF.set_title("Dataset atlas", loc="left", fontsize=10)
# Table-like grid: tissues x platforms
datasets = [
    ("GSE183456", "Human kidney", "Visium", "3010 spots"),
    ("GSE220442", "Human brain (VC)", "Visium", "4179 spots"),
    ("GSE237183", "Human glioma", "Visium", "4 sections"),
    ("GSE338525", "Human liver", "Visium", "4992 spots"),
    ("GSE169749", "Mouse colon", "Visium", "2715 spots"),
    ("GSE179572", "Mouse brain", "Visium", "4992 spots"),
    ("GSE263789", "Mouse AD brain", "Stereo-seq", "20.7M DNB"),
    ("MOB (spvAPA)", "Mouse OB", "ST", "260 spots"),
]
# render as 8 rows
ytop = 0.84
rh = 0.085
axF.text(0.03, 0.93, "8 tissues · 2 species · 2 platforms", fontsize=8, fontweight="bold", color=BLUE)
for i, (acc, tiss, plat, sc) in enumerate(datasets):
    yy = ytop - i*rh
    bg = "#F5F5F5" if i % 2 == 0 else "#FAFAFA"
    axF.add_patch(Rectangle((0.02, yy-rh*0.45), 0.96, rh*0.9, fc=bg, ec="none"))
    axF.text(0.04, yy, acc, fontsize=7, fontweight="bold", va="center", color="#222", family="monospace")
    axF.text(0.28, yy, tiss, fontsize=7, va="center")
    col = BLUE if plat == "Visium" else ORANGE
    axF.text(0.60, yy, plat, fontsize=7, va="center", color=col, fontweight="bold")
    axF.text(0.80, yy, sc, fontsize=6.8, va="center", color="#555")
axF.text(0.50, 0.07, "11 calibration samples · up to 100k spots",
         ha="center", fontsize=7.5, color=GREEN, fontweight="bold", style="italic")
panel_label(axF, "F")

save(fig, "fig1_framework_overview.png")
print("Figure 1 done")

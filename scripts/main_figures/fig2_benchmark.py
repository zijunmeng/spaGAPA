"""Figure 2: spaGAPA-GP vs mean / competitor benchmark (6 panels).

Honest framing (driven by the real benchmark numbers):
- Panels A-D: real per-dataset benchmark (transparent_comparison.csv). The
  per-gene *mean* is a very strong pointwise RMSE baseline and beats
  spaGAPA-GP on RMSE/Pearson/Spearman in aggregate. spaGAPA-GP wins clearly
  only on *spatial fidelity* (gradient recovery) -- it is the only method
  besides spatial-KNN that reconstructs a spatial gradient rather than
  collapsing to a per-gene constant.
- Panel E: a single representative high-spatial-signal gene (peak_94938,
  Moran's I = +0.50) from GSE183456, showing truth / mean / REAL spaGAPA-GP
  posterior / |GP-truth|. The per-spot posterior is produced by the actual
  `spagapa.SparseGPImputer` (NOT a proxy); see
  `fig2E_run_real_gp.py` and `_data/fig2E_gse183456_peak_94938_gp_posterior.csv`.
- Panel F: real per-gene stratification from 150 genes
  (`supplementary_figures/_cache/s4_stratification.csv`). The per-gene mean
  beats spaGAPA-GP on RMSE in ~85% of genes; spaGAPA-GP's value is NOT
  global RMSE minimisation but spatial reconstruction + calibrated
  uncertainty + scalable probabilistic inference.

All panel text >= 7 pt; figure width = PAGE_WIDTH_IN (~7 in). Explanatory
prose lives in this docstring / figure caption, not inside panels.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import (setup_rc, panel_label, save, PAGE_WIDTH_IN,
                    METHOD_COLORS, METHOD_ORDER, METHOD_LABELS,
                    BLUE, ORANGE, GREEN, SKYBLU, YELLOW, RED, GREY, BLACK)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import matplotlib.gridspec as gridspec

setup_rc()

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA = os.path.join(ROOT, "pipeline_output")

# ============================================================
# Load real data
# ============================================================
df = pd.read_csv(f"{DATA}/benchmark_mean_transparent/transparent_comparison.csv")
df = df[df.method.isin(METHOD_ORDER)].copy()

s4 = pd.read_csv(f"{DATA}/supplementary_figures/_cache/s4_stratification.csv")
s4 = s4.dropna(subset=["morans_i", "delta_rmse"])
s4 = s4[np.isfinite(s4["delta_rmse"])].copy()
s4["quintile"] = pd.qcut(s4["morans_i"], 5, labels=[1, 2, 3, 4, 5]).astype(int)

# Panel E: real GP posterior on peak_94938 (produced by fig2E_run_real_gp.py)
GENE2E = "peak_94938"
gp_path = f"{DATA}/main_figures/_data/fig2E_gse183456_{GENE2E}_gp_posterior.csv"
gp_df = pd.read_csv(gp_path)
import json
gp_meta = json.load(open(f"{DATA}/main_figures/_data/fig2E_gse183456_{GENE2E}_meta.json"))

# ============================================================
# Figure layout: PAGE_WIDTH_IN wide, 4 rows x 4 cols
#   row 0: A (schematic, span 2) | B1 RMSE | B2 corr
#   row 1: C spatial fidelity (span 2) | D accuracy-spatial (span 2)
#   row 2-3: E gene maps (2x4 mini inside span 2) | F boxplot (span 2)
# ============================================================
FIG_H = 9.4
fig = plt.figure(figsize=(PAGE_WIDTH_IN, FIG_H))
gs = fig.add_gridspec(
    3, 4, height_ratios=[1.05, 1.05, 1.35],
    hspace=0.62, wspace=0.42,
    left=0.07, right=0.975, top=0.965, bottom=0.065,
)
axA  = fig.add_subplot(gs[0, 0:2])
axB1 = fig.add_subplot(gs[0, 2])
axB2 = fig.add_subplot(gs[0, 3])
axC  = fig.add_subplot(gs[1, 0:2])
axD  = fig.add_subplot(gs[1, 2:])
# row 2 split: E (left 2 cols) and F (right 2 cols)
sub_gs_E = gridspec.GridSpecFromSubplotSpec(
    2, 4, subplot_spec=gs[2, 0:2], wspace=0.04, hspace=0.18)
axF = fig.add_subplot(gs[2, 2:])

# ---------------- Panel A: masking design schematic ----------------
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.set_title("Masking design: 20% per-gene holdout", loc="left", fontsize=9)
np.random.seed(42)
gx, gy = 0.06, 0.22
cw, ch = 0.019, 0.060
nrows, ncols = 8, 38
obs = np.random.rand(nrows, ncols) < 0.30
held = np.zeros_like(obs, dtype=bool)
for r in range(nrows):
    idx = np.where(obs[r])[0]
    np.random.shuffle(idx)
    held[r, idx[:int(len(idx) * 0.20)]] = True
for r in range(nrows):
    for c in range(ncols):
        col = RED if held[r, c] else (BLUE if obs[r, c] else "#E8E8E8")
        axA.add_patch(Rectangle((gx + c * cw, gy + (nrows - 1 - r) * ch),
                                cw * 0.9, ch * 0.9, fc=col, ec="white", lw=0.2))
axA.text(0.06, 0.80, "Gene", fontsize=7.5, fontweight="bold")
axA.text(0.06, 0.74, "↓", fontsize=8)
axA.text(0.52, 0.11, "Spot →", fontsize=7.5, fontweight="bold", ha="center")
leg = [Patch(fc=BLUE, label="Observed (train+cal)"),
       Patch(fc=RED, label="Held-out (test)"),
       Patch(fc="#E8E8E8", label="Unobserved")]
axA.legend(handles=leg, fontsize=6.6, loc="upper right",
           bbox_to_anchor=(0.99, 0.96), handlelength=1.1, borderpad=0.3)
panel_label(axA, "A", x=-0.06, y=1.10)

# ---------------- Panel B1: entry-wise RMSE ----------------
for i, m in enumerate(METHOD_ORDER):
    vals = df[df.method == m]["rmse"].dropna().values
    xs = np.full(len(vals), i) + np.random.RandomState(i).normal(0, 0.05, len(vals))
    axB1.scatter(xs, vals, s=22, color=METHOD_COLORS[m], edgecolor="white",
                 lw=0.5, zorder=3, alpha=0.9)
    axB1.plot([i - 0.13, i + 0.13], [vals.mean(), vals.mean()],
              color="black", lw=1.4, zorder=4)
axB1.set_xticks(range(len(METHOD_ORDER)))
axB1.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER],
                     rotation=35, ha="right", fontsize=6.5)
axB1.set_ylabel("RMSE (↓)", fontsize=8)
axB1.set_title("Entry-wise RMSE", fontsize=8.5)
axB1.tick_params(labelsize=7)
panel_label(axB1, "B", x=-0.30, y=1.14)

# ---------------- Panel B2: Pearson & Spearman ----------------
width = 0.35
xpos = np.arange(len(METHOD_ORDER))
for k, (met, col, lab) in enumerate([("pearson", BLUE, "Pearson r"),
                                     ("spearman", GREEN, "Spearman ρ")]):
    means = [df[df.method == m][met].dropna().mean() for m in METHOD_ORDER]
    stds  = [df[df.method == m][met].dropna().std()  for m in METHOD_ORDER]
    off = (k - 0.5) * width
    axB2.bar(xpos + off, means, width, color=col, yerr=stds, capsize=2,
             edgecolor="white", lw=0.5, label=lab)
    for i, m in enumerate(METHOD_ORDER):
        vals = df[df.method == m][met].dropna().values
        xs = np.full(len(vals), xpos[i] + off) + np.random.RandomState(i + 10).normal(0, 0.03, len(vals))
        axB2.scatter(xs, vals, s=10, color="black", alpha=0.6, zorder=4)
axB2.set_xticks(xpos)
axB2.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER],
                     rotation=35, ha="right", fontsize=6.5)
axB2.set_ylabel("Correlation", fontsize=8)
axB2.set_ylim(0.65, 1.0)
axB2.legend(loc="lower right", fontsize=6.3, handlelength=1.1)
axB2.set_title("Rank correlation", fontsize=8.5)
axB2.tick_params(labelsize=7)

# ---------------- Panel C: spatial fidelity ----------------
axC.set_title("Spatial fidelity (gradient recovery)", loc="left", fontsize=9)
means = [df[df.method == m]["spatial_fidelity"].dropna().mean() for m in METHOD_ORDER]
axC.bar(range(len(METHOD_ORDER)), means,
        color=[METHOD_COLORS[m] for m in METHOD_ORDER],
        edgecolor="white", lw=0.8, width=0.65)
for i, m in enumerate(METHOD_ORDER):
    vals = df[df.method == m]["spatial_fidelity"].dropna().values
    xs = np.full(len(vals), i) + np.random.RandomState(i).normal(0, 0.04, len(vals))
    axC.scatter(xs, vals, color="black", s=16, zorder=4)
axC.set_xticks(range(len(METHOD_ORDER)))
axC.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], fontsize=7)
axC.set_ylabel("Spatial fidelity", fontsize=8)
axC.set_ylim(-0.05, 1.0)
axC.axhline(0, color="#888", lw=0.6, ls="--")
axC.annotate("spaGAPA-GP 0.42 vs Mean 0.00\n(GP reconstructs gradient; mean is flat)",
             xy=(0, 0.42), xytext=(0.65, 0.70), fontsize=6.8,
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9),
             color=BLUE, fontweight="bold")
axC.tick_params(labelsize=7)
panel_label(axC, "C", x=-0.07, y=1.08)

# ---------------- Panel D: accuracy vs spatial 2D ----------------
axD.set_title("Accuracy vs spatial fidelity", loc="left", fontsize=9)
for m in METHOD_ORDER:
    sub = df[df.method == m]
    axD.errorbar(sub["rmse"].mean(), sub["spatial_fidelity"].mean(),
                 xerr=sub["rmse"].std(), yerr=sub["spatial_fidelity"].std(),
                 fmt="o", color=METHOD_COLORS[m], ms=8, capsize=2.5,
                 markeredgecolor="white", lw=1.0,
                 label=METHOD_LABELS[m], zorder=4)
axD.annotate("Mean: lowest RMSE,\nzero gradient", xy=(0.080, 0.0), xytext=(0.17, 0.13),
             fontsize=6.5, arrowprops=dict(arrowstyle="->", color=GREY, lw=0.8), color=GREY)
axD.annotate("spaGAPA-GP: moderate RMSE,\nstrong gradient", xy=(0.122, 0.42), xytext=(0.20, 0.60),
             fontsize=6.5, arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8), color=BLUE)
axD.set_xlabel("Entry-wise RMSE (↓)", fontsize=8)
axD.set_ylabel("Spatial fidelity", fontsize=8)
axD.legend(fontsize=6.4, loc="upper right", handlelength=1.1)
axD.set_xlim(0.05, 0.30)
axD.tick_params(labelsize=7)
panel_label(axD, "D", x=-0.07, y=1.08)

# ============================================================
# Panel E: representative APA spatial maps (REAL spaGAPA-GP)
#   one high-spatial-signal gene (peak_94938, Moran's I = +0.50)
#   2 rows: row0 = the 4 maps; row1 = blank spacer (kept for symmetry)
#   truth / mean / spaGAPA-GP posterior / |GP - truth|
# ============================================================
x_e = gp_df["x"].values
y_e = gp_df["y"].values
truth = gp_df["truth"].values
mean_pred = gp_df["pred_mean"].values
gp_pred = gp_df["pred_gp"].values
err = np.abs(gp_pred - truth)
# nan-safe colour limits from observed truth
obs_mask = np.isfinite(truth) & (truth > 0)
vmin = float(np.percentile(truth[obs_mask], 2)) if obs_mask.any() else 0.0
vmax = float(np.percentile(truth[obs_mask], 98)) if obs_mask.any() else 1.0
err_vmax = max(float(np.nanmax(err)), 1e-3)

panels_e = [
    ("Truth",            truth,     "viridis", vmin, vmax),
    ("Mean",             mean_pred, "viridis", vmin, vmax),
    ("spaGAPA-GP",       gp_pred,   "viridis", vmin, vmax),
    ("|GP − truth|",     err,       "Reds",    0.0,   err_vmax),
]
for mi, (lab, pdat, cmap, lo, hi) in enumerate(panels_e):
    axm = fig.add_subplot(sub_gs_E[0, mi])
    sc = axm.scatter(x_e, y_e, c=np.nan_to_num(pdat, nan=0.0),
                     s=2.4, cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
    axm.set_xticks([]); axm.set_yticks([])
    for s in axm.spines.values():
        s.set_linewidth(0.5)
    axm.set_title(lab, fontsize=7, fontweight="bold", pad=2)
    # compact shared colourbar-less; one bar under the error panel
# Hide the unused second sub-row (kept for visual breathing room / caption)
for mi in range(4):
    axhid = fig.add_subplot(sub_gs_E[1, mi]); axhid.axis("off")

# title + panel label on an invisible anchor over the E region
axE_anchor = fig.add_subplot(gs[2, 0:2]); axE_anchor.axis("off")
axE_anchor.set_title(
    f"Representative APA spatial map — {GENE2E}  (Moran's I = {gp_meta['morans_i']:+.2f})",
    loc="left", fontsize=8.5,
)
panel_label(axE_anchor, "E", x=-0.04, y=1.12)
axE_anchor.set_xlim(0, 1); axE_anchor.set_ylim(0, 1)

# ============================================================
# Panel F: per-gene ΔRMSE (GP − mean) by Moran's I quintile
#   REAL data from s4_stratification.csv (150 genes, GSE183456, 20% holdout).
#   Honest: mean beats GP on RMSE in ~85% of genes; GP's value is spatial
#   reconstruction + calibrated uncertainty, not global RMSE minimisation.
# ============================================================
quint_colors = [GREY, SKYBLU, GREEN, ORANGE, BLUE]
data_F = [s4[s4["quintile"] == q]["delta_rmse"].values for q in range(1, 6)]
pos_F = np.arange(1, 6)
bp = axF.boxplot(data_F, positions=pos_F, widths=0.6, patch_artist=True,
                 showfliers=False, medianprops={"color": "black", "lw": 1.1})
for patch, c in zip(bp["boxes"], quint_colors):
    patch.set_facecolor(c); patch.set_alpha(0.8); patch.set_edgecolor("#555")
for w in bp["whiskers"]:
    w.set_color("#555"); w.set_linewidth(0.9)
for cap in bp["caps"]:
    cap.set_color("#555"); cap.set_linewidth(0.9)
axF.axhline(0, color=RED, lw=1.0, ls="--",
            label="GP = mean (no RMSE advantage)")
# median Δ label per quintile
mi_med = s4.groupby("quintile")["morans_i"].median()
d_med = s4.groupby("quintile")["delta_rmse"].median()
ymax = max(np.nanmax(np.concatenate(data_F)), 0.0)
for q in range(1, 6):
    axF.annotate(f"med\n{d_med[q]:+.3f}", xy=(q, ymax * 0.96),
                 fontsize=6.2, ha="center", color="#333")
axF.set_xticks(pos_F)
axF.set_xticklabels(
    [f"Q{q}\nI={mi_med[q]:+.2f}" for q in range(1, 6)],
    fontsize=6.8,
)
axF.set_xlabel("Moran's I quintile  (Q1 = low → Q5 = high spatial signal)", fontsize=8)
axF.set_ylabel("Δ RMSE  (spaGAPA-GP − mean)", fontsize=8)
axF.set_title("Per-gene RMSE gap is positive (mean wins); gap widens at high spatial signal",
              loc="left", fontsize=8.2)
axF.legend(loc="upper left", fontsize=6.4, handlelength=1.2)
axF.tick_params(labelsize=7)
# mean Δ line for context
mean_delta = s4["delta_rmse"].mean()
frac_gp = (s4["delta_rmse"] < 0).mean() * 100
axF.text(0.5, -0.30,
         f"n = {len(s4)} genes (GSE183456, 20% hold-out).  Mean Δ = {mean_delta:+.3f};  "
         f"spaGAPA-GP beats mean in only {frac_gp:.0f}% of genes.\n"
         f"spaGAPA-GP's value is spatial reconstruction + calibrated uncertainty + scalable inference, not RMSE.",
         transform=axF.transAxes, fontsize=6.2, ha="center",
         style="italic", color="#555")
panel_label(axF, "F", x=-0.10, y=1.08)

save(fig, "fig2_spatial_vs_mean_benchmark.png")
print("Figure 2 done")

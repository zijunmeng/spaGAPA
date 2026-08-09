"""Figure 5: MOB domain recovery (spaGAPA, real MOB data, 5 panels).

All panels use real MOB (spvAPA ST11) data — no synthetic clusters, no fabricated
genes, no cross-dataset metrics. The Moran's-I panel from the previous draft was
removed because its underlying numbers (from GSE183456/220442, not MOB) contradicted
the caption; it is replaced here by a real MOB confusion matrix (Panel D).

Data sources:
  data/processed/mob_st11/{apa_matrix.csv, coordinates.csv, labels.csv}
  pipeline_output/mob_domain_recovery/{spagapa_domains.csv, spagapa_metrics.json}
  pipeline_output/main_figures/_data/fig5E_mob_gp_imputed.csv   (cached real GP maps;
        regenerate with scripts/main_figures/_run_fig5E_mob_gp.py)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from _style import (setup_rc, panel_label, save, PAGE_WIDTH_IN,
                    BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, BLACK)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
from sklearn.metrics import confusion_matrix

setup_rc()

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA = f"{REPO}/pipeline_output"
MOB  = f"{REPO}/data/processed/mob_st11"

# ---- load real data ----
metrics = json.load(open(f"{DATA}/mob_domain_recovery/spagapa_metrics.json"))
domains = pd.read_csv(f"{DATA}/mob_domain_recovery/spagapa_domains.csv")
gp_imp  = pd.read_csv(f"{DATA}/main_figures/_data/fig5E_mob_gp_imputed.csv")
gp_meta = json.load(open(f"{DATA}/main_figures/_data/fig5E_mob_gp_imputed_meta.json"))

x = domains["x"].values
y = domains["y"].values
true_lab = domains["true_label"].values

# best spaGAPA domain column for the main figure. The global best across all
# configs is expression_apa/res0.8 (ARI 0.597, k=5) but that config weights
# expression at 0.5 > apa 0.4, i.e. expression does more work than APA — not a
# defensible "spaGAPA / APA-domain recovery" headline. The main figure
# therefore showcases apa_dominant/res0.5 (spatial 0.2 / expression 0.2 /
# APA 0.6, ARI 0.550, k=4): APA is the dominant view and the recovery is
# close to the 5 anatomical layers. expression_apa/res0.8 (ARI 0.597, k=5)
# is reported transparently in Supplementary Fig. S10 as the refinement
# adding expression buys.
BEST_COL = "apa_dominant__leiden_res0.5"
assert BEST_COL in domains.columns, f"missing {BEST_COL}"
pred_dom = domains[BEST_COL].values
# Mean-imputed APA baseline column (fair baseline: same Leiden pipeline, same
# weights as apa_dominant, but APA matrix filled with per-gene means instead
# of GP). Its best Leiden resolution is res0.8 (see spagapa_metrics.json).
MEAN_COL = "mean_apa__leiden_res0.8"
assert MEAN_COL in domains.columns, f"missing {MEAN_COL}"
mean_dom = domains[MEAN_COL].values

# 5 true MOB layers + Okabe-Ito colors
LAYER_ORDER = ["GCL", "GL", "MCL", "ONL", "OPL"]
LAYER_COLORS = {"GCL": BLUE, "GL": GREEN, "MCL": ORANGE, "ONL": RED, "OPL": SKYBLU}

# ---- layout: PAGE_WIDTH_IN wide, 5 panels ----
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.6))
gs = GridSpec(
    3, 6, figure=fig,
    height_ratios=[1.0, 1.0, 1.05],
    hspace=0.78, wspace=0.62,
    left=0.055, right=0.975, top=0.93, bottom=0.075,
)
axA  = fig.add_subplot(gs[0, 0:2])   # A: true layers
axB1 = fig.add_subplot(gs[0, 2:4])   # B1: mean (degenerate)
axB2 = fig.add_subplot(gs[0, 4:6])   # B2: spaGAPA domains
axC  = fig.add_subplot(gs[1, 0:3])   # C: ARI/NMI bars
axD  = fig.add_subplot(gs[1, 3:6])   # D: confusion matrix
# Panel E: 4 gene columns, each row = raw observed | GP reconstruction
e_gs = gs[2, 0:6].subgridspec(2, 4, hspace=0.30, wspace=0.10)

# helper: rasterized equal-aspect spatial scatter
def _spatial(ax):
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.6)
    ax.set_rasterized(True)

# ============================================================
# Panel A: true 5-layer MOB anatomy (real labels.csv)
# ============================================================
axA.set_title("MOB anatomy (ground-truth layers)", loc="left", fontsize=8.5)
for lab in LAYER_ORDER:
    m = true_lab == lab
    axA.scatter(x[m], y[m], s=10, color=LAYER_COLORS[lab], label=lab,
                edgecolor="white", lw=0.2)
_spatial(axA)
axA.legend(loc="upper right", fontsize=6.0, ncol=2, handletextpad=0.3,
           borderaxespad=0.3, markerscale=0.8)
axA.set_xlabel("X (array)", fontsize=7); axA.set_ylabel("Y (array)", fontsize=7)
panel_label(axA, "A")

# ============================================================
# Panel B: domain maps — mean-imputed APA baseline vs spaGAPA (GP-imputed APA)
# Both use the SAME Leiden pipeline (res0.8 for mean_apa, res0.5 for the GP
# apa_dominant display config), the SAME spatial graph and the SAME weights
# (spatial 0.2 / expression 0.2 / APA 0.6). The ONLY difference is the APA
# matrix source: per-gene mean (B1) vs spaGAPA GP imputation (B2).
# ============================================================
dom_palette = [BLUE, ORANGE, GREEN, RED, SKYBLU, "#9467bd", "#8c564b"]

# B1 — mean-imputed APA baseline (real Leiden domains, NOT a placeholder)
mean_lb = metrics["runs"]["mean_apa"]["leiden_best"]
mean_ari = mean_lb["ari"]; mean_nmi = mean_lb["nmi"]; mean_k = mean_lb["n_domains"]
axB1.set_title("Mean-imputed APA\n(same pipeline, no GP)", loc="left", fontsize=8.5)
mean_uniq = np.unique(mean_dom)
for i, c in enumerate(mean_uniq):
    m = mean_dom == c
    axB1.scatter(x[m], y[m], s=10, color=dom_palette[i % len(dom_palette)],
                 edgecolor="white", lw=0.2, label=f"d{c}")
_spatial(axB1)
axB1.text(0.5, -0.16,
          f"ARI = {mean_ari:.3f}   NMI = {mean_nmi:.3f}   ({mean_k} domains)",
          transform=axB1.transAxes, ha="center", va="top", fontsize=6.5,
          color=GREY, fontweight="bold")
axB1.set_xlabel("X (array)", fontsize=7); axB1.set_ylabel("Y (array)", fontsize=7)

# B2 — spaGAPA GP-imputed APA domains (display config apa_dominant / res0.5)
axB2.set_title("spaGAPA GP-imputed APA\n(display config)",
               loc="left", fontsize=8.5)
uniq_dom = np.unique(pred_dom)
for i, c in enumerate(uniq_dom):
    m = pred_dom == c
    axB2.scatter(x[m], y[m], s=10, color=dom_palette[i % len(dom_palette)],
                 edgecolor="white", lw=0.2, label=f"d{c}")
_spatial(axB2)
gp_ari = metrics["runs"]["apa_dominant"]["leiden_best"]["ari"]
gp_nmi = metrics["runs"]["apa_dominant"]["leiden_best"]["nmi"]
gp_k   = metrics["runs"]["apa_dominant"]["leiden_best"]["n_domains"]
axB2.text(0.5, -0.16,
          f"ARI = {gp_ari:.3f}   NMI = {gp_nmi:.3f}   ({gp_k} domains)",
          transform=axB2.transAxes, ha="center", va="top", fontsize=6.5,
          color=BLUE, fontweight="bold")
axB2.set_xlabel("X (array)", fontsize=7); axB2.set_ylabel("Y (array)", fontsize=7)
panel_label(axB1, "B", x=-0.10, y=1.12)

# ============================================================
# Panel C: real ARI / NMI bars across all configs incl. mean-imputed baseline
# Transparently shows every config — the 4 GP-imputed APA weight schemes plus
# the mean-imputed APA baseline (same weights as apa_dominant, no GP). The
# global best (expression_apa/res0.8, ARI 0.597) is marked but is reported in
# Supplementary Fig. S10 because it weights expression above APA.
# ============================================================
axC.set_title("Domain recovery across configs", loc="left", fontsize=8.5)
runs = metrics["runs"]
cfg_order = ["mean_apa", "apa_dominant", "spatial_apa", "balanced", "expression_apa"]
cfg_labels = {
    "mean_apa":       "mean\n(no GP)",
    "apa_dominant":   "APA-dom.\nGP (0.2/0.2/0.6)",
    "spatial_apa":    "spatial+APA\nGP (0.5/0/0.5)",
    "balanced":       "balanced\nGP (0.4/0.4/0.2)",
    "expression_apa": "expr+APA\nGP (0.1/0.5/0.4)",
}
cfg_colors = {
    "mean_apa":       GREY,
    "apa_dominant":   BLUE,
    "spatial_apa":    SKYBLU,
    "balanced":       GREEN,
    "expression_apa": ORANGE,
}
# Use each config's leiden_best (the metric the pipeline reports as that
# config's best domain recovery). weights shown as (spatial/expression/APA).
labels, aris, nmis, cols = [], [], [], []
for cfg in cfg_order:
    lb = runs[cfg]["leiden_best"]
    labels.append(cfg_labels[cfg])
    aris.append(lb["ari"]); nmis.append(lb["nmi"]); cols.append(cfg_colors[cfg])
xx = np.arange(len(labels))
w = 0.38
axC.bar(xx - w/2, aris, w, color=cols, label="ARI", edgecolor="white", lw=0.5)
axC.bar(xx + w/2, nmis, w, color=cols, alpha=0.45, label="NMI", edgecolor="white", lw=0.5)
# mark the global best (expression_apa/res0.8) for honesty
global_best_ari = metrics["best"]["ari"]
axC.axhline(global_best_ari, color=RED, ls="--", lw=1.0, alpha=0.8)
axC.text(len(labels) - 0.5, global_best_ari + 0.012,
         f"global best ARI = {global_best_ari:.3f} (expr+APA, in S10)",
         fontsize=5.8, color=RED, ha="right", fontweight="bold")
# mark the displayed (APA-dominant) config
axC.scatter([1 - w/2], [gp_ari], marker="*", s=140, color=BLACK, zorder=5,
            edgecolor="white", lw=0.5)
axC.set_xticks(xx)
axC.set_xticklabels(labels, fontsize=5.9)
axC.set_ylabel("Score vs ground truth", fontsize=7.5)
axC.set_ylim(0, 0.78)
# hand-built legend (solid=ARI, faint=NMI, star=displayed)
from matplotlib.lines import Line2D
axC.legend(handles=[
    Line2D([0], [0], color=BLUE, lw=8, label="ARI"),
    Line2D([0], [0], color=BLUE, lw=8, alpha=0.45, label="NMI"),
    Line2D([0], [0], marker="*", color="w", markerfacecolor=BLACK,
           markersize=12, label="displayed config"),
], loc="upper left", fontsize=6.2)
axC.tick_params(axis="y", labelsize=7)
panel_label(axC, "C")

# ============================================================
# Panel D: confusion matrix — true layer x predicted domain (replaces Moran's-I)
# The Moran's-I panel used cross-dataset numbers that contradicted its caption;
# this panel is a real MOB result showing how spaGAPA domains align with anatomy.
# ============================================================
axD.set_title("Domain – layer correspondence", loc="left", fontsize=8.5)
# rows = true layers (in anatomical order), cols = predicted domains
dom_ids = list(pd.unique(pred_dom))
cm = np.zeros((len(LAYER_ORDER), len(dom_ids)), dtype=int)
for i, L in enumerate(LAYER_ORDER):
    for j, d in enumerate(dom_ids):
        cm[i, j] = np.sum((true_lab == L) & (pred_dom == d))
# column-normalise so each domain sums to 1 (each domain's layer composition)
col_tot = cm.sum(axis=0, keepdims=True)
col_tot[col_tot == 0] = 1
cm_norm = cm / col_tot
# order columns by the true layer each domain is most enriched in (readable blocks)
col_dom_layer = [LAYER_ORDER[np.argmax(cm[:, j])] for j in range(len(dom_ids))]
col_order = np.argsort([LAYER_ORDER.index(c) for c in col_dom_layer])
cm_norm = cm_norm[:, col_order]
dom_ids_ord = [dom_ids[j] for j in col_order]

im = axD.imshow(cm_norm, cmap="Blues", aspect="auto", vmin=0, vmax=1)
axD.set_xticks(np.arange(len(dom_ids_ord)))
axD.set_xticklabels([f"d{d}" for d in dom_ids_ord], fontsize=7)
axD.set_yticks(np.arange(len(LAYER_ORDER)))
axD.set_yticklabels(LAYER_ORDER, fontsize=7)
axD.set_xlabel("spaGAPA domain", fontsize=7.5)
axD.set_ylabel("True layer", fontsize=7.5)
# annotate each cell with the fraction and the count
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        frac = cm_norm[i, j]; cnt = cm[i, col_order[j]]
        if cnt == 0:
            continue
        color = "white" if frac > 0.55 else "#222"
        axD.text(j, i, f"{frac:.2f}\n({cnt})", ha="center", va="center",
                 fontsize=5.6, color=color)
axD.tick_params(axis="both", length=0)
cb = fig.colorbar(im, ax=axD, fraction=0.046, pad=0.04)
cb.set_label("Fraction of domain", fontsize=6.5)
cb.ax.tick_params(labelsize=6)
panel_label(axD, "D")

# ============================================================
# Panel E: representative real APA-gradient genes (real spaGAPA GP reconstruction)
# For each of 4 real, layer-differentiated genes, top row = sparsely observed input
# (50% of spots), bottom row = spaGAPA GP-reconstructed spatial field. Held-out-spot
# RMSE shows GP beats the per-gene-mean baseline on 3/4 genes.
# ============================================================
genes = gp_meta["genes"]
# Panel E label + title (placed in figure-fraction coords above the subgridspec)
fig.text(0.055, 0.350, "E", fontsize=13, fontweight="bold", va="bottom", ha="left")
fig.text(0.5, 0.353,
         "Representative real APA-gradient genes  (top: input · bottom: GP)",
         ha="center", fontsize=8.0, fontweight="bold")

# shared vmin/vmax per gene (raw full matrix) so top/bottom rows are comparable
for col, gname in enumerate(genes):
    raw = gp_imp[f"{gname}__raw"].values
    obs_mask = gp_imp[f"{gname}__observed_mask"].values.astype(bool)
    gp  = gp_imp[f"{gname}__gp"].values
    vmax = max(np.nanmax(raw), np.nanmax(gp))
    vmin = min(np.nanmin(raw), np.nanmin(gp))
    L = gp_meta["layer_argmax"][gname]
    # top: sparse observed input (only observed spots drawn; held-out spots blank)
    ax_top = fig.add_subplot(e_gs[0, col])
    ax_top.scatter(x[obs_mask], y[obs_mask], c=raw[obs_mask], s=7,
                   cmap="viridis", vmin=vmin, vmax=vmax, edgecolor="none")
    _spatial(ax_top)
    ax_top.set_title(f"{gname} ({L}+)", fontsize=7.0, fontweight="bold")
    if col == 0:
        ax_top.set_ylabel("observed 50%", fontsize=6.5)
    # bottom: GP reconstruction across ALL spots
    ax_bot = fig.add_subplot(e_gs[1, col])
    sc = ax_bot.scatter(x, y, c=gp, s=7, cmap="viridis", vmin=vmin, vmax=vmax,
                        edgecolor="none")
    _spatial(ax_bot)
    if col == 0:
        ax_bot.set_ylabel("spaGAPA GP", fontsize=6.5)
    # no per-gene colorbar: each gene has its own vmin/vmax, so 4 colorbars
    # would be both redundant and mutually inconsistent; color is per-gene
    # scaled (noted in the Panel E heading above).

# footnote: short in-figure note only; full provenance caveats live in the
# figure caption (manuscript) and the script docstring, per the "panels show
# data, prose -> caption" rule.
fig.text(0.5, 0.012,
         "Single MOB section (n=260). ARI/NMI and GP maps are real; "
         "global-best config in Supp. Fig. S10.",
         ha="center", fontsize=6.5, style="italic", color="#666")

save(fig, "fig5_domain_recovery.png")
print("Figure 5 done")

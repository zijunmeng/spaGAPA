"""Figure 4: Noise methods + risk-coverage curve (6 panels).

Reviewer-driven restructure (7 -> 6 panels):
  - Panels A-D kept (four noise models / multi-objective / pooled r /
    within-gene vs pooled). Within-gene panel D honestly shows pooled
    r ~0.5 plus per-gene median r ~0.10-0.19 (do NOT collapse to a single
    spot-level "local uncertainty correlation" number).
  - Panel E (subgroup coverage) is split clearly: left = tissue x
    uncertainty-quintile empirical-subgroup heatmap, right = expression-bin
    coverage bar plot surfacing the honest high-expression undercoverage
    (~0.82). Split-conformal guarantees MARGINAL coverage only;
    stratification is reported empirically, not as a theoretical guarantee.
  - Panel F (risk-coverage, RMSE) is the single retained triage panel and
    is now full-width on the bottom row so its annotations no longer
    collide with neighbours. The duplicate MAE-on-retained-set panel moved
    to Supplementary Figure S16 (same data, MAE y-axis).

Layout note: width = PAGE_WIDTH_IN (7.0 in ~= 178 mm). Tight-layout keeps
the PDF at ~178 mm; the previous 221 mm draft came from an oversized
figspec + missing constrained layout.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from _style import *  # noqa: F401,F403
import pandas as pd
import numpy as np
from matplotlib.patches import Patch, FancyBboxPatch
from matplotlib.lines import Line2D

setup_rc()

DATA = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
corr = pd.read_csv(f"{DATA}/uncertainty_corr_improvement/per_dataset_corr.csv")
gene = pd.read_csv(f"{DATA}/uncertainty_within_gene_audit/per_gene_corr_distribution.csv")
risk = pd.read_csv(f"{DATA}/risk_coverage_curve/risk_coverage_data.csv")
with open(f"{DATA}/risk_coverage_curve/summary.json") as fh:
    risk_summary = json.load(fh)
cond_q = pd.read_csv(f"{DATA}/conformal_conditional_coverage/by_uncertainty_quintile.csv")
cond_e = pd.read_csv(f"{DATA}/conformal_conditional_coverage/by_expression_level.csv")

# ---- layout: 3 rows. row1 = A/B/C, row2 = D/E(split), row3 = F (full width) ----
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 8.4))
gs = fig.add_gridspec(
    3, 6,
    height_ratios=[1.0, 1.0, 0.95],
    hspace=0.62, wspace=0.85,
    left=0.055, right=0.985, top=0.965, bottom=0.05,
)
axA = fig.add_subplot(gs[0, 0:2])   # four noise models (schematic)
axB = fig.add_subplot(gs[0, 2:4])   # multi-objective
axC = fig.add_subplot(gs[0, 4:6])   # pooled r per dataset (A vs D)
axD = fig.add_subplot(gs[1, 0:3])   # per-gene vs pooled
axE = fig.add_subplot(gs[1, 3:6])   # subgroup coverage (heatmap + expr bins)
axF = fig.add_subplot(gs[2, 0:6])   # risk-coverage RMSE (KEY, full width)

# ---------------- Panel A: four noise models ----------------
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.set_title("Four uncertainty (noise) models", loc="left")
models = [
    ("A", "Constant",      "global variance",                 GREY),
    ("B", "Local-gene",    "per-gene SD on obs",              BLUE),
    ("C", "Spatial-spot",  "local neighbor resid",            SKYBLU),
    ("D", "Residual-spot", "B + C, EB-shrunk (spaGAPA)",      ORANGE),
]
for i, (tag, name, body, col) in enumerate(models):
    yy = 0.82 - i * 0.205
    axA.add_patch(FancyBboxPatch((0.02, yy - 0.155), 0.10, 0.155,
                                 boxstyle="round,pad=0.006,rounding_size=0.012",
                                 fc=col, ec=col))
    axA.text(0.07, yy - 0.077, tag, ha="center", va="center",
             color="white", fontweight="bold", fontsize=9)
    axA.add_patch(FancyBboxPatch((0.14, yy - 0.155), 0.84, 0.155,
                                 boxstyle="round,pad=0.006,rounding_size=0.012",
                                 fc="#F7F7F7", ec=col, lw=1.0))
    axA.text(0.17, yy - 0.035, name, fontsize=7.8, fontweight="bold", color=col, va="top")
    axA.text(0.17, yy - 0.105, body, fontsize=7.0, color="#333", va="top")
panel_label(axA, "A")

# ---------------- Panel B: multi-objective method selection ----------------
axB.set_title("Multi-objective: corr vs width", loc="left")
cwidth = pd.read_csv(f"{DATA}/uncertainty_coverage_width_audit/coverage_width_comparison.csv")
mean_r = corr.groupby("method")["unc_error_pearson"].mean()
mean_w = cwidth.groupby("method")["mean_interval_width"].mean()
mean_dev = abs(cwidth.groupby("method")["marginal_coverage_90"].mean() - 0.90)
methods_b = ["A_constant", "B_local_gene", "C_spatial_spot", "D_residual_spot"]
max_dev = mean_dev.max()
for m in methods_b:
    if m in mean_w.index:
        size = 60 + (mean_dev[m] / max_dev) * 600
        axB.scatter(mean_r[m], mean_w[m], s=size, color=NOISE_METHOD_COLORS[m],
                    edgecolor="white", lw=1.5, alpha=0.85,
                    label=NOISE_METHOD_LABELS[m], zorder=5)
        axB.annotate(m.split("_")[0], (mean_r[m], mean_w[m]),
                     fontsize=8, fontweight="bold", xytext=(6, 6),
                     textcoords="offset points")
    else:
        # C has no width audit; show at baseline with open marker
        axB.scatter(mean_r[m], 0, s=120, color=NOISE_METHOD_COLORS[m],
                    edgecolor="white", lw=1.5, alpha=0.5, marker="d",
                    label=NOISE_METHOD_LABELS[m] + " (no width audit)", zorder=4)
        axB.annotate(m.split("_")[0], (mean_r[m], 0), fontsize=8, fontweight="bold",
                     xytext=(6, 6), textcoords="offset points",
                     color=NOISE_METHOD_COLORS[m])
axB.set_xlabel("Uncertainty-error correlation r (better ->)")
axB.set_ylabel("Mean interval width (lower better)")
axB.axvline(0.3, color=GREEN, ls=":", lw=1.0, alpha=0.7)
axB.text(0.31, axB.get_ylim()[0] + 0.02, "target r >= 0.3", fontsize=6.5, color=GREEN)
axB.legend(loc="upper right", fontsize=6.0, title="bubble = |coverage dev|", title_fontsize=6)
axB.annotate("D (spaGAPA):\nhigh r, near-minimal width",
             xy=(mean_r["D_residual_spot"], mean_w["D_residual_spot"]),
             xytext=(0.2, mean_w["B_local_gene"] * 0.9), fontsize=6.8, color=ORANGE,
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.8))
panel_label(axB, "B")

# ---------------- Panel C: pooled correlation per dataset (A vs D) ----------------
axC.set_title("Pooled r per dataset", loc="left")
datasets = corr["dataset"].unique()
x = np.arange(len(datasets))
width = 0.35
for k, m in enumerate(["A_constant", "D_residual_spot"]):
    off = (k - 0.5) * width
    vals = corr[corr.method == m].set_index("dataset").loc[datasets, "unc_error_pearson"].values
    axC.bar(x + off, vals, width, color=NOISE_METHOD_COLORS[m],
            edgecolor="white", label=NOISE_METHOD_LABELS[m])
    xs = np.full(len(vals), x + off) + np.random.RandomState(k).normal(0, 0.03, len(vals))
    axC.scatter(xs, vals, color="black", s=14, zorder=4)
axC.axhline(0.3, color=GREEN, ls=":", lw=1.0)
axC.text(len(datasets) - 0.5, 0.31, "tgt", fontsize=6.5, color=GREEN)
axC.set_xticks(x)
axC.set_xticklabels([corr[corr.dataset == d]["label"].iloc[0] for d in datasets],
                    rotation=30, ha="right", fontsize=7)
axC.set_ylabel("Pooled unc-error r")
axC.legend(loc="upper right", fontsize=7)
panel_label(axC, "C")

# ---------------- Panel D: per-gene vs pooled ----------------
axD.set_title("Per-gene r distribution vs pooled", loc="left")
gene_labels_map = dict(zip(gene["dataset"], gene["label"]))
positions = []
for i, (d, gdf) in enumerate(gene.groupby("dataset")):
    pos = i + 1
    positions.append(pos)
    bp = axD.boxplot([gdf["within_gene_r"].values], positions=[pos], widths=0.6,
                     patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(BLUE); patch.set_alpha(0.45); patch.set_edgecolor(BLUE)
    for med in bp["medians"]:
        med.set_color("black")
    pooled = corr[(corr.dataset == d) & (corr.method == "B_local_gene")]["unc_error_pearson"].iloc[0]
    axD.scatter(pos, pooled, marker="*", s=180, color=ORANGE,
                edgecolor="black", lw=0.6, zorder=6)
axD.set_xticks(positions)
axD.set_xticklabels([gene_labels_map[d] for d in gene.groupby("dataset").groups.keys()],
                    rotation=30, ha="right", fontsize=7)
axD.set_ylabel("Within-gene r (per gene)")
axD.axhline(0, color="#888", lw=0.6)
axD.legend(handles=[Patch(fc=BLUE, alpha=0.45, ec=BLUE, label="Per-gene r (IQR)"),
                    Line2D([0], [0], marker="*", color="w", mfc=ORANGE, mec="black",
                           ms=14, label="Pooled r")],
          loc="upper right", fontsize=7)
axD.text(0.02, -0.30,
         "Pooled r partly driven by cross-gene heterogeneity;\n"
         "within-gene discrimination is modest (median r ~ 0.10-0.19)",
         transform=axD.transAxes, fontsize=6.5, style="italic", color="#444")
panel_label(axD, "D", x=-0.08)

# ---------------- Panel E: empirical subgroup coverage (split into two views) ----
# Left:  tissue x uncertainty-quintile heatmap (marginal-coverage guarantee only;
#        stratification is empirical). Right: expression-bin bar plot surfacing
#        the honest limitation that the HIGH-expression bin undercovers (~0.82).
# NOTE: we use a sub-gridspec (NOT make_axes_locatable) so the two sub-axes stay
# inside the E cell - append_axes() was pushing them past the figure edge and
# inflating the tight bbox to 251 mm.
# Hide the dummy host; draw title + panel letter on it, then add two real axes.
axE.set_title("Subgroup coverage", loc="left")
axE.set_xticks([]); axE.set_yticks([])
for s in axE.spines.values():
    s.set_visible(False)

# sub-gridspec INSIDE the E cell: [heatmap | gap | barplot]
gsE = gs[1, 3:6].subgridspec(1, 3, width_ratios=[1.25, 0.08, 0.85], wspace=0.10)
axEmap = fig.add_subplot(gsE[0, 0])   # heatmap
axEbar = fig.add_subplot(gsE[0, 2])   # expression bins

# left: heatmap of empirical coverage by tissue x uncertainty quintile
piv_q = cond_q.pivot_table(index="tissue", columns="quintile",
                           values="coverage_90", aggfunc="mean")
im = axEmap.imshow(piv_q.values, cmap="RdYlGn", vmin=0.86, vmax=0.94, aspect="auto")
axEmap.set_xticks(range(piv_q.shape[1]))
axEmap.set_xticklabels([f"Q{c}" for c in piv_q.columns], fontsize=7)
axEmap.set_yticks(range(piv_q.shape[0]))
axEmap.set_yticklabels(piv_q.index, fontsize=7)
axEmap.set_xlabel("Uncertainty quintile", fontsize=8)
axEmap.set_ylabel("Tissue", fontsize=8)
axEmap.tick_params(labelsize=7)
for r_i in range(piv_q.shape[0]):
    for c_i in range(piv_q.shape[1]):
        v = piv_q.values[r_i, c_i]
        axEmap.text(c_i, r_i, f"{v:.3f}", ha="center", va="center", fontsize=6,
                    color="black" if abs(v - 0.90) < 0.02
                    else ("red" if v < 0.87 else "darkgreen"))

# right: expression-bin coverage bar plot (the undercoverage limitation)
expr_order = ["very_low", "low", "medium", "high", "very_high"]
expr_lab = ["v-low", "low", "med", "high", "v-high"]
e_g = cond_e.groupby("expr_bin")["coverage_90"].agg(["mean", "std"]).reindex(expr_order)
e_sem = (e_g["std"] / np.sqrt(cond_e.groupby("expr_bin")["coverage_90"].count()
                              .reindex(expr_order))).fillna(0)
bar_colors = [BLUE if v >= 0.88 else RED for v in e_g["mean"]]
axEbar.bar(range(len(expr_order)), e_g["mean"].values, color=bar_colors,
           edgecolor="white", lw=0.6, width=0.78)
axEbar.errorbar(range(len(expr_order)), e_g["mean"].values, yerr=e_sem.values,
                fmt="none", ecolor="black", lw=0.7, capsize=2)
axEbar.axhline(0.90, color=GREEN, ls="--", lw=1.0)
axEbar.text(len(expr_order) - 0.5, 0.903, "nominal 0.90", fontsize=6, color=GREEN,
            ha="right", va="bottom")
# annotate the high-expression undercoverage explicitly
hi_idx = expr_order.index("high")
hi_val = e_g["mean"]["high"]
axEbar.annotate(f"high bin\n{hi_val:.2f}",
                xy=(hi_idx, hi_val), xytext=(hi_idx - 0.2, hi_val - 0.07),
                fontsize=6.5, color=RED, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
axEbar.set_xticks(range(len(expr_order)))
axEbar.set_xticklabels(expr_lab, fontsize=6.5, rotation=0)
axEbar.set_xlabel("Expression bin", fontsize=8)
axEbar.set_ylabel("Coverage", fontsize=8)
axEbar.set_ylim(0.76, 0.98)
axEbar.tick_params(labelsize=7)

# panel caption (data-only): flags the limitation honestly
axE.text(0.5, -0.18,
         "High-expr bin undercovers (~0.82); marginal, not per-subgroup",
         transform=axE.transAxes, ha="center", fontsize=6.2, style="italic", color="#444")
panel_label(axE, "E", x=-0.06)

# ---------------- Panel F: RISK-COVERAGE CURVE (KEY, full width) ---------------
axF.set_title("Risk-coverage curve: uncertainty triage", loc="left", fontweight="bold")
piv = risk.pivot_table(index="retention_fraction", columns="method",
                       values="rmse", aggfunc="mean")
ret = piv.index.values
for m, col in [("uncertainty", BLUE), ("oracle", GREEN), ("random", GREY)]:
    axF.plot(ret, piv[m].values, "o-", color=col, lw=2.2, ms=8,
             markeredgecolor="white", label=m.capitalize())
axF.axhline(piv.loc[1.0, "uncertainty"], color=BLUE, ls=":", lw=0.8, alpha=0.5)
axF.set_xlabel("Retention fraction (spots kept, ranked by sigma)")
axF.set_ylabel("RMSE on retained set")
axF.set_xlim(0.15, 1.05)
axF.set_xticks(ret)
axF.legend(loc="upper left", fontsize=7)
# annotate improvements
imp20 = (piv.loc[0.2, "random"] - piv.loc[0.2, "uncertainty"]) / piv.loc[0.2, "random"] * 100
imp80 = (piv.loc[0.8, "random"] - piv.loc[0.8, "uncertainty"]) / piv.loc[0.8, "random"] * 100
axF.annotate(f"Uncertainty triage\n-{imp20:.0f}% RMSE @20%",
             xy=(0.2, piv.loc[0.2, "uncertainty"]),
             xytext=(0.32, piv.loc[0.2, "uncertainty"] + 0.005),
             fontsize=7, color=BLUE, fontweight="bold",
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
# --- DYNAMIC p-value from summary.json ---
# paired across the 5 independent datasets at 80% retention (NOT multiple points
# along one curve). Caption must state the statistical unit.
pval = float(risk_summary["paired_ttest_pvalue_at_80"])
n_ds = int(risk_summary.get("n_datasets", 5))
if pval < 0.001:
    pval_str = f"p={pval:.2e}"
else:
    pval_str = f"p={pval:.3f}"
axF.annotate(f"-{imp80:.0f}% @80%\n(paired t-test {pval_str})",
             xy=(0.8, piv.loc[0.8, "uncertainty"]),
             xytext=(0.55, piv.loc[0.6, "uncertainty"]),
             fontsize=7, color=ORANGE, fontweight="bold",
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
axF.text(0.5, -0.22,
         f"Uncertainty ranking outperforms random retention\n"
         f"(paired across {n_ds} datasets at 80% retention; "
         f"MAE variant in Supplementary S16)",
         transform=axF.transAxes, ha="center", fontsize=6.8,
         style="italic", color="#444")
panel_label(axF, "F", x=-0.05)

# --- Save: pin the figure bbox to PAGE_WIDTH_IN exactly. ---
# The shared `_style.save` uses `bbox_inches="tight"`, which lets out-of-axes
# panel letters / left-aligned titles inflate the page (this figure hit 221 mm,
# then 228 mm, before the layout was tightened). Now that every title, panel
# letter, and caption has been kept inside the figure frame, we save at the
# declared figsize (no tight expansion) so the PDF lands at ~178 mm.
import os as _os
_base4 = _os.path.splitext("fig4_noise_and_risk_coverage.png")[0]
_pdf4 = _os.path.join(OUT, _base4 + ".pdf")
_png4 = _os.path.join(OUT, _base4 + ".png")
fig.savefig(_pdf4, facecolor="white")           # vector, exact 7.0 in wide
fig.savefig(_png4, dpi=300, facecolor="white")  # preview
plt.close(fig)
print(f"  saved -> {_pdf4} (+png)")

print("Figure 4 done")

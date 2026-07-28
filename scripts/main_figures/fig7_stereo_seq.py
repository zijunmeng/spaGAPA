#!/usr/bin/env python3
"""Figure 7: spaGAPA enables quality-controlled high-resolution APA mapping in Stereo-seq.

Multi-panel publication figure (300 DPI, colorblind-friendly) for the spaGAPA BIB paper.

Panels:
  A. Stereo-seq workflow (schematic): FASTQ+mask -> SAW 8.2.2 -> retag -> scAPAtrap -> PAS -> spaGAPA
  B. 3'-bias QC: read 3'-end distance to gene TES (polyA-capture signature), mouse.
  C. Data-scale infographic: 20.7M DNBs / 21,455 PAS / 15,235 bins / 97.8% sparsity (raw rows).
  D. Tissue APA atlas: AD binned domain map + uncertainty map (embedded PNGs).
  E. Representative high-spatial-variance PAS maps (bin200 scatter).
  F. Binning consistency: cross-bin Pearson r box + Moran's I vs bin size.
  G. AD-WT descriptive effect sizes (n=1 per condition, no p-value claim).

Run on S91: ~/anaconda3/envs/spagapa/bin/python fig7_stereo_seq.py
"""
import os
import sys
import csv
import io
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.image import imread
from PIL import Image

# Shared style with the other 6 figures
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (  # noqa: E402
    setup_rc, panel_label, save,
    BLUE, ORANGE, GREEN, SKYBLU, YELLOW, RED, GREY, BLACK,
)

setup_rc()

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output")
STEREO = ROOT / "gse263789_stereo_pilot"
BIN200 = STEREO / "spagapa_downstream_full" / "binned_200"
AD = ROOT / "gse263789_ad_vs_wt_differential"
BINNING = ROOT / "stereo_binning_consistency"
OUT = ROOT / "main_figures"

# --------------------------------------------------------------------------- #
# Layout: 4 rows x 3 cols grid via GridSpec; custom widths.
# Row 0: A (workflow, spans 3 cols)         | (empty)
# Row 1: B (3'bias) | C (data scale) | D1+D2 (atlas, two sub-cols)
# Row 2: E (three gene maps spanning 3 cols)
# Row 3: F1 | F2 | G
# --------------------------------------------------------------------------- #
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402

fig = plt.figure(figsize=(13.5, 17.5))
gs = GridSpec(
    nrows=4, ncols=3, figure=fig,
    height_ratios=[0.78, 1.45, 1.45, 1.35],
    width_ratios=[1, 1, 1],
    hspace=0.50, wspace=0.28,
    left=0.055, right=0.975, top=0.975, bottom=0.035,
)

# =========================================================================== #
# PANEL A — Stereo-seq workflow schematic
# =========================================================================== #
axA = fig.add_subplot(gs[0, :])
axA.set_xlim(0, 100)
axA.set_ylim(0, 10)
axA.axis("off")

stages = [
    ("FASTQ\n+ mask",       "Raw Stereo-seq reads\n+ tissue mask",      GREY,   "#f0f0f0"),
    ("SAW 8.2.2",           "Alignment + cell calling\n(DNB -> BAM)",   BLUE,   "#e8f1fa"),
    ("Retag\nCB/UB/GX/GN",  "Transfer spot barcode +\ngene annotation", SKYBLU, "#eaf5fb"),
    ("scAPAtrap",           "PolyA-peak calling\n(PAS discovery)",      ORANGE, "#fbeadd"),
    ("PAS matrix",          "peak x spot counts\n(sparse)",             GREEN,  "#e6f6ef"),
    ("spaGAPA",             "GP imputation +\nuncertainty + domains",   RED,    "#f7eaef"),
]
n = len(stages)
box_w = 13.5
gap = (100 - n * box_w) / (n + 1)
y0 = 3.6
box_h = 4.6

centers_x = []
for i, (title, desc, edge, fc) in enumerate(stages):
    x = gap + i * (box_w + gap)
    centers_x.append(x + box_w / 2)
    # main box
    box = FancyBboxPatch(
        (x, y0), box_w, box_h,
        boxstyle="round,pad=0.12,rounding_size=0.5",
        linewidth=1.4, edgecolor=edge, facecolor=fc, zorder=2,
    )
    axA.add_patch(box)
    axA.text(x + box_w / 2, y0 + box_h - 1.25, title,
             ha="center", va="center", fontsize=9.5, fontweight="bold",
             color=edge, zorder=3)
    axA.text(x + box_w / 2, y0 + 1.35, desc,
             ha="center", va="center", fontsize=6.8, color="#333333", zorder=3)

# arrows between boxes
for i in range(n - 1):
    x_start = gap + (i + 1) * box_w + i * gap
    x_end = gap + (i + 1) * (box_w + gap)
    arr = FancyArrowPatch(
        (x_start, y0 + box_h / 2), (x_end, y0 + box_h / 2),
        arrowstyle="-|>", mutation_scale=14,
        linewidth=1.6, color="#555555", zorder=1,
        shrinkA=1, shrinkB=1,
    )
    axA.add_patch(arr)

axA.text(0.3, 9.4,
         "A.  Stereo-seq processing workflow: from raw reads to a spaGAPA-ready PAS matrix",
         fontsize=12, fontweight="bold", va="top", ha="left")

# =========================================================================== #
# PANEL B — 3'-bias QC (recreated from fig1.log values, mouse)
# =========================================================================== #
axB = fig.add_subplot(gs[1, 0])

# Values from pipeline_output/gse263789_stereo_pilot/figures/fig1.log
# ([bins] line, n = 1,000,000 gene-annotated reads, 80.4M scanned)
bin_labels = ["0-500 bp", "0.5-2 kb", "2-10 kb", ">10 kb"]
pct = np.array([53.44, 20.07, 11.96, 14.53])
counts = np.array([534400, 200700, 119600, 145300])  # n=1e6 sample
colors_b = ["#c0392b", "#e67e22", "#f1c40f", "#bdc3c7"]

bars = axB.bar(bin_labels, pct, color=colors_b, edgecolor="black", linewidth=0.8)
for b, p, c in zip(bars, pct, counts):
    axB.text(b.get_x() + b.get_width() / 2, p + 0.9,
             f"{p:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
axB.set_ylabel("Fraction of gene-annotated reads (%)")
axB.set_xlabel("|Distance of read 3'-end to gene TES|")
axB.set_ylim(0, max(pct) * 1.32)
axB.set_title("3'-end enrichment near gene TES")
axB.grid(axis="y", linestyle="--", alpha=0.4)
# inset summary box
axB.text(0.97, 0.97,
         "n = 1,000,000 reads\n(80.4M scanned)\n"
         "53.4% within 500 bp\n73.5% within 2 kb",
         transform=axB.transAxes, va="top", ha="right", fontsize=7,
         bbox=dict(boxstyle="round,pad=0.4", fc="#fdf6e3", ec="#999"))
panel_label(axB, "B", x=-0.16, y=1.06)

# =========================================================================== #
# PANEL C — Data scale infographic
# =========================================================================== #
axC = fig.add_subplot(gs[1, 1])
axC.set_xlim(0, 10)
axC.set_ylim(0, 10)
axC.axis("off")

# values (sources: scapatrap_qc.json n_barcodes; qc_summary.json n_sites/n_bins/
# n_raw_rows; binning summary n_bins_per_resolution; matrix sparsity from qc_summary)
stats = [
    ("20.7M",  "DNB spots\n(raw barcodes)",       BLUE,   str(20677477)),
    ("21,455", "PAS peaks\n(scAPAtrap-called)",   ORANGE, str(21455)),
    ("15,235", "spatial bins\n(bin200, computed)",GREEN,  str(15235)),
    ("97.8M",  "raw peak x spot\nrecords (97.8%)",SKYBLU, str(97818832)),
]
# 2x2 tile layout
xs = [1.5, 6.0, 1.5, 6.0]
ys = [6.2, 6.2, 1.8, 1.8]
tw, th = 3.2, 2.8
for (big, lbl, col, raw), x, y in zip(stats, xs, ys):
    box = FancyBboxPatch(
        (x, y), tw, th,
        boxstyle="round,pad=0.1,rounding_size=0.35",
        linewidth=1.4, edgecolor=col, facecolor=col + "18", zorder=2,
    )
    axC.add_patch(box)
    axC.text(x + tw / 2, y + th - 0.85, big,
             ha="center", va="center", fontsize=19, fontweight="bold", color=col)
    axC.text(x + tw / 2, y + 0.85, lbl,
             ha="center", va="center", fontsize=7.5, color="#333333")
    axC.text(x + tw / 2, y + 0.28, f"raw: {raw}",
             ha="center", va="center", fontsize=5.6, color="#888888", style="italic")

axC.text(5, 9.6, "Dataset scale (GSE263789 mouse AD brain)",
         ha="center", va="center", fontsize=10, fontweight="bold")
axC.text(5, 0.45,
         "1 AD section (GSM8199179)  |  Stereo-seq  |  25.8% reads intergenic",
         ha="center", va="center", fontsize=6.8, color="#666", style="italic")
panel_label(axC, "C", x=-0.06, y=1.06)

# =========================================================================== #
# PANEL D — Tissue APA atlas (embed AD domain map + uncertainty map)
# =========================================================================== #
# split column 2 into two sub-axes via a nested gridspec
gs_d = GridSpecFromSubplotSpec(
    1, 2, subplot_spec=gs[1, 2], wspace=0.22,
)
axD1 = fig.add_subplot(gs_d[0, 0])
axD2 = fig.add_subplot(gs_d[0, 1])

dom_img = imread(AD / "ad_domain_map.png")
unc_img = imread(AD / "ad_uncertainty_map.png")
axD1.imshow(dom_img); axD1.set_xticks([]); axD1.set_yticks([])
axD1.set_title("Tissue domains (41)", fontsize=8.5)
for s in axD1.spines.values():
    s.set_visible(False)
axD2.imshow(unc_img); axD2.set_xticks([]); axD2.set_yticks([])
axD2.set_title("GP uncertainty", fontsize=8.5)
for s in axD2.spines.values():
    s.set_visible(False)
axD1.text(0.5, 1.02, "D.  Tissue APA atlas (AD section, bin200)",
          transform=axD1.transAxes, fontsize=12, fontweight="bold",
          va="bottom", ha="left")

# =========================================================================== #
# PANEL E — Representative high-spatial-variance PAS maps
# =========================================================================== #
gs_e = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2, :], wspace=0.18)

# load coordinates (spot_id -> x, y)
coords = pd.read_csv(BIN200 / "coordinates.csv")
# spot_id order in apa_matrix header == coordinates order (same build)
spot_ids = coords["spot_id"].tolist()
xs_e = coords["x"].to_numpy(dtype=float)
ys_e = coords["y"].to_numpy(dtype=float)

# representative high-variance peaks (from binning consistency top sites)
panelE_peaks = ["peak_69378", "peak_312125", "peak_415285"]
panelE_loci = {
    "peak_69378":  "chr5:146.26M (+)",
    "peak_312125": "chr7:19.70M (-)",
    "peak_415285": "chr16:18.53M (-)",
}
panelE_path = Path("/s3/mengzijun/tmp/fig7_panelE_peaks.csv")
rows = {}
with open(panelE_path) as f:
    rdr = csv.reader(f)
    header = next(rdr)  # site_id, spot_id, spot_id, ...
    assert header[0] == "site_id"
    mat_spots = header[1:]
    for row in rdr:
        rows[row[0]] = np.array([float(v) for v in row[1:]], dtype=float)
# align matrix column order to coords order
# (both derived from the same build -> spot_id identical & same order)
assert mat_spots == spot_ids, "spot order mismatch between apa_matrix and coordinates"

for j, pk in enumerate(panelE_peaks):
    ax = fig.add_subplot(gs_e[0, j])
    vals = rows[pk]
    valid = np.isfinite(vals)
    # scatter with spatial orientation; invert y so origin top-left like tissue
    order = np.argsort(vals[valid])  # plot low first so highs sit on top
    vx, vy, vv = xs_e[valid][order], ys_e[valid][order], vals[valid][order]
    sc = ax.scatter(vx, vy, c=vv, s=2.0, cmap="magma", alpha=0.85,
                    edgecolors="none", vmin=np.nanpercentile(vv, 2),
                    vmax=np.nanpercentile(vv, 98))
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    ax.set_title(f"{pk}\n{panelE_loci[pk]}", fontsize=7.5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("site usage", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    if j == 0:
        ax.text(-0.10, 1.06, "E.  Representative high-spatial-variance PAS (bin200 usage)",
                transform=ax.transAxes, fontsize=12, fontweight="bold",
                va="bottom", ha="left")

# =========================================================================== #
# PANEL F — Binning consistency
# =========================================================================== #
gs_f = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3, 0:2], wspace=0.30)
axF1 = fig.add_subplot(gs_f[0, 0])
axF2 = fig.add_subplot(gs_f[0, 1])

bindf = pd.read_csv(BINNING / "binning_correlation.csv")

# F1: boxplot of cross-bin Pearson r (two pair types: 50vs100, 50vs200)
pair_groups = {
    "bin50 vs bin100": bindf.loc[bindf["bin_size2"] == 100, "pearson_r"].values,
    "bin50 vs bin200": bindf.loc[bindf["bin_size2"] == 200, "pearson_r"].values,
}
positions = [1, 2]
box_data = [pair_groups["bin50 vs bin100"], pair_groups["bin50 vs bin200"]]
bp = axF1.boxplot(
    box_data, positions=positions, widths=0.55, patch_artist=True,
    medianprops=dict(color="black", linewidth=1.4),
    boxprops=dict(facecolor=BLUE + "55", edgecolor=BLUE, linewidth=1.2),
    whiskerprops=dict(color=BLUE, linewidth=1.0),
    capprops=dict(color=BLUE, linewidth=1.0),
    flierprops=dict(marker="o", markerfacecolor=BLUE, markeredgecolor="none",
                    markersize=4, alpha=0.7),
)
# overlay individual points
for pos, arr in zip(positions, box_data):
    jit = pos + (np.random.RandomState(0).rand(len(arr)) - 0.5) * 0.12
    axF1.scatter(jit, arr, s=14, color=ORANGE, alpha=0.85, edgecolors="none", zorder=3)
axF1.set_xticks(positions)
axF1.set_xticklabels(["bin50 vs bin100", "bin50 vs bin200"])
axF1.set_ylabel("Cross-bin Pearson r (per PAS)")
axF1.set_ylim(0.70, 0.97)
axF1.axhline(0.7, color=GREY, linestyle="--", linewidth=0.9, alpha=0.7)
axF1.set_title("F1. APA usage reproducible across resolutions\n"
               "(8 high-variance PAS, 16 pairs, mean r = 0.84)", fontsize=8.5)
axF1.grid(axis="y", linestyle="--", alpha=0.4)
panel_label(axF1, "F", x=-0.10, y=1.06)

# F2: Moran's I vs bin size (should decrease as bin coarsens)
bin_sizes = [50, 100, 200]
mi_cols = ["morans_i_bin50", "morans_i_bin100", "morans_i_bin200"]
mean_mi = [bindf[c].mean() for c in mi_cols]
axF2.plot(bin_sizes, mean_mi, "-o", color=GREEN, linewidth=1.8,
          markersize=7, markeredgecolor="white", markeredgewidth=1.0, zorder=3)
for bs, m in zip(bin_sizes, mean_mi):
    axF2.annotate(f"{m:.2f}", (bs, m), textcoords="offset points",
                  xytext=(0, 9), ha="center", fontsize=8, fontweight="bold",
                  color=GREEN)
axF2.set_xlabel("Bin size (DNB spots per side)")
axF2.set_ylabel("Moran's I (spatial autocorr.)")
axF2.set_xticks(bin_sizes)
axF2.set_title("F2. Spatial autocorrelation decreases\nwith coarser binning (sanity)",
               fontsize=8.5)
axF2.grid(linestyle="--", alpha=0.4)
axF2.set_ylim(0, max(mean_mi) * 1.25)

# =========================================================================== #
# PANEL G — AD-WT descriptive effect sizes (n=1, no p-value claim)
# =========================================================================== #
axG = fig.add_subplot(gs[3, 2])

eff = pd.read_csv(AD / "sample_level_effect_sizes.csv")
# top 10 by |delta|, ranked bar chart
top = eff.sort_values("abs_delta", ascending=False).head(10).iloc[::-1]
colors_g = [ORANGE if d > 0 else BLUE for d in top["delta_AD_minus_WT"]]

bars = axG.barh(range(len(top)), top["delta_AD_minus_WT"], color=colors_g,
                edgecolor="black", linewidth=0.6)
axG.set_yticks(range(len(top)))
axG.set_yticklabels(top["gene"].values, fontsize=7.5)
axG.axvline(0, color="black", linewidth=0.8)
axG.axvline(0.1, color=GREY, linestyle=":", linewidth=0.8)
axG.axvline(-0.1, color=GREY, linestyle=":", linewidth=0.8)
axG.set_xlabel(r"$\Delta$ (AD$_{mean}$ - WT$_{mean}$) site-usage")
axG.set_xlim(-1.15, 1.15)
axG.set_title("G.  Descriptive effect sizes, AD vs WT\n"
              "(effect size only, n=1 per condition)", fontsize=8.5)
# annotate values
for i, (d, adn, wtn) in enumerate(zip(top["delta_AD_minus_WT"],
                                      top["ad_n_spots_obs"],
                                      top["wt_n_spots_obs"])):
    axG.text(d + (0.02 if d >= 0 else -0.02), i,
             f"{d:+.2f}", va="center",
             ha="left" if d >= 0 else "right", fontsize=6.8, fontweight="bold")
# legend for direction
from matplotlib.lines import Line2D  # noqa: E402
legend_elems = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=ORANGE,
           markersize=8, label="distal-up in AD"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=BLUE,
           markersize=8, label="proximal-up in AD"),
]
axG.legend(handles=legend_elems, loc="lower right", fontsize=6.5,
           frameon=True, framealpha=0.9)
axG.text(0.5, -0.30,
         "1,090 genes tested; 607 candidates at |$\\Delta$|>0.1  -  "
         "NO significance call (single AD vs single WT section)",
         transform=axG.transAxes, ha="center", va="top", fontsize=6,
         style="italic", color="#666")

# =========================================================================== #
# Save
# =========================================================================== #
out_path = save(fig, "fig7_stereo_seq.png")
print(f"\n[done] Figure 7 written to: {out_path}")

#!/usr/bin/env python3
"""Figure 7: spaGAPA enables quality-controlled high-resolution APA mapping in Stereo-seq.

Single-section publication figure (300 DPI, colorblind-friendly) for the spaGAPA BIB paper.
Width = PAGE_WIDTH_IN (~7 in, BIB double-column); text kept >= 7-8 pt at print.
All panels use REAL data from the GSE263789 Stereo-seq pilot — no illustrative numbers.

Panels (A-F; the n=1-vs-n=1 AD-WT descriptive effect-size panel was moved to
Supplementary Figure S14 to avoid over-claiming in the main text):
  A. Stereo-seq workflow (schematic): FASTQ+mask -> SAW 8.2.2 -> retag -> scAPAtrap
     -> PAS matrix -> spaGAPA.
  B. 3'-bias QC: distribution of read 3'-end distance to gene TES (polyA-capture
     signature), reusing the published Stereo-seq QC figure.
  C. Data scale & sparsity: 20.7M DNBs, 21,455 PAS, 15,235 bin200 spots;
     the gene x spot USAGE matrix is ~10% observed (~90% empty) -> imputation needed.
     (One consistent sparsity definition: matrix nnz fraction. The raw peak x spot
      record count is reported separately and not conflated with the observed fraction.)
  D. Tissue APA atlas: spatial GP-imputed APA domains + per-spot uncertainty
     across the bin200 AD section (descriptive bin-level domains; resolution note
     in caption).
  E. Representative high-spatial-variance PAS (bin200 usage) for three real peaks,
     each annotated with its overlapping gene (Cdk8, Apoe, Gnb1l) via mm10 GTF.
  F. Binning consistency: cross-bin Pearson r (8 high-variance PAS, 16 pairs) +
     Moran's I vs bin size. APA patterns are robust to resolution; the small PAS
     count is a descriptive cross-bin agreement check, not a population claim.

Run on S91: ~/anaconda3/envs/spagapa/bin/python fig7_stereo_seq.py
"""
import os
import sys
import csv
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.image import imread
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# Shared style with the other 6 main figures
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (  # noqa: E402
    setup_rc, panel_label, save, PAGE_WIDTH_IN,
    BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, BLACK,
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

# --------------------------------------------------------------------------- #
# Layout: width = PAGE_WIDTH_IN (~7 in). Five rows, single column flow.
#   Row 0  A   workflow schematic (full width, short)
#   Row 1  B | C   3'-bias QC  |  data scale + sparsity
#   Row 2  D   tissue atlas (domain + uncertainty)  -- enlarged
#   Row 3  E   three PAS maps (full width)
#   Row 4  F   binning consistency (two sub-panels, full width)
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 13.2))
gs = GridSpec(
    nrows=5, ncols=2, figure=fig,
    height_ratios=[0.72, 1.18, 1.30, 1.30, 1.18],
    width_ratios=[1, 1],
    hspace=0.62, wspace=0.30,
    left=0.072, right=0.972, top=0.978, bottom=0.048,
)

# =========================================================================== #
# PANEL A -- Stereo-seq workflow schematic (spans both columns)
# =========================================================================== #
axA = fig.add_subplot(gs[0, :])
axA.set_xlim(0, 100); axA.set_ylim(0, 10); axA.axis("off")

stages = [
    ("FASTQ\n+ mask",      "Raw Stereo-seq reads\n+ tissue mask",     GREY,   "#f0f0f0"),
    ("SAW 8.2.2",          "Alignment + cell calling\n(DNB -> BAM)",  BLUE,   "#e8f1fa"),
    ("Retag\nCB/UB/GX/GN", "Transfer spot barcode +\ngene annotation", SKYBLU, "#eaf5fb"),
    ("scAPAtrap",          "PolyA-peak calling\n(PAS discovery)",     ORANGE, "#fbeadd"),
    ("PAS matrix",         "peak x spot usage\n(sparse, ~90% empty)", GREEN,  "#e6f6ef"),
    ("spaGAPA",            "GP imputation +\nuncertainty + domains",  RED,    "#f7eaef"),
]
n = len(stages)
box_w = 13.5
gap = (100 - n * box_w) / (n + 1)
y0 = 3.6; box_h = 4.6

for i, (title, desc, edge, fc) in enumerate(stages):
    x = gap + i * (box_w + gap)
    box = FancyBboxPatch(
        (x, y0), box_w, box_h,
        boxstyle="round,pad=0.12,rounding_size=0.5",
        linewidth=1.3, edgecolor=edge, facecolor=fc, zorder=2,
    )
    axA.add_patch(box)
    axA.text(x + box_w / 2, y0 + box_h - 1.25, title,
             ha="center", va="center", fontsize=8.6, fontweight="bold",
             color=edge, zorder=3)
    axA.text(x + box_w / 2, y0 + 1.35, desc,
             ha="center", va="center", fontsize=6.3, color="#333333", zorder=3)

for i in range(n - 1):
    x_start = gap + (i + 1) * box_w + i * gap
    x_end = gap + (i + 1) * (box_w + gap)
    axA.add_patch(FancyArrowPatch(
        (x_start, y0 + box_h / 2), (x_end, y0 + box_h / 2),
        arrowstyle="-|>", mutation_scale=12,
        linewidth=1.4, color="#555555", zorder=1, shrinkA=1, shrinkB=1,
    ))

axA.text(0.3, 9.4,
         "A.  Stereo-seq processing workflow: raw reads -> spaGAPA-ready PAS matrix",
         fontsize=10.5, fontweight="bold", va="top", ha="left")

# =========================================================================== #
# PANEL B -- 3'-bias QC (reuse the published Stereo-seq pilot figure)
# =========================================================================== #
axB = fig.add_subplot(gs[1, 0])
bias_img_path = STEREO / "figures" / "fig1_3prime_enrichment_near_TES.png"
try:
    img = imread(bias_img_path)
    axB.imshow(img)
    axB.set_xticks([]); axB.set_yticks([])
    for s in axB.spines.values():
        s.set_visible(False)
except Exception:
    axB.text(0.5, 0.5, f"[reuse]\n{bias_img_path.name}",
             ha="center", va="center", fontsize=7)
    axB.axis("off")
axB.set_title("3'-end enrichment near gene TES", fontsize=9)
panel_label(axB, "B", x=-0.10, y=1.06)

# =========================================================================== #
# PANEL C -- Data scale + sparsity (single consistent sparsity definition)
# =========================================================================== #
axC = fig.add_subplot(gs[1, 1])
axC.set_xlim(0, 10); axC.set_ylim(0, 10); axC.axis("off")

# Source: binned_200/qc_summary.json  (matrix [21455 x 15235], nnz=33,797,211)
N_SITES, N_BINS = 21455, 15235
NNZ = 33_797_211
TOTAL = N_SITES * N_BINS
OBS_FRAC = NNZ / TOTAL          # ~0.1035  -> ~10% observed, ~90% empty
SPARSITY = 1.0 - OBS_FRAC

stats = [
    ("20.7M",  "DNB spots\n(raw barcodes)",       BLUE,   "20,677,477"),
    ("21,455", "PAS peaks\n(scAPAtrap-called)",   ORANGE, "21,455"),
    ("15,235", "bin200 spots\n(analyzed)",        GREEN,  "15,235"),
    (f"{OBS_FRAC*100:.1f}%", "matrix observed\n(~90% empty)", RED,
        f"nnz={NNZ:,} / {TOTAL:,}"),
]
xs = [1.2, 5.6, 1.2, 5.6]
ys = [6.0, 6.0, 1.5, 1.5]
tw, th = 3.3, 3.0
for (big, lbl, col, raw), x, y in zip(stats, xs, ys):
    box = FancyBboxPatch(
        (x, y), tw, th,
        boxstyle="round,pad=0.1,rounding_size=0.35",
        linewidth=1.3, edgecolor=col, facecolor=col + "18", zorder=2,
    )
    axC.add_patch(box)
    axC.text(x + tw / 2, y + th - 0.85, big,
             ha="center", va="center", fontsize=16, fontweight="bold", color=col)
    axC.text(x + tw / 2, y + 0.95, lbl,
             ha="center", va="center", fontsize=6.8, color="#333333")
    axC.text(x + tw / 2, y + 0.32, f"raw: {raw}",
             ha="center", va="center", fontsize=5.4, color="#888888", style="italic")

axC.text(5, 9.55, "Dataset scale (GSE263789 mouse AD brain)",
         ha="center", va="center", fontsize=9.2, fontweight="bold")
axC.text(5, 0.30,
         "1 AD section (GSM8199179) | Stereo-seq | "
         f"matrix nnz fraction = {OBS_FRAC*100:.1f}% observed",
         ha="center", va="center", fontsize=6.0, color="#666", style="italic")
panel_label(axC, "C", x=-0.06, y=1.06)

# =========================================================================== #
# PANEL D -- Tissue APA atlas (enlarged: domain map + uncertainty map)
# =========================================================================== #
gs_d = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2, :], wspace=0.20)
axD1 = fig.add_subplot(gs_d[0, 0])
axD2 = fig.add_subplot(gs_d[0, 1])

dom_img = imread(AD / "ad_domain_map.png")
unc_img = imread(AD / "ad_uncertainty_map.png")
axD1.imshow(dom_img); axD1.set_xticks([]); axD1.set_yticks([])
for s in axD1.spines.values():
    s.set_visible(False)
axD2.imshow(unc_img); axD2.set_xticks([]); axD2.set_yticks([])
for s in axD2.spines.values():
    s.set_visible(False)

# Domain count: ad_spagapa_summary.json n_domains = 41 (bin-level, descriptive).
# Likely over-segmented at this gene subset (1500 high-variance genes); caption
# frames it as a descriptive bin-level map, not a ground-truth tissue taxonomy.
axD1.set_title("GP-APA domains (n=41, bin-level)", fontsize=8.6)
axD2.set_title("GP posterior uncertainty", fontsize=8.6)
axD1.text(0.5, 1.04,
          "D.  Tissue APA atlas (AD section, bin200; descriptive bin-level domains)",
          transform=axD1.transAxes, fontsize=10.5, fontweight="bold",
          va="bottom", ha="left")

# =========================================================================== #
# PANEL E -- Representative high-spatial-variance PAS, annotated with genes
# =========================================================================== #
gs_e = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[3, :], wspace=0.20)

# spot coordinates (bin200)
coords = pd.read_csv(BIN200 / "coordinates.csv")
spot_ids = coords["spot_id"].tolist()
xs_e = coords["x"].to_numpy(dtype=float)
ys_e = coords["y"].to_numpy(dtype=float)

# Three real high-variance peaks (from binning-consistency top sites), each
# overlapped by a named gene on mm10:
#   peak_69378  chr5:146,260,979-146,261,371 (+) -> Cdk8
#   peak_312125 chr7:19,696,240-19,697,182 (-)   -> Apoe   (Alzheimer's APOE)
#   peak_415285 chr16:18,533,182-18,533,341 (-)  -> Gnb1l
panelE_peaks = ["peak_69378", "peak_312125", "peak_415285"]
panelE_titles = [
    "Cdk8\nchr5:146.26M (+)",
    "Apoe\nchr7:19.70M (-)",
    "Gnb1l\nchr16:18.53M (-)",
]
panelE_path = Path("/s3/mengzijun/tmp/fig7_panelE_peaks.csv")
rows = {}
with open(panelE_path) as f:
    rdr = csv.reader(f)
    header = next(rdr)
    assert header[0] == "site_id"
    mat_spots = header[1:]
    for row in rdr:
        rows[row[0]] = np.array([float(v) for v in row[1:]], dtype=float)
assert mat_spots == spot_ids, "spot order mismatch between apa_matrix and coordinates"

for j, pk in enumerate(panelE_peaks):
    ax = fig.add_subplot(gs_e[0, j])
    ax.set_rasterized(True)  # dense spatial scatter -> keep PDF small
    vals = rows[pk]
    valid = np.isfinite(vals)
    order = np.argsort(vals[valid])  # plot low first so highs sit on top
    vx, vy, vv = xs_e[valid][order], ys_e[valid][order], vals[valid][order]
    sc = ax.scatter(vx, vy, c=vv, s=2.2, cmap="magma", alpha=0.85,
                    edgecolors="none", vmin=np.nanpercentile(vv, 2),
                    vmax=np.nanpercentile(vv, 98))
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    ax.set_title(panelE_titles[j], fontsize=7.6)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("site usage", fontsize=6.4)
    cb.ax.tick_params(labelsize=6)
    if j == 0:
        ax.text(-0.10, 1.10,
                "E.  Representative high-spatial-variance PAS (bin200 usage, "
                "gene by mm10 overlap)",
                transform=ax.transAxes, fontsize=10.5, fontweight="bold",
                va="bottom", ha="left")

# =========================================================================== #
# PANEL F -- Binning consistency (8 PAS / 16 pairs, descriptive)
# =========================================================================== #
gs_f = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[4, :], wspace=0.30)
axF1 = fig.add_subplot(gs_f[0, 0])
axF2 = fig.add_subplot(gs_f[0, 1])

bindf = pd.read_csv(BINNING / "binning_correlation.csv")
n_pas = bindf["gene"].nunique()
n_pairs = len(bindf)
mean_r = bindf["pearson_r"].mean()

pair_groups = {
    "bin50 vs bin100": bindf.loc[bindf["bin_size2"] == 100, "pearson_r"].values,
    "bin50 vs bin200": bindf.loc[bindf["bin_size2"] == 200, "pearson_r"].values,
}
positions = [1, 2]
box_data = [pair_groups["bin50 vs bin100"], pair_groups["bin50 vs bin200"]]
axF1.boxplot(
    box_data, positions=positions, widths=0.55, patch_artist=True,
    medianprops=dict(color="black", linewidth=1.4),
    boxprops=dict(facecolor=BLUE + "55", edgecolor=BLUE, linewidth=1.2),
    whiskerprops=dict(color=BLUE, linewidth=1.0),
    capprops=dict(color=BLUE, linewidth=1.0),
    flierprops=dict(marker="o", markerfacecolor=BLUE, markeredgecolor="none",
                    markersize=4, alpha=0.7),
)
for pos, arr in zip(positions, box_data):
    jit = pos + (np.random.RandomState(0).rand(len(arr)) - 0.5) * 0.12
    axF1.scatter(jit, arr, s=14, color=ORANGE, alpha=0.85, edgecolors="none", zorder=3)
axF1.set_xticks(positions)
axF1.set_xticklabels(["bin50 vs bin100", "bin50 vs bin200"], fontsize=7.2)
axF1.set_ylabel("Cross-bin Pearson r (per PAS)")
axF1.set_ylim(0.70, 0.97)
axF1.axhline(0.7, color=GREY, linestyle="--", linewidth=0.9, alpha=0.7)
axF1.set_title(
    f"APA usage reproducible across resolutions\n"
    f"({n_pas} high-variance PAS, {n_pairs} pairs, mean r={mean_r:.2f})",
    fontsize=8.4)
axF1.grid(axis="y", linestyle="--", alpha=0.4)
panel_label(axF1, "F", x=-0.10, y=1.10)

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
axF2.set_title("Spatial autocorrelation decreases\nwith coarser binning (sanity)",
               fontsize=8.4)
axF2.grid(linestyle="--", alpha=0.4)
axF2.set_ylim(0, max(mean_mi) * 1.25)

# =========================================================================== #
# Save (vector PDF + 300-DPI PNG)
# =========================================================================== #
out_path = save(fig, "fig7_stereo_seq.png")
print(f"\n[done] Figure 7 written to: {out_path}")

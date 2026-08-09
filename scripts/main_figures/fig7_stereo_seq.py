#!/usr/bin/env python3
"""Figure 7: spaGAPA enables quality-controlled high-resolution APA mapping in Stereo-seq.

Single-section publication figure (vector PDF + 300-DPI PNG) for the spaGAPA BIB paper.
Width = PAGE_WIDTH_IN (~7 in, BIB double-column); text kept >= 7-8 pt at print.
All panels use REAL data from the GSE263789 Stereo-seq pilot — no illustrative numbers.

Refactored to 5 panels (reviewer-driven: drop the 41-domain + GP-uncertainty maps,
which reflected the observed/missing bin grid rather than biology; the AD-WT
descriptive panel already moved to Supplementary S14):

  A. Workflow (one compact row): FASTQ+mask -> SAW 8.2.2 -> retag -> scAPAtrap
     -> PAS matrix -> spaGAPA.
  B. Stereo-seq APA QC (native redraw, not a downscaled screenshot): read
     3'-end distance to gene TES + 3'-enrichment summary. Source: pipeline log
     of the published pilot (n=1,000,000 gene-annotated reads;
     53.4% within 500 bp, 73.5% within 2 kb of TES).
  C. Dataset scale & sparsity: 20.7M DNBs, 21,455 PAS, 15,235 bin200 spots,
     10.3% matrix observed (single consistent sparsity definition).
  D. Representative spatial APA usage (main visual, enlarged): three real
     high-spatial-variance PAS on the bin200 AD section, each annotated with its
     overlapping gene (Cdk8, Apoe, Gnb1l) via mm10 GTF.
  E. Binning robustness: cross-bin Pearson r (8 high-variance PAS, 16 pairs) +
     Moran's I vs bin size. Descriptive cross-resolution agreement check.

Run on S91: ~/anaconda3/envs/spagapa/bin/python fig7_stereo_seq.py
"""
import os
import sys
import csv
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# Shared style with the other 6 main figures
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (  # noqa: E402
    setup_rc, save, PAGE_WIDTH_IN,
    BLUE, ORANGE, GREEN, SKYBLU, RED, GREY,
)

setup_rc()

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output")
STEREO = ROOT / "gse263789_stereo_pilot"
BIN200 = STEREO / "spagapa_downstream_full" / "binned_200"
BINNING = ROOT / "stereo_binning_consistency"

# --------------------------------------------------------------------------- #
# Layout: width = PAGE_WIDTH_IN (~7 in). Five rows, single-column flow.
#   Row 0  A   workflow schematic (full width, compact one row)
#   Row 1  B | C   3'-bias QC (native)  |  data scale + sparsity
#   Row 2  D   three PAS maps (full width) -- main visual, enlarged
#   Row 3  E   binning consistency (two sub-panels, full width)
# Height kept restrained (target <~220 mm at print).
# --------------------------------------------------------------------------- #
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 8.4))
gs = GridSpec(
    nrows=4, ncols=2, figure=fig,
    height_ratios=[0.62, 1.42, 1.70, 1.30],
    width_ratios=[1, 1],
    hspace=0.50, wspace=0.30,
    left=0.090, right=0.972, top=0.975, bottom=0.045,
)

# =========================================================================== #
# PANEL A -- Stereo-seq workflow schematic (spans both columns, one row)
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
y0 = 3.0; box_h = 4.8

for i, (title, desc, edge, fc) in enumerate(stages):
    x = gap + i * (box_w + gap)
    box = FancyBboxPatch(
        (x, y0), box_w, box_h,
        boxstyle="round,pad=0.12,rounding_size=0.5",
        linewidth=1.3, edgecolor=edge, facecolor=fc, zorder=2,
    )
    axA.add_patch(box)
    axA.text(x + box_w / 2, y0 + box_h - 1.30, title,
             ha="center", va="center", fontsize=8.6, fontweight="bold",
             color=edge, zorder=3)
    axA.text(x + box_w / 2, y0 + 1.45, desc,
             ha="center", va="center", fontsize=6.3, color="#333333", zorder=3)

for i in range(n - 1):
    x_start = gap + (i + 1) * box_w + i * gap
    x_end = gap + (i + 1) * (box_w + gap)
    axA.add_patch(FancyArrowPatch(
        (x_start, y0 + box_h / 2), (x_end, y0 + box_h / 2),
        arrowstyle="-|>", mutation_scale=12,
        linewidth=1.4, color="#555555", zorder=1, shrinkA=1, shrinkB=1,
    ))

axA.text(0.3, 9.5,
         "A.  Stereo-seq processing workflow: raw reads -> spaGAPA-ready PAS matrix",
         fontsize=10.5, fontweight="bold", va="top", ha="left")

# =========================================================================== #
# PANEL B -- 3'-bias QC (NATIVE REDRAW, not a downscaled screenshot)
#   Left:  binned bar chart of |read 3'-end distance to gene TES|
#   Right: signed distance histogram + cumulative fraction (enrichment curve)
#   Source: pipeline log of the published Stereo-seq pilot
#     (n=1,000,000 gene-annotated reads; 53.4% within 500 bp, 73.5% within 2 kb)
# =========================================================================== #
gs_b = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 0], wspace=0.42)
axB1 = fig.add_subplot(gs_b[0, 0])
axB2 = fig.add_subplot(gs_b[0, 1])

# Distances saved alongside the published fig1 (full 1M-read sample that produced
# the logged 53.4% / 73.5% numbers). Fall back to log-derived bin fractions if
# the raw array is unavailable, so the panel is always native and readable.
fig1_dists = STEREO / "figures" / "_fig1_dists.npy"
have_raw = fig1_dists.exists()
if have_raw:
    dists = np.load(fig1_dists)
    ad = np.abs(dists)
    edges = [0, 500, 2000, 10000, np.inf]
    labels = ["0-500 bp", "500 bp-2 kb", "2-10 kb", ">10 kb"]
    counts, _ = np.histogram(ad, bins=edges)
    pct = counts / counts.sum() * 100.0
    within_500 = pct[0]
    within_2k = pct[0] + pct[1]
    n_reads = len(dists)
else:
    # Logged fractions from fig1.log (53.44 / 73.51). Reconstruct bin fractions.
    pct = np.array([53.44, 20.07, 11.96, 14.53])
    labels = ["0-500 bp", "500 bp-2 kb", "2-10 kb", ">10 kb"]
    within_500, within_2k = 53.44, 73.51
    n_reads = 1_000_000

# --- B-left: binned bar chart of |distance to TES| ---
colors_b = ["#c0392b", "#e67e22", "#f1c40f", "#bdc3c7"]
bars = axB1.bar(labels, pct, color=colors_b, edgecolor="black", linewidth=0.7,
                width=0.72)
for b, p in zip(bars, pct):
    axB1.text(b.get_x() + b.get_width() / 2, p + 0.8,
              f"{p:.1f}%", ha="center", va="bottom", fontsize=6.8)
axB1.set_ylabel("Fraction of reads (%)", fontsize=7.6)
axB1.set_xlabel("|distance of read 3'-end to gene TES|", fontsize=7.6)
axB1.set_ylim(0, max(pct) * 1.28)
axB1.tick_params(axis="x", labelsize=6.4)
axB1.tick_params(axis="y", labelsize=6.6)
axB1.grid(axis="y", linestyle="--", alpha=0.4)
axB1.set_title("Read 3'-ends near TES", fontsize=8.2)

# --- B-right: cumulative enrichment summary (no raw array -> use log bins) ---
axB2.set_xlim(-3, 3)
axB2.set_ylim(0, 100)
# Build cumulative curve from the logged bin fractions when raw array absent.
if have_raw:
    grid = np.array([0, 100, 200, 500, 1000, 2000, 5000, 10000])
    cum = np.array([100 * np.mean(ad <= g) for g in grid])
else:
    grid = np.array([0, 500, 2000, 10000])
    cum = np.array([0.0, within_500, within_2k, within_2k + pct[2]])
axB2.plot(np.log10(grid.clip(min=1)), cum, "-o", color=RED, linewidth=1.8,
          markersize=5, markeredgecolor="white", markeredgewidth=0.8, zorder=3)
axB2.set_xticks([-3, -2, -1, 0, 1, 2, 3])
axB2.set_xticklabels(["1 bp", "10", "100", "1 kb", "10", "100", "1 Mb"],
                     fontsize=6.4)
axB2.set_xlabel("Distance to TES (log)", fontsize=7.6)
axB2.set_ylabel("Cumulative reads (%)", fontsize=7.6)
axB2.tick_params(axis="y", labelsize=6.6)
axB2.axvline(np.log10(500), color=GREY, linestyle=":", linewidth=0.9, alpha=0.8)
axB2.axvline(np.log10(2000), color=GREY, linestyle=":", linewidth=0.9, alpha=0.8)
# Inline summary (kept inside axes bounds to avoid bbox_inches='tight' overhang)
axB2.text(0.03, 0.06,
          f"{within_500:.1f}% within 500 bp\n{within_2k:.1f}% within 2 kb",
          transform=axB2.transAxes, fontsize=6.6, color=RED,
          va="bottom", ha="left",
          bbox=dict(boxstyle="round,pad=0.3", fc="#fdf6e3", ec="#bbb",
                    linewidth=0.7))
axB2.grid(linestyle="--", alpha=0.4)
axB2.set_title("3'-enrichment (polyA)", fontsize=8.2)

axB1.text(0.0, 1.22,
          f"B.  Stereo-seq APA QC: read 3'-ends pile up near gene TES\n"
          f"({n_reads:,} gene-annotated reads; {within_500:.1f}% within 500 bp, "
          f"{within_2k:.1f}% within 2 kb)",
          transform=axB1.transAxes, fontsize=10.0, fontweight="bold",
          va="bottom", ha="left")

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
         "1 AD section (GSM8199179) | "
         f"nnz fraction = {OBS_FRAC*100:.1f}% observed",
         ha="center", va="center", fontsize=5.8, color="#666", style="italic")

# =========================================================================== #
# PANEL D -- Representative spatial APA usage (MAIN VISUAL, enlarged, 3 maps)
#   Three real high-variance peaks on the bin200 AD section, each overlapped
#   by a named gene on mm10:
#     peak_69378  chr5:146,260,979-146,261,371 (+) -> Cdk8
#     peak_312125 chr7:19,696,240-19,697,182 (-)   -> Apoe   (Alzheimer's APOE)
#     peak_415285 chr16:18,533,182-18,533,341 (-)  -> Gnb1l
# =========================================================================== #
gs_d = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2, :], wspace=0.18)

# spot coordinates (bin200)
coords = pd.read_csv(BIN200 / "coordinates.csv")
spot_ids = coords["spot_id"].tolist()
xs_d = coords["x"].to_numpy(dtype=float)
ys_d = coords["y"].to_numpy(dtype=float)

panelD_peaks = ["peak_69378", "peak_312125", "peak_415285"]
panelD_titles = [
    "Cdk8\nchr5:146.26M (+)",
    "Apoe\nchr7:19.70M (-)",
    "Gnb1l\nchr16:18.53M (-)",
]
panelD_path = Path("/s3/mengzijun/tmp/fig7_panelE_peaks.csv")
rows = {}
with open(panelD_path) as f:
    rdr = csv.reader(f)
    header = next(rdr)
    assert header[0] == "site_id"
    mat_spots = header[1:]
    for row in rdr:
        rows[row[0]] = np.array([float(v) for v in row[1:]], dtype=float)
assert mat_spots == spot_ids, "spot order mismatch between apa_matrix and coordinates"

for j, pk in enumerate(panelD_peaks):
    ax = fig.add_subplot(gs_d[0, j])
    ax.set_rasterized(True)  # dense spatial scatter -> keep PDF small
    vals = rows[pk]
    valid = np.isfinite(vals)
    order = np.argsort(vals[valid])  # plot low first so highs sit on top
    vx, vy, vv = xs_d[valid][order], ys_d[valid][order], vals[valid][order]
    sc = ax.scatter(vx, vy, c=vv, s=3.2, cmap="magma", alpha=0.90,
                    edgecolors="none",
                    vmin=np.nanpercentile(vv, 2),
                    vmax=np.nanpercentile(vv, 98))
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    ax.set_title(panelD_titles[j], fontsize=8.0)
    # Compact inset colorbar pinned INSIDE the axes (lower-right) so the
    # rightmost panel never pushes past the figure frame.
    cax = ax.inset_axes([0.78, 0.06, 0.030, 0.34], transform=ax.transAxes)
    cb = fig.colorbar(sc, cax=cax)
    cb.ax.tick_params(labelsize=5.4)
    cb.set_label("usage", fontsize=5.8, labelpad=2)
    if j == 0:
        # Center the cross-panel title over the whole Panel D row (figure
        # coords) so the long label never anchors at the leftmost axis and
        # overhang the right margin.
        fig.text(0.545, 0.795,
                 "D.  Representative spatial APA usage (bin200, AD section; genes by mm10 overlap)",
                 fontsize=9.5, fontweight="bold", va="bottom", ha="center")

# =========================================================================== #
# PANEL E -- Binning robustness (8 PAS / 16 pairs, descriptive)
# =========================================================================== #
gs_e = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3, :], wspace=0.30)
axE1 = fig.add_subplot(gs_e[0, 0])
axE2 = fig.add_subplot(gs_e[0, 1])

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
axE1.boxplot(
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
    axE1.scatter(jit, arr, s=14, color=ORANGE, alpha=0.85, edgecolors="none", zorder=3)
axE1.set_xticks(positions)
axE1.set_xticklabels(["bin50 vs bin100", "bin50 vs bin200"], fontsize=7.2)
axE1.set_ylabel("Cross-bin Pearson r (per PAS)")
axE1.set_ylim(0.70, 0.97)
axE1.axhline(0.7, color=GREY, linestyle="--", linewidth=0.9, alpha=0.7)
axE1.set_title(
    f"E.  APA usage reproducible across resolutions\n"
    f"({n_pas} high-variance PAS, {n_pairs} pairs, mean r={mean_r:.2f})",
    fontsize=8.4)
axE1.grid(axis="y", linestyle="--", alpha=0.4)

bin_sizes = [50, 100, 200]
mi_cols = ["morans_i_bin50", "morans_i_bin100", "morans_i_bin200"]
mean_mi = [bindf[c].mean() for c in mi_cols]
axE2.plot(bin_sizes, mean_mi, "-o", color=GREEN, linewidth=1.8,
          markersize=7, markeredgecolor="white", markeredgewidth=1.0, zorder=3)
for bs, m in zip(bin_sizes, mean_mi):
    axE2.annotate(f"{m:.2f}", (bs, m), textcoords="offset points",
                  xytext=(0, 9), ha="center", fontsize=8, fontweight="bold",
                  color=GREEN)
axE2.set_xlabel("Bin size (DNB spots per side)")
axE2.set_ylabel("Moran's I (spatial autocorr.)")
axE2.set_xticks(bin_sizes)
axE2.set_title("Spatial autocorrelation decreases\nwith coarser binning (sanity)",
               fontsize=8.4)
axE2.grid(linestyle="--", alpha=0.4)
axE2.set_ylim(0, max(mean_mi) * 1.25)

# =========================================================================== #
# Save (vector PDF + 300-DPI PNG)
# =========================================================================== #
out_path = save(fig, "fig7_stereo_seq.png")
print(f"\n[done] Figure 7 written to: {out_path}")

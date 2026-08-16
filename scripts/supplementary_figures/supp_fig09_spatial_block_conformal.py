#!/usr/bin/env python
"""Supplementary Figure S9: Spatial-block vs random-split conformal coverage.

Grouped bar chart of empirical coverage under random splitting vs spatial-block
splitting, per dataset, at 80/90/95%. Confirms split-conformal intervals remain
valid (within ~2 pp of nominal) even under the much harder spatially-disjoint
block split — i.e. coverage is not an artefact of spatial autocorrelation
leaking adjacent spots across the train/test boundary.

Data:
  pipeline_output/spatial_block_conformal/block_split_coverage.csv
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, GREY, save_supp

# PAGE_WIDTH_IN is defined in main_figures/_style.py (the shared BIB page
# geometry constant, ~7 in); mirror it here without modifying _style.py.
PAGE_WIDTH_IN = 7.0

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    df = pd.read_csv(os.path.join(ROOT, "pipeline_output/spatial_block_conformal/block_split_coverage.csv"))
    datasets = df["dataset"].unique()
    print(f"[S9] {len(datasets)} datasets", flush=True)

    # Print-width layout: panels stacked vertically so each keeps the full
    # PAGE_WIDTH_IN (side-by-side squeezed the 8-pt tick labels).
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 6.6))
    gs = fig.add_gridspec(2, 1, hspace=0.52)

    # ---------- Panel A: coverage random vs block, 3 levels ----------
    axA = fig.add_subplot(gs[0, 0])  # row 0 (stacked layout)
    x = np.arange(len(datasets))
    w = 0.13
    levels = [("80%", "coverage_80", 0.80), ("90%", "coverage_90", 0.90), ("95%", "coverage_95", 0.95)]
    split_pairs = [("random", BLUE), ("block", ORANGE)]
    handles_made = False
    for li, (lname, ccol, nom) in enumerate(levels):
        for si, (split, col) in enumerate(split_pairs):
            vals = []
            for ds in datasets:
                sub = df[(df["dataset"] == ds) & (df["split_type"] == split)]
                vals.append(sub[ccol].values[0] if len(sub) else np.nan)
            offset = (li - 1) * (2 * w) + (si - 0.5) * w
            label = f"{lname} {split}" if not handles_made else None
            axA.bar(x + offset, vals, width=w, color=col, edgecolor="white", linewidth=0.35,
                    alpha=0.55 + 0.15 * li)
            handles_made = handles_made or (li == 0 and si == 1)
    for nom in [0.80, 0.90, 0.95]:
        axA.axhline(nom, color="#555", lw=0.7, ls=":")
    # Custom legend: split type.
    from matplotlib.patches import Patch
    axA.legend(handles=[Patch(facecolor=BLUE, label="Random split"),
                        Patch(facecolor=ORANGE, label="Block split"),
                        *([Patch(facecolor=GREY, alpha=a, label=l) for a, l in [(0.55,"80%"),(0.7,"90%"),(0.85,"95%")]])],
              loc="lower right", ncol=2, fontsize=7)
    axA.set_xticks(x); axA.set_xticklabels([d.replace("GSE","") for d in datasets], fontsize=8)
    axA.set_xlabel("Dataset"); axA.set_ylabel("Empirical coverage")
    axA.set_ylim(0.88, 1.005)
    axA.set_title("Random vs spatial-block split: coverage stays valid", loc="left", fontsize=9.5)
    panel_label(axA, "A", x=-0.06, y=1.05)

    # ---------- Panel B: coverage gap (block - random) at 90% ----------
    axB = fig.add_subplot(gs[1, 0])  # row 1 (stacked layout)
    gaps = []
    for ds in datasets:
        rnd = df[(df["dataset"]==ds) & (df["split_type"]=="random")]["coverage_90"].values[0]
        blk = df[(df["dataset"]==ds) & (df["split_type"]=="block")]["coverage_90"].values[0]
        gaps.append(blk - rnd)
    colors = [GREEN if abs(g) < 0.01 else (ORANGE if g > 0 else BLUE) for g in gaps]
    axB.bar(x, [g*100 for g in gaps], color=colors, edgecolor="white", linewidth=0.5)
    axB.axhline(0, color="#333", lw=0.8)
    axB.set_xticks(x); axB.set_xticklabels([d.replace("GSE","") for d in datasets], fontsize=8)
    axB.set_xlabel("Dataset")
    axB.set_ylabel("Coverage gap  (block − random, pp)  @ 90%")
    axB.set_title("Block vs random coverage gap at 90%", loc="left", fontsize=9.5)
    for i, g in enumerate(gaps):
        axB.annotate(f"{g*100:+.2f}", xy=(x[i], g*100), xytext=(0, 4 if g>=0 else -8),
                     textcoords="offset points", ha="center", fontsize=7, fontweight="bold")
    axB.text(0.5, -0.20,
             "Gap near 0 => coverage is not driven by spatial leakage across the\n"
             "train/test boundary. Positive gap = block split is more conservative.",
             transform=axB.transAxes, fontsize=6.8, ha="center", style="italic", color="#444")
    panel_label(axB, "B", x=-0.06, y=1.05)

    fig.suptitle("Supplementary Figure S9 — Spatial-block conformal coverage",
                 fontsize=11, fontweight="bold", y=1.01)
    save_supp(fig, "supp_fig09_spatial_block_conformal.png")
    print("[S9] done", flush=True)


if __name__ == "__main__":
    main()

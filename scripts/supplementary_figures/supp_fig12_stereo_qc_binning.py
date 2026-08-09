#!/usr/bin/env python
"""Supplementary Figure S12: Stereo-seq QC + binning sensitivity.

Panel A: PAS count distribution by chromosome (from scAPAtrap peak calls on the
         Stereo-seq mouse-brain sample).
Panel B: Read 3'-end-to-TES signed-distance histogram (500k gene-annotated
         reads), showing the 3'-enrichment (polyA-capture) signature of
         Stereo-seq.
Panel C: Binning sensitivity — cross-bin Pearson r between bin sizes (50/100,
         50/200) and Moran's I decay with bin size, over 8 high-coverage peaks.

Data:
  pipeline_output/gse263789_stereo_pilot/gsm8199179_full/scapatrap_raw/peaks_meta.csv.gz
    (21,455 PAS — the full Stereo-seq run; matches Fig 7. Previously the
     8,659-PAS pilot subset was used, which conflicted with Fig 7.)
  pipeline_output/supplementary_figures/_cache/s12_tes_distances.npy  (500k reads)
  pipeline_output/stereo_binning_consistency/binning_correlation.csv

Note: binning-correlation.csv only contains the 50→100 and 50→200 comparisons
(no bin500 was run); the panel reports the available pairs with a caption note.
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, save_supp

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    peaks = pd.read_csv(os.path.join(ROOT, "pipeline_output/gse263789_stereo_pilot/gsm8199179_full/scapatrap_raw/peaks_meta.csv.gz"))
    dists = np.load(os.path.join(ROOT, "pipeline_output/supplementary_figures/_cache/s12_tes_distances.npy"))
    binning = pd.read_csv(os.path.join(ROOT, "pipeline_output/stereo_binning_consistency/binning_correlation.csv"))

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.28, height_ratios=[1.0, 1.0])

    # ---------- Panel A: PAS count by chromosome ----------
    axA = fig.add_subplot(gs[0, 0])
    chr_counts = peaks["chr"].value_counts()
    # drop very small scaffolds; keep standard chromosomes.
    chr_counts = chr_counts[chr_counts.index.isin([f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"])]
    order = sorted(chr_counts.index, key=lambda c: (len(c), c))
    order = sorted(order, key=lambda c: int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 99)
    vals = [chr_counts[c] for c in order]
    axA.bar(range(len(order)), vals, color=BLUE, edgecolor="white", lw=0.4)
    axA.set_xticks(range(len(order)))
    axA.set_xticklabels([c.replace("chr", "") for c in order], fontsize=7.5)
    axA.set_xlabel("Chromosome")
    axA.set_ylabel("PAS count (scAPAtrap peaks)")
    axA.set_title(f"PAS count by chromosome (n = {len(peaks):,} peaks)", loc="left", fontsize=9.5)
    axA.yaxis.set_major_locator(MaxNLocator(5))
    panel_label(axA, "A", x=-0.08, y=1.05)

    # ---------- Panel B: 3'-end distance to TES histogram ----------
    axB = fig.add_subplot(gs[0, 1])
    # Clip to a sensible window for visualisation; report stats in caption.
    clip = dists[(dists >= -3000) & (dists <= 3000)]
    bins = np.arange(-3000, 3001, 40)
    axB.hist(clip, bins=bins, color=ORANGE, edgecolor="white", lw=0.2, alpha=0.9)
    axB.axvline(0, color="black", lw=1.2, ls="--", label="TES (gene 3' end)")
    med = float(np.median(dists))
    frac200 = float(np.mean(np.abs(dists) <= 200))
    axB.axvline(med, color=RED, lw=1.4, label=f"median = {med:.0f} bp")
    axB.set_xlabel("Read 3'-end − TES  (signed, bp)")
    axB.set_ylabel("# reads")
    axB.set_title("3'-end enrichment near gene TES (polyA-capture signature)", loc="left", fontsize=9.5)
    axB.legend(loc="upper left", fontsize=7.5)
    axB.text(0.98, 0.92,
             f"n = {dists.size:,} gene-annotated reads\n"
             f"{frac200*100:.0f}% within ±200 bp of TES\n"
             f"median = {med:.0f} bp",
             transform=axB.transAxes, ha="right", va="top", fontsize=7.2,
             bbox=dict(facecolor="white", edgecolor="#bbb", boxstyle="round,pad=0.3"))
    panel_label(axB, "B", x=-0.08, y=1.05)

    # ---------- Panel C: binning sensitivity ----------
    axC = fig.add_subplot(gs[1, :])
    # Cross-bin Pearson r per peak (grouped by pair).
    w = 0.36
    # Two pairs: (50,100) and (50,200).
    p50100 = binning[binning["bin_size2"] == 100].reset_index(drop=True)
    p50200 = binning[binning["bin_size2"] == 200].reset_index(drop=True)
    labels = [g.replace("peak_", "") for g in p50100["gene"]]
    xp = np.arange(len(p50100))
    axC.bar(xp - w/2, p50100["pearson_r"], width=w, color=BLUE, edgecolor="white", lw=0.4,
            label="bin 50 ↔ 100")
    axC.bar(xp + w/2, p50200["pearson_r"], width=w, color=ORANGE, edgecolor="white", lw=0.4,
            label="bin 50 ↔ 200")
    axC.axhline(p50100["pearson_r"].mean(), color=BLUE, lw=0.9, ls=":", alpha=0.7)
    axC.axhline(p50200["pearson_r"].mean(), color=ORANGE, lw=0.9, ls=":", alpha=0.7)
    axC.set_xticks(xp); axC.set_xticklabels(labels, fontsize=7, rotation=35, ha="right")
    axC.set_xlabel("PAS peak")
    axC.set_ylabel("Cross-bin Pearson r (per peak)")
    axC.set_ylim(0, 1.02)
    axC.set_title("Binning sensitivity: cross-bin consistency (8 high-coverage peaks)", loc="left", fontsize=9.5)
    axC.legend(loc="lower right", fontsize=8)
    # annotation: means.
    axC.text(0.012, 0.04,
             f"mean r:  50↔100 = {p50100['pearson_r'].mean():.3f}   "
             f"50↔200 = {p50200['pearson_r'].mean():.3f}",
             transform=axC.transAxes, fontsize=7.5, color="#333",
             bbox=dict(facecolor="white", edgecolor="#bbb", boxstyle="round,pad=0.3"))
    # Moran's I decay (right inset).
    inset = axC.inset_axes([0.72, 0.22, 0.26, 0.55])
    mi_means = [binning["morans_i_bin50"].mean(),
                binning["morans_i_bin100"].mean(),
                binning["morans_i_bin200"].mean()]
    inset.plot([50, 100, 200], mi_means, "-o", color=GREEN, ms=7, lw=2.0)
    inset.set_xlabel("bin size", fontsize=7.5)
    inset.set_ylabel("mean Moran's I", fontsize=7.5)
    inset.set_title("Spatial signal vs bin size", fontsize=8)
    inset.tick_params(labelsize=7)
    panel_label(axC, "C", x=-0.04, y=1.05)

    fig.suptitle("Supplementary Figure S12 — Stereo-seq QC + binning sensitivity (GSE263789 mouse brain)",
                 fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.005,
             "Panel C reports the available bin pairs (50↔100, 50↔200); a bin-500 comparison was not run for this pilot. "
             "Moran's I decays with coarser binning as expected.",
             ha="center", fontsize=6.8, style="italic", color="#555")
    save_supp(fig, "supp_fig12_stereo_qc_binning.png")
    print(f"[S12] done (peaks={len(peaks)}, reads={dists.size})", flush=True)


if __name__ == "__main__":
    main()

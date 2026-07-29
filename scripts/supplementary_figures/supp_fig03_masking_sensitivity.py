#!/usr/bin/env python
"""Supplementary Figure S3: Masking-level sensitivity.

Re-runs the imputation benchmark at 10/20/30/50% masking (not just 20%) on
GSE183456, for the three Python-native methods (spaGAPA-GP, mean, spatial-KNN).
Grouped bar chart: RMSE (median per-gene, robust) by method × masking level.

stAPAminer and spvAPA (R-based) were not re-run at non-20% levels in this
supplementary sweep; the 20% values for all 5 methods are already reported in
the main benchmark and in Fig 2. The three Python methods are shown here to
demonstrate the masking-level trend is stable.

Data: pipeline_output/supplementary_figures/_cache/s3_masking_sweep.csv
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, METHOD_COLORS, METHOD_LABELS, save_supp

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
CACHE = os.path.join(ROOT, "pipeline_output/supplementary_figures/_cache", "s3_masking_sweep.csv")
setup_rc()


def main():
    df = pd.read_csv(CACHE)
    print(f"[S3] cache: {len(df)} rows", flush=True)

    methods = ["spaGAPA-GP", "mean", "spatial-KNN"]
    fracs = sorted(df["mask_fraction"].unique())

    fig = plt.figure(figsize=(12, 5.5))
    gs = fig.add_gridspec(1, 2, wspace=0.25, width_ratios=[1.3, 1.0])

    # ---------- Panel A: median per-gene RMSE by method × masking level ----------
    axA = fig.add_subplot(gs[0, 0])
    x = np.arange(len(fracs))
    w = 0.26
    for i, mth in enumerate(methods):
        ys = [df[(df["mask_fraction"] == f) & (df["method"] == mth)]["rmse_median_gene"].values[0] for f in fracs]
        axA.bar(x + (i - 1) * w, ys, width=w, color=METHOD_COLORS[mth], edgecolor="white", lw=0.5,
                label=METHOD_LABELS[mth])
        for xi, yi in zip(x + (i - 1) * w, ys):
            axA.annotate(f"{yi:.3f}", xy=(xi, yi), xytext=(0, 2), textcoords="offset points",
                         ha="center", fontsize=6.2, rotation=90)
    axA.set_xticks(x); axA.set_xticklabels([f"{int(f*100)}%" for f in fracs], fontsize=8)
    axA.set_xlabel("Masking level (fraction of observed entries held out)")
    axA.set_ylabel("Median per-gene held-out RMSE")
    axA.set_title("RMSE by method × masking level", loc="left", fontsize=9.5)
    axA.legend(loc="upper left", fontsize=8)
    panel_label(axA, "A", x=-0.06, y=1.05)

    # ---------- Panel B: RMSE growth rate (slope vs masking level) ----------
    axB = fig.add_subplot(gs[0, 1])
    slopes = []
    for mth in methods:
        ys = [df[(df["mask_fraction"] == f) & (df["method"] == mth)]["rmse_median_gene"].values[0] for f in fracs]
        # linear slope through the points.
        sl = np.polyfit(fracs, ys, 1)[0]
        slopes.append(sl)
    xb = np.arange(len(methods))
    bars = axB.bar(xb, slopes, color=[METHOD_COLORS[m] for m in methods], edgecolor="white", lw=0.5, width=0.55)
    for b, s in zip(bars, slopes):
        axB.annotate(f"{s:.3f}\nRMSE/mask%", xy=(b.get_x()+b.get_width()/2, s),
                     xytext=(0, 4), textcoords="offset points", ha="center", fontsize=7.5, fontweight="bold")
    axB.set_xticks(xb); axB.set_xticklabels([METHOD_LABELS[m] for m in methods], fontsize=8)
    axB.set_ylabel("RMSE growth rate  (Δ RMSE per +10% masking)")
    axB.set_title("Methods degrade similarly with masking", loc="left", fontsize=9.5)
    panel_label(axB, "B", x=-0.10, y=1.05)

    fig.suptitle("Supplementary Figure S3 — Masking-level sensitivity (GSE183456)",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.text(0.5, -0.04,
             "n = 120 genes (obs-fraction 0.10–0.50). Median per-gene RMSE shown (robust to sparse-gene outliers). "
             "stAPAminer/spvAPA (R-based) are reported at 20% in Fig 2; only the Python methods are swept here.",
             ha="center", fontsize=6.6, style="italic", color="#555")
    save_supp(fig, "supp_fig03_masking_sensitivity.png")
    print("[S3] done", flush=True)


if __name__ == "__main__":
    main()

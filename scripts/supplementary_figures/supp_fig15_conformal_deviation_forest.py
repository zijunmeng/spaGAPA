#!/usr/bin/env python
"""Supplementary Figure S15: Per-sample conformal coverage deviation forest.

Promoted from the former main Figure 3 Panel D so the main figure (now 5 panels,
~178 mm) stays at print width. Shows, per sample, the empirical − nominal
coverage deviation at all three nominal levels (80/90/95%) as a forest plot,
ordered by 90% deviation. This is the per-sample complement to the pooled
calibration curve (main Fig 3C): it makes explicit that no individual sample
breaches the binomial-CI envelope.

Data:
  pipeline_output/conformal_validation/all_samples_coverage.csv  (11 samples)
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binom
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import (setup_rc, panel_label, BLUE, ORANGE, GREEN, GREY, BLACK,
                    save_supp)

# PAGE_WIDTH_IN is defined in main_figures/_style.py (the shared BIB page
# geometry constant, ~7 in); mirror it here without modifying _style.py.
PAGE_WIDTH_IN = 7.0

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    cov = pd.read_csv(os.path.join(ROOT, "pipeline_output/conformal_validation/all_samples_coverage.csv"))
    cov_sorted = cov.sort_values("coverage_90").reset_index(drop=True)
    n = len(cov_sorted)
    print(f"[S15] {n} samples", flush=True)

    # Print-width layout: PAGE_WIDTH_IN canvas; caption wrapped so the
    # bbox_inches="tight" PDF stays inside the ~178 mm print width.
    fig, ax = plt.subplots(figsize=(PAGE_WIDTH_IN, 4.4))
    fig.subplots_adjust(left=0.20, right=0.96, top=0.90, bottom=0.18)

    ys = np.arange(n)
    dev80 = cov_sorted["coverage_80"].values - 0.80
    dev90 = cov_sorted["coverage_90"].values - 0.90
    dev95 = cov_sorted["coverage_95"].values - 0.95

    # faint ±0.5% guide band
    ax.axvspan(-0.005, 0.005, color="#EEE", lw=0)
    ax.vlines([-0.005, 0.005], -0.6, n - 0.4, color="#DDD", lw=0.7, ls=":")

    # binomial-CI envelope (per-sample, from its own n_test) at 90%, as ±around 0
    hi = binom.ppf(0.975, cov_sorted["n_test"].values, 0.90) / cov_sorted["n_test"].values - 0.90
    lo = binom.ppf(0.025, cov_sorted["n_test"].values, 0.90) / cov_sorted["n_test"].values - 0.90
    ax.fill_betweenx(ys - 0.18, lo, hi, color=ORANGE, alpha=0.10, zorder=0)

    ax.hlines(ys, -0.012, 0.012, color="#EEE", lw=0.8)
    ax.errorbar(dev90, ys,        fmt="o", color=ORANGE, ms=7, label="90%", zorder=5)
    ax.errorbar(dev80, ys - 0.22, fmt="s", color=BLUE,   ms=6, label="80%", zorder=4)
    ax.errorbar(dev95, ys + 0.22, fmt="^", color=GREEN,  ms=6, label="95%", zorder=4)

    ax.axvline(0, color=BLACK, lw=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{d}\n{s.split('gsm')[1]}"
                        for d, s in zip(cov_sorted["dataset"], cov_sorted["sample"])],
                       fontsize=7)
    ax.set_xlabel("Empirical − nominal coverage")
    ax.set_xlim(-0.012, 0.012)
    ax.set_ylim(-0.7, n - 0.3)
    # format x as signed percentage points
    xt = np.array([-0.010, -0.005, 0, 0.005, 0.010])
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{v*100:+.1f}%" for v in xt], fontsize=8)
    ax.set_title("Coverage deviation from nominal (per sample)", loc="left")
    ax.legend(loc="lower right", fontsize=8, ncol=3)
    panel_label(ax, "A", x=-0.04, y=1.06)

    fig.suptitle("Supplementary Figure S15 — Per-sample coverage deviation forest",
                 fontsize=11, fontweight="bold", y=0.985)
    fig.text(0.5, 0.025,
             "Each row = one sample (n=11). Orange band = 95% binomial CI for the 90% level (per-sample n_test).\n"
             "Max |dev|: 0.49% (80%), 0.47% (90%), 0.24% (95%).",
             ha="center", va="bottom", fontsize=7, style="italic", color="#555")

    save_supp(fig, "supp_fig15_conformal_deviation_forest.png")
    print("[S15] done", flush=True)


if __name__ == "__main__":
    main()

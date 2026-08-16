#!/usr/bin/env python
"""Supplementary Figure S8: Per-gene uncertainty–error correlation.

Faceted histogram of per-gene Pearson r (uncertainty ↔ absolute error) for each
of the 5 datasets, with the pooled (across-gene) r marked as a vertical line.
The per-gene r is weak on average (consistent with the marginal r ~0.07), but
the distribution is centred above zero — uncertainty is informative within
genes, not just globally.

Data:
  pipeline_output/uncertainty_within_gene_audit/per_gene_corr_distribution.csv
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, save_supp

# PAGE_WIDTH_IN is defined in main_figures/_style.py (the shared BIB page
# geometry constant, ~7 in); mirror it here without modifying _style.py.
PAGE_WIDTH_IN = 7.0

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    df = pd.read_csv(os.path.join(ROOT, "pipeline_output/uncertainty_within_gene_audit/per_gene_corr_distribution.csv"))
    datasets = list(df.groupby(["dataset", "label"]).groups.keys())
    print(f"[S8] {len(df)} genes across {len(datasets)} datasets", flush=True)

    DS_COLOR = {d: c for d, c in zip([k[0] for k in datasets], [BLUE, ORANGE, GREEN, SKYBLU, RED][: len(datasets)])}

    # Pooled r across all genes (per dataset) — weighted by n_test.
    pooled = (df.groupby("dataset")
              .apply(lambda g: np.average(g["within_gene_r"], weights=g["n_test"]))
              .to_dict())

    # Print-width layout: 2 columns (not 3) so each panel keeps room for the
    # 8-pt x-axis label at PAGE_WIDTH_IN.
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 8.2))
    n = len(datasets)
    ncol = 2 if n >= 2 else n
    nrow = int(np.ceil(n / ncol))
    gs = fig.add_gridspec(nrow, ncol, hspace=0.45, wspace=0.28)

    for i, (ds, lab) in enumerate(datasets):
        ax = fig.add_subplot(gs[i // ncol, i % ncol])
        sub = df[df["dataset"] == ds]
        rvals = sub["within_gene_r"].values
        # Clip extreme tails for visual clarity, keep all data in stats.
        clip = rvals[np.isfinite(rvals)]
        bins = np.linspace(-0.5, 1.0, 40)
        ax.hist(clip, bins=bins, color=DS_COLOR[ds], alpha=0.85, edgecolor="white", linewidth=0.3)
        med = float(np.median(clip))
        p = pooled[ds]
        ax.axvline(med, color="black", lw=1.2, ls="--", label=f"median {med:+.2f}")
        ax.axvline(p, color=RED, lw=1.4, ls="-", label=f"pooled r {p:+.2f}")
        ax.axvline(0, color="#888", lw=0.6)
        ax.set_title(f"{ds}  ({lab})", fontsize=9)
        ax.set_xlabel("Per-gene Pearson r (uncertainty ↔ |error|)", fontsize=8)
        ax.set_ylabel("# genes", fontsize=8)
        ax.legend(loc="upper left", fontsize=7)
        ax.set_xlim(-0.5, 1.0)
        if i == 0:
            panel_label(ax, "A", x=-0.14, y=1.08)
    # Hide unused subplots.
    for j in range(n, nrow * ncol):
        fig.add_subplot(gs[j // ncol, j % ncol]).axis("off")

    fig.suptitle("Supplementary Figure S8 — Per-gene uncertainty–error\ncorrelation distribution",
                 fontsize=11, fontweight="bold", y=0.995, va="top")
    fig.text(0.5, 0.005,
             "Each panel = one dataset; genes weighted by #test points. Pooled r (red) is the gene-aggregated\n"
             "correlation; median (black dashed) is the per-gene typical value.",
             ha="center", va="bottom", fontsize=7, style="italic", color="#555")
    save_supp(fig, "supp_fig08_per_gene_uncertainty.png")
    print("[S8] done", flush=True)


if __name__ == "__main__":
    main()

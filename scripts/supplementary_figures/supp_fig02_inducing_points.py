#!/usr/bin/env python
"""Supplementary Figure S2: Sparse GP inducing-point sensitivity.

Two panels:
  A: Held-out RMSE vs number of inducing points m (50/100/200/500). Shows
     overall and median-per-gene RMSE — both are flat across m, i.e. the GP
     posterior mean is robust to the inducing-point count.
  B: Wall time vs m — scales near-linearly (the only practical cost of larger m).

Data: pipeline_output/supplementary_figures/_cache/s2_inducing_sweep.csv
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

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
CACHE = os.path.join(ROOT, "pipeline_output/supplementary_figures/_cache", "s2_inducing_sweep.csv")
setup_rc()


def main():
    df = pd.read_csv(CACHE)
    print(f"[S2] cache: {len(df)} rows\n{df}", flush=True)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 3.4))

    m = df["n_inducing"].values

    # ---------- Panel A: RMSE vs m ----------
    axA.plot(m, df["rmse_overall"], "-o", color=BLUE, ms=7, lw=2.0, label="Overall RMSE")
    axA.plot(m, df["rmse_median_gene"], "-s", color=ORANGE, ms=7, lw=2.0, label="Median per-gene RMSE")
    axA.fill_between(m, df["rmse_median_gene"], df["rmse_overall"], color=BLUE, alpha=0.08)
    # Head-room so the top-right annotation box clears the flat RMSE line.
    axA.set_ylim(0.0, 1.3 * df["rmse_overall"].max())
    axA.set_xscale("log")
    axA.set_xticks(m); axA.set_xticklabels([str(v) for v in m])
    axA.set_xlabel("# inducing points  m  (log scale)")
    axA.set_ylabel("Held-out RMSE")
    axA.set_title("RMSE is flat across m (robust)", loc="left", fontsize=9.5)
    # Annotate the spread.
    span = df["rmse_overall"].max() - df["rmse_overall"].min()
    axA.annotate(f"RMSE span across m = {span:.2e}\n(essentially unchanged)",
                 xy=(0.98, 0.97), xycoords="axes fraction", ha="right", va="top",
                 fontsize=7.5, color=GREEN, fontweight="bold",
                 bbox=dict(facecolor="white", edgecolor="#bbb", boxstyle="round,pad=0.3"))
    axA.legend(loc="center right", fontsize=8)
    panel_label(axA, "A", x=-0.07, y=1.05)

    # ---------- Panel B: wall time vs m ----------
    axB.plot(m, df["wall_s"], "-o", color=GREEN, ms=7, lw=2.2)
    for mv, tv in zip(m, df["wall_s"]):
        axB.annotate(f"{tv:.1f}s", xy=(mv, tv), xytext=(0, 8), textcoords="offset points",
                     ha="center", fontsize=8, fontweight="bold")
    # linear fit reference.
    coef = np.polyfit(m, df["wall_s"], 1)
    xs = np.array([m.min(), m.max()])
    axB.plot(xs, np.polyval(coef, xs), "--", color=GREY, lw=1.0, label=f"linear fit ({coef[0]*1000:.1f} ms/m)")
    axB.set_xscale("log")
    axB.set_xticks(m); axB.set_xticklabels([str(v) for v in m])
    axB.set_xlabel("# inducing points  m  (log scale)")
    axB.set_ylabel("Wall time (s)")
    axB.set_title("Wall time scales ~linearly with m", loc="left", fontsize=9.5)
    axB.legend(loc="upper left", fontsize=8)
    panel_label(axB, "B", x=-0.07, y=1.05)

    fig.suptitle("Supplementary Figure S2 — Sparse GP inducing-point sensitivity\n(GSE183456)",
                 fontsize=11, fontweight="bold", y=1.06)
    fig.text(0.5, -0.04,
             f"n = {int(df['n_genes'].iloc[0])} genes, 20% held-out, mask seed 42. "
             "Larger m buys no accuracy but costs time —\nthe default m=100 is the efficient choice.",
             ha="center", fontsize=7, style="italic", color="#555")
    save_supp(fig, "supp_fig02_inducing_points.png")
    print("[S2] done", flush=True)


if __name__ == "__main__":
    main()

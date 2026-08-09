#!/usr/bin/env python
"""Supplementary Figure S13: Batch correction — QN vs Harmony (APA matrix).

Panel A: mean pairwise PCC (biological signal recovery) vs batch signal
         (residual batch artefact) for each correction method. On the APA matrix
         spaGAPA-QN recovers most of the PCC gain of Harmony while keeping batch
         signal near zero; Harmony collapses per-gene variance.
Panel B: variance preservation — per-gene variance ratio (corrected / before)
         distribution. QN is tightly centred at 1 (variance preserved);
         Harmony drives the ratio toward 0.

CAVEAT: Harmony was designed for expression counts, not APA matrices, so the
comparison here is APA-specific. Batch correction is Pillar 2 and is still under
evaluation; the Harmony "over-correction" reading should be treated as a
preliminary APA-matrix observation, not a general conclusion about Harmony.

Data:
  pipeline_output/bias_correction_v2/harmony_comparison.csv
  pipeline_output/bias_correction_v2/multi_scenario_results.csv
  pipeline_output/bias_correction_v2/biology_preservation.csv
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

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    h = pd.read_csv(os.path.join(ROOT, "pipeline_output/bias_correction_v2/harmony_comparison.csv"))
    bio = pd.read_csv(os.path.join(ROOT, "pipeline_output/bias_correction_v2/biology_preservation.csv"))

    METHOD_COLOR = {
        "before":             GREY,
        "spaGAPA_QN":         BLUE,
        "spaGAPA_linear":     SKYBLU,
        "spaGAPA_linear_preserve": GREEN,
        "Harmony":            ORANGE,
    }
    METHOD_LABEL = {
        "before":             "Before",
        "spaGAPA_QN":         "spaGAPA-QN",
        "spaGAPA_linear":     "spaGAPA-linear",
        "spaGAPA_linear_preserve": "spaGAPA-linear (preserve)",
        "Harmony":            "Harmony",
    }
    METHOD_ORDER = ["before", "spaGAPA_QN", "spaGAPA_linear", "Harmony"]

    h = h.set_index("method").loc[METHOD_ORDER].reset_index()

    fig = plt.figure(figsize=(12, 5.5))
    gs = fig.add_gridspec(1, 2, wspace=0.28, width_ratios=[1.2, 1.0])

    # ---------- Panel A: PCC vs batch signal scatter ----------
    axA = fig.add_subplot(gs[0, 0])
    for _, r in h.iterrows():
        c = METHOD_COLOR[r["method"]]
        axA.scatter(r["batch_signal"], r["mean_pairwise_pcc"], s=180, color=c,
                    edgecolor="black", lw=0.9, zorder=3, label=METHOD_LABEL[r["method"]])
        axA.annotate(METHOD_LABEL[r["method"]], xy=(r["batch_signal"], r["mean_pairwise_pcc"]),
                     xytext=(8, 6), textcoords="offset points", fontsize=7.5, fontweight="bold")
    axA.axvline(0, color="#888", lw=0.7, ls=":")
    axA.set_xlabel("Batch signal  (residual artefact; lower = better)")
    axA.set_ylabel("Mean pairwise PCC  (biological signal; higher = better)")
    axA.set_title("QN recovers biological signal without Harmony's variance collapse (APA matrix)", loc="left", fontsize=9.5)
    # Ideal region shading.
    axA.axvspan(-0.01, 0.01, alpha=0.08, color=GREEN, zorder=0)
    axA.text(0.02, axA.get_ylim()[0] + 0.02, "low batch\nsignal", fontsize=6.5, color=GREEN, style="italic")
    panel_label(axA, "A", x=-0.07, y=1.05)

    # ---------- Panel B: variance preservation + condition-prediction accuracy ----------
    axB = fig.add_subplot(gs[0, 1])
    x = np.arange(len(METHOD_ORDER))
    w = 0.38
    var_med = h.set_index("method").loc[METHOD_ORDER]["median_var_ratio"].values
    bio_map = bio.set_index("method")
    bio_acc = bio_map.loc[METHOD_ORDER]["condition_prediction_accuracy"].values
    axB.bar(x - w/2, var_med, width=w, color=[METHOD_COLOR[m] for m in METHOD_ORDER],
            edgecolor="white", lw=0.5, label="Median var ratio")
    axB.axhline(1.0, color="#555", lw=0.8, ls="--", label="var ratio = 1 (preserved)")
    axB2 = axB.twinx()
    axB2.bar(x + w/2, bio_acc, width=w, color=[METHOD_COLOR[m] for m in METHOD_ORDER], alpha=0.45,
             edgecolor="white", lw=0.5, label="Condition-pred. accuracy")
    axB2.set_ylim(0, 0.65)
    axB2.set_ylabel("Condition-prediction accuracy", fontsize=8.5, color="#555")
    axB2.tick_params(axis="y", labelsize=8)
    axB2.spines["right"].set_visible(True)
    axB.set_xticks(x); axB.set_xticklabels([METHOD_LABEL[m] for m in METHOD_ORDER], fontsize=7.5, rotation=15, ha="right")
    axB.set_ylabel("Median variance ratio (corrected / before)", fontsize=8.5)
    axB.set_title("Variance preservation + downstream-biology recovery", loc="left", fontsize=9.5)
    # custom legend.
    from matplotlib.patches import Patch
    axB.legend(handles=[Patch(facecolor="#888", label="Median var ratio"),
                        Patch(facecolor="#888", alpha=0.45, label="Condition-pred. accuracy"),
                        plt.Line2D([0], [0], color="#555", ls="--", label="var ratio = 1")],
               loc="upper left", fontsize=6.8)
    panel_label(axB, "B", x=-0.10, y=1.05)

    fig.suptitle("Supplementary Figure S13 — Batch correction on the APA matrix: spaGAPA-QN vs Harmony (GSE237183)",
                 fontsize=10.5, fontweight="bold", y=1.01)
    n_genes_txt = f"n = {h.set_index('method').loc['before','n_genes']}" if 'n_genes' in h.columns else 'n = 1,164 genes'
    fig.text(0.5, -0.06,
             f"{n_genes_txt}. On the APA matrix Harmony maximises pairwise PCC but collapses per-gene variance, "
             "while spaGAPA-QN keeps the variance ratio near 1.0 and recovers most of the PCC gain. "
             "Caveat: Harmony was designed for expression counts, not APA matrices, so this is an APA-specific "
             "comparison; batch correction (Pillar 2) is still under evaluation and the variance-collapse reading "
             "is preliminary, not a general statement about Harmony.",
             ha="center", fontsize=6.5, style="italic", color="#555")
    save_supp(fig, "supp_fig13_batch_correction.png")
    print("[S13] done", flush=True)


if __name__ == "__main__":
    main()

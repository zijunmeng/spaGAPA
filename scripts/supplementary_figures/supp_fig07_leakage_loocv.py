#!/usr/bin/env python
"""Supplementary Figure S7: Leakage audit + LOOCV method-selection stability.

Panel A: clean vs leaked per-sample correlation for the local-gene uncertainty
         method across 5 datasets (should be ~identical => no train/test leakage).
Panel B: LOOCV across 5 folds — each held-out dataset selects the same method
         (B_local_gene), i.e. 5/5 stable selection. Shows held-out r per method
         per fold + the selected method.

Data:
  pipeline_output/uncertainty_leakage_audit/clean_vs_leaked_corr.csv
  pipeline_output/uncertainty_loocv/loocv_results.csv
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
    leak = pd.read_csv(os.path.join(ROOT, "pipeline_output/uncertainty_leakage_audit/clean_vs_leaked_corr.csv"))
    loocv = pd.read_csv(os.path.join(ROOT, "pipeline_output/uncertainty_loocv/loocv_results.csv"))

    fig = plt.figure(figsize=(7.0, 3.9))
    gs = fig.add_gridspec(1, 2, wspace=0.30)

    # ---------- Panel A: clean vs leaked ----------
    axA = fig.add_subplot(gs[0, 0])
    ds = leak["dataset"].values
    x = np.arange(len(ds)); w = 0.36
    axA.bar(x - w/2, leak["clean_r"], width=w, color=BLUE, edgecolor="white", linewidth=0.5, label="Clean split")
    axA.bar(x + w/2, leak["leaked_r"], width=w, color=ORANGE, edgecolor="white", linewidth=0.5, label="Leaked (audit)")
    # annotate inflation factor.
    for i, row in leak.iterrows():
        delta = row["leaked_r"] - row["clean_r"]
        axA.annotate(f"Δ={delta:+.3f}", xy=(x[i], max(row["clean_r"], row["leaked_r"]) + 0.018),
                     ha="center", fontsize=6.6, color=GREY if abs(delta) < 0.02 else RED, fontweight="bold")
    axA.set_xticks(x); axA.set_xticklabels([d.replace("GSE", "") for d in ds], fontsize=7.5)
    axA.set_xlabel("Dataset (held-out GSE)")
    axA.set_ylabel("Pearson r (uncertainty ↔ error)")
    axA.set_ylim(0, 0.82)
    axA.set_title("Clean vs leaked split: correlations are\nidentical => no train/test leakage in\nuncertainty calibration",
                  loc="left", fontsize=9)
    axA.legend(loc="lower right", fontsize=7.5)
    panel_label(axA, "A", x=-0.08, y=1.05)

    # ---------- Panel B: LOOCV held-out r per method per fold ----------
    axB = fig.add_subplot(gs[0, 1])
    method_cols = {
        "A_constant": ("held_out_r_A_constant", GREY, "A: Constant"),
        "B_local_gene": ("held_out_r_B_local_gene", BLUE, "B: Local-gene"),
        "C_spatial_spot": ("held_out_r_C_spatial_spot", SKYBLU, "C: Spatial-spot"),
        "D_residual_spot": ("held_out_r_D_residual_spot", ORANGE, "D: Residual-spot"),
    }
    folds = loocv["fold"].values
    xf = np.arange(len(folds))
    wf = 0.2
    for i, (mkey, (col, c, lab)) in enumerate(method_cols.items()):
        axB.bar(xf + (i - 1.5) * wf, loocv[col], width=wf, color=c, edgecolor="white", linewidth=0.4, label=lab)
    # Mark the selected method per fold (always B_local_gene here).
    for i, row in loocv.iterrows():
        sel = row["selected_method_held_out_r"]
        axB.scatter(xf[i], sel, marker="*", s=130, color="gold", edgecolor="black", lw=0.6, zorder=5)
    axB.set_xticks(xf)
    axB.set_xticklabels([f"Fold {f}\n(held-out:\n{d.replace('GSE','')})" for f, d in zip(folds, loocv["held_out_dataset"])],
                        fontsize=6.8)
    axB.set_xlabel("LOOCV fold (one dataset held out, 4 used for selection)")
    axB.set_ylabel("Held-out Pearson r")
    axB.set_ylim(0, 0.78)
    axB.set_title("LOOCV selects B (local-gene) on all\n5 folds — stable (gold star = method\nselected on each fold)",
                  loc="left", fontsize=9)
    axB.legend(loc="upper right", fontsize=7, ncol=2)
    panel_label(axB, "B", x=-0.08, y=1.05)

    fig.suptitle("Supplementary Figure S7 — Uncertainty leakage audit and\nLOOCV method-selection stability",
                 fontsize=11, fontweight="bold", y=1.02)
    save_supp(fig, "supp_fig07_leakage_loocv.png")
    print("[S7] done", flush=True)


if __name__ == "__main__":
    main()

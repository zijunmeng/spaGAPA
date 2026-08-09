#!/usr/bin/env python
"""Supplementary Figure S16: MAE-on-retained-set risk-coverage curve.

Moved out of main Figure 4 per reviewer request (7 -> 6 panels): the main
figure's risk-coverage panel keeps RMSE; this supplementary panel is the
duplicate triage view expressed as MAE. Both answer "does uncertainty
ranking beat random retention?" and the two metrics agree, so only RMSE
is shown in the main figure to avoid redundancy.

Design:
  - Left:  MAE vs retention (uncertainty / oracle / random), pooled across
           the 5 datasets, with the same triage annotations as main Fig 4F
           (expressed as MAE here).
  - Right: per-dataset MAE improvement at 80% retention (paired dots +
           mean bar), making the paired statistical unit explicit: each
           dot is one of the 5 independent datasets, paired
           (uncertainty vs random) at 80% retention.

Data:
  pipeline_output/risk_coverage_curve/risk_coverage_data.csv
  pipeline_output/risk_coverage_curve/summary.json
"""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import (setup_rc, panel_label, save_supp,
                    BLUE, ORANGE, GREEN, GREY, RED)

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
PAGE_WIDTH_IN = 7.0   # mirror main_figures/_style.py (do NOT import/edit it)
setup_rc()

DATA = os.path.join(ROOT, "pipeline_output/risk_coverage_curve")
risk = pd.read_csv(os.path.join(DATA, "risk_coverage_data.csv"))
with open(os.path.join(DATA, "summary.json")) as fh:
    risk_summary = json.load(fh)

# dataset short labels (kidney / brain / colon / liver / ...)
LABELS = {
    "GSE183456": "kidney", "GSE220442": "brain", "GSE169749": "colon",
    "GSE338525": "liver", "GSE237183": "brain-2",
}


def main():
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 3.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.32,
                          left=0.08, right=0.985, top=0.92, bottom=0.16)
    axL = fig.add_subplot(gs[0, 0])   # MAE vs retention (pooled)
    axR = fig.add_subplot(gs[0, 1])   # per-dataset improvement at 80%

    # ---- Left: MAE vs retention, pooled across datasets ----
    mae_piv = risk.pivot_table(index="retention_fraction", columns="method",
                               values="mae", aggfunc="mean")
    ret = mae_piv.index.values
    for m, col, lab in [("uncertainty", BLUE, "Uncertainty"),
                        ("oracle", GREEN, "Oracle"),
                        ("random", GREY, "Random")]:
        axL.plot(ret, mae_piv[m].values, "o-", color=col, lw=2.0, ms=6,
                 markeredgecolor="white", label=lab)
    axL.set_xlabel("Retention fraction (spots kept, ranked by sigma)")
    axL.set_ylabel("MAE on retained set")
    axL.set_title("MAE risk-coverage (pooled, 5 datasets)", loc="left", fontsize=9)
    axL.set_xlim(0.15, 1.05)
    axL.set_xticks(ret)
    axL.legend(loc="upper left", fontsize=7)
    # annotate the triage gap at 80% retention (MAE basis)
    imp80_mae = ((mae_piv.loc[0.8, "random"] - mae_piv.loc[0.8, "uncertainty"])
                 / mae_piv.loc[0.8, "random"] * 100)
    axL.annotate(f"-{imp80_mae:.0f}% MAE @80%",
                 xy=(0.8, mae_piv.loc[0.8, "uncertainty"]),
                 xytext=(0.45, mae_piv.loc[0.6, "uncertainty"]),
                 fontsize=7, color=BLUE, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
    panel_label(axL, "A")

    # ---- Right: per-dataset MAE improvement at 80% retention (paired) ----
    # The statistical unit is ONE dataset; this makes the paired t-test unit
    # explicit (n=5 independent datasets), matching main Fig 4F's caption.
    imp80 = {}
    for ds in risk["dataset"].unique():
        sub = risk[(risk.dataset == ds) & (risk.retention_fraction == 0.8)]
        u = sub.loc[sub.method == "uncertainty", "mae"].iloc[0]
        r_ = sub.loc[sub.method == "random", "mae"].iloc[0]
        imp80[ds] = (r_ - u) / r_ * 100
    order = list(imp80.keys())
    vals = [imp80[d] for d in order]
    labs = [LABELS.get(d, d) for d in order]
    x = np.arange(len(order))
    axR.bar(x, vals, color=BLUE, edgecolor="white", width=0.62,
            label="per dataset")
    axR.scatter(x, vals, color="black", s=18, zorder=4)
    mean_v = float(np.mean(vals))
    axR.axhline(mean_v, color=ORANGE, ls="--", lw=1.2,
                label=f"mean {mean_v:.0f}%")
    axR.set_xticks(x)
    axR.set_xticklabels(labs, rotation=30, ha="right", fontsize=7)
    axR.set_ylabel("MAE improvement @80%\n(random - uncertainty) / random")
    axR.set_title("Per-dataset MAE gain @80%", loc="left", fontsize=9)
    axR.legend(loc="upper right", fontsize=7)
    axR.set_ylim(0, max(vals) * 1.18)
    # paired t-test p (same pairing unit as main Fig 4F)
    pval = float(risk_summary["paired_ttest_pvalue_at_80"])
    pval_str = f"p={pval:.3f}" if pval >= 0.001 else f"p={pval:.2e}"
    n_ds = int(risk_summary.get("n_datasets", 5))
    axR.text(0.5, -0.30,
             f"paired t-test {pval_str}, n={n_ds} datasets\n"
             f"(same pairing unit as main Fig 4F; RMSE there, MAE here)",
             transform=axR.transAxes, ha="center", fontsize=6.5,
             style="italic", color="#444")
    panel_label(axR, "B")

    fig.suptitle("Supplementary Figure S16 - MAE risk-coverage (companion to Fig 4F)",
                 fontsize=10, fontweight="bold", y=0.99)
    save_supp(fig, "supp_fig16_mae_risk_coverage.png")
    print("[S16] done", flush=True)


if __name__ == "__main__":
    main()

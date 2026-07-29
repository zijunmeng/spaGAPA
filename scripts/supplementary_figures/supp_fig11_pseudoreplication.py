#!/usr/bin/env python
"""Supplementary Figure S11: Pseudoreplication analysis.

Panel A: spot-level vs donor-level (sample-level) significant gene counts.
         Spot-level (pseudoreplicated) calls 91 genes significant; the correct
         donor-level test (n=3 ctrl vs n=3 AD) calls 0 — classic pseudoreplication
         inflation. Direction-agreement genes = 29.
Panel B: effect-size agreement for the 29 direction-agreement genes: scatter of
         spot-level Δ vs donor-level Δ (the directions agree, magnitudes shrink
         at donor level).
Panel C: uncertainty-guided filter on the spot-level calls — global spot filter
         (robust, 84/91 survive halving) vs gene-specific filter (37% of calls
         dropped as depending on unreliable spots).

Data:
  pipeline_output/gse220442_sample_level_differential/{spot_vs_sample_comparison.csv,summary.json}
  pipeline_output/uncertainty_guided_differential/{results.json,gene_lists_comparison.csv}
"""
import os, sys, json
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
    comp = pd.read_csv(os.path.join(ROOT, "pipeline_output/gse220442_sample_level_differential/spot_vs_sample_comparison.csv"))
    summ = json.load(open(os.path.join(ROOT, "pipeline_output/gse220442_sample_level_differential/summary.json")))
    ug = json.load(open(os.path.join(ROOT, "pipeline_output/uncertainty_guided_differential/results.json")))

    n_spot_sig = int(comp["spot_sig"].sum())
    n_sample_sig = int(comp["sample_sig"].sum())
    n_dir = int(comp["direction_agree"].sum())
    n_testable = summ["n_genes_testable_sample_level"]

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28, height_ratios=[1.0, 1.0])

    # ---------- Panel A: pseudoreplication inflation ----------
    axA = fig.add_subplot(gs[0, 0])
    cats = ["Donor-level\n(correct replicate:\nn=3 vs n=3)",
            "Spot-level\n(pseudoreplication:\n~29k spots)",
            "Direction\nagreement"]
    vals = [n_sample_sig, n_spot_sig, n_dir]
    colors = [GREEN, RED, ORANGE]
    bars = axA.bar(cats, vals, color=colors, edgecolor="white", lw=0.6, width=0.6)
    for b, v in zip(bars, vals):
        axA.annotate(f"{v}", xy=(b.get_x() + b.get_width()/2, v),
                     xytext=(0, 3), textcoords="offset points", ha="center",
                     fontsize=11, fontweight="bold")
    axA.set_ylabel(f"# significant genes  (of {n_testable} testable)")
    axA.set_title("Pseudoreplication inflates the call set", loc="left", fontsize=9.5)
    axA.text(0.5, -0.30,
             "Spot-level t-test treats spots as independent replicates => 91 sig.\n"
             "Donor-level test (correct unit) => 0 sig; 29 spot-calls agree in direction.",
             transform=axA.transAxes, ha="center", fontsize=7, style="italic", color="#444")
    panel_label(axA, "A", x=-0.10, y=1.05)

    # ---------- Panel B: effect-size agreement (29 direction-agree genes) ----------
    axB = fig.add_subplot(gs[0, 1])
    da = comp[comp["direction_agree"] == 1].copy()
    axB.scatter(da["sample_delta"], da["spot_delta"], s=34, color=BLUE,
                edgecolor="white", lw=0.4, alpha=0.85, zorder=3)
    # unity line (donor vs spot effect).
    lim = max(np.abs(da["sample_delta"]).max(), np.abs(da["spot_delta"]).max()) * 1.1
    axB.plot([-lim, lim], [-lim, lim], color=GREY, lw=1.0, ls="--", label="y = x (equal effect)")
    axB.axhline(0, color="#aaa", lw=0.6); axB.axvline(0, color="#aaa", lw=0.6)
    # label the strongest-effect genes.
    da["abs_delta"] = np.abs(da["spot_delta"])
    for _, row in da.nlargest(5, "abs_delta").iterrows():
        axB.annotate(row["gene"], xy=(row["sample_delta"], row["spot_delta"]),
                     xytext=(4, 4), textcoords="offset points", fontsize=6.5, fontweight="bold")
    axB.set_xlim(-lim, lim); axB.set_ylim(-lim, lim)
    axB.set_xlabel("Donor-level Δ (distal usage, AD − ctrl)")
    axB.set_ylabel("Spot-level Δ (AD − ctrl)")
    axB.set_title(f"Effect sizes agree in direction (n={len(da)} genes)", loc="left", fontsize=9.5)
    axB.legend(loc="lower right", fontsize=7.5)
    panel_label(axB, "B", x=-0.10, y=1.05)

    # ---------- Panel C: uncertainty-guided filtering ----------
    axC = fig.add_subplot(gs[1, :])
    # Two sub-stories: (left) global spot filter retention curve; (right) gene-specific filter.
    # Use a 2-axis inset.
    sub = axC.inset_axes([0.07, 0.14, 0.40, 0.74])
    fc = ug["fine_retention_curve"]
    drops = fc["drop_pct"]
    nsig = fc["n_sig"]
    pct = fc["pct_baseline_retained"]
    sub.plot(drops, nsig, "-o", color=BLUE, ms=6, lw=2.0, label="# sig genes")
    sub2 = sub.twinx()
    sub2.plot(drops, pct, "--s", color=ORANGE, ms=5, lw=1.6, label="% baseline retained")
    sub.set_xlabel("Most-uncertain spots dropped (%)", fontsize=8)
    sub.set_ylabel("# significant genes", fontsize=8, color=BLUE)
    sub2.set_ylabel("% baseline retained", fontsize=8, color=ORANGE)
    sub.tick_params(axis="y", labelcolor=BLUE, labelsize=7.5); sub.tick_params(axis="x", labelsize=7.5)
    sub2.tick_params(axis="y", labelcolor=ORANGE, labelsize=7.5)
    sub.set_title("(i) Global spot filter — robust", loc="left", fontsize=8.5, color=BLUE)
    sub2.spines["right"].set_visible(True)
    # annotate endpoints.
    sub.annotate(f"{nsig[-1]} sig\n({pct[-1]:.0f}% kept)", xy=(drops[-1], nsig[-1]),
                 xytext=(-6, 6), textcoords="offset points", fontsize=7, ha="right",
                 color=BLUE, fontweight="bold")

    subR = axC.inset_axes([0.58, 0.14, 0.40, 0.74])
    compD = ug["comparison_D_gene_specific_weighted"]
    cats2 = ["Baseline\n(all spots)", "Survive\ngene-specific\nfilter", "Dropped by\nfilter\n(false-pos.)"]
    n_base = ug["comparison_A"]["n_sig_baseline"]
    n_survive = round(n_base * compD["retention_of_baseline_pct"] / 100.0)
    n_dropped = compD["n_dropped_by_weighting"]
    v2 = [n_base, n_survive, n_dropped]
    colors2 = [GREY, GREEN, RED]
    xpos = np.arange(len(v2))
    bars2 = subR.bar(xpos, v2, color=colors2, edgecolor="white", lw=0.5, width=0.62)
    subR.set_xticks(xpos)
    subR.set_xticklabels(cats2, fontsize=7)
    for b, v in zip(bars2, v2):
        subR.annotate(f"{v}", xy=(b.get_x()+b.get_width()/2, v), xytext=(0, 3),
                      textcoords="offset points", ha="center", fontsize=10, fontweight="bold")
    subR.set_ylabel("# baseline-sig genes", fontsize=8)
    subR.set_title(f"(ii) Gene-specific uncertainty filter — {compD['false_positive_reduction_pct']:.0f}% flagged",
                   loc="left", fontsize=8.5, color=RED)
    subR.tick_params(axis="x", labelsize=7.5); subR.tick_params(axis="y", labelsize=7.5)

    axC.axis("off")
    axC.set_title("Uncertainty-guided filtering audits the spot-level call set",
                  loc="left", fontsize=9.5, pad=6)
    axC.text(0.5, 0.03,
             "(i) Dropping the most-uncertain spots barely changes the call count (robust). "
             "(ii) A gene-specific filter flags 37% of baseline calls as resting on that gene's unreliable spots.",
             ha="center", fontsize=7, style="italic", color="#444", transform=axC.transAxes)
    panel_label(axC, "C", x=0.0, y=0.98)

    fig.suptitle("Supplementary Figure S11 — Pseudoreplication + uncertainty-guided differential APA (GSE220442)",
                 fontsize=11, fontweight="bold", y=0.995)
    save_supp(fig, "supp_fig11_pseudoreplication.png")
    print(f"[S11] done (spot_sig={n_spot_sig}, sample_sig={n_sample_sig}, dir_agree={n_dir})", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Supplementary Figure S6: Full conformal coverage across all 11 samples.

Panels:
  A: Empirical vs nominal coverage at 80/90/95% (grouped bar, per sample).
  B: Interval width = 2*qhat at each level (per sample).
  C: Interval score (Gneiting & Raftery 2007, "Winkler score") per sample,
     computed exactly per observation from the persisted (y_true, lo, hi)
     bounds (calibrate_uncertainty_all_datasets.py now writes
     per_observation_bounds.npz + winkler_{80,90,95} columns). Lower = better.
     If the per-observation bounds are unavailable it falls back to the
     qhat+RMSE approximation and labels the panel accordingly.
  D: Empirical subgroup coverage by uncertainty quintile (deviation from nominal
     90%), summarised across samples — a marginal, not conditional, diagnostic.

Data:
  pipeline_output/conformal_validation/all_samples_coverage.csv  (11 samples, incl. winkler_{80,90,95})
  pipeline_output/conformal_validation/per_observation_bounds.npz (per-obs y/lo/hi)
  pipeline_output/conformal_conditional_coverage/by_uncertainty_quintile.csv
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, SKYBLU, RED, GREY
import _style as _S

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    cov = pd.read_csv(os.path.join(ROOT, "pipeline_output/conformal_validation/all_samples_coverage.csv"))
    # The conformal pipeline reports 11 samples; the CSV has 11 rows (12 incl. header).
    cov = cov.sort_values(["dataset", "sample"]).reset_index(drop=True)
    n = len(cov)
    print(f"[S6] {n} samples", flush=True)
    labels = [f"{d.replace('GSE','')}\n{s.replace('gsm','')}" for d, s in zip(cov["dataset"], cov["sample"])]

    levels = [(0.80, "coverage_80", "qhat_80", BLUE),
              (0.90, "coverage_90", "qhat_90", ORANGE),
              (0.95, "coverage_95", "qhat_95", GREEN)]

    # Panel C uses the exact per-observation Winkler (Gneiting–Raftery) score
    # from the persisted bounds when available; see the Panel C block below.
    fig = plt.figure(figsize=(7.0, 7.5))
    fig.subplots_adjust(left=0.12, right=0.97, top=0.90, bottom=0.09)
    gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.34)

    # ---------- Panel A: coverage ----------
    axA = fig.add_subplot(gs[0, 0])
    x = np.arange(n); w = 0.26
    for i, (nom, ccol, _, col) in enumerate(levels):
        axA.bar(x + (i - 1) * w, cov[ccol], width=w, color=col, edgecolor="white", linewidth=0.4,
                label=f"{int(nom*100)}% target")
    for nom in [0.80, 0.90, 0.95]:
        axA.axhline(nom, color="#555", lw=0.7, ls=":", zorder=0)
    axA.set_xticks(x); axA.set_xticklabels(labels, fontsize=6.3, rotation=90, ha="center")
    axA.set_ylim(0.76, 0.985)
    axA.set_ylabel("Empirical coverage")
    axA.set_xlabel("Sample (GSE / GSM)")
    axA.set_title("Marginal coverage matches nominal level", loc="left")
    axA.legend(title="Target", loc="lower right", fontsize=7)
    panel_label(axA, "A", x=-0.06, y=1.05)

    # ---------- Panel B: interval width (2*qhat) ----------
    axB = fig.add_subplot(gs[0, 1])
    for i, (nom, _, qcol, col) in enumerate(levels):
        axB.bar(x + (i - 1) * w, 2 * cov[qcol], width=w, color=col, edgecolor="white", linewidth=0.4,
                label=f"{int(nom*100)}%")
    axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=6.3, rotation=90, ha="center")
    axB.set_ylabel("Interval width  (2 $\\cdot$ $\\hat{q}$)")
    axB.set_xlabel("Sample")
    axB.set_title("Interval width grows with target level", loc="left")
    axB.legend(title="Target", loc="upper left", fontsize=7)
    panel_label(axB, "B", x=-0.06, y=1.05)

    # ---------- Panel C: interval score (Winkler) ----------
    # Exact Gneiting–Raftery (2007) interval score, averaged per sample over the
    # held-out test set from the persisted per-observation (y, lo, hi) bounds.
    # S = (hi-lo) + (2/alpha)*(lo-y)*1[y<lo] + (2/alpha)*(y-hi)*1[y>hi].
    # Falls back to the qhat+RMSE approximation only if the winkler_* columns
    # are absent (e.g. an older cached CSV) and labels the panel accordingly.
    axC = fig.add_subplot(gs[1, 0])
    has_winkler = all(f"winkler_{lvl}" in cov.columns for _, _, lvl, _ in
                      [(0.80, None, 80, None), (0.90, None, 90, None), (0.95, None, 95, None)])
    if has_winkler:
        for i, (nom, _ccol, qcol, col) in enumerate(levels):
            lvl = int(round(nom * 100))
            axC.bar(x + (i - 1) * w, cov[f"winkler_{lvl}"], width=w, color=col,
                    edgecolor="white", linewidth=0.4, label=f"{int(nom*100)}%")
        score_ylabel = "Interval score\n(Winkler; lower = better)"
        score_title = "Exact Winkler score\nper held-out spot"
        score_formal = True
    else:
        for i, (nom, ccol, qcol, col) in enumerate(levels):
            alpha = 1 - nom
            width = 2 * cov[qcol]
            outside_pen = (2.0 / alpha) * cov["rmse"] * (1.0 - cov[ccol])
            score = width + outside_pen
            axC.bar(x + (i - 1) * w, score, width=w, color=col, edgecolor="white",
                    linewidth=0.4, label=f"{int(nom*100)}%")
        score_ylabel = "Interval score (approx. from $\\hat{q}$ + RMSE; lower = better)"
        score_title = "Interval score (approximation — per-obs bounds unavailable)"
        score_formal = False
    axC.set_xticks(x); axC.set_xticklabels(labels, fontsize=6.3, rotation=90, ha="center")
    axC.set_ylabel(score_ylabel)
    axC.set_xlabel("Sample")
    axC.set_title(score_title, loc="left")
    axC.legend(title="Target", loc="upper left", fontsize=7)
    panel_label(axC, "C", x=-0.06, y=1.05)

    # ---------- Panel D: empirical subgroup coverage by uncertainty quintile ----------
    axD = fig.add_subplot(gs[1, 1])
    cc = pd.read_csv(os.path.join(ROOT, "pipeline_output/conformal_conditional_coverage/by_uncertainty_quintile.csv"))
    # Mean deviation across samples, per quintile, at 90% target.
    summary = cc.groupby("quintile")["deviation"].agg(["mean", "std", "count"]).reset_index()
    qx = np.arange(len(summary))
    bars = axD.bar(qx, summary["mean"], yerr=summary["std"], capsize=3,
                   color=[SKYBLU, GREEN, GREY, ORANGE, RED][: len(summary)],
                   edgecolor="white", linewidth=0.5, error_kw={"lw": 0.8, "ecolor": "#555"})
    axD.axhline(0, color="#333", lw=0.8)
    axD.set_xticks(qx)
    axD.set_xticklabels([f"Q{int(q)}\n({'low unc' if q==1 else 'high unc' if q==len(summary) else ''})".strip("\n")
                         for q in summary["quintile"]], fontsize=7.5)
    axD.set_ylabel("Coverage − 90% target  (deviation)")
    axD.set_xlabel("Uncertainty quintile (per-sample GP posterior std)")
    axD.set_title("Empirical subgroup coverage by uncertainty quintile", loc="left")
    axD.text(0.5, -0.32,
             "Positive deviation = conservative (over-cover); empirical subgroup check,\n"
             "not a formal conditional guarantee. Highest-uncertainty quintile over-covers.",
             transform=axD.transAxes, fontsize=6.8, ha="center", style="italic", color="#444")
    panel_label(axD, "D", x=-0.06, y=1.05)

    fig.suptitle("Supplementary Figure S6 — Full conformal coverage (11 samples × 80/90/95%)",
                 fontsize=10, fontweight="bold", y=0.995)
    # Caption note on interval score (formal per-observation Winkler, or the
    # approximation if the per-obs bounds were unavailable).
    if score_formal:
        score_note = (
            "Panel C: exact Gneiting–Raftery (2007) Winkler score per held-out spot, "
            "S = (hi−lo) + (2/α)·(lo−y)·1[y<lo] + (2/α)·(y−hi)·1[y>hi],\n"
            "averaged per sample over the test set (bounds in per_observation_bounds.npz). "
            "Coverage (A) and width (B) are exact; Panel D is an empirical subgroup check,\n"
            "not a formal conditional-coverage guarantee.")
    else:
        score_note = (
            "Panel C interval score is the Gneiting–Raftery score approximated from per-sample qhat and RMSE "
            "(width + (2/α)·RMSE·(1−coverage)); per-observation bounds were unavailable,\n"
            "so only this aggregate proxy could be computed. "
            "Coverage (A) and width (B) are exact; Panel D is an empirical subgroup check, "
            "not a formal conditional-coverage guarantee.")
    fig.text(0.5, 0.012, score_note,
             ha="center", fontsize=6.8, style="italic", color="#555")
    # Explicit save (no bbox_inches='tight'): rotated tick labels make the tight
    # bbox expand far beyond the canvas (a matplotlib quirk), so save at the
    # exact figsize — all artists are placed within the canvas margins above.
    import matplotlib.pyplot as _plt
    _plt.ioff()
    fig.savefig(os.path.join(_S.OUT, "supp_fig06_conformal_coverage.pdf"), facecolor="white")
    fig.savefig(os.path.join(_S.OUT, "supp_fig06_conformal_coverage.png"), dpi=300, facecolor="white")
    _plt.close(fig)
    print(f"[S6] done -> {_S.OUT}/supp_fig06_conformal_coverage.pdf (+png)", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Supplementary Figure S14: AD vs WT descriptive APA effect sizes (Stereo-seq).

Houses the n=1-vs-n=1 exploratory effect-size analysis that was moved OUT of
main Figure 7 to avoid over-claiming in the main text. All data are REAL
(per-gene sample-level effect sizes from the GSE263789 AD vs WT Stereo-seq
sections); the panel is explicitly framed as exploratory / effect-size based
with NO p-value, NO FDR (biological n = 1 per condition -> no valid
between-individual test exists; see sample_level_note.md).

Panels:
  A. Top-10 genes by |Delta = AD_mean - WT_mean| (ranked horizontal bars,
     coloured by APA direction). A prominent caveat banner flags n=1.
  B. Direction split of the 607 candidate genes (|Delta| > 0.1): distal-up
     vs proximal-up in AD; plus a small annotation of the 1,090-gene universe
     tested.

Data:
  pipeline_output/gse263789_ad_vs_wt_differential/sample_level_effect_sizes.csv
  (gene, ad_mean, wt_mean, delta_AD_minus_WT, abs_delta, ad_n_spots_obs,
   wt_n_spots_obs, direction, candidate_differential)

Run on S91: ~/anaconda3/envs/spagapa/bin/python supp_fig14_stereo_adwt_descriptive.py
"""
import os
import sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import (  # noqa: E402
    setup_rc, panel_label, save_supp,
    BLUE, ORANGE, GREEN, RED, GREY, BLACK,
)

# PAGE_WIDTH_IN is defined in main_figures/_style.py (the shared BIB page
# geometry constant, ~7 in); mirror it here without modifying _style.py.
PAGE_WIDTH_IN = 7.0

setup_rc()

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
AD = f"{ROOT}/gse263789_ad_vs_wt_differential"
eff = pd.read_csv(f"{AD}/sample_level_effect_sizes.csv")

n_tested = len(eff)
n_cand = int(eff["candidate_differential"].sum())
n_distal = int(((eff["candidate_differential"] == 1) &
                (eff["delta_AD_minus_WT"] > 0)).sum())
n_prox = int(((eff["candidate_differential"] == 1) &
              (eff["delta_AD_minus_WT"] < 0)).sum())

fig = plt.figure(figsize=(PAGE_WIDTH_IN, 5.0))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.55, 1.0], wspace=0.32,
              left=0.115, right=0.965, top=0.80, bottom=0.20)

# --------------------------------------------------------------------------- #
# Panel A -- top-10 by |Delta|
# --------------------------------------------------------------------------- #
axA = fig.add_subplot(gs[0, 0])
top = eff.sort_values("abs_delta", ascending=False).head(10).iloc[::-1]
colors_a = [ORANGE if d > 0 else BLUE for d in top["delta_AD_minus_WT"]]
axA.barh(range(len(top)), top["delta_AD_minus_WT"], color=colors_a,
         edgecolor="black", linewidth=0.6)
axA.set_yticks(range(len(top)))
axA.set_yticklabels(top["gene"].values, fontsize=7.5)
axA.axvline(0, color="black", linewidth=0.8)
axA.axvline(0.1, color=GREY, linestyle=":", linewidth=0.8)
axA.axvline(-0.1, color=GREY, linestyle=":", linewidth=0.8)
axA.set_xlabel(r"$\Delta$ (AD$_{mean}$ - WT$_{mean}$) site usage")
axA.set_xlim(-1.20, 1.20)
axA.set_title("Top-10 genes by |effect size|", fontsize=9.5)
for i, d in enumerate(top["delta_AD_minus_WT"]):
    axA.text(d + (0.02 if d >= 0 else -0.02), i, f"{d:+.2f}", va="center",
             ha="left" if d >= 0 else "right", fontsize=6.8, fontweight="bold")
legend_elems = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=ORANGE,
           markersize=8, label="distal-up in AD"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=BLUE,
           markersize=8, label="proximal-up in AD"),
]
axA.legend(handles=legend_elems, loc="lower right", fontsize=7,
           frameon=True, framealpha=0.9)
panel_label(axA, "A", x=-0.18, y=1.06)

# --------------------------------------------------------------------------- #
# Panel B -- candidate direction split
# --------------------------------------------------------------------------- #
axB = fig.add_subplot(gs[0, 1])
cats = ["Distal-up\n(AD)", "Proximal-up\n(AD)"]
vals = [n_distal, n_prox]
bars = axB.bar(cats, vals, color=[ORANGE, BLUE], edgecolor="black", linewidth=0.7,
               width=0.55)
for b, v in zip(bars, vals):
    axB.text(b.get_x() + b.get_width() / 2, v + max(vals) * 0.02, str(v),
             ha="center", va="bottom", fontsize=9, fontweight="bold")
axB.set_ylabel(f"Candidate genes (|$\\Delta$| > 0.1)\nout of {n_tested} tested")
axB.set_ylim(0, max(vals) * 1.30)
axB.set_title(f"{n_cand} candidates of {n_tested} genes\n"
              "(effect size only, no p-value)", fontsize=9.5)
axB.grid(axis="y", linestyle="--", alpha=0.4)
panel_label(axB, "B", x=-0.22, y=1.06)

# --------------------------------------------------------------------------- #
# Prominent n=1 caveat banner across the top of the figure
# --------------------------------------------------------------------------- #
fig.text(0.5, 0.945,
         "Supplementary Figure S14.  AD vs WT Stereo-seq APA effect sizes "
         "(GSE263789)",
         ha="center", va="center", fontsize=11, fontweight="bold")
fig.text(0.5, 0.905,
         r"Biological n = 1 per condition (1 AD vs 1 WT section) "
         r"$\Rightarrow$ NO statistical test / NO FDR. "
         r"Per-gene $\Delta$=AD$_{mean}$-WT$_{mean}$ is reported for "
         r"hypothesis generation only.",
         ha="center", va="center", fontsize=7.6, color=RED, fontweight="bold")
fig.text(0.5, 0.045,
         "Spot-level (pseudoreplicated) calls are invalid here and are NOT "
         "shown; see sample_level_note.md. "
         "Candidates require validation in an independent cohort with "
         "biological replicates.",
         ha="center", va="center", fontsize=6.6, style="italic", color="#666")

out = save_supp(fig, "supp_fig14_stereo_adwt_descriptive.png")
print(f"\n[done] Supplementary Figure S14 written to: {out}")

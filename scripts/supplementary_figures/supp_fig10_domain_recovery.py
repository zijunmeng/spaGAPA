#!/usr/bin/env python
"""Supplementary Figure S10: Domain recovery details on MOB.

Panel A: Leiden resolution sweep (ARI/NMI vs resolution) for every config:
         the 4 GP-imputed APA weight schemes (apa_dominant, balanced,
         spatial_apa, expression_apa) plus the mean-imputed APA baseline
         (mean_apa, same weights as apa_dominant, no GP). The best resolution
         per config is marked.
Panel B: Best-resolution ARI/NMI per config. This is where the headline
         expression_apa / res0.8 result (ARI 0.597, k=5) is reported: adding
         expression (0.5 > apa 0.4) refines the recovery to all 5 anatomical
         layers, and is shown here rather than in the main figure because the
         main figure restricts itself to APA-dominant configs. The mean-impute
         baseline (ARI 0.576, k=4) is also reported here for honesty: per-gene
         mean filling already captures substantial domain structure on MOB,
         and spaGAPA's GP imputation provides a modest improvement at this
         config (apa_dominant 0.550) while expression refinement (0.597) does
         the rest.

Data:
  pipeline_output/mob_domain_recovery/spagapa_metrics.json

Note: The metrics file stores a single run per (configuration, resolution), so
random-seed stability is NOT assessed here; the robustness axis is the set of
weight configurations + the mean-impute baseline. This is stated in the figure
caption.
"""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
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
    with open(os.path.join(ROOT, "pipeline_output/mob_domain_recovery/spagapa_metrics.json")) as fh:
        m = json.load(fh)

    n_true = m["n_mob_layers"]
    runs = m["runs"]
    CONFIG_COLOR = {
        "apa_dominant":   BLUE,
        "balanced":       ORANGE,
        "spatial_apa":    GREEN,
        "expression_apa": RED,
        "mean_apa":       GREY,
    }
    CONFIG_LABEL = {
        "apa_dominant":   "APA-dominant GP (0.2/0.2/0.6)",
        "balanced":       "Balanced GP (0.4/0.4/0.2)",
        "spatial_apa":    "Spatial+APA GP (0.5/0/0.5)",
        "expression_apa": "Expr+APA GP (0.1/0.5/0.4)",
        "mean_apa":       "Mean-imputed APA (0.2/0.2/0.6, no GP)",
    }
    # display order: mean baseline first, then APA-centric GP configs, then the
    # expression-heavy refinement last (so the headline expression_apa result
    # reads as an additive improvement on top of the APA-dominant base).
    CONFIG_ORDER = ["mean_apa", "apa_dominant", "spatial_apa",
                    "balanced", "expression_apa"]

    # Print-width layout: panels stacked vertically (side-by-side left no room
    # for the 3-line per-config annotations in Panel B at PAGE_WIDTH_IN).
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.6))
    gs = fig.add_gridspec(2, 1, hspace=0.55)

    # ---------- Panel A: resolution sweep ----------
    axA = fig.add_subplot(gs[0, 0])  # row 0 (stacked layout)
    for cfg in CONFIG_ORDER:
        run = runs[cfg]
        res = sorted(run["leiden"].keys(), key=lambda s: float(s.split("_")[1]))
        xs = [float(r.split("_")[1]) for r in res]
        aris = [run["leiden"][r]["ari"] for r in res]
        nmis = [run["leiden"][r]["nmi"] for r in res]
        c = CONFIG_COLOR[cfg]
        axA.plot(xs, aris, "-o", color=c, ms=5, lw=1.8, label=f"{cfg} (ARI)")
        # mark best.
        bx = float(run["leiden_best"]["resolution"])
        ba = run["leiden_best"]["ari"]
        axA.scatter([bx], [ba], s=130, marker="*", color=c, edgecolor="black", lw=0.8, zorder=5)
    axA.set_xlabel("Leiden resolution")
    axA.set_ylabel("ARI  (vs MOB anatomical layers)")
    axA.set_title("Leiden resolution sweep (all configs incl. mean baseline)", loc="left", fontsize=9.5)
    axA.legend(loc="lower right", fontsize=6.8, ncol=2)
    axA.text(0.5, -0.22, "Star = best resolution per config. Truth = 5 MOB anatomical layers.",
             transform=axA.transAxes, ha="center", fontsize=6.8, style="italic", color="#444")
    panel_label(axA, "A", x=-0.05, y=1.05)

    # ---------- Panel B: best-of-config ARI/NMI bar ----------
    axB = fig.add_subplot(gs[1, 0])  # row 1 (stacked layout)
    cfgs = CONFIG_ORDER
    aris = [runs[c]["leiden_best"]["ari"] for c in cfgs]
    nmis = [runs[c]["leiden_best"]["nmi"] for c in cfgs]
    ndom = [runs[c]["leiden_best"]["n_domains"] for c in cfgs]
    x = np.arange(len(cfgs)); w = 0.36
    axB.bar(x - w/2, aris, width=w, color=[CONFIG_COLOR[c] for c in cfgs], edgecolor="white", lw=0.5, label="ARI")
    axB.bar(x + w/2, nmis, width=w, color=[CONFIG_COLOR[c] for c in cfgs], alpha=0.5, edgecolor="white", lw=0.5, label="NMI")
    for i, (a, n, d) in enumerate(zip(aris, nmis, ndom)):
        axB.annotate(f"ARI={a:.2f}\nNMI={n:.2f}\n(k={d})", xy=(x[i], max(a, n)),
                     xytext=(0, 4), textcoords="offset points", ha="center", fontsize=5.8)
    axB.set_xticks(x)
    axB.set_xticklabels(["Mean\n(no GP)", "APA-dom\nGP", "Spat+APA\nGP",
                         "Bal\nGP", "Expr+APA\nGP"], fontsize=6.8)
    axB.set_ylim(0, 0.82)
    axB.set_ylabel("Best-resolution score")
    axB.set_xlabel("Weight configuration (spatial / expr / APA)")
    axB.set_title("Best-resolution ARI/NMI per config\n(incl. mean-impute baseline)", loc="left", fontsize=9.5)
    axB.legend(loc="upper right", fontsize=7.5)
    panel_label(axB, "B", x=-0.06, y=1.05)

    fig.suptitle("Supplementary Figure S10 — MOB domain recovery:\nLeiden sweep + config comparison",
                 fontsize=11, fontweight="bold", y=0.995, va="top")
    fig.text(0.5, 0.005,
             f"n = {m['n_spots']} spots, {n_true} true MOB layers. "
             "Bars show 4 GP-imputed APA weight schemes plus the mean-imputed APA baseline\n"
             "(same weights as APA-dominant, no GP). The headline expr+APA config "
             "(ARI 0.597, k=5) refines recovery\n"
             "to all 5 layers; the main figure showcases the APA-dominant config "
             "(expression 0.2 < APA 0.6). Random-seed stability not assessed.",
             ha="center", va="bottom", fontsize=6.5, style="italic", color="#555")
    save_supp(fig, "supp_fig10_domain_recovery.png")
    print("[S10] done", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Supplementary Figure S10: Domain recovery details on MOB.

Panel A: Leiden resolution sweep (ARI/NMI vs resolution) for each of the 4
         weight configurations (apa_dominant, balanced, spatial_apa,
         expression_apa). The chosen "best" resolution per config is marked.
Panel B: Weight-configuration robustness: ARI/NMI under the 4 weight regimes
         (each at its best Leiden resolution). The result is stable across
         weight choices, i.e. the figure is a weight-configuration robustness
         test (NOT a random-seed stability test).

Data:
  pipeline_output/mob_domain_recovery/spagapa_metrics.json

Note: The metrics file stores a single run per (configuration, resolution), so
random-seed stability is NOT assessed here; the robustness axis is the 4 weight
configurations (spatial/expression/APA mixing). This is stated in the figure
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
    }
    CONFIG_LABEL = {
        "apa_dominant":   "APA-dominant (0.2/0.2/0.6)",
        "balanced":       "Balanced (0.4/0.4/0.2)",
        "spatial_apa":    "Spatial+APA (0.4/0/0.6)",
        "expression_apa": "Expr+APA (0/0.4/0.6)",
    }

    fig = plt.figure(figsize=(12, 5.3))
    gs = fig.add_gridspec(1, 2, wspace=0.22, width_ratios=[1.5, 1.0])

    # ---------- Panel A: resolution sweep ----------
    axA = fig.add_subplot(gs[0, 0])
    for cfg, run in runs.items():
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
    axA.set_title("Leiden resolution sweep (4 weight configs)", loc="left", fontsize=9.5)
    axA.legend(loc="lower right", fontsize=6.8, ncol=2)
    axA.text(0.5, -0.22, "Star = best resolution per config. Truth = 5 MOB anatomical layers.",
             transform=axA.transAxes, ha="center", fontsize=6.8, style="italic", color="#444")
    panel_label(axA, "A", x=-0.05, y=1.05)

    # ---------- Panel B: best-of-config ARI/NMI bar ----------
    axB = fig.add_subplot(gs[0, 1])
    cfgs = list(runs.keys())
    aris = [runs[c]["leiden_best"]["ari"] for c in cfgs]
    nmis = [runs[c]["leiden_best"]["nmi"] for c in cfgs]
    ndom = [runs[c]["leiden_best"]["n_domains"] for c in cfgs]
    x = np.arange(len(cfgs)); w = 0.36
    axB.bar(x - w/2, aris, width=w, color=[CONFIG_COLOR[c] for c in cfgs], edgecolor="white", lw=0.5, label="ARI")
    axB.bar(x + w/2, nmis, width=w, color=[CONFIG_COLOR[c] for c in cfgs], alpha=0.5, edgecolor="white", lw=0.5, label="NMI")
    for i, (a, n, d) in enumerate(zip(aris, nmis, ndom)):
        axB.annotate(f"ARI={a:.2f}\nNMI={n:.2f}\n(k={d})", xy=(x[i], max(a, n)),
                     xytext=(0, 4), textcoords="offset points", ha="center", fontsize=6.3)
    axB.set_xticks(x)
    axB.set_xticklabels(["APA-\ndom", "Bal", "Spat+\nAPA", "Expr+\nAPA"], fontsize=7.5)
    axB.set_ylim(0, 0.82)
    axB.set_ylabel("Best-resolution score")
    axB.set_xlabel("Weight configuration")
    axB.set_title("Weight-configuration robustness\n(each config at its best resolution)", loc="left", fontsize=9.5)
    axB.legend(loc="upper right", fontsize=7.5)
    panel_label(axB, "B", x=-0.10, y=1.05)

    fig.suptitle("Supplementary Figure S10 — MOB domain recovery: Leiden sweep + weight-config robustness",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.text(0.5, -0.05,
             f"n = {m['n_spots']} spots, {n_true} true MOB layers. "
             "This is a weight-configuration robustness test (4 spatial/expression/APA weightings\n"
             "× Leiden resolution); random-seed stability across runs was not assessed here.",
             ha="center", fontsize=6.8, style="italic", color="#555")
    save_supp(fig, "supp_fig10_domain_recovery.png")
    print("[S10] done", flush=True)


if __name__ == "__main__":
    main()

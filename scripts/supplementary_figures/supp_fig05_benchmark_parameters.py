#!/usr/bin/env python
"""Supplementary Figure S5: Benchmark parameter table.

Renders benchmark_fairness/parameter_table.csv as a publication-style figure
table (5 methods: language, neighbour basis, key params, # free params, seed,
single-run wall time). Wall time is taken from the cached transparent head-to-
head benchmark (GSE183456); note is added that values are single-run on the
same hardware.

Data:
  pipeline_output/benchmark_fairness/parameter_table.csv
  pipeline_output/benchmark_mean_transparent/transparent_comparison.csv (wall time)
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, METHOD_ORDER, save_supp

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def main():
    params = pd.read_csv(os.path.join(ROOT, "pipeline_output/benchmark_fairness/parameter_table.csv"))
    wall = pd.read_csv(os.path.join(ROOT, "pipeline_output/benchmark_mean_transparent/transparent_comparison.csv"))
    # mean wall time across the two transparent datasets, per method.
    wall_mean = wall.groupby("method")["wall_time_s"].mean().round(1).to_dict()

    SEED_NOTE = "42 (single run; deterministic for mean/spatial-KNN; stAPAminer/spvAPA use their own RNG init)"

    # Order rows to match the canonical method order.
    order = [m for m in METHOD_ORDER if m in set(params["method"])]
    params = params.set_index("method").loc[order].reset_index()

    fig, ax = plt.subplots(figsize=(13.5, 5.2))
    ax.axis("off")

    # Build a compact, human-readable cell text from the raw CSV fields.
    def compact(s, n):
        s = str(s)
        return s if len(s) <= n else s[: n - 1] + "..."

    cols = ["Method", "Language", "Neighbour basis", "Key parameters", "# Free params",
            "Wall time (s)\n(GSE183456, mean of 2)", "Seed / repeats"]
    cell_text = []
    for _, r in params.iterrows():
        m = r["method"]
        # Combine k / n_inducing / length-scale / noise into a short key-param string.
        kp_parts = []
        if str(r.get("k", "")) not in ("", "nan"):
            kp_parts.append(f"k={r['k']}")
        if str(r.get("n_inducing", "")) not in ("", "nan") and "n/a" not in str(r.get("n_inducing", "")).lower():
            kp_parts.append(f"m={r['n_inducing']}")
        if str(r.get("length_scale", "")) not in ("", "nan") and "n/a" not in str(r.get("length_scale", "")).lower():
            kp_parts.append(f"ls={compact(r['length_scale'], 18)}")
        if str(r.get("noise", "")) not in ("", "nan") and "n/a" not in str(r.get("noise", "")).lower():
            kp_parts.append(f"noise={compact(r['noise'], 12)}")
        if str(r.get("normalization", "")) not in ("", "nan") and "none" not in str(r.get("normalization", "")).lower():
            kp_parts.append(f"norm={compact(r['normalization'], 14)}")
        kp = "; ".join(kp_parts) if kp_parts else "—"
        wt = f"{wall_mean.get(m, np.nan):.1f}" if m in wall_mean else "—"
        cell_text.append([m, compact(r["language"], 22), compact(r["neighbour_basis"], 30),
                          compact(kp, 46), str(r["n_free_params"]), wt, compact(SEED_NOTE, 40)])

    tbl = ax.table(cellText=cell_text, colLabels=cols, loc="center",
                   cellLoc="left", colLoc="center", colWidths=[0.09, 0.10, 0.16, 0.24, 0.07, 0.10, 0.24])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.4)
    tbl.scale(1.0, 1.7)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor("#0072B2")
            cell.set_text_props(color="white", fontweight="bold", ha="center")
        else:
            cell.get_text().set_horizontalalignment("left")
            if r % 2 == 0:
                cell.set_facecolor("#f5f7fa")
            # method-name col bold.
            if c == 0:
                cell.set_text_props(fontweight="bold")
            # # free params col centre.
            if c in (4, 5):
                cell.get_text().set_horizontalalignment("center")

    ax.set_title("Supplementary Figure S5 — Benchmark parameters, free-parameter count, and runtime\n"
                 "All five methods run on identical hardware (S91); same held-out mask per dataset. "
                 "Seed = 42; single run (deterministic for mean / spatial-KNN).",
                 loc="left", fontsize=9.5, pad=14)
    panel_label(ax, "A", x=-0.02, y=1.0)
    save_supp(fig, "supp_fig05_benchmark_parameters.png")
    print("[S5] done", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Supplementary Figure S1: Dataset overview.

Panel A: bar chart of PAS count per sample (32 samples across 8 GSE series),
        coloured by tissue type.
Panel B: summary table (species / platform / tissue / spot count / sample count)
        per GSE series.

Data: data/processed/*_scapatrap/qc_summary.json

Note: the supp brief stated "9 GSE x 32 samples"; the processed directory
contains 8 GSE series x 32 samples. We report the data as it exists.
"""
import os, sys, glob, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, BLACK, EXTRA_COLORS, save_supp

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
setup_rc()


def load():
    rows = []
    for p in sorted(glob.glob(os.path.join(ROOT, "data/processed/*_scapatrap/qc_summary.json"))):
        d = json.load(open(p))
        gse = d["dataset"].split("_gsm")[0].upper()
        # Normalise tissue to a short label for colouring.
        tiss_raw = d.get("tissue", "unknown")
        tiss = tiss_raw.lower()
        if "brain" in tiss or "metast" in tiss or "glioma" in tiss or "nf1" in tiss:
            tissue = "brain"
        elif "kidney" in tiss:
            tissue = "kidney"
        elif "colon" in tiss:
            tissue = "colon"
        elif "liver" in tiss:
            tissue = "liver"
        elif "skin" in tiss or "psoriasis" in tiss:
            tissue = "skin"
        else:
            tissue = "other"
        rows.append({
            "dataset": d["dataset"], "gse": gse,
            "species": d.get("species", "?"),
            "platform": d.get("platform", "?"),
            "tissue_raw": tiss_raw, "tissue": tissue,
            "n_spots": d.get("n_spots", 0),
            "n_pas": d.get("n_called_sites", 0),
            "n_apa_sites": d.get("n_apa_usage_sites", 0),
            "gsm": d["dataset"].split("_gsm")[1].split("_")[0] if "_gsm" in d["dataset"] else d["dataset"],
        })
    return pd.DataFrame(rows)


def main():
    df = load()
    print(f"[S1] {len(df)} samples across {df['gse'].nunique()} GSE series", flush=True)

    TISSUE_CMAP = {
        "brain":  ORANGE,
        "kidney": BLUE,
        "colon":  GREEN,
        "liver":  SKYBLU,
        "skin":   RED,
        "other":  GREY,
    }

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.34)

    # ---------- Panel A: PAS count per sample ----------
    axA = fig.add_subplot(gs[0])
    order = np.arange(len(df))
    colors = [TISSUE_CMAP[t] for t in df["tissue"]]
    axA.bar(order, df["n_pas"], color=colors, edgecolor="white", linewidth=0.4, width=0.82)
    # GSE group separators + labels.
    gse_groups = df.groupby("gse")
    gse_starts = []
    gse_labels = []
    for gse, sub in gse_groups:
        idx = sub.index.tolist()
        gse_starts.append((min(idx) + max(idx)) / 2.0)
        gse_labels.append(gse)
        if min(idx) > 0:
            axA.axvline(min(idx) - 0.5, color="#bbb", lw=0.6, ls="--")
    axA.set_xticks(gse_starts)
    axA.set_xticklabels(gse_labels, rotation=0, fontsize=7.5)
    axA.set_ylabel("PAS count (scAPAtrap called sites)")
    axA.set_xlabel("Sample (grouped by GSE series)")
    axA.set_title("PAS count per sample across the 32-sample benchmark", loc="left")
    axA.set_xlim(-0.7, len(df) - 0.3)
    handles = [Patch(facecolor=c, label=t.capitalize()) for t, c in TISSUE_CMAP.items() if t in df["tissue"].unique()]
    axA.legend(handles=handles, title="Tissue", loc="upper right", ncol=2, title_fontsize=7.5)
    panel_label(axA, "A", x=-0.04, y=1.04)

    # ---------- Panel B: per-GSE summary table ----------
    axB = fig.add_subplot(gs[1])
    axB.axis("off")
    # Aggregate per GSE.
    agg = df.groupby("gse").agg(
        n_samples=("dataset", "count"),
        species=("species", lambda s: ", ".join(sorted(set(s)))),
        tissue=("tissue_raw", lambda s: ", ".join(sorted(set(s)))),
        platform=("platform", lambda s: ", ".join(sorted(set(s)))),
        spots_mean=("n_spots", "mean"),
        spots_total=("n_spots", "sum"),
        pas_mean=("n_pas", "mean"),
    ).reset_index()
    agg["spots_mean"] = agg["spots_mean"].round(0).astype(int)
    agg["pas_mean"] = agg["pas_mean"].round(0).astype(int)

    cols = ["GSE", "# Samples", "Species", "Tissue", "Platform", "Mean spots", "Total spots", "Mean PAS"]
    cell_text = []
    for _, r in agg.iterrows():
        tissue_short = r["tissue"] if len(r["tissue"]) < 38 else r["tissue"][:35] + "..."
        cell_text.append([r["gse"], str(r["n_samples"]), r["species"].capitalize(),
                          tissue_short, r["platform"], f"{r['spots_mean']:,}",
                          f"{r['spots_total']:,}", f"{r['pas_mean']:,}"])
    # Totals row.
    cell_text.append(["TOTAL", str(len(df)), "—", "—", "—",
                      f"{int(df['n_spots'].mean()):,}",
                      f"{int(df['n_spots'].sum()):,}",
                      f"{int(df['n_pas'].mean()):,}"])

    tbl = axB.table(cellText=cell_text, colLabels=cols, loc="upper center",
                    cellLoc="center", colLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.6)
    tbl.scale(1.0, 1.32)
    # Style header + totals row.
    n_rows = len(cell_text) + 1
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        cell.set_linewidth(0.5)
        if r == 0:
            cell.set_facecolor("#0072B2")
            cell.set_text_props(color="white", fontweight="bold")
        elif r == n_rows - 1:
            cell.set_facecolor("#e8eef4")
            cell.set_text_props(fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f5f7fa")
        # left-align text cols.
        if c in (2, 3, 4):
            cell._loc = "left"
            cell.get_text().set_horizontalalignment("left")
            cell.PAD = 0.04
    axB.set_title("Per-series summary (species / platform / tissue / spot & PAS counts)", loc="center", pad=12)
    panel_label(axB, "B", x=-0.04, y=1.04)

    fig.suptitle("Supplementary Figure S1 — Dataset overview (8 GSE series, 32 samples)",
                 fontsize=11, fontweight="bold", y=0.995)
    save_supp(fig, "supp_fig01_dataset_overview.png")
    print("[S1] done", flush=True)


if __name__ == "__main__":
    main()

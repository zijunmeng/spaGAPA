"""Figure 6: Scalability (6 panels).

Panel F = ablation, moved to supplementary per spec.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import *
import pandas as pd
import numpy as np
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch, FancyArrowPatch

setup_rc()

DATA = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
rt = pd.read_csv(f"{DATA}/benchmark_runtime/runtime_table.csv")

# method display grouping
RUN_COLORS = {
    "spaGAPA-highres_fast": BLUE,
    "spaGAPA-highres_accuracy": SKYBLU,
    "stAPAminer": ORANGE,
    "spvAPA": RED,
}
RUN_LABELS = {
    "spaGAPA-highres_fast": "spaGAPA (fast)",
    "spaGAPA-highres_accuracy": "spaGAPA (accuracy)",
    "stAPAminer": "stAPAminer",
    "spvAPA": "spvAPA",
}
RUN_ORDER = ["spaGAPA-highres_fast", "spaGAPA-highres_accuracy", "stAPAminer", "spvAPA"]

fig = plt.figure(figsize=(13, 8))
gs = fig.add_gridspec(3, 4, hspace=0.60, wspace=0.40,
                       left=0.07, right=0.97, top=0.92, bottom=0.08)
axA = fig.add_subplot(gs[0, 0:2])
axB = fig.add_subplot(gs[0, 2:4])
axC = fig.add_subplot(gs[1, 0:2])
axD = fig.add_subplot(gs[1, 2:4])
axE = fig.add_subplot(gs[2, 0:4])

# ---------------- Panel A: complexity schematic ----------------
axA.set_title("Algorithmic complexity", loc="left")
n = np.linspace(500, 100000, 200)
# O(n^2) KNN, O(n M^2) sparse GP, O(n) mean
o2 = (n/1000)**2 * 0.001
ogp = (n/1000) * (0.15**2) * 8
omean = (n/1000) * 0.5
axA.plot(n/1000, o2, color=RED, lw=2.2, label="O(N²) KNN\n(stAPAminer, spvAPA)")
axA.plot(n/1000, ogp, color=BLUE, lw=2.2, label="O(N·M²) sparse GP\n(spaGAPA, M=150)")
axA.plot(n/1000, omean, color=GREY, lw=1.8, ls="--", label="O(N) mean")
axA.set_yscale("log")
axA.set_xlabel("Spots N (×1000)")
axA.set_ylabel("Relative cost (log)")
axA.legend(loc="upper left", fontsize=6.8)
axA.axvspan(42, 100, alpha=0.06, color=ORANGE)
yl = axA.get_ylim()
axA.text(71, yl[0]*2.0, "42k-100k\n(fails for O(N²))", fontsize=6.5, ha="center", color=ORANGE, fontweight="bold")
panel_label(axA, "A")

# ---------------- Panel B: runtime scaling (log-log) ----------------
axB.set_title("Runtime scaling (log-log)", loc="left")
for m in RUN_ORDER:
    sub = rt[(rt.method==m) & (rt.status=="COMPLETED")].sort_values("n_spots")
    if len(sub) == 0:
        continue
    axB.plot(sub["n_spots"], sub["time_s"], "o-", color=RUN_COLORS[m], lw=2,
             ms=7, markeredgecolor="white", label=RUN_LABELS[m])
# mark failures
for m in RUN_ORDER:
    fail = rt[(rt.method==m) & (rt.status!="COMPLETED") & (rt.status!="SKIPPED")]
    for _, row in fail.iterrows():
        marker = "x" if row["status"]=="TIMEOUT" else "v"
        axB.scatter([row["n_spots"]], [row["time_s"]], marker=marker, s=80,
                    color=RUN_COLORS[m], edgecolor="black", lw=0.7, zorder=6)
axB.set_xscale("log"); axB.set_yscale("log")
axB.set_xlabel("Spots N")
axB.set_ylabel("Wall time (s)")
axB.axhline(1200, color=RED, ls=":", lw=1.0, alpha=0.7)
axB.text(1100, 1300, "1200s wall cap", fontsize=6.5, color=RED)
fail_legend = [axB.scatter([],[], marker="x", color="black", s=60, label="TIMEOUT"),
               axB.scatter([],[], marker="v", color="black", s=60, label="FAILED")]
handles, labels = axB.get_legend_handles_labels()
axB.legend(handles=handles+fail_legend, loc="upper left", fontsize=6.3)
panel_label(axB, "B")

# ---------------- Panel C: peak memory scaling ----------------
axC.set_title("Peak memory scaling", loc="left")
for m in RUN_ORDER:
    sub = rt[(rt.method==m) & (rt.status=="COMPLETED")].sort_values("n_spots")
    if len(sub)==0: continue
    axC.plot(sub["n_spots"], sub["peak_mem_mb"]/1024, "s-", color=RUN_COLORS[m], lw=2,
             ms=7, markeredgecolor="white", label=RUN_LABELS[m])
axC.set_xscale("log")
axC.set_xlabel("Spots N")
axC.set_ylabel("Peak memory (GB)")
axC.legend(loc="upper left", fontsize=6.3)
panel_label(axC, "C")

# ---------------- Panel D: success/failure matrix ----------------
axD.set_title("Completion matrix", loc="left")
scales = sorted(rt["n_spots"].unique())
mat = np.full((len(RUN_ORDER), len(scales)), np.nan)
status_mat = np.empty(mat.shape, dtype=object)
for i, m in enumerate(RUN_ORDER):
    for j, sc in enumerate(scales):
        row = rt[(rt.method==m) & (rt.n_spots==sc)]
        if len(row):
            s = row["status"].iloc[0]
            status_mat[i,j] = s
            if s == "COMPLETED": mat[i,j] = 1
            elif s == "SKIPPED": mat[i,j] = 0.5
            else: mat[i,j] = 0  # timeout/failed
# custom colormap: red (fail) -> yellow (skip) -> green (complete)
from matplotlib.colors import ListedColormap
cmap = ListedColormap([RED, "#F0E442", GREEN])
im = axD.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
axD.set_yticks(range(len(RUN_ORDER)))
axD.set_yticklabels([RUN_LABELS[m] for m in RUN_ORDER], fontsize=7)
axD.set_xticks(range(len(scales)))
axD.set_xticklabels([f"{s//1000}k" for s in scales], fontsize=7)
axD.set_xlabel("Scale")
# annotate each cell
label_map = {"COMPLETED":"✓", "SKIPPED":"—", "TIMEOUT":"✗", "FAILED(rc=1)":"✗"}
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        s = status_mat[i,j]
        if s:
            txt = "✓" if s=="COMPLETED" else ("—" if s=="SKIPPED" else "✗")
            axD.text(j, i, txt, ha="center", va="center", fontsize=9, fontweight="bold",
                     color="white" if s=="COMPLETED" else "black")
leg = [Patch(fc=GREEN, label="Completed"),
       Patch(fc="#F0E442", label="Skipped (out of scope)"),
       Patch(fc=RED, label="Timeout / OOM")]
axD.legend(handles=leg, loc="upper right", fontsize=6, bbox_to_anchor=(1.0, -0.12), ncol=3)
panel_label(axD, "D")

# ---------------- Panel E: accuracy-runtime Pareto ----------------
axE.set_title("Accuracy-runtime Pareto frontier", loc="left")
# Use benchmark_mean_transparent for accuracy, runtime_table for time at ~4k spots
tc = pd.read_csv(f"{DATA}/benchmark_mean_transparent/transparent_comparison.csv")
tc = tc[tc.method.isin(["spaGAPA-GP","mean","stAPAminer","spvAPA","spatial-KNN"])].dropna(subset=["spatial_fidelity","wall_time_s"])
# average across the 2 datasets
agg = tc.groupby("method").agg({"spatial_fidelity":"mean","wall_time_s":"mean","pearson":"mean"}).reset_index()
# map to runtime table method names for scaling note
pareto_colors = {"spaGAPA-GP":BLUE, "mean":GREY, "stAPAminer":ORANGE, "spvAPA":RED, "spatial-KNN":GREEN}
for _, row in agg.iterrows():
    m = row["method"]
    axE.scatter(row["wall_time_s"], row["spatial_fidelity"], s=140, color=pareto_colors.get(m, BLACK),
                edgecolor="black", lw=0.8, label={"spaGAPA-GP":"spaGAPA-GP","mean":"Mean","stAPAminer":"stAPAminer","spvAPA":"spvAPA","spatial-KNN":"Spatial-KNN"}[m], zorder=5)
    offy = 0.03
    axE.annotate({"spaGAPA-GP":"spaGAPA-GP","mean":"Mean","stAPAminer":"stAPAminer","spvAPA":"spvAPA","spatial-KNN":"Spatial-KNN"}[m],
                 (row["wall_time_s"], row["spatial_fidelity"]),
                 xytext=(5, 5 if m!="spvAPA" else -12), textcoords="offset points", fontsize=7, fontweight="bold")
# shade Pareto region (frontier = mean fast+accurate-low-spatial, spaGAPA, spatial-KNN)
axE.set_xlabel("Wall time (s, 2 datasets avg)")
axE.set_ylabel("Spatial fidelity")
axE.set_xscale("log")
axE.set_xlim(1, 300)
axE.axhline(0, color="#888", lw=0.5, ls="--")
axE.legend(loc="center right", fontsize=7)
axE.text(0.5, -0.18, "At matched accuracy, spaGAPA-GP is 5-7× faster than stAPAminer/spvAPA;\nmean/spatial-KNN are fast but recover no spatial gradient",
         transform=axE.transAxes, ha="center", fontsize=6.8, style="italic", color="#444")
panel_label(axE, "E")

# Panel F note (ablation moved to supp)
fig.text(0.5, 0.005, "Panel F (inducing-point ablation) moved to Supplementary per reviewer spec.",
         ha="center", fontsize=6.8, style="italic", color="#888")

save(fig, "fig6_scalability.png")
print("Figure 6 done")

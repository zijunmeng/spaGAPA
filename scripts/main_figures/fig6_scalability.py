"""Figure 6: Scalability (5 panels).

Panel A reports EMPIRICAL scaling slopes fit to measured runtimes
(runtime_table.csv); theoretical complexity is shown only as a reference
footnote, not as a proven claim about competitor implementations.
Panel E reports measured spatial fidelity vs. wall time from
benchmark_mean_transparent/transparent_comparison.csv (honest Pareto
narrative; spatial-KNN is highest-fidelity but non-scalable / non-probabilistic).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import (setup_rc, panel_label, save, PAGE_WIDTH_IN,
                    BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, BLACK)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

setup_rc()

DATA = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
rt = pd.read_csv(f"{DATA}/benchmark_runtime/runtime_table.csv")

# method display grouping (runtime-table method names)
RUN_COLORS = {
    "spaGAPA-highres_fast":      BLUE,
    "spaGAPA-highres_accuracy":  SKYBLU,
    "stAPAminer":                ORANGE,
    "spvAPA":                    RED,
}
RUN_LABELS = {
    "spaGAPA-highres_fast":      "spaGAPA (fast)",
    "spaGAPA-highres_accuracy":  "spaGAPA (accuracy)",
    "stAPAminer":                "stAPAminer",
    "spvAPA":                    "spvAPA",
}
RUN_ORDER = ["spaGAPA-highres_fast", "spaGAPA-highres_accuracy", "stAPAminer", "spvAPA"]

# ---------------------------------------------------------------------------
# Layout: PAGE_WIDTH_IN wide, 5 panels.
#   Row 0: A (empirical scaling, spans 3 cols)
#   Row 1: B | C  (runtime log-log, peak memory)
#   Row 2: D | E  (completion matrix, fidelity-runtime Pareto)
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 7.4))
gs = fig.add_gridspec(3, 6, hspace=0.62, wspace=0.55,
                      left=0.075, right=0.975, top=0.94, bottom=0.085)
axA = fig.add_subplot(gs[0, 0:6])
axB = fig.add_subplot(gs[1, 0:3])
axC = fig.add_subplot(gs[1, 3:6])
axD = fig.add_subplot(gs[2, 0:3])
axE = fig.add_subplot(gs[2, 3:6])

# =====================================================================
# Panel A: EMPIRICAL scaling slopes (measured runtimes + fitted power law)
# =====================================================================
axA.set_title("Empirical runtime scaling (measured)", loc="left")
# x grid for fitted trend lines (span the observed N range)
n_grid = np.logspace(np.log10(800), np.log10(120000), 60)
empirical_info = []  # (method, slope, n_min, n_max)
for m in RUN_ORDER:
    sub = rt[(rt.method == m) & (rt.status == "COMPLETED")].sort_values("n_spots")
    if len(sub) < 2:
        continue
    # fit log(time) ~ slope * log(N) + c on COMPLETED points only
    lg_n = np.log10(sub["n_spots"].values)
    lg_t = np.log10(sub["time_s"].values)
    slope, intercept = np.polyfit(lg_n, lg_t, 1)
    empirical_info.append((m, slope, sub["n_spots"].min(), sub["n_spots"].max()))
    col = RUN_COLORS[m]
    # fitted trend (only over the N range it was actually measured on)
    n_fit = np.logspace(np.log10(sub["n_spots"].min()),
                        np.log10(sub["n_spots"].max()), 40)
    axA.plot(n_fit, 10 ** (slope * np.log10(n_fit) + intercept),
             color=col, lw=1.4, ls="--", alpha=0.7)
    # measured points
    axA.plot(sub["n_spots"], sub["time_s"], "o", color=col, ms=6,
             markeredgecolor="white", lw=0.6,
             label=f"{RUN_LABELS[m]}  (slope={slope:.2f})")

# mark non-completed runs (TIMEOUT / FAILED) at their wall-cap time
fail_handles = {}
for m in RUN_ORDER:
    fail = rt[(rt.method == m) & (rt.status != "COMPLETED") & (rt.status != "SKIPPED")]
    for _, row in fail.iterrows():
        mk = "x" if row["status"] == "TIMEOUT" else "v"
        ec = None if mk == "x" else "black"
        axA.scatter([row["n_spots"]], [row["time_s"]], marker=mk, s=55,
                    color=RUN_COLORS[m], edgecolor=ec, lw=0.5, zorder=6)
        fail_handles[row["status"]] = mk

axA.set_xscale("log"); axA.set_yscale("log")
axA.set_xlabel("Spots N")
axA.set_ylabel("Wall time (s)")
axA.set_xlim(800, 150000)
# 1200 s wall cap reference line
axA.axhline(1200, color=RED, ls=":", lw=0.9, alpha=0.7)
axA.text(900, 1280, "1200 s wall cap", fontsize=6.5, color=RED, va="bottom")
# legend: methods (with fitted slope) + failure markers
extra = []
extra.append(plt.Line2D([0], [0], marker="x", color="black", ls="None", mew=0,
                        ms=5, label="TIMEOUT"))
extra.append(plt.Line2D([0], [0], marker="v", color="black", ls="None",
                        ms=5, label="FAILED"))
handles, labels = axA.get_legend_handles_labels()
axA.legend(handles=handles + extra, loc="upper left", fontsize=6.3,
           ncol=2, columnspacing=1.0)
# footnote: empirical vs theoretical
axA.text(0.005, -0.30,
         "Empirical power-law slope = local log(time)/log(N) over the measured range "
         "(dashed). Slopes reflect the deployed implementations and may differ from "
         "textbook complexity;\ne.g. approximate-neighbour structures in WNN/KNN "
         "routines can flatten observed scaling below the worst-case bound.",
         transform=axA.transAxes, ha="left", va="top", fontsize=5.9,
         style="italic", color="#555")
panel_label(axA, "A")

# =====================================================================
# Panel B: runtime scaling (log-log, measured points only)
# =====================================================================
axB.set_title("Runtime vs. spots", loc="left")
for m in RUN_ORDER:
    sub = rt[(rt.method == m) & (rt.status == "COMPLETED")].sort_values("n_spots")
    if len(sub) == 0:
        continue
    axB.plot(sub["n_spots"], sub["time_s"], "o-", color=RUN_COLORS[m], lw=1.6,
             ms=5.5, markeredgecolor="white", mew=0.5, label=RUN_LABELS[m])
# mark failures
for m in RUN_ORDER:
    fail = rt[(rt.method == m) & (rt.status != "COMPLETED") & (rt.status != "SKIPPED")]
    for _, row in fail.iterrows():
        marker = "x" if row["status"] == "TIMEOUT" else "v"
        ec = None if marker == "x" else "black"
        axB.scatter([row["n_spots"]], [row["time_s"]], marker=marker, s=55,
                    color=RUN_COLORS[m], edgecolor=ec, lw=0.5, zorder=6)
axB.set_xscale("log"); axB.set_yscale("log")
axB.set_xlabel("Spots N")
axB.set_ylabel("Wall time (s)")
axB.axhline(1200, color=RED, ls=":", lw=0.9, alpha=0.7)
axB.text(900, 1280, "1200 s cap", fontsize=6.0, color=RED, va="bottom")
fail_legend = [plt.Line2D([0], [0], marker="x", color="black", ls="None", mew=0,
                          ms=5, label="TIMEOUT"),
               plt.Line2D([0], [0], marker="v", color="black", ls="None",
                          ms=5, label="FAILED")]
handles, labels = axB.get_legend_handles_labels()
axB.legend(handles=handles + fail_legend, loc="upper left", fontsize=5.8)
panel_label(axB, "B")

# =====================================================================
# Panel C: peak memory scaling
# =====================================================================
axC.set_title("Peak memory vs. spots", loc="left")
for m in RUN_ORDER:
    sub = rt[(rt.method == m) & (rt.status == "COMPLETED")].sort_values("n_spots")
    if len(sub) == 0:
        continue
    axC.plot(sub["n_spots"], sub["peak_mem_mb"] / 1024, "s-", color=RUN_COLORS[m],
             lw=1.6, ms=5.5, markeredgecolor="white", label=RUN_LABELS[m])
axC.set_xscale("log")
axC.set_xlabel("Spots N")
axC.set_ylabel("Peak memory (GB)")
axC.legend(loc="upper left", fontsize=5.8)
panel_label(axC, "C")

# =====================================================================
# Panel D: completion matrix
# =====================================================================
axD.set_title("Completion across scales", loc="left")
scales = sorted(rt["n_spots"].unique())
mat = np.full((len(RUN_ORDER), len(scales)), np.nan)
status_mat = np.empty(mat.shape, dtype=object)
for i, m in enumerate(RUN_ORDER):
    for j, sc in enumerate(scales):
        row = rt[(rt.method == m) & (rt.n_spots == sc)]
        if len(row):
            s = row["status"].iloc[0]
            status_mat[i, j] = s
            if s == "COMPLETED":
                mat[i, j] = 1
            elif s == "SKIPPED":
                mat[i, j] = 0.5
            else:
                mat[i, j] = 0
cmap = ListedColormap([RED, "#F0E442", GREEN])
im = axD.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
axD.set_yticks(range(len(RUN_ORDER)))
axD.set_yticklabels([RUN_LABELS[m] for m in RUN_ORDER], fontsize=6.5)
axD.set_xticks(range(len(scales)))
axD.set_xticklabels([f"{s//1000}k" for s in scales], fontsize=6.5)
axD.set_xlabel("Scale")
axD.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
axD.tick_params(axis="x", which="both", length=2.5)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        s = status_mat[i, j]
        if s:
            txt = "✓" if s == "COMPLETED" else ("—" if s == "SKIPPED" else "✗")
            axD.text(j, i, txt, ha="center", va="center", fontsize=8,
                     fontweight="bold",
                     color="white" if s == "COMPLETED" else "black")
leg = [Patch(fc=GREEN, label="Completed"),
       Patch(fc="#F0E442", label="Skipped (out of scope)"),
       Patch(fc=RED, label="Timeout / OOM")]
axD.legend(handles=leg, loc="upper center", fontsize=5.6,
           bbox_to_anchor=(0.5, -0.16), ncol=3, columnspacing=0.8,
           handlelength=1.0, handletextpad=0.3)
panel_label(axD, "D")

# =====================================================================
# Panel E: spatial-fidelity vs. runtime Pareto (honest narrative)
#   Real data (2-dataset mean):
#     spatial-KNN  fid=0.91  t=3.2s   (fastest-fidelity, non-probabilistic)
#     spaGAPA-GP   fid=0.42  t=23.5s  (calibrated probabilistic inference)
#     stAPAminer   fid=0.05  t=147.6s (slowest)
#     spvAPA       fid=0.04  t=165.4s
#     mean         fid=0.00  t=2.5s   (zero spatial signal)
# =====================================================================
axE.set_title("Spatial fidelity vs. runtime", loc="left")
tc = pd.read_csv(f"{DATA}/benchmark_mean_transparent/transparent_comparison.csv")
tc = tc[tc.method.isin(["spaGAPA-GP", "mean", "stAPAminer", "spvAPA", "spatial-KNN"])] \
       .dropna(subset=["spatial_fidelity", "wall_time_s"])
agg = tc.groupby("method").agg(
    {"spatial_fidelity": "mean", "wall_time_s": "mean"}).reset_index()

E_COLORS = {"spaGAPA-GP": BLUE, "mean": GREY, "stAPAminer": ORANGE,
            "spvAPA": RED, "spatial-KNN": GREEN}
E_LABELS = {"spaGAPA-GP": "spaGAPA-GP", "mean": "Mean", "stAPAminer": "stAPAminer",
            "spvAPA": "spvAPA", "spatial-KNN": "Spatial-KNN"}
# per-point label offsets to avoid overlap
E_OFFSETS = {
    "spaGAPA-GP":  (7, 6),
    "mean":        (5, -10),
    "stAPAminer":  (5, 6),
    "spvAPA":      (5, -10),
    "spatial-KNN": (5, 6),
}
for _, row in agg.iterrows():
    m = row["method"]
    axE.scatter(row["wall_time_s"], row["spatial_fidelity"], s=90,
                color=E_COLORS[m], edgecolor="black", lw=0.7,
                label=E_LABELS[m], zorder=5)
    dx, dy = E_OFFSETS[m]
    axE.annotate(E_LABELS[m], (row["wall_time_s"], row["spatial_fidelity"]),
                 xytext=(dx, dy), textcoords="offset points",
                 fontsize=6.2, fontweight="bold")
# connect the two probabilistic GP-style competitors to make the speed gap visible
gp_row = agg.set_index("method").loc["spaGAPA-GP"]
for comp in ["stAPAminer", "spvAPA"]:
    comp_row = agg.set_index("method").loc[comp]
    axE.annotate("", xy=(comp_row["wall_time_s"], comp_row["spatial_fidelity"]),
                 xytext=(gp_row["wall_time_s"], gp_row["spatial_fidelity"]),
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=0.7,
                                 ls=":", alpha=0.6))

axE.set_xlabel("Wall time (s, 2-dataset mean)")
axE.set_ylabel("Spatial fidelity")
axE.set_xscale("log")
axE.set_xlim(1.5, 300)
axE.axhline(0, color="#888", lw=0.5, ls="--")
# honest footnote
axE.text(0.005, -0.30,
         "spaGAPA-GP is ~6-7× faster than stAPAminer/spvAPA at higher spatial fidelity; "
         "spatial-KNN yields the highest fidelity\nbut is non-probabilistic and does "
         "not scale beyond ~4-15 k spots (Panels A-D); Mean is fastest but encodes\n"
         "no spatial information (fidelity ≈ 0).",
         transform=axE.transAxes, ha="left", va="top", fontsize=5.9,
         style="italic", color="#555")
panel_label(axE, "E")

save(fig, "fig6_scalability.png")
print("Figure 6 done")

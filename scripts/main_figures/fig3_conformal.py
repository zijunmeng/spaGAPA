"""Figure 3: Conformal marginal coverage (5 panels, ~178 mm wide).

11 samples, 523,174 test points, 80/90/95% nominal. Pooled empirical coverage
0.8005 / 0.8996 / 0.9497; mean abs deviation 0.21% / 0.16% / 0.10% and max
deviation 0.49% / 0.47% / 0.24% respectively (so the 80% level is NOT strictly
within 0.2%).

Layout (Plan B, cleaner 5-panel main figure):
  A: split-conformal flow (full-width strip)
  B: per-sample coverage at 80/90/95% (grouped bars + binomial CI band)
  C: nominal-empirical calibration curve (per-sample thin + pooled bold)
  D: interval width — Constant (A) vs Residual-spot (D) boxplot
  E: per-spot spatial instance — posterior mu / sigma / |error| / covered-vs-miss
The per-sample coverage-deviation forest (former Panel D) moved to Supplementary
Figure S15 (`supp_fig15_conformal_deviation_forest.py`).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import *
import pandas as pd
import numpy as np
from scipy.stats import binom
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.gridspec as gridspec

setup_rc()

DATA = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
cov = pd.read_csv(f"{DATA}/conformal_validation/all_samples_coverage.csv")
cwidth = pd.read_csv(f"{DATA}/uncertainty_coverage_width_audit/coverage_width_comparison.csv")

# ---- Figure geometry: 7.0 in wide (~178 mm) x 8.6 in tall, 3 rows ----
# Row 0: Panel A (flow) full width, short. Row 1: B/C/D side by side.
# Row 2: Panel E (spatial) full width.
fig = plt.figure(figsize=(PAGE_WIDTH_IN, 8.6))
gs = fig.add_gridspec(
    3, 6,
    height_ratios=[1.05, 2.55, 2.55],
    hspace=0.62, wspace=0.95,
    left=0.115, right=0.955, top=0.955, bottom=0.055,
)
axA = fig.add_subplot(gs[0, 0:6])              # full-width flow strip
axB = fig.add_subplot(gs[1, 0:2])              # per-sample coverage
axC = fig.add_subplot(gs[1, 2:4])              # calibration curve
axD = fig.add_subplot(gs[1, 4:6])              # interval width
axE = fig.add_subplot(gs[2, 0:6])              # spatial instance (host)

# ---------------- Panel A: conformal flow (horizontal strip) ----------------
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.set_title("Split-conformal prediction", loc="left")

def box(ax, xy, w, h, text, fc, tc="white", fs=7.5):
    b = FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.008,rounding_size=0.02",
                       fc=fc, ec=fc, lw=1.0)
    ax.add_patch(b)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center",
            color=tc, fontsize=fs, fontweight="bold")

def arr(ax, p0, p1, c=BLACK, lw=1.3):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="->", mutation_scale=11,
                                 color=c, lw=lw, shrinkA=1, shrinkB=1))

# four boxes across the strip + an annotation line beneath
box(axA, (0.03, 0.42), 0.205, 0.34, "Calibration\nset\n(n_cal)", BLUE, fs=7.2)
box(axA, (0.285, 0.42), 0.205, 0.34, "Nonconformity\n|yᵢ − ŷᵢ|", ORANGE, fs=7)
box(axA, (0.545, 0.42), 0.175, 0.34, "Quantile\nqₜ", GREEN, fs=7.2)
box(axA, (0.765, 0.42), 0.205, 0.34, "Test\ncoverage", "#555", fs=7.2)
arr(axA, (0.235, 0.59), (0.285, 0.59))
arr(axA, (0.49, 0.59), (0.545, 0.59))
arr(axA, (0.72, 0.59), (0.765, 0.59))
axA.text(0.5, 0.20, "Interval:  ŷ ± qₜ   →   P(y ∈ interval) ≥ 1−α",
         ha="center", fontsize=8.2, color="#222", fontweight="bold")
axA.text(0.5, 0.07, "11 samples · 523,174 test points",
         ha="center", fontsize=7.5, color=GREEN)
panel_label(axA, "A", x=-0.005, y=1.06)

# ---------------- Panel B: per-sample coverage at 80/90/95 ----------------
axB.set_title("Marginal coverage per sample", loc="left")
targets = [(0.80, "coverage_80", BLUE), (0.90, "coverage_90", ORANGE),
           (0.95, "coverage_95", GREEN)]
xs = np.arange(len(cov))
width = 0.26
mean_n = cov["n_test"].mean()
for k, (nom, col_name, col) in enumerate(targets):
    off = (k - 1) * width
    vals = cov[col_name].values
    axB.bar(xs + off, vals - nom, width, bottom=nom, color=col,
            edgecolor="white", lw=0.4, label=f"{int(nom*100)}% nominal")
    axB.scatter(xs + off, vals, color="black", s=11, zorder=4)
axB.axhline(0.80, color=BLUE, lw=0.8, ls=":", alpha=0.6)
axB.axhline(0.90, color=ORANGE, lw=0.8, ls=":", alpha=0.6)
axB.axhline(0.95, color=GREEN, lw=0.8, ls=":", alpha=0.6)
axB.set_xticks(xs)
axB.set_xticklabels([s.split("gsm")[1] if "gsm" in s else s
                     for s in cov["sample"].values],
                    rotation=45, ha="right", fontsize=6.6)
axB.set_ylabel("Empirical coverage")
axB.set_xlabel("Sample")
axB.set_ylim(0.77, 0.965)
axB.legend(loc="lower right", fontsize=6.8, title="nominal", title_fontsize=6.5)
# binomial CI band for 90% (use mean n_test)
ci = binom.ppf([0.025, 0.975], mean_n, 0.90) / mean_n
axB.axhspan(ci[0], ci[1], color=ORANGE, alpha=0.10)
axB.text(10.3, 0.902, "90% binomial CI", fontsize=5.9, color=ORANGE, va="center")
panel_label(axB, "B", x=-0.16)

# ---------------- Panel C: nominal-empirical calibration curve ----------------
axC.set_title("Calibration curve: nominal vs empirical", loc="left")
# per-sample thin lines (we only have 80/90/95; plot those + connect)
for _, row in cov.iterrows():
    pts_x = [0.80, 0.90, 0.95]
    pts_y = [row["coverage_80"], row["coverage_90"], row["coverage_95"]]
    axC.plot(pts_x, pts_y, color=BLUE, alpha=0.18, lw=0.7, marker=".", ms=2)
# pooled means
pool_x = [0.80, 0.90, 0.95]
pool_y = [cov["coverage_80"].mean(), cov["coverage_90"].mean(),
          cov["coverage_95"].mean()]
axC.plot(pool_x, pool_y, color=BLUE, lw=2.6, marker="o", ms=8, label="Pooled mean",
         markeredgecolor="white", zorder=5)
axC.plot([0, 1], [0, 1], color=GREY, ls="--", lw=1.0, label="Perfect (y=x)")
axC.set_xlim(0.74, 0.97); axC.set_ylim(0.78, 0.965)
axC.set_xlabel("Nominal coverage (1−α)")
axC.set_ylabel("Empirical coverage")
axC.legend(loc="upper left", fontsize=7)
# annotate: report MAD and max deviation per level, precisely (from summary.json).
# 80%: MAD 0.21%, max 0.49%; 90%: MAD 0.16%, max 0.47%; 95%: MAD 0.10%, max 0.24%.
# NB: 80% MAD = 0.21%, so a blanket "within 0.2%" claim does NOT hold for 80%.
axC.annotate(f"MAD / max dev vs nominal\n"
             f"80%: 0.21% / 0.49%\n"
             f"90%: 0.16% / 0.47%\n"
             f"95%: 0.10% / 0.24%",
             xy=(0.90, pool_y[1]), xytext=(0.775, 0.825), fontsize=6.6, color=BLUE,
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8))
panel_label(axC, "C", x=-0.16)

# ---------------- Panel D: interval width A vs D ----------------
axD.set_title("Interval width (A vs D)", loc="left")
methods_d = ["A_constant", "D_residual_spot"]
data_d = [cwidth[cwidth.method == m]["mean_interval_width"].values for m in methods_d]
bp = axD.boxplot(data_d, labels=["A: Constant", "D: Residual-spot"],
                 patch_artist=True, widths=0.5, showfliers=False)
for patch, m in zip(bp["boxes"], methods_d):
    patch.set_facecolor(NOISE_METHOD_COLORS[m]); patch.set_alpha(0.7)
    patch.set_edgecolor("white")
for med in bp["medians"]:
    med.set_color("black")
# overlay per-dataset dots
for i, m in enumerate(methods_d):
    vals = cwidth[cwidth.method == m]["mean_interval_width"].values
    xs_d = np.full(len(vals), i + 1) + np.random.RandomState(i).normal(0, 0.04, len(vals))
    axD.scatter(xs_d, vals, color="black", s=18, zorder=4, alpha=0.7)
    axD.set_rasterized(True)
axD.set_ylabel("Mean interval width\n(APA usage units)")
axD.set_xlabel("")
axD.tick_params(axis="x", labelsize=7)
panel_label(axD, "D", x=-0.16)

# ---------------- Panel E: spatial instance (REAL spaGAPA posterior) ----------------
# Per-spot data produced by scripts/main_figures/_fig3F_dump_perspot.py: runs the
# actual SparseGPImputer on one gene (peak_20919; Moran's I=+0.081, obs-frac=0.22)
# then split-conformal (20% held out, 50/50 calibrate/test, qhat = |y-yhat|
# quantile on calibration half). mu/sigma are the real spaGAPA posterior; the
# interval [mu-qhat, mu+qhat] is an honest split-conformal interval at 90%.
GENE_E = "peak_20919"
perspot_csv = f"{DATA}/main_figures/_data/fig3F_{GENE_E}_perspot.csv"
psd = pd.read_csv(perspot_csv)
x = psd["x"].values
y = psd["y"].values
mu = psd["posterior_mean"].values
sigma = psd["sigma"].values
abs_err = psd["abs_error"].values            # NaN where unobserved
observed = psd["observed"].values.astype(bool)
covered = psd["covered_90"].values.astype(bool)   # 1 only where observed & in interval
qhat = float(psd["qhat"].iloc[0])
n_obs = int(observed.sum())
n_cov = int(covered.sum())

axE.axis("off")
axE.set_title(f"Per-spot posterior + 90% conformal interval  "
              f"({GENE_E}, GSE183456)",
              loc="left", fontsize=8)
panel_label(axE, "E", x=-0.005, y=1.18)
pos = axE.get_position()
sub = gridspec.GridSpec(1, 4, left=pos.x0, right=pos.x1,
                        bottom=pos.y0 + 0.10, top=pos.y0 + pos.height - 0.16,
                        wspace=0.55)

# tile order: posterior mean | local uncertainty sigma | |error| (observed) | covered/miss
def _spatial_ax(i, title):
    axm = fig.add_subplot(sub[i])
    axm.set_xticks([]); axm.set_yticks([])
    axm.set_aspect("equal")
    for s in axm.spines.values():
        s.set_linewidth(0.5)
    axm.set_title(title, fontsize=7, fontweight="bold")
    axm.set_rasterized(True)
    return axm

# 1) posterior mean (ALL spots, incl. imputed)
vm_mu = np.percentile(mu, 99)
axm = _spatial_ax(0, "Posterior mean μ")
sc = axm.scatter(x, y, c=mu, s=3.0, cmap="viridis", vmin=0, vmax=vm_mu, rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.4)

# 2) local uncertainty sigma (per-spot posterior SD, ALL spots)
axm = _spatial_ax(1, "Posterior σ")
sc = axm.scatter(x, y, c=sigma, s=3.0, cmap="magma", vmin=0,
                 vmax=np.percentile(sigma, 99), rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.4)

# 3) |error| only on observed spots (truth known); unobserved = grey background
axm = _spatial_ax(2, "|observed − μ|")
axm.scatter(x, y, c=GREY, s=3.0, alpha=0.25, rasterized=True)   # all spots faint
vm_e = np.nanmax(abs_err) if np.isfinite(np.nanmax(abs_err)) else 0.001
sc = axm.scatter(x[observed], y[observed], c=abs_err[observed], s=3.2,
                 cmap="Reds", vmin=0, vmax=vm_e, rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=6); cb.outline.set_linewidth(0.4)

# 4) covered (green) vs miss (red) among observed; unobserved = grey
axm = _spatial_ax(3, f"Covered @90% (q̂={qhat:.2f})")
axm.scatter(x, y, c=GREY, s=3.0, alpha=0.20, rasterized=True)
miss = observed & (~covered)
axm.scatter(x[covered], y[covered], c=GREEN, s=3.6, rasterized=True, label="covered")
axm.scatter(x[miss],     y[miss],     c=RED,   s=6.0, rasterized=True, label="miss",
            marker="X")
axm.legend(loc="lower left", fontsize=6, markerscale=0.6, handletextpad=0.1)

axE.text(0.5, 0.02,
         f"Interval = μ ± q̂, q̂={qhat:.3f} (90% split-conformal)   ·   "
         f"observed coverage {n_cov}/{n_obs} = {n_cov/n_obs:.2f}",
         transform=axE.transAxes, ha="center", fontsize=6.6, style="italic", color="#444")

# NOTE on save(): _style.save() calls savefig with bbox_inches="tight", which
# expands the canvas to grab rotated tick labels / panel letters that poke past
# the axes — on this figure that inflated the PDF to ~219 mm. We instead save
# at the exact PAGE_WIDTH_IN so the PDF is exactly ~178 mm (BIB double column),
# having tuned the gridspec margins above so nothing is clipped.
_BASE = "fig3_conformal_marginal_coverage"
_OUTDIR = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/main_figures"
fig.savefig(f"{_OUTDIR}/{_BASE}.pdf", facecolor="white")            # vector, exact 7.0in
fig.savefig(f"{_OUTDIR}/{_BASE}.png", dpi=300, facecolor="white")   # preview
plt.close(fig)
print(f"  saved -> {_OUTDIR}/{_BASE}.pdf (+png)")
print("Figure 3 done")

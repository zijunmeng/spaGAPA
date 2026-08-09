"""Figure 3: Conformal marginal coverage (6 panels).

11 samples, 80/90/95% nominal. Pooled empirical coverage 0.8005 / 0.8996 /
0.9497; mean abs deviation 0.21% / 0.16% / 0.10% and max deviation 0.49% /
0.47% / 0.24% respectively (so the 80% level is NOT strictly within 0.2%).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import *
import pandas as pd
import numpy as np
from scipy.stats import binom
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch, FancyArrowPatch

setup_rc()

DATA = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output"
cov = pd.read_csv(f"{DATA}/conformal_validation/all_samples_coverage.csv")
cond_q = pd.read_csv(f"{DATA}/conformal_conditional_coverage/by_uncertainty_quintile.csv")
cond_r = pd.read_csv(f"{DATA}/conformal_conditional_coverage/by_spatial_region.csv")
cond_e = pd.read_csv(f"{DATA}/conformal_conditional_coverage/by_expression_level.csv")
cwidth = pd.read_csv(f"{DATA}/uncertainty_coverage_width_audit/coverage_width_comparison.csv")

fig = plt.figure(figsize=(13, 8.5))
gs = fig.add_gridspec(3, 6, hspace=0.60, wspace=0.85,
                       left=0.06, right=0.97, top=0.92, bottom=0.07)
axA = fig.add_subplot(gs[0, 0:3])
axB = fig.add_subplot(gs[0, 3:6])
axC = fig.add_subplot(gs[1, 0:3])
axD = fig.add_subplot(gs[1, 3:6])
axE = fig.add_subplot(gs[2, 0:3])
axF = fig.add_subplot(gs[2, 3:6])

# ---------------- Panel A: conformal flow ----------------
axA.set_xlim(0,1); axA.set_ylim(0,1); axA.axis("off")
axA.set_title("Split-conformal prediction", loc="left")
def box(ax, xy, w, h, text, fc, tc="white", fs=7.5):
    b = FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.008,rounding_size=0.02",
                       fc=fc, ec=fc, lw=1.0)
    ax.add_patch(b)
    ax.text(xy[0]+w/2, xy[1]+h/2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="bold")
def arr(ax, p0, p1, c=BLACK, lw=1.3):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="->", mutation_scale=11,
                                 color=c, lw=lw, shrinkA=1, shrinkB=1))
box(axA, (0.03, 0.62), 0.22, 0.20, "Calibration\nset\n(n cal)", BLUE, fs=7.2)
box(axA, (0.30, 0.62), 0.22, 0.20, "Nonconformity\n|yᵢ − ŷᵢ|", ORANGE, fs=7)
box(axA, (0.57, 0.62), 0.18, 0.20, "Quantile\nqₜ", GREEN, fs=7.2)
box(axA, (0.79, 0.62), 0.17, 0.20, "Test\ncoverage", "#555", fs=7.2)
arr(axA, (0.25, 0.72), (0.30, 0.72))
arr(axA, (0.52, 0.72), (0.57, 0.72))
arr(axA, (0.75, 0.72), (0.79, 0.72))
axA.text(0.5, 0.50, "Interval:  ŷ ± qₜ   →   P(y ∈ interval) ≥ 1−α",
         ha="center", fontsize=8, color="#222", fontweight="bold")
axA.text(0.5, 0.40, "11 samples · 523,174 test points", ha="center", fontsize=7.5, color=GREEN)
panel_label(axA, "A")

# ---------------- Panel B: per-sample coverage at 80/90/95 ----------------
axB.set_title("Marginal coverage per sample", loc="left")
targets = [(0.80, "coverage_80", BLUE), (0.90, "coverage_90", ORANGE), (0.95, "coverage_95", GREEN)]
xs = np.arange(len(cov))
width = 0.26
# pooled binomial CI (use mean n_test)
mean_n = cov["n_test"].mean()
for k, (nom, col_name, col) in enumerate(targets):
    off = (k-1)*width
    vals = cov[col_name].values
    axB.bar(xs+off, vals-nom, width, bottom=nom, color=col, edgecolor="white", lw=0.4,
            label=f"{int(nom*100)}% nominal")
    axB.scatter(xs+off, vals, color="black", s=10, zorder=4)
axB.axhline(0.80, color=BLUE, lw=0.8, ls=":", alpha=0.6)
axB.axhline(0.90, color=ORANGE, lw=0.8, ls=":", alpha=0.6)
axB.axhline(0.95, color=GREEN, lw=0.8, ls=":", alpha=0.6)
axB.set_xticks(xs)
axB.set_xticklabels([s.split("gsm")[1] if "gsm" in s else s for s in cov["sample"].values],
                     rotation=45, ha="right", fontsize=6.2)
axB.set_ylabel("Empirical coverage")
axB.set_xlabel("Sample")
axB.set_ylim(0.77, 0.965)
axB.legend(loc="lower right", fontsize=6.8, title="nominal", title_fontsize=6.5)
# binomial CI band for 90%
ci = binom.ppf([0.025, 0.975], mean_n, 0.90)/mean_n
axB.axhspan(ci[0], ci[1], color=ORANGE, alpha=0.10)
axB.text(10.3, 0.902, "90% binomial CI", fontsize=5.8, color=ORANGE, va="center")
panel_label(axB, "B")

# ---------------- Panel C: nominal-empirical calibration curve ----------------
axC.set_title("Calibration curve: nominal vs empirical", loc="left")
nominals = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
# we only have 80/90/95 from samples; use them plus interpolate a monotone curve
# per-sample thin lines then pooled bold
cols = {0.80:"coverage_80", 0.90:"coverage_90", 0.95:"coverage_95"}
for _, row in cov.iterrows():
    pts_x = [0.80, 0.90, 0.95]
    pts_y = [row["coverage_80"], row["coverage_90"], row["coverage_95"]]
    axC.plot(pts_x, pts_y, color=BLUE, alpha=0.18, lw=0.7, marker=".", ms=2)
# pooled means
pool_x = [0.80, 0.90, 0.95]
pool_y = [cov["coverage_80"].mean(), cov["coverage_90"].mean(), cov["coverage_95"].mean()]
axC.plot(pool_x, pool_y, color=BLUE, lw=2.6, marker="o", ms=8, label="Pooled mean",
         markeredgecolor="white", zorder=5)
axC.plot([0,1],[0,1], color=GREY, ls="--", lw=1.0, label="Perfect (y=x)")
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
             xy=(0.90, pool_y[1]), xytext=(0.775, 0.825), fontsize=6.4, color=BLUE,
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.8))
panel_label(axC, "C")

# ---------------- Panel D: coverage deviation forest ----------------
axD.set_title("Coverage deviation from nominal", loc="left")
# show 90% deviations per sample as forest plot
cov_sorted = cov.sort_values("coverage_90").reset_index(drop=True)
ys = np.arange(len(cov_sorted))
dev90 = cov_sorted["coverage_90"].values - 0.90
dev80 = cov_sorted["coverage_80"].values - 0.80
dev95 = cov_sorted["coverage_95"].values - 0.95
axD.hlines(ys, -0.05, 0.05, color="#EEE", lw=0.8)
axD.errorbar(dev90, ys, xerr=None, fmt="o", color=ORANGE, ms=6, label="90%")
axD.errorbar(dev80, ys-0.18, xerr=None, fmt="s", color=BLUE, ms=5, label="80%")
axD.errorbar(dev95, ys+0.18, xerr=None, fmt="^", color=GREEN, ms=5, label="95%")
axD.axvline(0, color=BLACK, lw=1.0)
axD.set_yticks(ys)
axD.set_yticklabels([f"{d}\n{s.split('gsm')[1]}" for d, s in zip(cov_sorted["dataset"], cov_sorted["sample"])],
                    fontsize=5.8)
axD.set_xlabel("Empirical − nominal coverage")
axD.set_xlim(-0.012, 0.012)
axD.legend(loc="lower right", fontsize=6.5, ncol=3)
panel_label(axD, "D")

# ---------------- Panel E: interval width A vs D ----------------
axE.set_title("Interval width: Constant vs Residual-spot", loc="left")
# boxplot of width per method across datasets (pooled by method)
methods_e = ["A_constant", "D_residual_spot"]
data_e = [cwidth[cwidth.method==m]["mean_interval_width"].values for m in methods_e]
bp = axE.boxplot(data_e, labels=["A: Constant\n(baseline)", "D: Residual-spot\n(spaGAPA)"],
                 patch_artist=True, widths=0.5, showfliers=False)
for patch, m in zip(bp["boxes"], methods_e):
    patch.set_facecolor(NOISE_METHOD_COLORS[m]); patch.set_alpha(0.7)
    patch.set_edgecolor("white")
for med in bp["medians"]: med.set_color("black")
# overlay per-dataset dots
for i, m in enumerate(methods_e):
    vals = cwidth[cwidth.method==m]["mean_interval_width"].values
    xs = np.full(len(vals), i+1) + np.random.RandomState(i).normal(0,0.04,len(vals))
    axE.scatter(xs, vals, color="black", s=18, zorder=4, alpha=0.7)
axE.set_ylabel("Mean interval width\n(APA usage units)")
axE.set_xlabel("")
axE.text(0.5, -0.22, "D ≈ same width as A → tighter,\nadaptive intervals without inflation",
         transform=axE.transAxes, ha="center", fontsize=6.8, style="italic", color="#444")
panel_label(axE, "E")

# ---------------- Panel F: spatial instance (REAL spaGAPA posterior) ----------------
# Per-spot data produced by scripts/main_figures/_fig3F_dump_perspot.py: runs the
# actual SparseGPImputer on one gene (peak_20919; Moran's I=+0.081, obs-frac=0.22)
# then split-conformal (20% held out, 50/50 calibrate/test, qhat = |y-yhat|
# quantile on calibration half). mu/sigma are the real spaGAPA posterior; the
# interval [mu-qhat, mu+qhat] is an honest split-conformal interval at 90%.
import matplotlib.gridspec as gridspec
GENE_F = "peak_20919"
perspot_csv = f"{DATA}/main_figures/_data/fig3F_{GENE_F}_perspot.csv"
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

axF.axis("off")
axF.set_title(f"Per-spot spaGAPA posterior + 90% conformal interval\n"
              f"{GENE_F}, GSE183456 (n={n_obs} observed spots)", loc="left", fontsize=8)
panel_label(axF, "F", x=-0.06, y=1.10)
pos = axF.get_position()
sub = gridspec.GridSpec(1, 4, left=pos.x0, right=pos.x1,
                        bottom=pos.y0+0.04, top=pos.y0+pos.height-0.08, wspace=0.05)

# tile order: posterior mean | local uncertainty sigma | |error| (observed) | covered/miss
def _spatial_ax(i, title):
    axm = fig.add_subplot(sub[i])
    axm.set_xticks([]); axm.set_yticks([])
    axm.set_aspect("equal")
    for s in axm.spines.values(): s.set_linewidth(0.5)
    axm.set_title(title, fontsize=6.8, fontweight="bold")
    axm.set_rasterized(True)
    return axm

# 1) posterior mean (ALL spots, incl. imputed)
vm_mu = np.percentile(mu, 99)
axm = _spatial_ax(0, "Posterior mean μ")
sc = axm.scatter(x, y, c=mu, s=2.6, cmap="viridis", vmin=0, vmax=vm_mu, rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=5.5); cb.outline.set_linewidth(0.4)

# 2) local uncertainty sigma (per-spot posterior SD, ALL spots)
axm = _spatial_ax(1, "Posterior σ")
sc = axm.scatter(x, y, c=sigma, s=2.6, cmap="magma", vmin=0,
                 vmax=np.percentile(sigma, 99), rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=5.5); cb.outline.set_linewidth(0.4)

# 3) |error| only on observed spots (truth known); unobserved = grey background
axm = _spatial_ax(2, "|observed − μ|")
axm.scatter(x, y, c=GREY, s=2.6, alpha=0.25, rasterized=True)   # all spots faint
vm_e = np.nanmax(abs_err) if np.isfinite(np.nanmax(abs_err)) else 0.001
sc = axm.scatter(x[observed], y[observed], c=abs_err[observed], s=2.8,
                 cmap="Reds", vmin=0, vmax=vm_e, rasterized=True)
cb = fig.colorbar(sc, ax=axm, fraction=0.046, pad=0.04, shrink=0.85)
cb.ax.tick_params(labelsize=5.5); cb.outline.set_linewidth(0.4)

# 4) covered (green) vs miss (red) among observed; unobserved = grey
axm = _spatial_ax(3, f"Covered @90% (q̂={qhat:.2f})")
axm.scatter(x, y, c=GREY, s=2.6, alpha=0.20, rasterized=True)
miss = observed & (~covered)
axm.scatter(x[covered], y[covered], c=GREEN, s=3.2, rasterized=True, label="covered")
axm.scatter(x[miss],     y[miss],     c=RED,   s=5.0, rasterized=True, label="miss",
            marker="X")
axm.legend(loc="lower left", fontsize=5.5, markerscale=0.6, handletextpad=0.1)

axF.text(0.5, 0.005,
         f"Interval = μ ± q̂, q̂={qhat:.3f} (90% split-conformal)   ·   "
         f"observed coverage {n_cov}/{n_obs} = {n_cov/n_obs:.2f}",
         transform=axF.transAxes, ha="center", fontsize=6.2, style="italic", color="#444")

save(fig, "fig3_conformal_marginal_coverage.png")
print("Figure 3 done")

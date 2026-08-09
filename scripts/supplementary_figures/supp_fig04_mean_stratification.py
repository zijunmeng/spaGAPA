#!/usr/bin/env python
"""Supplementary Figure S4: Mean-baseline stratification.

Box plot of (GP RMSE − mean RMSE) per gene, stratified by the gene's spatial
signal (Moran's I quintile). On these sparse Visium APA data the per-gene mean
is a strong RMSE baseline: GP beats it in only a minority of genes (~15%), and
the GP−mean RMSE gap does not shrink with spatial signal — if anything the
high-spatial-signal quintile shows the *largest* positive median Δ (GP worse),
so the GP advantage does NOT grow with spatial signal.

Pipeline (cached to _cache/s4_stratification.csv):
  - Load GSE183456 apa matrix + coordinates.
  - Sample N_TEST genes with enough observations; mask 20% of observed entries.
  - Compute Moran's I for each gene (on observed entries, spatial weights from
    k-nearest spot neighbours).
  - Run GP and mean imputation; compute per-gene held-out RMSE.
  - Stratify Δ = RMSE(GP) − RMSE(mean) by Moran's I quintile.
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import setup_rc, panel_label, BLUE, ORANGE, GREEN, SKYBLU, RED, GREY, save_supp

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
CACHE = os.path.join(ROOT, "pipeline_output/supplementary_figures/_cache", "s4_stratification.csv")
setup_rc()


def morans_i(values, W):
    """Global Moran's I for a 1D array of values with precomputed sparse weight
    matrix W (scipy.sparse). Uses only finite entries."""
    finite = np.isfinite(values)
    if finite.sum() < 5:
        return np.nan
    x = values[finite]
    xbar = x.mean()
    z = x - xbar
    denom = (z ** 2).sum()
    if denom == 0:
        return 0.0
    Wf = W[finite][:, finite]
    s0 = Wf.sum()
    if s0 == 0:
        return np.nan
    z2 = z ** 2
    # moran = (N / S0) * (z' W z) / (z' z)
    n = finite.sum()
    cross = (Wf @ z) @ z
    return (n / s0) * (cross / denom)


def main():
    if os.path.exists(CACHE):
        print(f"[S4] cache hit -> {CACHE}", flush=True)
        res = pd.read_csv(CACHE)
    else:
        from scipy.sparse import csr_matrix
        from sklearn.neighbors import NearestNeighbors
        from spagapa import SparseGPImputer

        apa = pd.read_csv(os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap/apa_matrix.csv"), index_col=0)
        coords = pd.read_csv(os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap/coordinates.csv"))
        X = coords[["x", "y"]].values.astype(float)
        nn = NearestNeighbors(n_neighbors=2).fit(X); d, _ = nn.kneighbors(X)
        median_nn = float(np.median(d[:, 1]))

        # Spatial weights: 8 nearest neighbours, row-normalised, sparse (computed once).
        K = 8
        knn = NearestNeighbors(n_neighbors=K + 1).fit(X)
        _, idx = knn.kneighbors(X)
        rows = np.repeat(np.arange(len(X)), K)
        cols = idx[:, 1:].ravel()
        Wraw = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(len(X), len(X)))
        # symmetric + row-normalised.
        Wsym = (Wraw + Wraw.T) * 0.5
        deg = np.asarray(Wsym.sum(axis=1)).ravel()
        deg[deg == 0] = 1.0
        W = Wsym.multiply(1.0 / deg[:, None]).tocsr()

        # Pick genes with enough observations for a stable per-gene RMSE.
        obs_frac = apa.notna().mean(axis=1)
        keep = apa[obs_frac.between(0.05, 0.55)].sample(n=min(150, (obs_frac.between(0.05, 0.55)).sum()), random_state=42)
        V = keep.values.astype(float)
        rng = np.random.RandomState(42)
        mask = np.zeros_like(V, dtype=bool)
        for g in range(V.shape[0]):
            obs = np.where(np.isfinite(V[g]))[0]
            n = max(1, int(round(0.2 * len(obs))))
            mask[g, rng.choice(obs, size=n, replace=False)] = True

        imp = SparseGPImputer(n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1)
        recs = []
        for g in range(V.shape[0]):
            vals = V[g].astype(float)
            tm = np.isfinite(vals) & (~mask[g])
            # Moran's I on observed entries (truth).
            mi = morans_i(vals, W)
            # mean baseline.
            pred_mean = np.full(len(vals), np.nanmean(vals[tm]))
            # GP.
            try:
                pred_gp, _ = imp.impute(X, vals, mask=tm, return_uncertainty=True)
            except Exception:
                pred_gp = pred_mean.copy()
            hg = mask[g] & np.isfinite(vals)
            rmse_gp = np.sqrt(np.nanmean((vals[hg] - pred_gp[hg]) ** 2))
            rmse_mean = np.sqrt(np.nanmean((vals[hg] - pred_mean[hg]) ** 2))
            recs.append({"gene_idx": g, "morans_i": mi, "rmse_gp": rmse_gp, "rmse_mean": rmse_mean,
                         "delta_rmse": rmse_gp - rmse_mean, "n_obs": int(tm.sum())})
            if (g + 1) % 25 == 0:
                print(f"[S4] {g+1}/{V.shape[0]} genes done", flush=True)
        res = pd.DataFrame(recs)
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        res.to_csv(CACHE, index=False)
        print(f"[S4] wrote {CACHE}", flush=True)

    # ---------- Stratify by Moran's I quintile ----------
    res = res.dropna(subset=["morans_i", "delta_rmse"])
    res = res[np.isfinite(res["delta_rmse"])]
    res["quintile"] = pd.qcut(res["morans_i"], 5, labels=[1, 2, 3, 4, 5])
    res["quintile"] = res["quintile"].astype(int)

    # ---------- Plot ----------
    fig = plt.figure(figsize=(11, 5.5))
    gs = fig.add_gridspec(1, 2, wspace=0.32, width_ratios=[1.1, 1.0])

    # Panel A: box of Δ RMSE by Moran's I quintile.
    axA = fig.add_subplot(gs[0, 0])
    data = [res[res["quintile"] == q]["delta_rmse"].values for q in range(1, 6)]
    positions = np.arange(1, 6)
    quint_colors = [GREY, SKYBLU, GREEN, ORANGE, BLUE]
    bp = axA.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                     showfliers=False, medianprops={"color": "black", "lw": 1.2})
    for patch, c in zip(bp["boxes"], quint_colors):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("#555")
    axA.axhline(0, color=RED, lw=1.0, ls="--", label="GP = mean (no advantage)")
    # label n per quintile and median MI.
    mi_med = res.groupby("quintile")["morans_i"].median()
    for q in range(1, 6):
        n = (res["quintile"] == q).sum()
        axA.annotate(f"n={n}\nI={mi_med[q]:+.2f}", xy=(q, axA.get_ylim()[0]),
                     xytext=(0, 4), textcoords="offset points", ha="center", fontsize=6.5, color="#444")
    axA.set_xticks(positions)
    axA.set_xticklabels([f"Q{q}\n({'low spatial' if q==1 else 'high spatial' if q==5 else ''})".strip() for q in range(1, 6)],
                        fontsize=7.5)
    axA.set_xlabel("Moran's I quintile (gene-level spatial autocorrelation)")
    axA.set_ylabel("Δ RMSE  (GP − mean)")
    axA.set_title("Mean is a strong baseline; GP advantage does not grow with spatial signal", loc="left", fontsize=9.5)
    axA.legend(loc="upper left", fontsize=7.5)
    panel_label(axA, "A", x=-0.07, y=1.05)

    # Panel B: scatter Δ RMSE vs Moran's I (raw, with lowess-ish trend).
    axB = fig.add_subplot(gs[0, 1])
    axB.scatter(res["morans_i"], res["delta_rmse"], s=14, alpha=0.55, color=BLUE, edgecolor="white", lw=0.3)
    axB.axhline(0, color=RED, lw=1.0, ls="--")
    # binned mean trend.
    bins = np.linspace(res["morans_i"].min(), res["morans_i"].max(), 9)
    res["mb"] = pd.cut(res["morans_i"], bins, labels=False, include_lowest=True)
    trend = res.groupby("mb").agg(mi=("morans_i", "median"), d=("delta_rmse", "median")).dropna()
    axB.plot(trend["mi"], trend["d"], color=ORANGE, lw=2.2, marker="o", ms=4, label="median Δ RMSE (binned)")
    axB.set_xlabel("Moran's I (per gene)")
    axB.set_ylabel("Δ RMSE  (GP − mean)")
    axB.set_title("Δ RMSE does not fall as Moran's I grows", loc="left", fontsize=9.5)
    axB.legend(loc="upper left", fontsize=7.5)
    panel_label(axB, "B", x=-0.07, y=1.05)

    fig.suptitle("Supplementary Figure S4 — Mean-baseline stratification by gene-level spatial signal",
                 fontsize=11, fontweight="bold", y=1.04)
    # Honest caption: report what the data actually shows.
    frac_gp_beats = (res["delta_rmse"] < 0).mean() * 100
    q_hi = res[res["quintile"] == 5]
    q_lo = res[res["quintile"] == 1]
    fig.text(0.5, -0.02,
             f"n = {len(res)} genes (GSE183456, 20% held-out). Negative Δ = GP beats mean. "
             f"GP beats mean in {frac_gp_beats:.0f}% of genes overall; the per-quintile spread "
             f"(Q5 median Δ={q_hi['delta_rmse'].median():+.3f} vs Q1 {q_lo['delta_rmse'].median():+.3f}) "
             f"is small — for these sparse Visium APA data the per-gene mean is a strong baseline.",
             ha="center", fontsize=7, style="italic", color="#555")
    save_supp(fig, "supp_fig04_mean_stratification.png")
    print(f"[S4] done (n_genes={len(res)}, mean Δ={res['delta_rmse'].mean():+.4f})", flush=True)


if __name__ == "__main__":
    main()

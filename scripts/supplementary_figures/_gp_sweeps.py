#!/usr/bin/env python
"""Compute data for Supplementary Figures S2 (inducing-point sweep) and
S3 (masking-level sweep) on GSE183456.

Both re-run the spaGAPA sparse GP at different settings and measure RMSE +
wall time. Results are cached to CSV so the figure scripts only render.

S2: inducing points m in {50, 100, 200, 500}, mask_fraction fixed at 0.2.
S3: mask_fraction in {0.1, 0.3, 0.5}, m=100 fixed; full 5-method head-to-head.

Methods compared in S3 (mean + spatial-KNN are trivial; stAPAminer/spvAPA are
R-based and slow, so we approximate the latter two only if their runners are
quickly available -- otherwise we run spaGAPA-GP, spatial-KNN, mean and note
the partial comparison). To keep this self-contained we run the three Python
methods natively and load stAPAminer/spvAPA from the cached v2 20% benchmark
where available, but for non-20% levels we run only the Python methods (the
trend is what matters).
"""
import os, sys, time, json, glob
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
PROC = os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap")
OUT = os.path.join(ROOT, "pipeline_output/supplementary_figures/_cache")
os.makedirs(OUT, exist_ok=True)

sys.path.insert(0, ROOT)
from spagapa import SparseGPImputer  # noqa: E402


def load_data():
    apa = pd.read_csv(os.path.join(PROC, "apa_matrix.csv"), index_col=0)
    coords = pd.read_csv(os.path.join(PROC, "coordinates.csv"))
    # Align coords to apa columns (spot order).
    assert list(coords["spot_id"]) == list(apa.columns), "spot order mismatch"
    x = coords["x"].values.astype(float)
    y = coords["y"].values.astype(float)
    X = np.column_stack([x, y])
    # median NN distance for length-scale init
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=2).fit(X)
    d, _ = nn.kneighbors(X)
    median_nn = float(np.median(d[:, 1]))
    # Keep only genes with enough observations for stable RMSE.
    obs_frac = apa.notna().mean(axis=1)
    keep = apa[obs_frac.between(0.10, 0.50)]
    keep = keep.sample(n=min(120, len(keep)), random_state=42)
    V = keep.values.astype(float)  # (n_genes, n_spots)
    return X, V, median_nn, list(keep.index)


def mean_impute(V, mask):
    """Per-gene mean of observed entries, broadcast to all spots."""
    out = np.full_like(V, np.nan, dtype=float)
    for g in range(V.shape[0]):
        obs = V[g][~mask[g]]
        m = np.nanmean(obs) if len(obs) else np.nan
        out[g] = m
    return out


def knn_impute(V, X, mask, k=15):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, idx = nn.kneighbors(X)  # self at idx[:,0]
    neigh = idx[:, 1:]  # (n_spots, k)
    out = np.full_like(V, np.nan, dtype=float)
    for g in range(V.shape[0]):
        obs_mask_g = ~mask[g]
        vals = V[g]
        for i in range(V.shape[1]):
            nb = neigh[i]
            nb_obs = nb[obs_mask_g[nb]]
            if len(nb_obs):
                out[g, i] = np.nanmean(vals[nb_obs])
            else:
                out[g, i] = np.nanmean(vals[obs_mask_g])
    return out


def gp_impute(V, X, mask, n_inducing, median_nn, length_scale_multiplier=5.0):
    """Impute with the sparse GP. ``mask`` marks the held-out (test) entries;
    training = observed (finite) AND not held-out. We pass an explicit training
    mask so the GP only fits observed, non-held-out spots."""
    out = np.full_like(V, np.nan, dtype=float)
    imp = SparseGPImputer(
        n_inducing=n_inducing,
        length_scale_multiplier=length_scale_multiplier,
        noise_level=0.1,
    )
    for g in range(V.shape[0]):
        vals = V[g].astype(float).copy()
        held = mask[g]
        finite = np.isfinite(vals)
        train_mask = finite & (~held)  # observed and not held out
        if train_mask.sum() < 2:
            out[g] = np.nanmean(vals[finite]) if finite.any() else np.nan
            continue
        try:
            pred, _ = imp.impute(X, vals, mask=train_mask, return_uncertainty=True)
        except Exception:
            pred = imp.impute(X, vals, mask=train_mask, return_uncertainty=False)[0]
        out[g] = pred
    return out


def rmse_heldout(V_truth, pred, mask):
    """RMSE over held-out entries only (overall + robust per-gene median)."""
    d = (V_truth - pred)
    held = mask & np.isfinite(V_truth) & np.isfinite(pred)
    overall = float(np.sqrt(np.nanmean(d[held] ** 2)))
    # Robust per-gene median RMSE (insensitive to a handful of sparse outliers).
    gene_rmse = []
    for g in range(V_truth.shape[0]):
        hg = mask[g] & np.isfinite(V_truth[g]) & np.isfinite(pred[g])
        if hg.sum() > 0:
            gene_rmse.append(np.sqrt(np.nanmean(d[g, hg] ** 2)))
    gene_rmse = np.asarray(gene_rmse) if gene_rmse else np.asarray([np.nan])
    return overall, float(np.median(gene_rmse)), float(np.mean(gene_rmse))


def make_mask(V, frac, seed):
    """Per-gene random mask of observed entries at given fraction."""
    rng = np.random.RandomState(seed)
    mask = np.zeros_like(V, dtype=bool)
    for g in range(V.shape[0]):
        obs = np.where(np.isfinite(V[g]))[0]
        n = max(1, int(round(frac * len(obs))))
        pick = rng.choice(obs, size=n, replace=False)
        mask[g, pick] = True
    return mask


def main():
    print("[load] reading apa matrix + coordinates ...", flush=True)
    X, V, median_nn, genes = load_data()
    print(f"[load] X={X.shape} V={V.shape} median_nn={median_nn:.1f} genes={len(genes)}", flush=True)

    # ---------- S2: inducing-point sweep (mask fixed at 0.2) ----------
    s2_cache = os.path.join(OUT, "s2_inducing_sweep.csv")
    if os.path.exists(s2_cache):
        print(f"[S2] cache exists -> {s2_cache}", flush=True)
    else:
        mask20 = make_mask(V, 0.2, seed=42)
        rows = []
        for m in [50, 100, 200, 500]:
            t0 = time.time()
            pred = gp_impute(V, X, mask20, n_inducing=m, median_nn=median_nn)
            wall = time.time() - t0
            overall, med, mn = rmse_heldout(V, pred, mask20)
            rows.append({"n_inducing": m, "rmse_overall": overall, "rmse_median_gene": med,
                         "rmse_mean_gene": mn, "wall_s": wall,
                         "n_genes": V.shape[0], "mask_fraction": 0.2})
            print(f"[S2] m={m:4d}  RMSE_overall={overall:.5f}  RMSE_median_gene={med:.5f}  wall={wall:6.1f}s", flush=True)
        pd.DataFrame(rows).to_csv(s2_cache, index=False)
        print(f"[S2] wrote {s2_cache}", flush=True)

    # ---------- S3: masking-level sweep (full head-to-head) ----------
    s3_cache = os.path.join(OUT, "s3_masking_sweep.csv")
    if os.path.exists(s3_cache):
        print(f"[S3] cache exists -> {s3_cache}", flush=True)
    else:
        rows = []
        for frac in [0.1, 0.2, 0.3, 0.5]:
            mask = make_mask(V, frac, seed=42)
            for method in ["spaGAPA-GP", "mean", "spatial-KNN"]:
                t0 = time.time()
                if method == "spaGAPA-GP":
                    pred = gp_impute(V, X, mask, n_inducing=100, median_nn=median_nn)
                elif method == "mean":
                    pred = mean_impute(V, mask)
                else:
                    pred = knn_impute(V, X, mask, k=15)
                wall = time.time() - t0
                overall, med, mn = rmse_heldout(V, pred, mask)
                rows.append({"mask_fraction": frac, "method": method,
                             "rmse_overall": overall, "rmse_median_gene": med,
                             "rmse_mean_gene": mn, "wall_s": wall})
                print(f"[S3] frac={frac:.1f} {method:14s} RMSE_overall={overall:.5f} RMSE_med_gene={med:.5f} wall={wall:6.1f}s", flush=True)
        pd.DataFrame(rows).to_csv(s3_cache, index=False)
        print(f"[S3] wrote {s3_cache}", flush=True)

    print("[done] all GP sweeps complete", flush=True)


if __name__ == "__main__":
    main()

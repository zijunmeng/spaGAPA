#!/usr/bin/env python3
"""Count-native GP vs Gaussian GP vs mean baseline — REAL data head-to-head.

The decisive experiment for the count-native thesis: the manuscript
discloses that the per-gene MEAN wins entry-wise RMSE over the Gaussian-GP
imputation of usage fractions (0.080 vs 0.122).  If the misspecified
Gaussian likelihood is a driver, the binomial-native GP should close (or
flip) that gap on the same masking protocol.

Protocol mirrors the published head-to-head (20% per-gene mask, seed 42):
on GSE183456 (kidney) + GSE220442 (brain), per gene with >= 2 sites:
  k = counts of the gene's most-distal (max-coordinate) site
  n = total site counts of the gene at that spot
Methods compared on masked entries:
  gauss_frac : SparseGPImputer on k/n fractions (current spaGAPA)
  count_gp   : BinomialUsageGP on (k, n) — empirical-logit heteroscedastic
  mean       : per-gene mean of observed fractions

Usage: OPENBLAS_NUM_THREADS=8 python scripts/count_gp_real_validation.py
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from scipy.spatial import cKDTree

from spagapa.imputation import SparseGPImputer
from spagapa.imputation.count_gp import BinomialUsageGP

DATA_ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/data/processed"
SAMPLES = [
    ("GSE183456", "gsm6047774", "gse183456_gsm6047774_scapatrap"),
    ("GSE220442", "gsm6801751", "gse220442_gsm6801751_scapatrap"),
]
MASK_FRACTION = 0.2
MAX_GENES = 400
OUT_DIR = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/count_gp_validation"

def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def load_gene_counts(dirname):
    """Wide-format loader (vectorized): K (distal counts), N (totals), xy."""
    counts = pd.read_csv(os.path.join(DATA_ROOT, dirname, "apa_site_counts.csv.gz"),
                         index_col=0)
    sites = pd.read_csv(os.path.join(DATA_ROOT, dirname, "apa_sites.csv.gz"))
    coords = pd.read_csv(os.path.join(DATA_ROOT, dirname, "coordinates.csv"),
                         index_col=0)
    xy_df = coords[["x", "y"]] if {"x", "y"} <= set(coords.columns) else coords.iloc[:, :2]
    spot_order = counts.columns.astype(str)
    xy = (xy_df.assign(_k=xy_df.index.astype(str)).set_index("_k")
          .reindex(spot_order).values.astype(float))

    si = sites.set_index("site_id").reindex(counts.index.astype(str))
    gene = si["gene_name"].astype(str).fillna("na").values
    coordv = si["coord"].astype(float).fillna(0.0).values
    C = counts.values.astype(np.float32)

    order = np.lexsort((-coordv, gene))          # gene asc, coord desc
    gene_s, coord_s, C_s = gene[order], coordv[order], C[order]
    boundaries = np.flatnonzero(np.r_[True, gene_s[1:] != gene_s[:-1]])
    gene_names = gene_s[boundaries]
    sizes = np.diff(np.r_[boundaries, len(gene_s)])

    multi = sizes >= 2                            # genes with >= 2 sites
    K = C_s[boundaries[multi]].astype(float)      # first row = max coord (desc)
    N = np.add.reduceat(C_s, boundaries, axis=0).astype(float)[multi]
    names = gene_names[multi].tolist()

    obs = (N > 0).sum(axis=1)
    keep = np.where(obs >= 10)[0]
    keep = keep[np.argsort(-obs[keep])][:MAX_GENES]
    return K[keep], N[keep], xy, [names[i] for i in keep]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = []
    for dataset, sample, dirname in SAMPLES:
        log(f"=== {dataset}/{sample} ===")
        K, N, xy, genes = load_gene_counts(dirname)
        log(f"  {K.shape[0]} genes x {K.shape[1]} spots "
            f"(observed {100.0*(N>0).mean():.1f}%)")
        rng = np.random.default_rng(42)

        frac = np.where(N > 0, K / np.maximum(N, 1), 0.0)
        k_m, n_m, frac_m = K.copy(), N.copy(), frac.copy()
        for g in range(K.shape[0]):
            obs_idx = np.where(N[g] > 0)[0]
            m = rng.choice(obs_idx, size=max(1, int(len(obs_idx) * MASK_FRACTION)),
                           replace=False)
            k_m[g, m], n_m[g, m], frac_m[g, m] = 0, 0, 0.0
        held = (N > 0) & (n_m == 0)
        g_idx, s_idx = np.where(held)
        truth = K[g_idx, s_idx] / np.maximum(N[g_idx, s_idx], 1)

        nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
        ls = float(np.median(nn)) * 5
        n_ind = min(500, max(100, xy.shape[0] // 100))

        # 1) Gaussian GP on fractions (current)
        t0 = time.time()
        b1 = SparseGPImputer(n_inducing=n_ind, length_scale=ls,
                             noise_level=0.1).fit_batch(xy, frac_m, mask=frac_m > 0,
                                                        verbose=False)
        p1_all = b1.impute(return_uncertainty=False)
        p1 = (p1_all[0] if isinstance(p1_all, tuple) else p1_all)[g_idx, s_idx]
        t_gauss = time.time() - t0

        # 2) count-native GP
        t0 = time.time()
        gp2 = BinomialUsageGP(n_inducing=n_ind, length_scale=ls, total_floor=2)
        gp2.fit(xy, k_m, n_m)
        p2 = gp2.predict_usage()[g_idx, s_idx]
        t_count = time.time() - t0

        # 3) per-gene mean
        mean_hat = np.array([frac_m[g][frac_m[g] > 0].mean() for g in g_idx])

        for name, pred, tt in [("gauss_frac", p1, t_gauss), ("count_gp", p2, t_count),
                               ("mean", mean_hat, 0.0)]:
            rmse = float(np.sqrt(np.mean((pred - truth) ** 2)))
            r = float(np.corrcoef(pred, truth)[0, 1])
            rows.append({"dataset": dataset, "sample": sample, "method": name,
                         "n_test": len(truth), "rmse": rmse, "pearson_r": r,
                         "seconds": round(tt, 1)})
            log(f"  {name:10s} RMSE={rmse:.4f} r={r:.3f} ({tt:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "count_vs_gauss_vs_mean.csv"), index=False)
    log(f"done -> {OUT_DIR}/count_vs_gauss_vs_mean.csv")


if __name__ == "__main__":
    main()

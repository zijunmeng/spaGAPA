#!/usr/bin/env python3
"""v2 statistical upgrade — full 11-sample validation sweep.

Extends the published conformal validation (11 samples, 523k test points)
with the two v2 methods, same masking protocol (20% per-gene, seed 42):

  A gauss_global  : Gaussian GP on fractions + global split conformal
                    (reproduces the published pipeline)
  B count_global  : BinomialUsageGP + global split conformal
  C count_mondrian: BinomialUsageGP + Mondrian (5 gene-depth bins)
  mean            : per-gene mean baseline (RMSE reference)

Per sample: RMSE (+mean baseline), marginal coverage at 80/90/95,
per-bin coverage at 90% (B vs C), uncertainty-error Pearson r.

Usage: OPENBLAS_NUM_THREADS=8 python scripts/v2_calibration_all_datasets.py
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from scipy.spatial import cKDTree

from spagapa.imputation import SparseGPImputer
from spagapa.imputation.count_gp import BinomialUsageGP
from spagapa.imputation.calibration import (
    ConformalCalibrator,
    MondrianConformalCalibrator,
    bin_by_quantiles,
)

DATA_ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/data/processed"
SAMPLES = [
    ("GSE237183", "gsm7596587"), ("GSE237183", "gsm7596595"),
    ("GSE237183", "gsm7596604"), ("GSE183456", "gsm6047774"),
    ("GSE179572", "gsm5420751"), ("GSE220442", "gsm6801751"),
    ("GSE220442", "gsm6801753"), ("GSE169749", "gsm5213483"),
    ("GSE263303", "gsm8189356"), ("GSE338525", "gsm9876373"),
    ("GSE338525", "gsm9876374"),
]
MASK_FRACTION = 0.2
CAL_FRACTION = 0.5
ALPHAS = [(0.2, 80), (0.1, 90), (0.05, 95)]
MAX_GENES = 300
OUT_DIR = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/v2_calibration"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_gene_counts(dirname):
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

    order = np.lexsort((-coordv, gene))
    gene_s, C_s = gene[order], C[order]
    boundaries = np.flatnonzero(np.r_[True, gene_s[1:] != gene_s[:-1]])
    sizes = np.diff(np.r_[boundaries, len(gene_s)])
    multi = sizes >= 2
    K = C_s[boundaries[multi]].astype(float)
    N = np.add.reduceat(C_s, boundaries, axis=0).astype(float)[multi]

    obs = (N > 0).sum(axis=1)
    keep = np.where(obs >= 5)[0]
    keep = keep[np.argsort(-obs[keep])][:MAX_GENES]
    return K[keep], N[keep], xy


def run_sample(dataset, sample):
    dirname = f"{dataset.lower()}_{sample}_scapatrap"
    K, N, xy = load_gene_counts(dirname)
    rng = np.random.default_rng(42)

    frac = np.where(N > 0, K / np.maximum(N, 1), 0.0)
    k_m, n_m, frac_m = K.copy(), N.copy(), frac.copy()
    for g in range(K.shape[0]):
        obs_idx = np.where(N[g] > 0)[0]
        if len(obs_idx) < 5:
            continue
        m = rng.choice(obs_idx, size=max(1, int(len(obs_idx) * MASK_FRACTION)),
                       replace=False)
        k_m[g, m], n_m[g, m], frac_m[g, m] = 0, 0, 0.0
    held = (N > 0) & (n_m == 0)
    g_idx, s_idx = np.where(held)
    truth = K[g_idx, s_idx] / np.maximum(N[g_idx, s_idx], 1)

    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    ls = float(np.median(nn)) * 5
    n_ind = min(500, max(100, xy.shape[0] // 100))

    # --- A: gauss_global (published pipeline) ---
    t0 = time.time()
    bA = SparseGPImputer(n_inducing=n_ind, length_scale=ls,
                         noise_level=0.1).fit_batch(xy, frac_m, mask=frac_m > 0,
                                                    verbose=False)
    pA_all, sA_all = bA.impute(return_uncertainty=True)
    pA, sA = pA_all[g_idx, s_idx], sA_all[g_idx, s_idx]
    errA = np.abs(truth - pA)
    tA = time.time() - t0

    # --- B/C: count-native GP ---
    t0 = time.time()
    gp = BinomialUsageGP(n_inducing=n_ind, length_scale=ls, total_floor=2)
    gp.fit(xy, k_m, n_m)
    pB, sB = gp.predict_usage(return_std=True)
    pB, sB = pB[g_idx, s_idx], sB[g_idx, s_idx]
    errB = np.abs(truth - pB)
    tB = time.time() - t0

    # mean baseline
    mean_hat = np.array([frac_m[g][frac_m[g] > 0].mean() for g in g_idx])

    # cal/test split (shared)
    perm = np.random.default_rng(7).permutation(len(truth))
    n_cal = int(len(truth) * CAL_FRACTION)
    cal_i, te_i = perm[:n_cal], perm[n_cal:]

    # deployable conditioning: gene mean observed depth
    gene_depth = np.array([n_m[g][n_m[g] > 0].mean() for g in g_idx])
    bins = bin_by_quantiles(gene_depth, n_bins=5)

    rows = []
    # calibration objects per alpha
    calA = {a: ConformalCalibrator(alpha=a, mode="locally_adaptive")
            for a, _ in ALPHAS}
    calB = {a: ConformalCalibrator(alpha=a, mode="locally_adaptive")
            for a, _ in ALPHAS}
    monB = {a: MondrianConformalCalibrator(alpha=a, mode="locally_adaptive",
                                           min_group_size=200)
            for a, _ in ALPHAS}
    for a, lvl in ALPHAS:
        calA[a].fit(errA[cal_i], sA[cal_i])
        calB[a].fit(errB[cal_i], sB[cal_i])
        monB[a].fit(errB[cal_i], sB[cal_i], bins[cal_i])

    def cov(lo, hi):
        return float(np.mean((truth[te_i] >= lo[te_i]) & (truth[te_i] <= hi[te_i])))

    for a, lvl in ALPHAS:
        loA, hiA = calA[a].predict(pA, sA)
        loB, hiB = calB[a].predict(pB, sB)
        loC, hiC = monB[a].predict(pB, sB, bins)
        rows.append(dict(dataset=dataset, sample=sample, n_test=len(te_i),
                         alpha=a, level=lvl,
                         cov_A_gauss=float(cov(loA, hiA)),
                         cov_B_count=float(cov(loB, hiB)),
                         cov_C_mondrian=float(cov(loC, hiC))))
        if a == 0.1:  # per-bin at 90%
            for b in range(5):
                m = te_i[bins[te_i] == b]
                if m.sum() < 30:
                    continue
                rows.append(dict(
                    dataset=dataset, sample=sample, n_test=int(m.sum()),
                    alpha=a, level=f"90_bin{b+1}",
                    cov_A_gauss=float(np.mean((truth[m] >= loA[m]) & (truth[m] <= hiA[m]))),
                    cov_B_count=float(np.mean((truth[m] >= loB[m]) & (truth[m] <= hiB[m]))),
                    cov_C_mondrian=float(np.mean((truth[m] >= loC[m]) & (truth[m] <= hiC[m]))),
                ))

    rmse = lambda p: float(np.sqrt(np.mean((p[te_i] - truth[te_i]) ** 2)))
    summary = dict(
        dataset=dataset, sample=sample, n_genes=K.shape[0], n_spots=xy.shape[0],
        n_test=len(te_i),
        rmse_A=rmse(pA), rmse_B=rmse(pB), rmse_mean=rmse(mean_hat),
        unc_err_r_A=float(np.corrcoef(sA, errA)[0, 1]),
        unc_err_r_B=float(np.corrcoef(sB, errB)[0, 1]),
        wall_A=round(tA, 1), wall_B=round(tB, 1),
    )
    return rows, summary


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_rows, summaries = [], []
    for dataset, sample in SAMPLES:
        log(f"=== {dataset}/{sample} ===")
        try:
            rows, s = run_sample(dataset, sample)
            all_rows.extend(rows)
            summaries.append(s)
            log(f"  rmse A={s['rmse_A']:.4f} B={s['rmse_B']:.4f} "
                f"mean={s['rmse_mean']:.4f} | unc-err r A={s['unc_err_r_A']:.2f} "
                f"B={s['unc_err_r_B']:.2f}")
        except Exception as e:
            log(f"  FAIL {type(e).__name__}: {e}")
    pd.DataFrame(all_rows).to_csv(os.path.join(OUT_DIR, "coverage_3config.csv"),
                                  index=False)
    pd.DataFrame(summaries).to_csv(os.path.join(OUT_DIR, "sample_summary.csv"),
                                   index=False)
    cov = pd.DataFrame(all_rows)
    marg = cov[cov.level.isin([80, 90, 95])]
    print("\n=== marginal coverage (mean over 11 samples) ===")
    print(marg.groupby("level")[["cov_A_gauss", "cov_B_count",
                                 "cov_C_mondrian"]].mean().round(4))
    bins90 = cov[cov.level.astype(str).str.startswith("90_bin")]
    print("\n=== 90% per-depth-bin coverage (mean over samples) ===")
    print(bins90.groupby("level")[["cov_A_gauss", "cov_B_count",
                                   "cov_C_mondrian"]].mean().round(4))
    log("done")


if __name__ == "__main__":
    main()

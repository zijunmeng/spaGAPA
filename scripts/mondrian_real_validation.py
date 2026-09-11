#!/usr/bin/env python3
"""Mondrian vs global conformal on REAL data — the subgroup-coverage fix.

Reproduces the published conditional-coverage failure mode (high-expression
bins undercut the 90% marginal target) with the SAME pipeline as
scripts/conditional_coverage_validation.py (20% per-gene mask, SparseGPImputer
config, 50/50 cal/test split), then shows MondrianConformalCalibrator
(group-conditional split conformal) restores per-bin coverage.

Conditioning variable: per-gene mean OBSERVED usage (deployable — known at
prediction time from the gene's observed spots), 5 quantile bins.

Usage:
  OPENBLAS_NUM_THREADS=8 python scripts/mondrian_real_validation.py
"""
from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from scipy.spatial import cKDTree

from spagapa.imputation import SparseGPImputer
from spagapa.imputation.calibration import (
    ConformalCalibrator,
    MondrianConformalCalibrator,
    bin_by_quantiles,
)

DATA_ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/data/processed"
SAMPLES = [
    ("GSE183456", "gsm6047774", "gse183456_gsm6047774_scapatrap"),
    ("GSE220442", "gsm6801751", "gse220442_gsm6801751_scapatrap"),
]
MASK_FRACTION = 0.2
CAL_FRACTION = 0.5
ALPHA = 0.1
MAX_GENES = 400
MIN_SPOTS = 500
MIN_GENES_WITH_OBS = 50
OUT_DIR = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/mondrian_validation"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_sample(dir_basename: str):
    apa = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "apa_matrix.csv"),
                      index_col=0)
    coords = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "coordinates.csv"),
                         index_col=0)
    xy_df = coords[["x", "y"]] if {"x", "y"} <= set(coords.columns) else coords.iloc[:, :2]
    xy_df.columns = ["x", "y"]
    common = apa.columns.intersection(xy_df.index)
    apa, xy_df = apa.loc[:, common], xy_df.loc[common]
    return apa.values.astype(float), xy_df.values.astype(float)


def fit_sample(dir_basename: str, rng):
    """Mirror conditional_coverage_validation.fit_sample; return cal/test arrays
    plus per-point gene-level mean observed usage (deployable conditioning)."""
    values, xy = load_sample(dir_basename)
    n_spots = values.shape[1]
    obs_count = (values > 0).sum(axis=1)
    keep_rows = np.where(obs_count >= 5)[0]
    keep = keep_rows[np.argsort(-obs_count[keep_rows])][:MAX_GENES]
    sub = values[keep]

    masked = sub.copy()
    for g in range(sub.shape[0]):
        obs_idx = np.where(sub[g] > 0)[0]
        if len(obs_idx) < 5:
            continue
        m = rng.choice(obs_idx, size=max(1, int(len(obs_idx) * MASK_FRACTION)),
                       replace=False)
        masked[g, m] = 0.0

    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    base = SparseGPImputer(n_inducing=min(500, max(100, n_spots // 100)),
                           length_scale=float(np.median(nn)) * 5,
                           noise_level=0.1)
    tg = time.time()
    batch = base.fit_batch(xy, masked, mask=masked > 0, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    log(f"  GP fit+impute {time.time()-tg:.1f}s")

    held = (sub > 0) & (masked <= 0)
    g_idx, s_idx = np.where(held)
    gene_expr = np.array([sub[g][sub[g] > 0].mean() for g in range(sub.shape[0])])
    truths, preds, unc = values_keep(sub, g_idx, s_idx), gp_pred[g_idx, s_idx], gp_unc[g_idx, s_idx]
    errors = np.abs(truths - preds)
    cond = gene_expr[g_idx]  # per-point conditioning variable

    perm = rng.permutation(len(truths))
    n_cal = int(len(truths) * CAL_FRACTION)
    return dict(cal_i=perm[:n_cal], test_i=perm[n_cal:], truths=truths,
                preds=preds, unc=unc, errors=errors, cond=cond)


def values_keep(sub, g_idx, s_idx):
    return sub[g_idx, s_idx]


def per_bin_coverage(err, half, bins, labels):
    rows = []
    for b in sorted(np.unique(bins[~np.isnan(bins)])):
        m = bins == b
        rows.append({"bin": labels[int(b)], "n": int(m.sum()),
                     "coverage": float(np.mean(err[m] <= half[m])),
                     "mean_half_width": float(np.mean(half[m]))})
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_rows = []
    for dataset, sample, dirname in SAMPLES:
        log(f"=== {dataset}/{sample} ===")
        d = fit_sample(dirname, np.random.default_rng(42))
        cal_i, te_i = d["cal_i"], d["test_i"]
        bins = bin_by_quantiles(d["cond"], n_bins=5)
        labels = ["q1_low", "q2", "q3", "q4", "q5_high"]

        glob = ConformalCalibrator(alpha=ALPHA, mode="locally_adaptive")
        glob.fit(d["errors"][cal_i], d["unc"][cal_i])
        lo, hi = glob.predict(d["preds"], d["unc"])
        half_g = (hi - lo) / 2

        mon = MondrianConformalCalibrator(alpha=ALPHA, mode="locally_adaptive",
                                          min_group_size=200)
        mon.fit(d["errors"][cal_i], d["unc"][cal_i], bins[cal_i])
        lo2, hi2 = mon.predict(d["preds"], d["unc"], bins)
        half_m = (hi2 - lo2) / 2

        for r in per_bin_coverage(d["errors"][te_i], half_g[te_i], bins[te_i], labels):
            all_rows.append({"dataset": dataset, "sample": sample,
                             "method": "global", **r})
        for r in per_bin_coverage(d["errors"][te_i], half_m[te_i], bins[te_i], labels):
            all_rows.append({"dataset": dataset, "sample": sample,
                             "method": "mondrian", **r})
        all_rows.append({"dataset": dataset, "sample": sample, "method": "global",
                         "bin": "OVERALL", "n": len(te_i),
                         "coverage": float(np.mean(d["errors"][te_i] <= half_g[te_i])),
                         "mean_half_width": float(np.mean(half_g[te_i]))})
        all_rows.append({"dataset": dataset, "sample": sample, "method": "mondrian",
                         "bin": "OVERALL", "n": len(te_i),
                         "coverage": float(np.mean(d["errors"][te_i] <= half_m[te_i])),
                         "mean_half_width": float(np.mean(half_m[te_i]))})
    df = pd.DataFrame(all_rows)
    df.to_csv(os.path.join(OUT_DIR, "mondrian_vs_global.csv"), index=False)
    print(df.pivot_table(index=["dataset", "bin"], columns="method",
                         values="coverage").round(4).to_string())
    log(f"done -> {OUT_DIR}/mondrian_vs_global.csv")


if __name__ == "__main__":
    main()

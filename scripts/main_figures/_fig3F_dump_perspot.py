#!/usr/bin/env python
"""One-shot: generate REAL per-spot spaGAPA posterior + split-conformal
intervals for ONE gene in GSE183456, dumped as a CSV for Figure 3 Panel F.

This replaces the previous Panel F, which painted an RBF-smoothing / local-SD
proxy and mislabeled it as "spaGAPA posterior mean / interval width".  Here we
actually run ``SparseGPImputer`` (the model used everywhere else in the paper)
on a single gene with adequate observations and clear spatial structure, then
run a split-conformal calibration (hold out 20%, split 50/50 into calibrate /
test, take the |y - yhat| quantile on the calibration half) to obtain an honest
per-spot interval [yhat - qhat, yhat + qhat].

Output
------
pipeline_output/main_figures/_data/fig3F_<gene>_perspot.csv
    columns: spot, x, y, observed, posterior_mean, sigma, qhat,
             abs_error, covered  (covered uses the 90% qhat)

Gene: peak_20919  (Moran's I = +0.081, obs-frac = 0.218; chosen as the gene
with the strongest spatial autocorrelation among the 0.15-0.30 obs-frac band,
so the spatial maps are informative rather than near-flat).
"""
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402

DATA_DIR = os.path.join(REPO, "data/processed/gse183456_gsm6047774_scapatrap")
OUT_DIR = os.path.join(REPO, "pipeline_output/main_figures/_data")
os.makedirs(OUT_DIR, exist_ok=True)

GENE = "peak_20919"
MASK_FRAC = 0.20          # hold out 20% of observed entries
CAL_FRAC = 0.50           # split held-out 50/50 calibrate / test
ALPHA = 0.10              # 90% nominal
SEED = 42


def main():
    apa = pd.read_csv(os.path.join(DATA_DIR, "apa_matrix.csv"), index_col=0)
    coord = pd.read_csv(os.path.join(DATA_DIR, "coordinates.csv"), index_col=0)
    # align spots
    apa = apa[coord.index]
    X = coord[["x", "y"]].values.astype(float)
    n_spots = X.shape[0]

    assert GENE in apa.index, f"{GENE} not in apa matrix"
    y = apa.loc[GENE].values.astype(float)   # observed APA usage per spot (0 = unobserved)
    obs_mask = y > 0
    n_obs = int(obs_mask.sum())
    print(f"[F] gene={GENE}  spots={n_spots}  observed={n_obs} ({n_obs/n_spots:.3f})")

    # ---- hold out 20% of OBSERVED entries (truth we will recover) ----
    rng = np.random.default_rng(SEED)
    obs_idx = np.where(obs_mask)[0]
    n_hold = max(1, int(round(MASK_FRAC * len(obs_idx))))
    hold_i = rng.choice(obs_idx, size=n_hold, replace=False)
    hold_mask = np.zeros(n_spots, dtype=bool)
    hold_mask[hold_i] = True
    # train mask = observed AND not held-out  (this is what the GP sees)
    train_mask = obs_mask & (~hold_mask)

    y_train = y.copy()
    y_train[~train_mask] = 0.0   # SparseGPImputer treats 0 as "missing"

    # ---- GP fit (mirrors calibrate_uncertainty_all_datasets.py) ----
    kdt = cKDTree(X)
    nn = kdt.query(X, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    length_scale = nn_dist * 5
    n_inducing = min(500, max(100, n_spots // 100))
    print(f"[F] GP: n_inducing={n_inducing}  length_scale={length_scale:.1f}  nn_dist={nn_dist:.2f}")
    imp = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale,
                          noise_level=0.1)
    mu, sigma = imp.impute(X, y_train, mask=train_mask, return_uncertainty=True)
    mu = np.asarray(mu, dtype=float).ravel()
    sigma = np.asarray(sigma, dtype=float).ravel()
    print(f"[F] posterior mean range [{mu.min():.3f}, {mu.max():.3f}], "
          f"sigma range [{sigma.min():.4f}, {sigma.max():.4f}]")

    # ---- split-conformal on the held-out observed entries ----
    truths = y[hold_mask]
    preds = mu[hold_mask]
    errs = np.abs(truths - preds)
    perm = rng.permutation(len(truths))
    n_cal = int(len(truths) * CAL_FRAC)
    cal_i = perm[:n_cal]
    test_i = perm[n_cal:]
    cal_err = errs[cal_i]

    # split-conformal qhat at 90%:  ceil((n_cal+1)*(1-alpha))/n_cal quantile
    n_cal_n = len(cal_err)
    q_level = min(1.0, np.ceil((n_cal_n + 1) * (1 - ALPHA)) / n_cal_n)
    qhat = float(np.quantile(cal_err, q_level))
    test_cov = float(np.mean(np.abs(truths[test_i] - preds[test_i]) <= qhat))
    print(f"[F] n_cal={n_cal_n}  n_test={len(test_i)}  qhat90={qhat:.4f}  "
          f"test coverage@90={test_cov:.4f}")

    # ---- per-spot table ----
    covered = np.abs(y - mu) <= qhat          # only meaningful where observed
    covered_obs = np.where(obs_mask, covered, False)
    abs_err = np.where(obs_mask, np.abs(y - mu), np.nan)

    df = pd.DataFrame({
        "spot": coord.index.values,
        "x": X[:, 0],
        "y": X[:, 1],
        "observed": obs_mask.astype(int),
        "observed_value": np.where(obs_mask, y, np.nan),
        "posterior_mean": mu,
        "sigma": sigma,
        "qhat": qhat,                       # same scalar for every spot (90% split-conformal)
        "interval_width_2qhat": 2 * qhat,
        "abs_error": abs_err,
        "covered_90": covered_obs.astype(int),
    })
    out = os.path.join(OUT_DIR, f"fig3F_{GENE}_perspot.csv")
    df.to_csv(out, index=False)
    print(f"[F] wrote {out}  ({len(df)} rows)")
    print(f"[F] observed spots: {int(obs_mask.sum())}  "
          f"covered@90 among observed: {int(covered_obs.sum())} "
          f"({covered_obs.sum()/obs_mask.sum():.4f})")


if __name__ == "__main__":
    main()

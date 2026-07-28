#!/usr/bin/env python3
"""Expert-reviewer must-add results #1 and #2 (shared GP+Method-D infrastructure).

RESULT 1 (Fig 4G): Risk-coverage curve
    On 5 representative datasets, mask 20%, fit GP with Method D (residual_spot),
    conformal-calibrate on a 50% held-out split, then on the test set rank spots
    by calibrated uncertainty and retain the top X% most confident predictions
    (X = 20/40/60/80/100%). Compare uncertainty-ranked vs random (avg of 10
    shuffles) vs oracle (rank by true error). Report RMSE + MAE at each threshold.

RESULT 2 (Suppl Fig S9): Spatial-block conformal split
    On the same 5 datasets, K-means-cluster coordinates into K=5 spatial blocks.
    Assign 3 blocks -> calibration, 2 blocks -> test (block-level split). Fit GP
    on the 80% training spots, conformal-calibrate on the calibration-block
    spots, evaluate coverage on the test-block spots. Compare against the random
    spot split.

Both results are produced in one run because they share the expensive GP fit
(Method D residual_spot, 2-pass).

Protocol hyperparameters mirror calibrate_uncertainty.py / test_uncertainty_noise_methods.py.
Method D = residual-based per-spot noise (2-pass empirical Bayes), the
recommended final noise model.

Outputs
-------
pipeline_output/risk_coverage_curve/
    risk_coverage_data.csv   dataset, retention_fraction, method, rmse, mae, n_retained
    figure.png               multi-panel risk-coverage curve (publication quality)
    summary.json             mean RMSE improvement @80%/@50%, paired t-test p-values
pipeline_output/spatial_block_conformal/
    block_split_coverage.csv dataset, split_type, coverage_90, coverage_80, coverage_95, n_test
    figure.png               side-by-side bar chart random vs block coverage
    summary.json             mean coverage, max deviation, exchangeability verdict
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, ttest_rel
from sklearn.cluster import KMeans

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402

DATA_ROOT = os.path.join(REPO, "data", "processed")
RC_OUT = os.path.join(REPO, "pipeline_output", "risk_coverage_curve")
SB_OUT = os.path.join(REPO, "pipeline_output", "spatial_block_conformal")

# 5 representative datasets (kidney, brain, mouse colon, liver, glioma).
SAMPLES = [
    ("GSE183456",  "gse183456_gsm6047774_scapatrap"),
    ("GSE220442",  "gse220442_gsm6801751_scapatrap"),
    ("GSE169749",  "gse169749_gsm5213483_scapatrap"),
    ("GSE338525",  "gse338525_gsm9876373_scapatrap"),
    ("GSE237183",  "gse237183_gsm7596587_scapatrap"),
]

# ---- protocol hyperparameters ----
MASK_FRACTION = 0.20
CAL_FRACTION = 0.50           # of held-out -> calibration (rest = test)
SEED = 42
MAX_GENES = 300               # top-K most-observed sites per sample
MIN_SPOTS = 500
MIN_GENES_WITH_OBS = 50
RESID_K = 8                   # kNN neighbourhood for residual-variance binning
NOISE_SCALE = 1000.0          # kernel-amplitude scaling (mirrors local_noise default)
N_BLOCKS = 5                  # K-means spatial blocks for block-conformal
N_CAL_BLOCKS = 3              # blocks -> calibration (rest -> test)
N_RANDOM_SHUFFLES = 10        # for random retention baseline in risk-coverage
RETENTION_FRACS = [0.20, 0.40, 0.60, 0.80, 1.00]
ALPHA_LEVELS = [(0.20, 80), (0.10, 90), (0.05, 95)]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ===========================================================================
# Data loading + masking (mirrors test_uncertainty_noise_methods.prep_sample)
# ===========================================================================
def load_sample(dir_basename: str):
    apa = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "apa_matrix.csv"),
                      index_col=0)
    coords = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "coordinates.csv"),
                         index_col=0)
    if {"x", "y"} <= set(coords.columns):
        xy_df = coords[["x", "y"]]
    else:
        xy_df = coords.iloc[:, :2].astype(float)
        xy_df.columns = ["x", "y"]
    common = apa.columns.intersection(xy_df.index)
    if len(common) > 0:
        apa = apa.loc[:, common]
        xy_df = xy_df.loc[common]
    elif apa.shape[1] != xy_df.shape[0]:
        raise ValueError(f"spot mismatch: apa={apa.shape}, coords={xy_df.shape}")
    xy = xy_df.values.astype(float)
    values = apa.values.astype(float)
    return values, xy


def prep_sample(dir_basename: str, rng: np.random.Generator):
    values, xy = load_sample(dir_basename)
    n_genes_full, n_spots = values.shape
    obs_count = (values > 0).sum(axis=1)
    keep_rows = np.where(obs_count >= 5)[0]
    if len(keep_rows) == 0:
        raise RuntimeError("no genes with >=5 observations")
    order = keep_rows[np.argsort(-obs_count[keep_rows])]
    keep = order[:MAX_GENES]
    values_sub = values[keep]
    n_genes = values_sub.shape[0]
    if n_spots < MIN_SPOTS:
        raise RuntimeError(f"too few spots ({n_spots} < {MIN_SPOTS})")
    if n_genes < MIN_GENES_WITH_OBS:
        raise RuntimeError(f"too few genes ({n_genes} < {MIN_GENES_WITH_OBS})")

    # mask 20% per row
    masked = values_sub.copy()
    for g in range(n_genes):
        obs_idx = np.where(values_sub[g] > 0)[0]
        if len(obs_idx) < 5:
            continue
        n_mask = max(1, int(len(obs_idx) * MASK_FRACTION))
        m = rng.choice(obs_idx, size=n_mask, replace=False)
        masked[g, m] = 0.0

    kdt = cKDTree(xy)
    nn = kdt.query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    length_scale = nn_dist * 5
    n_inducing = min(500, max(100, n_spots // 100))
    return {
        "values": values_sub, "masked": masked, "xy": xy,
        "n_genes": n_genes, "n_spots": n_spots, "nn_dist": nn_dist,
        "length_scale": length_scale, "n_inducing": n_inducing, "tree": kdt,
    }


# ===========================================================================
# Method D: residual-based per-spot noise (2-pass empirical Bayes)
# (verbatim logic from test_uncertainty_noise_methods.estimate_residual_noise)
# ===========================================================================
def estimate_residual_noise(values, masked, xy, residual_pred, tree,
                            k=RESID_K, scale=NOISE_SCALE):
    n_genes, n_spots = masked.shape
    k_use = min(k, n_spots - 1)
    _, idx = tree.query(tree.data, k=k_use + 1)
    nbr_idx = idx[:, 1:]
    spot_noise = np.empty((n_genes, n_spots), dtype=float)
    obs_mask = masked > 0
    for g in range(n_genes):
        v = masked[g]
        obs = obs_mask[g]
        resid = np.where(obs, v - residual_pred[g], np.nan)
        nbr_resid = resid[nbr_idx]
        with np.errstate(all="ignore"):
            with np.testing.suppress_warnings() as sup:
                sup.filter(RuntimeWarning, "Degrees of freedom <= 0 for slice.")
                var = np.nanvar(nbr_resid, axis=1)
        cnt = np.sum(np.isfinite(nbr_resid), axis=1)
        glob = np.nanvar(resid)
        if not np.isfinite(glob) or glob <= 0:
            glob = float(np.var(v[obs])) * 0.5 if obs.any() else 0.01
        bad = (cnt < 2) | ~np.isfinite(var) | (var <= 0)
        if np.any(bad):
            var[bad] = max(glob, 1e-6)
        var = np.clip(var, 1e-6, None)
        spot_noise[g] = var * scale
    return spot_noise


def fit_method_d(p: dict):
    """Return (gp_pred[n_genes,n_spots], gp_unc[n_genes,n_spots]) via Method D."""
    # pass 1: constant noise -> residuals on training spots
    base1 = SparseGPImputer(n_inducing=p["n_inducing"],
                            length_scale=p["length_scale"], noise_level=0.1)
    batch1 = base1.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                             verbose=False)
    pred1, _ = batch1.impute(return_uncertainty=False)
    spot_noise = estimate_residual_noise(p["values"], p["masked"], p["xy"],
                                         pred1, p["tree"])
    # pass 2: re-fit with per-spot noise
    base2 = SparseGPImputer(n_inducing=p["n_inducing"],
                            length_scale=p["length_scale"])
    batch2 = base2.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                             verbose=False, spot_noise=spot_noise)
    return batch2.impute(return_uncertainty=True)


def gather_held(values, masked, gp_pred, gp_unc):
    """(gene_idx, spot_idx, truth, pred, unc, abs_err) for held-out entries."""
    held = (values > 0) & (masked <= 0)
    g_idx, s_idx = np.where(held)
    truth = values[g_idx, s_idx]
    pred = gp_pred[g_idx, s_idx]
    unc = gp_unc[g_idx, s_idx]
    err = np.abs(truth - pred)
    return g_idx, s_idx, truth, pred, unc, err


# ===========================================================================
# RESULT 1: Risk-coverage curve
# ===========================================================================
def conformal_qhat(cal_err, alpha):
    """Split-conformal quantile of absolute errors at level 1-alpha."""
    n = len(cal_err)
    if n < 2:
        return float(np.quantile(cal_err, 1 - alpha)) if n else np.inf
    q = np.quantile(cal_err, np.ceil((1 - alpha) * (n + 1)) / n,
                    method="higher")
    return float(q)


def _rmse(err_subset):
    err_subset = np.asarray(err_subset, dtype=float)
    return float(np.sqrt(np.mean(err_subset ** 2))) if err_subset.size else float("nan")


def _mae(err_subset):
    err_subset = np.asarray(err_subset, dtype=float)
    return float(np.mean(err_subset)) if err_subset.size else float("nan")


def risk_coverage_for_dataset(test_err, test_score, rng, dataset):
    """Compute RMSE/MAE at each retention fraction for 3 ranking methods.

    test_err   : abs error of every test point
    test_score : uncertainty score (lower = more confident) for every test point.
                 For conformal-calibrated uncertainty we use the locally-adaptive
                 interval half-width q_hat(alpha=0.1) * std; if that is degenerate
                 we fall back to raw GP std.
    """
    n = len(test_err)
    order_unc = np.argsort(test_score)          # ascending uncertainty
    order_oracle = np.argsort(test_err)         # ascending true error
    rows = []
    for frac in RETENTION_FRACS:
        n_keep = max(1, int(round(n * frac)))
        # uncertainty-ranked (most-confident first)
        k_unc = order_unc[:n_keep]
        rows.append(dict(dataset=dataset, retention_fraction=frac,
                         method="uncertainty",
                         rmse=_rmse(test_err[k_unc]), mae=_mae(test_err[k_unc]),
                         n_retained=int(n_keep)))
        # oracle (upper bound)
        k_or = order_oracle[:n_keep]
        rows.append(dict(dataset=dataset, retention_fraction=frac,
                         method="oracle",
                         rmse=_rmse(test_err[k_or]), mae=_mae(test_err[k_or]),
                         n_retained=int(n_keep)))
        # random (avg of 10 shuffles)
        rmses, maes = [], []
        for _ in range(N_RANDOM_SHUFFLES):
            perm = rng.permutation(n)[:n_keep]
            rmses.append(_rmse(test_err[perm]))
            maes.append(_mae(test_err[perm]))
        rows.append(dict(dataset=dataset, retention_fraction=frac,
                         method="random",
                         rmse=float(np.mean(rmses)), mae=float(np.mean(maes)),
                         n_retained=int(n_keep)))
    return rows


# ===========================================================================
# RESULT 2: Spatial-block conformal coverage
# ===========================================================================
def spatial_block_labels(xy, seed):
    """K-means cluster coordinates into N_BLOCKS spatial blocks."""
    km = KMeans(n_clusters=N_BLOCKS, random_state=seed, n_init=10)
    return km.fit_predict(xy.astype(float))


def block_conformal_coverage(p, gp_pred, gp_unc, block_labels, rng):
    """Block-level conformal split: 3 blocks calibrate, 2 blocks test.

    Coverage computed on held-out test-block spots. GP already fit on the
    80% observed training spots (masked matrix). Returns dict of coverage
    at 80/90/95% + n_test.
    """
    g_idx, s_idx, truth, pred, unc, err = gather_held(p["values"], p["masked"],
                                                      gp_pred, gp_unc)
    blocks_at_held = block_labels[s_idx]

    # pick N_CAL_BLOCKS distinct blocks for calibration, rest for test
    uniq_blocks = np.unique(block_labels)
    if len(uniq_blocks) < 2:
        return None
    cal_blocks = rng.choice(uniq_blocks, size=min(N_CAL_BLOCKS, len(uniq_blocks) - 1),
                            replace=False)
    cal_mask = np.isin(blocks_at_held, cal_blocks)
    test_mask = ~cal_mask
    if cal_mask.sum() < 10 or test_mask.sum() < 10:
        # fall back: flip assignment
        cal_mask, test_mask = test_mask, cal_mask
    cal_err = err[cal_mask]
    test_err = err[test_mask]
    test_unc = unc[test_mask]
    test_truth = truth[test_mask]
    test_pred = pred[test_mask]

    out = {"n_test": int(test_mask.sum()),
           "n_calibration": int(cal_mask.sum())}
    for alpha, level in ALPHA_LEVELS:
        q_hat = conformal_qhat(cal_err, alpha)
        # use locally-adaptive interval half-width q_hat * std (calibrated
        # heteroscedastic uncertainty). Coverage target = 1-alpha.
        half = q_hat * test_unc
        lo, hi = test_pred - half, test_pred + half
        cov = float(np.mean((test_truth >= lo) & (test_truth <= hi)))
        out[f"qhat_{level}"] = q_hat
        out[f"coverage_{level}"] = cov
    return out


def random_conformal_coverage(p, gp_pred, gp_unc, rng):
    """Random spot split: 50% calibrate, 50% test (the protocol's reference)."""
    g_idx, s_idx, truth, pred, unc, err = gather_held(p["values"], p["masked"],
                                                      gp_pred, gp_unc)
    perm = rng.permutation(len(err))
    n_cal = int(len(err) * CAL_FRACTION)
    cal_idx, test_idx = perm[:n_cal], perm[n_cal:]
    cal_err = err[cal_idx]
    test_err = err[test_idx]
    test_unc = unc[test_idx]
    test_truth = truth[test_idx]
    test_pred = pred[test_idx]
    out = {"n_test": int(len(test_idx)),
           "n_calibration": int(n_cal)}
    for alpha, level in ALPHA_LEVELS:
        q_hat = conformal_qhat(cal_err, alpha)
        half = q_hat * test_unc
        lo, hi = test_pred - half, test_pred + half
        cov = float(np.mean((test_truth >= lo) & (test_truth <= hi)))
        out[f"qhat_{level}"] = q_hat
        out[f"coverage_{level}"] = cov
    return out


# ===========================================================================
# Plotting (publication quality)
# ===========================================================================
def plot_risk_coverage(df_rc, summary, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    datasets = list(dict.fromkeys(df_rc["dataset"]))
    nd = len(datasets)
    ncol = min(3, nd)
    nrow = int(np.ceil(nd / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow),
                             sharex=True, squeeze=False)
    colors = {"uncertainty": "#0072B2", "random": "#999999", "oracle": "#D55E00"}
    styles = {"uncertainty": "-", "random": "--", "oracle": ":"}
    labels = {"uncertainty": "Uncertainty-ranked",
              "random": "Random (mean of 10)",
              "oracle": "Oracle (upper bound)"}
    for ax, ds in zip(axes.ravel(), datasets):
        sub = df_rc[df_rc["dataset"] == ds]
        for m in ("uncertainty", "random", "oracle"):
            s = sub[sub["method"] == m].sort_values("retention_fraction")
            ax.plot(s["retention_fraction"], s["rmse"], styles[m],
                    color=colors[m], lw=2.0, label=labels[m])
        ax.set_title(ds, fontsize=10)
        ax.set_xlim(0.15, 1.05)
        ax.grid(alpha=0.25, lw=0.5)
        ax.tick_params(labelsize=8)
    for ax in axes.ravel()[nd:]:
        ax.axis("off")
    # shared labels + legend
    handles, lbls = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, lbls, loc="lower center", ncol=3, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.supxlabel("Fraction of predictions retained (most-confident first)",
                  fontsize=11)
    fig.supylabel("RMSE of retained predictions", fontsize=11)
    imp80 = summary["mean_rmse_improvement_pct_at_80"]
    imp50 = summary["mean_rmse_improvement_pct_at_50"]
    p80 = summary["paired_ttest_pvalue_at_80"]
    fig.suptitle(
        f"Risk-coverage: uncertainty ranking is actionable\n"
        f"mean RMSE improvement vs random: +{imp80:.1f}% @80% retained "
        f"(p={p80:.2g}), +{imp50:.1f}% @50%", fontsize=11)
    fig.tight_layout(rect=(0, 0.03, 1, 0.93))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_block_coverage(df_cov, summary, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    datasets = list(dict.fromkeys(df_cov["dataset"]))
    x = np.arange(len(datasets))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8, 4.5))
    cov_rand = np.array([summary["_per_dataset"][d]["random"]["coverage_90"]
                         for d in datasets])
    cov_block = np.array([summary["_per_dataset"][d]["block"]["coverage_90"]
                          for d in datasets])
    b1 = ax.bar(x - w / 2, cov_rand, w, label="Random spot split",
                color="#0072B2", edgecolor="black", lw=0.4)
    b2 = ax.bar(x + w / 2, cov_block, w, label="Spatial-block split",
                color="#D55E00", edgecolor="black", lw=0.4)
    ax.axhline(0.90, color="black", ls="--", lw=1.0, alpha=0.7,
               label="Nominal 90% target")
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Empirical coverage at 90%", fontsize=11)
    ax.set_xlabel("Dataset", fontsize=11)
    ax.set_ylim(0, 1.05)
    md = summary["max_abs_deviation_block_from_nominal"]
    verdict = summary["exchangeability_concern_supported"]
    ax.set_title(
        f"Random vs spatial-block conformal coverage (90% target)\n"
        f"mean random={summary['mean_coverage_random_90']:.3f}  "
        f"mean block={summary['mean_coverage_block_90']:.3f}  "
        f"max |block-0.90|={md:.3f}  "
        f"exchangeability concern={'SUPPORTED' if verdict else 'NOT supported (coverage holds)'}",
        fontsize=10)
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    ax.grid(axis="y", alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ===========================================================================
# Driver
# ===========================================================================
def main():
    os.makedirs(RC_OUT, exist_ok=True)
    os.makedirs(SB_OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rc_rows = []
    sb_rows = []
    per_ds_block = {}

    for dataset, dir_base in SAMPLES:
        log(f"=== {dataset} ({dir_base}) ===")
        t0 = time.time()
        try:
            p = prep_sample(dir_base, rng)
        except Exception as e:
            log(f"  SKIP prep failed: {e}")
            continue
        log(f"  loaded {p['n_genes']} genes x {p['n_spots']} spots, "
            f"nn_dist={p['nn_dist']:.1f}, ls={p['length_scale']:.1f}, "
            f"n_inducing={p['n_inducing']}")
        t1 = time.time()
        gp_pred, gp_unc = fit_method_d(p)
        log(f"  Method D fit done ({time.time()-t1:.1f}s)")

        # ---- held-out tuples ----
        g_idx, s_idx, truth, pred, unc, err = gather_held(
            p["values"], p["masked"], gp_pred, gp_unc)
        log(f"  held-out points: {len(err)}")

        # =================== RESULT 1: risk-coverage ====================
        # Split calibration / test (50/50), calibrate conformal q_hat at alpha=0.1
        # to obtain a calibrated per-point uncertainty score = q_hat * std.
        perm = rng.permutation(len(err))
        n_cal = int(len(err) * CAL_FRACTION)
        cal_idx_, test_idx_ = perm[:n_cal], perm[n_cal:]
        cal_err = err[cal_idx_]
        q_hat_90 = conformal_qhat(cal_err, alpha=0.10)
        test_err_rc = err[test_idx_]
        test_unc_rc = unc[test_idx_]
        test_score_rc = q_hat_90 * test_unc_rc   # calibrated half-width (lower = more confident)
        log(f"  [RC] n_test={len(test_err_rc)}  "
            f"qhat_90={q_hat_90:.4f}  "
            f"RMSE_all={_rmse(test_err_rc):.4f}")
        rc_rows.extend(risk_coverage_for_dataset(
            test_err_rc, test_score_rc, rng, dataset))

        # =================== RESULT 2: spatial-block conformal ===========
        block_labels = spatial_block_labels(p["xy"], SEED)
        # random split (reference) — already have held-out, just need cal/test
        rc_cov = random_conformal_coverage(p, gp_pred, gp_unc, rng)
        blk_cov = block_conformal_coverage(p, gp_pred, gp_unc, block_labels, rng)
        log(f"  [SB] random cov90={rc_cov['coverage_90']:.4f}  "
            f"block cov90={blk_cov['coverage_90']:.4f}  "
            f"(n_test_random={rc_cov['n_test']}, n_test_block={blk_cov['n_test']})")
        for split_type, cov in (("random", rc_cov), ("block", blk_cov)):
            sb_rows.append(dict(
                dataset=dataset, split_type=split_type,
                coverage_80=cov["coverage_80"], coverage_90=cov["coverage_90"],
                coverage_95=cov["coverage_95"], n_test=cov["n_test"],
                n_calibration=cov["n_calibration"]))
        per_ds_block[dataset] = {"random": rc_cov, "block": blk_cov}

        log(f"  total {time.time()-t0:.1f}s\n")

    # ---------------------- write RESULT 1 ----------------------------
    df_rc = pd.DataFrame(rc_rows)
    df_rc.to_csv(os.path.join(RC_OUT, "risk_coverage_data.csv"), index=False)

    # summary: RMSE improvement (%) uncertainty vs random at 80% and 50%
    def improvement_at(frac):
        imps, rand_rmse_all = [], []
        for ds in df_rc["dataset"].unique():
            s = df_rc[(df_rc["dataset"] == ds) & (df_rc["retention_fraction"] == frac)]
            u = s.loc[s["method"] == "uncertainty", "rmse"].values
            r = s.loc[s["method"] == "random", "rmse"].values
            if len(u) and len(r) and r[0] > 0:
                imps.append(float((r[0] - u[0]) / r[0] * 100.0))
                rand_rmse_all.append(r[0])
        return imps

    imp80 = improvement_at(0.80)
    imp50 = improvement_at(0.50)
    # paired t-test: per-dataset RMSE(uncertainty) vs RMSE(random) at 80%
    def paired_at(frac):
        u_list, r_list = [], []
        for ds in df_rc["dataset"].unique():
            s = df_rc[(df_rc["dataset"] == ds) & (df_rc["retention_fraction"] == frac)]
            u = s.loc[s["method"] == "uncertainty", "rmse"].values
            r = s.loc[s["method"] == "random", "rmse"].values
            if len(u) and len(r):
                u_list.append(u[0]); r_list.append(r[0])
        if len(u_list) >= 2:
            t, p = ttest_rel(u_list, r_list)
            return float(t), float(p)
        return float("nan"), float("nan")

    _, p80 = paired_at(0.80)

    rc_summary = {
        "n_datasets": int(df_rc["dataset"].nunique()),
        "datasets": list(df_rc["dataset"].unique()),
        "retention_fractions": RETENTION_FRACS,
        "mean_rmse_improvement_pct_at_80": float(np.mean(imp80)) if imp80 else float("nan"),
        "mean_rmse_improvement_pct_at_50": float(np.mean(imp50)) if imp50 else float("nan"),
        "per_dataset_improvement_pct_at_80": dict(zip(
            list(df_rc["dataset"].unique()), imp80)),
        "paired_ttest_pvalue_at_80": p80,
        "random_shuffles": N_RANDOM_SHUFFLES,
        "method": "Method D (residual_spot, 2-pass EB) + split-conformal (alpha=0.10) "
                  "half-width as uncertainty score",
        "interpretation": (
            "Positive improvement means uncertainty-ranked retention beats random: "
            "uncertainty is actionable for triage." if np.nanmean(imp80) > 0
            else "Uncertainty ranking does NOT beat random: uncertainty is not "
                 "informative for triage."),
    }
    with open(os.path.join(RC_OUT, "summary.json"), "w") as f:
        json.dump(rc_summary, f, indent=2)
    plot_risk_coverage(df_rc, rc_summary, os.path.join(RC_OUT, "figure.png"))
    log(f"[RC] wrote risk_coverage_data.csv, summary.json, figure.png")

    # ---------------------- write RESULT 2 ----------------------------
    df_sb = pd.DataFrame(sb_rows)
    df_sb.to_csv(os.path.join(SB_OUT, "block_split_coverage.csv"), index=False)
    rand90 = [per_ds_block[d]["random"]["coverage_90"] for d in per_ds_block]
    block90 = [per_ds_block[d]["block"]["coverage_90"] for d in per_ds_block]
    max_dev = float(np.max(np.abs(np.array(block90) - 0.90)))
    mean_block = float(np.mean(block90))
    mean_rand = float(np.mean(rand90))
    # concern supported if block coverage deviates from nominal by > 5pp on
    # average (i.e. drops well below 0.85).
    concern = bool(mean_block < 0.85 or max_dev > 0.10)
    sb_summary = {
        "n_datasets": int(len(per_ds_block)),
        "datasets": list(per_ds_block.keys()),
        "n_blocks": N_BLOCKS,
        "n_cal_blocks": N_CAL_BLOCKS,
        "mean_coverage_random_90": mean_rand,
        "mean_coverage_block_90": mean_block,
        "max_abs_deviation_block_from_nominal": max_dev,
        "mean_abs_deviation_random_from_nominal": float(
            np.mean(np.abs(np.array(rand90) - 0.90))),
        "mean_abs_deviation_block_from_nominal": float(
            np.mean(np.abs(np.array(block90) - 0.90))),
        "exchangeability_concern_supported": concern,
        "verdict": (
            "Spatial autocorrelation IS a concern: block-split coverage drops "
            "materially below nominal." if concern
            else "Exchangeability HOLDS: block-split coverage stays near nominal "
                 "(strong defense against the spatial-autocorrelation attack)."),
        "_per_dataset": {d: {"random": {k: v for k, vv in per_ds_block[d]["random"].items()
                                         for k2, v in [(k, vv)]},
                              "block": {k: v for k, vv in per_ds_block[d]["block"].items()
                                         for k2, v in [(k, vv)]}}
                         for d in per_ds_block},
    }
    # fix the nested dict construction above (it was malformed)
    sb_summary["_per_dataset"] = {
        d: {"random": dict(per_ds_block[d]["random"]),
            "block": dict(per_ds_block[d]["block"])}
        for d in per_ds_block
    }
    with open(os.path.join(SB_OUT, "summary.json"), "w") as f:
        json.dump(sb_summary, f, indent=2)
    plot_block_coverage(df_sb, sb_summary, os.path.join(SB_OUT, "figure.png"))
    log(f"[SB] wrote block_split_coverage.csv, summary.json, figure.png")

    log("\n=== DONE ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Conformal uncertainty-calibration validation across ALL processed Visium APA
datasets.

Mirrors scripts/calibrate_uncertainty.py (single-dataset split-conformal
benchmark) but loops over a curated, representative set of samples spanning
all 9 GSE accessions.  The goal (BIB Pillar 1) is to show that conformal
calibration achieves target marginal coverage (80/90/95%) *universally* --
i.e. on every tissue / organism / library-prep combo we have, not just on one
or two hand-picked datasets.

Protocol (per sample)
---------------------
1. Load apa_matrix.csv (site x spot) + coordinates.csv (x,y).
2. Subsample to the top-K most-observed sites (default 300) -- the GP fits one
   model per row, so we keep per-sample runtime < ~3 min.  Coverage statistics
   only need a large pool of held-out points, which the remaining rows provide.
3. Mask 20% of observed entries per site.
4. Fit SparseGPImputerBatch (n_inducing=min(500,max(100,n_spots//100)),
   length_scale=cKDTree NN dist x 5, noise_level=0.1).
5. Split held-out into calibration (50%) + test (50%).
6. Fit conformal calibrator (locally_adaptive) on calibration errors.
7. Report on test: raw unc-error Pearson, coverage at 80/90/95%, RMSE, RMSE
   after dropping top-20% uncertainty.

Outputs
-------
pipeline_output/conformal_validation/all_samples_coverage.csv
pipeline_output/conformal_validation/summary.json
pipeline_output/conformal_validation/coverage_figure.png
pipeline_output/conformal_validation/run.log
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

# --- repo-local imports (must run after sys.path has repo root) ---
REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402
from spagapa.imputation.calibration import (  # noqa: E402
    ConformalCalibrator,
    evaluate_coverage,
)

DATA_ROOT = os.path.join(REPO, "data", "processed")
OUT_DIR = os.path.join(REPO, "pipeline_output", "conformal_validation")

# Curated representative sample list (1-2 per GSE; for GSE237183 pick a spread
# of spot counts).  Each tuple: (dataset_label, sample_label, dir_basename)
SAMPLES = [
    ("GSE237183", "gsm7596587", "gse237183_gsm7596587_scapatrap"),   # 2541 spots
    ("GSE237183", "gsm7596595", "gse237183_gsm7596595_scapatrap"),   # 1899 spots (small)
    ("GSE237183", "gsm7596604", "gse237183_gsm7596604_scapatrap"),   # 4881 spots (large)
    ("GSE183456", "gsm6047774", "gse183456_gsm6047774_scapatrap"),   # 3010
    ("GSE179572", "gsm5420751", "gse179572_gsm5420751_scapatrap"),   # 4992
    ("GSE220442", "gsm6801751", "gse220442_gsm6801751_scapatrap"),   # 4179
    ("GSE220442", "gsm6801753", "gse220442_gsm6801753_scapatrap"),   # 4976
    ("GSE169749", "gsm5213483", "gse169749_gsm5213483_scapatrap"),   # mouse colon 2715
    ("GSE206391", "gsm6252925", "gse206391_gsm6252925_scapatrap"),   # skin 796 (LOW genes)
    ("GSE263303", "gsm8189356", "gse263303_gsm8189356_scapatrap"),   # mouse brain 1554
    ("GSE338525", "gsm9876373", "gse338525_gsm9876373_scapatrap"),   # liver 4992
    ("GSE338525", "gsm9876374", "gse338525_gsm9876374_scapatrap"),   # liver 4992
]

# ---- hyperparameters mirroring calibrate_uncertainty.py ----
MASK_FRACTION = 0.20
CAL_FRACTION = 0.50
FILTER_FRAC = 0.20
SEED = 42
MAX_GENES = 300          # subsample top-K most-observed sites per sample
MIN_SPOTS = 500
MIN_GENES_WITH_OBS = 50


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


def _safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) >= 2 and np.std(a) > 0 and np.std(b) > 0:
        try:
            r, _ = pearsonr(a, b)
            return float(r)
        except Exception:
            return float("nan")
    return float("nan")


def _rmse_filter(err: np.ndarray, score: np.ndarray, frac: float) -> float:
    """RMSE after dropping the top `frac` points by `score` (uncertainty)."""
    err = np.asarray(err, dtype=float)
    score = np.asarray(score, dtype=float)
    n = len(err)
    if n == 0:
        return float("nan")
    n_keep = int(np.ceil(n * (1.0 - frac)))
    order = np.argsort(score)
    keep = order[:n_keep]
    return float(np.sqrt(np.mean(err[keep] ** 2)))


def load_sample(dir_basename: str):
    """Load apa_matrix + coordinates, return (values, xy) aligned."""
    apa = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "apa_matrix.csv"),
                      index_col=0)
    coords = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "coordinates.csv"),
                         index_col=0)
    if {"x", "y"} <= set(coords.columns):
        xy_df = coords[["x", "y"]]
    else:
        xy_df = coords.iloc[:, :2].astype(float)
        xy_df.columns = ["x", "y"]
    # align columns <-> index
    common = apa.columns.intersection(xy_df.index)
    if len(common) > 0:
        apa = apa.loc[:, common]
        xy_df = xy_df.loc[common]
    else:
        # fall back to positional if names mismatch but counts match
        if apa.shape[1] != xy_df.shape[0]:
            raise ValueError(
                f"spot mismatch and no common names: apa={apa.shape}, "
                f"coords={xy_df.shape}")
    xy = xy_df.values.astype(float)
    values = apa.values.astype(float)
    return values, xy


def run_one(dataset: str, sample: str, dir_basename: str, rng: np.random.Generator):
    t0 = time.time()
    values, xy = load_sample(dir_basename)
    n_genes_full, n_spots = values.shape

    # subsample top-K most-observed rows (sites) for tractable GP fit
    obs_count = (values > 0).sum(axis=1)
    keep_rows = np.where(obs_count >= 5)[0]
    if len(keep_rows) == 0:
        return None, "no genes with >=5 observations"
    # rank by observation count desc, take top MAX_GENES
    order = keep_rows[np.argsort(-obs_count[keep_rows])]
    keep = order[:MAX_GENES]
    values_sub = values[keep]
    n_genes = values_sub.shape[0]

    if n_spots < MIN_SPOTS:
        return None, f"too few spots ({n_spots} < {MIN_SPOTS})"
    if n_genes < MIN_GENES_WITH_OBS:
        return None, f"too few genes with obs ({n_genes} < {MIN_GENES_WITH_OBS})"

    log(f"  {dataset}/{sample}: {n_genes} sites (of {n_genes_full}), "
        f"{n_spots} spots, observed frac={(values_sub>0).mean():.3f}")

    # ---- masking 20% per row ----
    # Held-out entries are identified downstream as (truth>0) & (masked==0),
    # so we just zero out the chosen fraction per row here.
    masked = values_sub.copy()
    for g in range(n_genes):
        obs_idx = np.where(values_sub[g] > 0)[0]
        if len(obs_idx) < 5:
            continue
        n_mask = max(1, int(len(obs_idx) * MASK_FRACTION))
        m = rng.choice(obs_idx, size=n_mask, replace=False)
        masked[g, m] = 0.0

    # ---- GP fit ----
    kdt = cKDTree(xy)
    nn = kdt.query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    length_scale = nn_dist * 5
    n_inducing = min(500, max(100, n_spots // 100))
    base = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale,
                           noise_level=0.1)
    train_mask = masked > 0
    tg = time.time()
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    log(f"  GP fit+impute: {time.time()-tg:.1f}s "
        f"(n_inducing={n_inducing}, length_scale={length_scale:.1f})")

    return _finish(dataset, sample, n_genes, n_spots, nn_dist, length_scale,
                   n_inducing, masked, values_sub, gp_pred, gp_unc, rng, t0)


def _finish(dataset, sample, n_genes, n_spots, nn_dist, length_scale,
            n_inducing, masked, values_sub, gp_pred, gp_unc, rng, t0):
    """Re-gather held-out predictions cleanly (recompute masks deterministically)."""
    # Recompute the held-out mask the SAME way run_one did, but using the same
    # rng state is not possible here.  Instead we identify held-out entries as
    # those that were observed (truth > 0) but are now zero in `masked`.
    truth_obs = values_sub > 0
    masked_zero = masked <= 0
    held_mask = truth_obs & masked_zero   # (n_genes, n_spots) bool
    n_held_total = int(held_mask.sum())
    if n_held_total < 20:
        return None, f"too few held-out points ({n_held_total})"

    g_idx, s_idx = np.where(held_mask)
    truths = values_sub[g_idx, s_idx]
    preds = gp_pred[g_idx, s_idx]
    unc = gp_unc[g_idx, s_idx]
    errors = np.abs(truths - preds)

    # split calibration / test (50/50)
    perm = rng.permutation(len(truths))
    n_cal = int(len(truths) * CAL_FRACTION)
    cal_i = perm[:n_cal]
    test_i = perm[n_cal:]
    cal_err = errors[cal_i]
    cal_std = unc[cal_i]
    test_truth = truths[test_i]
    test_pred = preds[test_i]
    test_std = unc[test_i]
    test_err = errors[test_i]

    raw_corr = _safe_corr(test_std, test_err)

    rmse_full = float(np.sqrt(np.mean(test_err ** 2)))
    rmse_after = _rmse_filter(test_err, test_std, FILTER_FRAC)
    rmse_improvement = (rmse_full - rmse_after) / rmse_full * 100.0 \
        if rmse_full > 0 else float("nan")

    # conformal locally-adaptive at 80/90/95
    row = {
        "dataset": dataset,
        "sample": sample,
        "n_genes": int(n_genes),
        "n_spots": int(n_spots),
        "n_held_out": int(len(truths)),
        "n_test": int(len(test_i)),
        "raw_unc_error_corr": raw_corr,
        "rmse": rmse_full,
        "rmse_after_filter": rmse_after,
        "rmse_improvement_pct": rmse_improvement,
    }
    for alpha, level in [(0.20, 80), (0.10, 90), (0.05, 95)]:
        cal = ConformalCalibrator(alpha=alpha, mode="locally_adaptive")
        cal.fit(cal_err, cal_std)
        lo, hi = cal.predict(test_pred, test_std)
        cov = evaluate_coverage(lo, hi, test_truth)
        row[f"coverage_{level}"] = cov
        row[f"qhat_{level}"] = float(cal.interval_.q_hat)

    row["gp_length_scale"] = float(length_scale)
    row["gp_n_inducing"] = int(n_inducing)
    row["gp_nn_dist"] = float(nn_dist)
    row["wall_s"] = float(time.time() - t0)
    return row, None


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    rng = np.random.default_rng(SEED)

    log_path = os.path.join(OUT_DIR, "run.log")
    log_fh = open(log_path, "w")

    def both(msg):
        log(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    both(f"Conformal validation over {len(SAMPLES)} candidate samples")
    both(f"MASK_FRACTION={MASK_FRACTION} CAL_FRACTION={CAL_FRACTION} "
         f"FILTER_FRAC={FILTER_FRAC} SEED={SEED} MAX_GENES={MAX_GENES}")
    both("=" * 70)

    rows = []
    skipped = []
    for dataset, sample, dir_basename in SAMPLES:
        both("")
        both(f">>> {dataset}/{sample} ({dir_basename})")
        try:
            row, err = run_one(dataset, sample, dir_basename, rng)
            if row is None:
                both(f"  SKIP: {err}")
                skipped.append({"dataset": dataset, "sample": sample,
                                "reason": err})
            else:
                rows.append(row)
                both(f"  corr={row['raw_unc_error_corr']:.4f}  "
                     f"cov80={row['coverage_80']:.3f} "
                     f"cov90={row['coverage_90']:.3f} "
                     f"cov95={row['coverage_95']:.3f}  "
                     f"rmse={row['rmse']:.4f}->"
                     f"{row['rmse_after_filter']:.4f} "
                     f"({row['rmse_improvement_pct']:+.1f}%)  "
                     f"[{row['wall_s']:.0f}s]")
        except Exception as e:
            tb = traceback.format_exc()
            both(f"  ERROR: {e}\n{tb}")
            skipped.append({"dataset": dataset, "sample": sample,
                            "reason": f"error: {e}"})

    if not rows:
        both("No samples completed successfully.")
        log_fh.close()
        return 1

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "all_samples_coverage.csv")
    df.to_csv(csv_path, index=False)
    both(f"\nWrote {csv_path}")

    # ---- summary ----
    summary = {
        "n_samples_validated": int(len(rows)),
        "n_samples_skipped": int(len(skipped)),
        "skipped": skipped,
        "targets": {"coverage_80": 0.80, "coverage_90": 0.90,
                    "coverage_95": 0.95},
        "mean_coverage_80": float(df["coverage_80"].mean()),
        "mean_coverage_90": float(df["coverage_90"].mean()),
        "mean_coverage_95": float(df["coverage_95"].mean()),
        "median_coverage_80": float(df["coverage_80"].median()),
        "median_coverage_90": float(df["coverage_90"].median()),
        "median_coverage_95": float(df["coverage_95"].median()),
        "std_coverage_80": float(df["coverage_80"].std()),
        "std_coverage_90": float(df["coverage_90"].std()),
        "std_coverage_95": float(df["coverage_95"].std()),
        "min_coverage_80": float(df["coverage_80"].min()),
        "min_coverage_90": float(df["coverage_90"].min()),
        "min_coverage_95": float(df["coverage_95"].min()),
        "max_coverage_95": float(df["coverage_95"].max()),
        "mean_raw_unc_error_corr": float(df["raw_unc_error_corr"].mean()),
        "median_raw_unc_error_corr": float(df["raw_unc_error_corr"].median()),
        "mean_rmse_improvement_pct": float(df["rmse_improvement_pct"].mean()),
        "median_rmse_improvement_pct": float(df["rmse_improvement_pct"].median()),
        "mean_rmse": float(df["rmse"].mean()),
        "total_genes_tested": int(df["n_genes"].sum()),
        "total_test_points": int(df["n_test"].sum()),
        "datasets_covered": sorted(df["dataset"].unique().tolist()),
    }
    # absolute deviation from target
    summary["mean_abs_dev_80"] = float((df["coverage_80"] - 0.80).abs().mean())
    summary["mean_abs_dev_90"] = float((df["coverage_90"] - 0.90).abs().mean())
    summary["mean_abs_dev_95"] = float((df["coverage_95"] - 0.95).abs().mean())
    # fraction within +/- 0.05 of target
    for lvl, tgt in [("80", 0.80), ("90", 0.90), ("95", 0.95)]:
        within = ((df[f"coverage_{lvl}"] - tgt).abs() <= 0.05).mean()
        summary[f"frac_within_5pct_{lvl}"] = float(within)

    json_path = os.path.join(OUT_DIR, "summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    both(f"Wrote {json_path}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(df)
        x = np.arange(n)
        w = 0.27
        fig, ax = plt.subplots(figsize=(max(10, n * 0.7), 6))
        ax.bar(x - w, df["coverage_80"], w, label="80%", color="#4C72B0")
        ax.bar(x,     df["coverage_90"], w, label="90%", color="#55A868")
        ax.bar(x + w, df["coverage_95"], w, label="95%", color="#C44E52")
        # target lines
        for lvl, tgt, off in [("80", 0.80, -w), ("90", 0.90, 0.0),
                              ("95", 0.95, w)]:
            ax.hlines(tgt, x[0] - 0.5, x[-1] + 0.5, colors="k",
                      linestyles="--", linewidths=0.8, alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r['dataset'][:8]}\n{r['sample']}"
                            for _, r in df.iterrows()],
                           rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Empirical coverage (test set)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("Conformal calibration coverage across Visium APA datasets\n"
                     "(dashed = nominal target; locally-adaptive split conformal)")
        ax.legend(loc="lower right", ncol=3)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        png_path = os.path.join(OUT_DIR, "coverage_figure.png")
        fig.savefig(png_path, dpi=140)
        both(f"Wrote {png_path}")
    except Exception as e:
        both(f"figure error: {e}")

    # ---- console summary ----
    both("")
    both("=" * 70)
    both("SUMMARY")
    both("=" * 70)
    both(f"samples validated : {summary['n_samples_validated']}  "
         f"(skipped {summary['n_samples_skipped']})")
    both(f"datasets covered  : {len(summary['datasets_covered'])} "
         f"{summary['datasets_covered']}")
    both(f"mean coverage 80/90/95: "
         f"{summary['mean_coverage_80']:.3f} / "
         f"{summary['mean_coverage_90']:.3f} / "
         f"{summary['mean_coverage_95']:.3f}  "
         f"(target 0.80/0.90/0.95)")
    both(f"mean |dev| from target: "
         f"{summary['mean_abs_dev_80']:.3f} / "
         f"{summary['mean_abs_dev_90']:.3f} / "
         f"{summary['mean_abs_dev_95']:.3f}")
    both(f"frac within +/-0.05 of target: "
         f"{summary['frac_within_5pct_80']:.2f} / "
         f"{summary['frac_within_5pct_90']:.2f} / "
         f"{summary['frac_within_5pct_95']:.2f}")
    both(f"mean raw unc-error corr : {summary['mean_raw_unc_error_corr']:.4f}")
    both(f"mean RMSE improvement   : "
         f"{summary['mean_rmse_improvement_pct']:+.2f}%")

    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

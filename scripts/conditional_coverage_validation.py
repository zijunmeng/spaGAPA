#!/usr/bin/env python3
"""
Conditional conformal-coverage validation.

Motivation
----------
The existing `calibrate_uncertainty_all_datasets.py` shows that split conformal
calibration hits its TARGET MARGINAL coverage (0.80/0.90/0.95) on every
tissue.  A reviewer can still object:

    "Coverage is only marginal (averaged over all points).  Spatial
    autocorrelation makes nearby points correlated, so exchangeability is
    violated, and the 90% number could be the average of over-coverage in one
    region and under-coverage in another."

Split conformal only *guarantees* marginal coverage under exchangeability; it
offers NO finite-sample guarantee conditional on a subgroup.  The honest way
to defend the calibration is therefore EMPIRICAL: show that coverage also
holds *conditionally* across the subgroups a reviewer would worry about:

    A. GP posterior-std quintile     (does the interval scale track error?)
    B. spatial region (k-means)      (does spatial autocorrelation break it?)
    C. gene expression level         (does it hold at low vs high counts?)
    D. tissue / dataset type         (already in marginal validation)

Protocol
--------
Per sample (5 representative, diverse tissues), replicate the EXACT pipeline
of `calibrate_uncertainty_all_datasets.py`:
    load -> top-300 sites -> mask 20% per row -> SparseGP fit ->
    split held-out into calibration (50%) + test (50%) ->
    fit locally-adaptive split conformal at alpha=0.1 ->
    predict 90% intervals on the TEST set.

Then, ON THE TEST SET ONLY (calibration set is never used for evaluation),
bin the test points by:
    - GP std quintile            (Q1=lowest ... Q5=highest uncertainty)
    - k-means spatial region     (~6 clusters on (x,y))
    - binned truth (expression)  (5 quantile bins of the observed truth)

and report empirical 90% coverage within each bin.

If every bin is ~0.90, conditional coverage is empirically reasonable ->
strong defence against the exchangeability / spatial-autocorrelation critique.
If some bin deviates, we report it honestly.

Outputs (all under pipeline_output/conformal_conditional_coverage/)
    by_uncertainty_quintile.csv
    by_spatial_region.csv
    by_expression_level.csv
    by_tissue.csv
    summary.json
    conditional_coverage_figure.png
    run.log
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
from sklearn.cluster import KMeans

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402
from spagapa.imputation.calibration import (  # noqa: E402
    ConformalCalibrator,
    evaluate_coverage,
)

DATA_ROOT = os.path.join(REPO, "data", "processed")
OUT_DIR = os.path.join(REPO, "pipeline_output", "conformal_conditional_coverage")

# 5 representative, diverse samples (Visium brain, kidney, skin; mouse colon,
# mouse brain; liver).  These are a subset of the 12 in
# calibrate_uncertainty_all_datasets.SAMPLES so the marginal numbers agree.
SAMPLES = [
    ("GSE237183", "gsm7596587", "gse237183_gsm7596587_scapatrap", "Visium-brain"),
    ("GSE183456", "gsm6047774", "gse183456_gsm6047774_scapatrap", "Visium-kidney"),
    ("GSE206391", "gsm6252925", "gse206391_gsm6252925_scapatrap", "Visium-skin"),
    ("GSE169749", "gsm5213483", "gse169749_gsm5213483_scapatrap", "Mouse-colon"),
    ("GSE338525", "gsm9876373", "gse338525_gsm9876373_scapatrap", "Visium-liver"),
]

# ---- hyperparameters mirroring calibrate_uncertainty_all_datasets.py ----
MASK_FRACTION = 0.20
CAL_FRACTION = 0.50
SEED = 42
MAX_GENES = 300
MIN_SPOTS = 500
MIN_GENES_WITH_OBS = 50
ALPHA = 0.10          # 90% target interval
TARGET = 1 - ALPHA    # 0.90
N_UNC_BINS = 5
N_SPATIAL_CLUSTERS = 6
N_EXPR_BINS = 5


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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
    else:
        if apa.shape[1] != xy_df.shape[0]:
            raise ValueError("spot mismatch and no common names")
    xy = xy_df.values.astype(float)
    values = apa.values.astype(float)
    return values, xy


def fit_sample(dataset, sample, dir_basename, rng):
    """Replicate the GP + conformal pipeline; return test-set arrays.

    Returns dict with the TEST-set (only) quantities needed for conditional
    evaluation: truth, pred, gp_std, abs_err, interval bounds, spot xy.
    """
    values, xy = load_sample(dir_basename)
    n_genes_full, n_spots = values.shape

    obs_count = (values > 0).sum(axis=1)
    keep_rows = np.where(obs_count >= 5)[0]
    order = keep_rows[np.argsort(-obs_count[keep_rows])]
    keep = order[:MAX_GENES]
    values_sub = values[keep]
    n_genes = values_sub.shape[0]
    if n_spots < MIN_SPOTS or n_genes < MIN_GENES_WITH_OBS:
        return None, "too few spots/genes"

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
    base = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale,
                           noise_level=0.1)
    train_mask = masked > 0
    tg = time.time()
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    log(f"  GP fit+impute: {time.time()-tg:.1f}s "
        f"(n_inducing={n_inducing}, length_scale={length_scale:.1f})")

    # held-out = truth>0 & masked<=0
    truth_obs = values_sub > 0
    held_mask = truth_obs & (masked <= 0)
    g_idx, s_idx = np.where(held_mask)
    truths = values_sub[g_idx, s_idx]
    preds = gp_pred[g_idx, s_idx]
    unc = gp_unc[g_idx, s_idx]
    errors = np.abs(truths - preds)
    held_xy = xy[s_idx]                      # spot coordinate per held-out point

    perm = rng.permutation(len(truths))
    n_cal = int(len(truths) * CAL_FRACTION)
    test_i = perm[n_cal:]
    cal_i = perm[:n_cal]

    cal = ConformalCalibrator(alpha=ALPHA, mode="locally_adaptive")
    cal.fit(errors[cal_i], unc[cal_i])
    lo_all, hi_all = cal.predict(preds, unc)

    return {
        "dataset": dataset, "sample": sample,
        "n_spots": int(n_spots), "n_genes": int(n_genes),
        "n_held": int(len(truths)),
        "n_cal": int(n_cal), "n_test": int(len(test_i)),
        "test_truth": truths[test_i],
        "test_pred": preds[test_i],
        "test_std": unc[test_i],
        "test_err": errors[test_i],
        "test_lo": lo_all[test_i],
        "test_hi": hi_all[test_i],
        "test_xy": held_xy[test_i],
        "qhat90": float(cal.interval_.q_hat),
    }, None


def coverage(lo, hi, y):
    return float(np.mean((y >= lo) & (y <= hi)))


def bin_by_uncertainty(d):
    """Sort test points by GP std into N_UNC_BINS quintiles."""
    std = d["test_std"]
    y = d["test_truth"]
    lo, hi = d["test_lo"], d["test_hi"]
    # use quantile cuts; labels=False -> integer bin 0..N-1
    # bins computed on this sample's own test_std distribution
    qs = np.quantile(std, np.linspace(0, 1, N_UNC_BINS + 1))
    # ensure strictly increasing edges (collapse dups)
    edges = np.unique(qs)
    if len(edges) < 2:
        return [(1, std.size, coverage(lo, hi, y))]
    bins = np.clip(np.digitize(std, edges[1:-1], right=False), 0, len(edges) - 2)
    out = []
    for b in range(len(edges) - 1):
        m = bins == b
        out.append((b + 1, int(m.sum()),
                    coverage(lo[m], hi[m], y[m]) if m.any() else float("nan")))
    return out


def bin_by_region(d, rng):
    xy = d["test_xy"]
    y = d["test_truth"]
    lo, hi = d["test_lo"], d["test_hi"]
    n_clusters = min(N_SPATIAL_CLUSTERS, max(2, len(xy) // 100))
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=SEED)
    labels = km.fit_predict(xy)
    out = []
    for b in range(n_clusters):
        m = labels == b
        out.append((b, int(m.sum()),
                    coverage(lo[m], hi[m], y[m]) if m.any() else float("nan")))
    return out


def bin_by_expression(d):
    """Bin by the observed truth (gene expression) into 5 quantile bins."""
    y = d["test_truth"]
    lo, hi = d["test_lo"], d["test_hi"]
    qs = np.quantile(y, np.linspace(0, 1, N_EXPR_BINS + 1))
    edges = np.unique(qs)
    if len(edges) < 2:
        return [(1, y.size, coverage(lo, hi, y))]
    bins = np.clip(np.digitize(y, edges[1:-1], right=False), 0, len(edges) - 2)
    labels = ["very_low", "low", "medium", "high", "very_high"]
    out = []
    for b in range(len(edges) - 1):
        m = bins == b
        out.append((labels[b] if b < len(labels) else f"bin{b}",
                    int(m.sum()),
                    coverage(lo[m], hi[m], y[m]) if m.any() else float("nan")))
    return out


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

    both(f"Conditional conformal coverage validation over {len(SAMPLES)} "
         f"representative samples (alpha={ALPHA}, target={TARGET})")
    both("=" * 70)

    unc_rows, reg_rows, expr_rows, tissue_rows = [], [], [], []
    all_test_cov = []   # for histogram panel D
    fitted = []

    for dataset, sample, dir_basename, tissue in SAMPLES:
        both("")
        both(f">>> {dataset}/{sample} ({tissue})")
        try:
            d, err = fit_sample(dataset, sample, dir_basename, rng)
            if d is None:
                both(f"  SKIP: {err}")
                continue
        except Exception as e:
            both(f"  ERROR: {e}\n{traceback.format_exc()}")
            continue

        overall = coverage(d["test_lo"], d["test_hi"], d["test_truth"])
        both(f"  n_test={d['n_test']}  overall 90% coverage={overall:.4f}  "
             f"qhat90={d['qhat90']:.3f}")
        all_test_cov.append(overall)

        # A. uncertainty quintile
        for q, n, cov in bin_by_uncertainty(d):
            unc_rows.append({
                "dataset": dataset, "sample": sample, "tissue": tissue,
                "quintile": q, "n_points": n, "coverage_90": cov,
                "target_90": TARGET,
                "deviation": cov - TARGET if not np.isnan(cov) else float("nan"),
            })
        both("  uncertainty-quintile coverage: " +
             " ".join(f"Q{r['quintile']}={r['coverage_90']:.3f}"
                      for r in unc_rows if r["sample"] == sample))

        # B. spatial region
        for rid, n, cov in bin_by_region(d, rng):
            reg_rows.append({
                "dataset": dataset, "sample": sample, "tissue": tissue,
                "region_id": rid, "n_points": n, "coverage_90": cov,
                "deviation": cov - TARGET if not np.isnan(cov) else float("nan"),
            })
        both("  spatial-region coverage: " +
             " ".join(f"R{r['region_id']}={r['coverage_90']:.3f}"
                      for r in reg_rows if r["sample"] == sample))

        # C. expression level
        for lab, n, cov in bin_by_expression(d):
            expr_rows.append({
                "dataset": dataset, "sample": sample, "tissue": tissue,
                "expr_bin": lab, "n_points": n, "coverage_90": cov,
                "deviation": cov - TARGET if not np.isnan(cov) else float("nan"),
            })
        both("  expression-bin coverage: " +
             " ".join(f"{r['expr_bin']}={r['coverage_90']:.3f}"
                      for r in expr_rows if r["sample"] == sample))

        # D. tissue (single row per sample)
        tissue_rows.append({
            "dataset": dataset, "sample": sample, "tissue": tissue,
            "n_points": d["n_test"], "coverage_90": overall,
            "deviation": overall - TARGET,
        })

        fitted.append(d)

    log_fh.flush()
    if not fitted:
        both("No samples completed.")
        log_fh.close()
        return 1

    # ---- write CSVs ----
    df_unc = pd.DataFrame(unc_rows)
    df_reg = pd.DataFrame(reg_rows)
    df_expr = pd.DataFrame(expr_rows)
    df_tissue = pd.DataFrame(tissue_rows)
    df_unc.to_csv(os.path.join(OUT_DIR, "by_uncertainty_quintile.csv"), index=False)
    df_reg.to_csv(os.path.join(OUT_DIR, "by_spatial_region.csv"), index=False)
    df_expr.to_csv(os.path.join(OUT_DIR, "by_expression_level.csv"), index=False)
    df_tissue.to_csv(os.path.join(OUT_DIR, "by_tissue.csv"), index=False)
    both("Wrote by_uncertainty_quintile.csv / by_spatial_region.csv / "
         "by_expression_level.csv / by_tissue.csv")

    # ---- summary stats ----
    def dev_stats(df, key):
        dev = df["deviation"].abs().dropna()
        return {
            f"mean_abs_dev_{key}": float(dev.mean()) if len(dev) else float("nan"),
            f"max_abs_dev_{key}": float(dev.max()) if len(dev) else float("nan"),
            f"median_abs_dev_{key}": float(dev.median()) if len(dev) else float("nan"),
            f"frac_within_0.05_{key}": float((dev <= 0.05).mean()) if len(dev) else float("nan"),
            f"frac_within_0.10_{key}": float((dev <= 0.10).mean()) if len(dev) else float("nan"),
            f"n_bins_{key}": int(len(dev)),
        }

    summary = {
        "n_samples": len(fitted),
        "target_90": TARGET,
        "samples": [{"dataset": d["dataset"], "sample": d["sample"],
                     "n_test": d["n_test"]} for d in fitted],
        "mean_overall_coverage": float(np.mean(all_test_cov)),
        "min_overall_coverage": float(np.min(all_test_cov)),
        "max_overall_coverage": float(np.max(all_test_cov)),
    }
    summary.update(dev_stats(df_unc, "by_quintile"))
    summary.update(dev_stats(df_reg, "by_region"))
    summary.update(dev_stats(df_expr, "by_expression"))
    summary.update(dev_stats(df_tissue, "by_tissue"))

    # combined across all three subgroup types
    all_dev = pd.concat([df_unc["deviation"], df_reg["deviation"],
                         df_expr["deviation"]]).abs().dropna()
    summary["mean_abs_dev_all_subgroups"] = float(all_dev.mean())
    summary["max_abs_dev_all_subgroups"] = float(all_dev.max())
    summary["frac_within_0.05_all_subgroups"] = float((all_dev <= 0.05).mean())
    summary["frac_within_0.10_all_subgroups"] = float((all_dev <= 0.10).mean())
    summary["n_all_subgroups"] = int(len(all_dev))

    # worst-case subgroup
    combined = pd.concat([
        df_unc.assign(subgroup_type="uncertainty_quintile",
                      subgroup_id=df_unc["quintile"].astype(str)),
        df_reg.assign(subgroup_type="spatial_region",
                      subgroup_id=df_reg["region_id"].astype(str)),
        df_expr.assign(subgroup_type="expression_level",
                       subgroup_id=df_expr["expr_bin"].astype(str)),
    ], ignore_index=True)
    worst = combined.reindex(combined["deviation"].abs().fillna(-1).sort_values().index)
    worst_row = worst.iloc[-1]
    summary["worst_subgroup"] = {
        "sample": str(worst_row["sample"]),
        "type": str(worst_row["subgroup_type"]),
        "id": str(worst_row["subgroup_id"]),
        "n_points": int(worst_row["n_points"]),
        "coverage_90": float(worst_row["coverage_90"]),
        "deviation": float(worst_row["deviation"]),
    }
    # top 5 worst
    summary["top5_worst_subgroups"] = [
        {"sample": str(r["sample"]), "type": str(r["subgroup_type"]),
         "id": str(r["subgroup_id"]), "n": int(r["n_points"]),
         "coverage_90": float(r["coverage_90"]),
         "deviation": float(r["deviation"])}
        for _, r in worst.tail(5).iterrows()
    ]

    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    both(f"Wrote summary.json")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        green, red = "#2ca02c", "#d62728"

        def add_target(ax):
            ax.axhline(TARGET, color="k", ls="--", lw=1.0, alpha=0.6,
                       label=f"target {TARGET:.2f}")

        # Panel A: uncertainty quintile boxplot (one row per sample-quintile)
        ax = axes[0, 0]
        df_unc["quintile_str"] = "Q" + df_unc["quintile"].astype(str)
        data_a = [df_unc.loc[df_unc["quintile"] == q, "coverage_90"].dropna().values
                  for q in range(1, N_UNC_BINS + 1)]
        pos_a = np.arange(1, N_UNC_BINS + 1)
        ax.boxplot(data_a, positions=pos_a, widths=0.5, showfliers=True,
                   patch_artist=True,
                   boxprops=dict(facecolor="#9ecae1", alpha=0.6))
        for q in range(1, N_UNC_BINS + 1):
            vals = df_unc.loc[df_unc["quintile"] == q, "coverage_90"].dropna()
            ax.scatter(np.full_like(vals, q, dtype=float) +
                       np.random.default_rng(0).uniform(-0.08, 0.08, len(vals)),
                       vals, color="#08519c", s=18, alpha=0.8, zorder=3)
        add_target(ax)
        ax.set_xticks(pos_a)
        ax.set_xticklabels([f"Q{q}\n({'low' if q==1 else 'high' if q==N_UNC_BINS else ''})"
                            for q in range(1, N_UNC_BINS + 1)])
        ax.set_ylim(0.70, 1.02)
        ax.set_ylabel("Empirical 90% coverage")
        ax.set_xlabel("GP posterior std quintile")
        ax.set_title(f"A. Coverage by GP-uncertainty quintile\n"
                     f"mean |dev|={summary['mean_abs_dev_by_quintile']:.3f}, "
                     f"max={summary['max_abs_dev_by_quintile']:.3f}")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)

        # Panel B: spatial region
        ax = axes[0, 1]
        max_reg = int(df_reg["region_id"].max()) + 1
        data_b = [df_reg.loc[df_reg["region_id"] == r, "coverage_90"].dropna().values
                  for r in range(max_reg)]
        ax.boxplot(data_b, positions=np.arange(max_reg), widths=0.5,
                   patch_artist=True,
                   boxprops=dict(facecolor="#a1d99b", alpha=0.6))
        for r in range(max_reg):
            vals = df_reg.loc[df_reg["region_id"] == r, "coverage_90"].dropna()
            ax.scatter(np.full_like(vals, r, dtype=float) +
                       np.random.default_rng(1).uniform(-0.08, 0.08, len(vals)),
                       vals, color="#31a354", s=18, alpha=0.8, zorder=3)
        add_target(ax)
        ax.set_xticks(np.arange(max_reg))
        ax.set_xticklabels([f"R{r}" for r in range(max_reg)])
        ax.set_ylim(0.70, 1.02)
        ax.set_ylabel("Empirical 90% coverage")
        ax.set_xlabel("Spatial region (k-means)")
        ax.set_title(f"B. Coverage by spatial region\n"
                     f"mean |dev|={summary['mean_abs_dev_by_region']:.3f}, "
                     f"max={summary['max_abs_dev_by_region']:.3f}")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)

        # Panel C: expression level
        ax = axes[1, 0]
        order_c = ["very_low", "low", "medium", "high", "very_high"]
        present = [b for b in order_c if b in df_expr["expr_bin"].unique()]
        data_c = [df_expr.loc[df_expr["expr_bin"] == b, "coverage_90"].dropna().values
                  for b in present]
        ax.boxplot(data_c, positions=np.arange(len(present)), widths=0.5,
                   patch_artist=True,
                   boxprops=dict(facecolor="#fdae6b", alpha=0.6))
        for i, b in enumerate(present):
            vals = df_expr.loc[df_expr["expr_bin"] == b, "coverage_90"].dropna()
            ax.scatter(np.full_like(vals, i, dtype=float) +
                       np.random.default_rng(2).uniform(-0.08, 0.08, len(vals)),
                       vals, color="#e6550d", s=18, alpha=0.8, zorder=3)
        add_target(ax)
        ax.set_xticks(np.arange(len(present)))
        ax.set_xticklabels(present, rotation=20)
        ax.set_ylim(0.70, 1.02)
        ax.set_ylabel("Empirical 90% coverage")
        ax.set_xlabel("Observed expression bin")
        ax.set_title(f"C. Coverage by expression level\n"
                     f"mean |dev|={summary['mean_abs_dev_by_expression']:.3f}, "
                     f"max={summary['max_abs_dev_by_expression']:.3f}")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)

        # Panel D: overall coverage distribution (per-sample overall)
        ax = axes[1, 1]
        ax.hist(all_test_cov, bins=np.linspace(0.80, 1.0, 21),
                color="#807dba", edgecolor="white")
        ax.axvline(TARGET, color="k", ls="--", lw=1.0,
                   label=f"target {TARGET:.2f}")
        ax.axvline(np.mean(all_test_cov), color=red, ls="-", lw=1.2,
                   label=f"mean={np.mean(all_test_cov):.3f}")
        ax.set_xlabel("Overall 90% coverage (per sample)")
        ax.set_ylabel("# samples")
        ax.set_title(f"D. Overall marginal coverage (n={len(all_test_cov)} samples)\n"
                     f"mean={summary['mean_overall_coverage']:.4f}, "
                     f"range=[{summary['min_overall_coverage']:.3f},"
                     f"{summary['max_overall_coverage']:.3f}]")
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(axis="y", alpha=0.3)

        fig.suptitle("Conditional conformal coverage (90% target) — does coverage "
                     "hold within subgroups, not just marginally?",
                     fontsize=12, y=1.00)
        fig.tight_layout()
        png = os.path.join(OUT_DIR, "conditional_coverage_figure.png")
        fig.savefig(png, dpi=140, bbox_inches="tight")
        both(f"Wrote {png}")
    except Exception as e:
        both(f"figure error: {e}\n{traceback.format_exc()}")

    # ---- console summary ----
    both("")
    both("=" * 70)
    both("CONDITIONAL COVERAGE SUMMARY")
    both("=" * 70)
    both(f"samples: {summary['n_samples']}  target_90: {TARGET}")
    both(f"overall marginal coverage: mean={summary['mean_overall_coverage']:.4f} "
         f"min={summary['min_overall_coverage']:.4f} "
         f"max={summary['max_overall_coverage']:.4f}")
    both(f"by uncertainty quintile : mean|dev|={summary['mean_abs_dev_by_quintile']:.4f} "
         f"max={summary['max_abs_dev_by_quintile']:.4f} "
         f"frac<=0.05={summary['frac_within_0.05_by_quintile']:.2f}")
    both(f"by spatial region       : mean|dev|={summary['mean_abs_dev_by_region']:.4f} "
         f"max={summary['max_abs_dev_by_region']:.4f} "
         f"frac<=0.05={summary['frac_within_0.05_by_region']:.2f}")
    both(f"by expression level     : mean|dev|={summary['mean_abs_dev_by_expression']:.4f} "
         f"max={summary['max_abs_dev_by_expression']:.4f} "
         f"frac<=0.05={summary['frac_within_0.05_by_expression']:.2f}")
    both(f"ALL subgroups combined  : mean|dev|={summary['mean_abs_dev_all_subgroups']:.4f} "
         f"max={summary['max_abs_dev_all_subgroups']:.4f} "
         f"frac<=0.05={summary['frac_within_0.05_all_subgroups']:.2f} "
         f"(n={summary['n_all_subgroups']})")
    both(f"worst subgroup: {summary['worst_subgroup']}")

    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

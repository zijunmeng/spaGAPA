#!/usr/bin/env python3
"""
AUDIT 4: Coverage + interval width under the new local_noise model.

EXPERT ATTACK (does r=0.55 come at a coverage / width cost?)
-------------------------------------------------------------
"You improved the raw GP uncertainty-error correlation from ~0.07 to ~0.55 by
switching from constant noise to local_noise.  But does that improvement
translate to BETTER calibrated, SHARPER conformal intervals -- or does it
break marginal coverage (via non-exchangeability of heteroscedastic noise),
or inflate interval width so much that the intervals become uninformative?"

PROTOCOL (per dataset x method)
-------------------------------
For each of the 5 datasets (test_uncertainty_noise_methods.py SAMPLES) and
each of three methods:

    A. constant noise (0.1)            -- baseline
    B. local_noise (per-gene kNN MAD)  -- the r=0.55 method
    D. residual_spot (2-pass EB)       -- the runner-up (r=0.51)

we replicate the conformal pipeline:
    load -> top-300 sites -> mask 20% per row -> SparseGP fit (per method) ->
    split held-out 50/50 into cal/test -> fit locally-adaptive split
    conformal at alpha=0.1 -> predict 90% intervals on the TEST set.

On the TEST set we then report:
    - marginal 90% coverage          (should be ~0.90 for all -- conformal)
    - mean interval WIDTH            (hi - lo; should NOT blow up for B/D)
    - interval score (Gneiting)      (coverage-weighted width; lower=better)
    - high-expr coverage             (coverage on the top-20% expression bin;
                                       Method A is ~0.81 here -- does B fix it?)
    - low-expr coverage              (bottom-20% bin; symmetry check)

VERDICT RULES (per the audit spec)
----------------------------------
    coverage_B ~= 0.90 AND width_B/A < 2  -> GOOD (sharper, calibrated)
    coverage_B ~= 0.90 AND width_B/A > 5  -> BAD  (intervals too wide)
    coverage_B != 0.90                    -> conformal broken by
                                             heteroscedastic non-exchangeability

This script reuses the GP / noise-estimator helpers from
scripts/test_uncertainty_noise_methods.py so the noise models are IDENTICAL
to the ones that produced r=0.55.

Outputs -> pipeline_output/uncertainty_coverage_width_audit/
    coverage_width_comparison.csv : dataset, method, marginal_coverage_90,
                                    mean_interval_width, interval_score,
                                    high_expr_coverage, low_expr_coverage,
                                    width_ratio_vs_A, n_test
    summary.json                  : per-method means, width ratios,
                                    high-expr improvement, overall verdict
    figure.png                    : (a) coverage by method
                                    (b) interval width by method
                                    (c) high-expr coverage by method
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

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

# Reuse the EXACT helpers (load, prep, noise estimators, method runners)
# from test_uncertainty_noise_methods.py so the noise models match the
# r=0.55 result bit-for-bit.
from scripts import test_uncertainty_noise_methods as tum  # noqa: E402

from spagapa.imputation.calibration import ConformalCalibrator  # noqa: E402

OUT_DIR = os.path.join(REPO, "pipeline_output",
                       "uncertainty_coverage_width_audit")

# Methods under audit (A baseline, B the r=0.55 winner, D the runner-up).
# C is omitted per the audit spec (which only asks for A/B/D), but the
# runner map in tum supports all four if needed.
AUDIT_METHODS = ["A_constant", "B_local_gene", "D_residual_spot"]

# Conformal hyperparameters (mirror conditional_coverage_validation.py).
ALPHA = 0.10
TARGET = 1 - ALPHA            # 0.90
CAL_FRACTION = 0.50
HIGH_EXPR_QUANTILE = 0.80     # top-20% expression bin
LOW_EXPR_QUANTILE = 0.20      # bottom-20% expression bin


def log(msg: str, fh) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def interval_score(y_true, lo, hi, alpha=ALPHA):
    """Gneiting & Raftery (2007) mean interval score (lower = better).

    IS = (hi - lo) + (2/alpha)*(lo - y)*(y < lo) + (2/alpha)*(y - hi)*(y > hi)
    """
    y_true = np.asarray(y_true, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    width = hi - lo
    penalty_lo = (2.0 / alpha) * (lo - y_true) * (y_true < lo)
    penalty_hi = (2.0 / alpha) * (y_true - hi) * (y_true > hi)
    return float(np.mean(width + penalty_lo + penalty_hi))


def coverage(y, lo, hi):
    y = np.asarray(y, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    if y.size == 0:
        return float("nan")
    return float(np.mean((y >= lo) & (y <= hi)))


def run_conformal_for_method(method, p, rng):
    """Fit GP with the chosen noise method, then run split conformal.

    Returns a dict of TEST-set arrays (truth, pred, std, lo, hi) plus
    calibration-set diagnostic corr.  Mirrors the conformal split in
    conditional_coverage_validation.fit_sample.

    IMPORTANT: conformal operates on the FLATTENED held-out arrays returned
    by tum.gather_held (1-D vectors of length n_held), NOT the full
    (n_genes, n_spots) matrices -- the cal/test split indexes into these
    flattened vectors, exactly as conditional_coverage_validation.fit_sample
    does (it indexes truths[test_i] etc. from the flattened held-out arrays).
    """
    pred_full, unc_full = tum.METHOD_RUNNERS[method](p, rng)
    truth, pp, u, err = tum.gather_held(p["values"], p["masked"],
                                        pred_full, unc_full)

    perm = rng.permutation(len(truth))
    n_cal = int(len(truth) * CAL_FRACTION)
    cal_i = perm[:n_cal]
    test_i = perm[n_cal:]

    cal = ConformalCalibrator(alpha=ALPHA, mode="locally_adaptive")
    cal.fit(err[cal_i], u[cal_i])
    # predict on the FULL flattened held-out set, then split cal/test
    lo_all, hi_all = cal.predict(pp, u)

    # cal-set unc-error corr (diagnostic -- relates to the r=0.55 finding)
    cal_corr = float(cal.interval_.calibration_unc_error_corr)

    return {
        "truth": truth, "pred": pp, "std": u, "err": err,
        "lo": lo_all, "hi": hi_all,
        "cal_i": cal_i, "test_i": test_i,
        "qhat90": float(cal.interval_.q_hat),
        "cal_unc_error_corr": cal_corr,
        "n_cal": int(n_cal), "n_test": int(len(test_i)),
    }


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    rng = np.random.default_rng(tum.SEED)

    log_path = os.path.join(OUT_DIR, "run.log")
    fh = open(log_path, "w")

    log(f"AUDIT 4: coverage + interval width under new noise model "
        f"(alpha={ALPHA}, target={TARGET})", fh)
    log(f"Methods under audit: {AUDIT_METHODS}", fh)
    log(f"Datasets: {[s[0] for s in tum.SAMPLES]}", fh)
    log("=" * 70, fh)

    rows = []
    per_dataset_test = {}   # for the figure

    for dataset, label, dir_basename in tum.SAMPLES:
        log("", fh)
        log(f">>> {dataset}/{label} ({dir_basename})", fh)
        try:
            t0 = time.time()
            p = tum.prep_sample(dir_basename, rng)
            log(f"  loaded: {p['n_genes']} genes x {p['n_spots']} spots, "
                f"obs_frac={(p['values']>0).mean():.3f}", fh)
        except Exception as e:
            log(f"  DATASET ERROR: {e}\n{traceback.format_exc()}", fh)
            continue

        # Pre-compute expression-bin edges ONCE on the full held-out truth
        # pool so high/low bins are comparable across methods on this dataset.
        # (We will re-derive truth-specific masks per method because each
        # method's held-out set is the SAME set of entries -- truth is
        # identical across methods -- so the bin edges are shared.)
        # Build the union held-out truth to define expression quantiles.
        truth_pool = p["values"][np.where((p["values"] > 0) & (p["masked"] <= 0))]
        if truth_pool.size == 0:
            log("  no held-out entries, SKIP", fh)
            continue
        hi_q = float(np.quantile(truth_pool, HIGH_EXPR_QUANTILE))
        lo_q = float(np.quantile(truth_pool, LOW_EXPR_QUANTILE))

        per_dataset_test[dataset] = {}

        width_A = None  # for ratio computation

        for method in AUDIT_METHODS:
            tm0 = time.time()
            try:
                # Each method gets its OWN fresh rng draw so the cal/test
                # split is identical across methods (deterministic given
                # SEED) -- this is what makes width comparisons fair.
                rng_method = np.random.default_rng(tum.SEED)
                r = run_conformal_for_method(method, p, rng_method)

                ti = r["test_i"]
                y = r["truth"][ti]
                lo_t = r["lo"][ti]
                hi_t = r["hi"][ti]
                std_t = r["std"][ti]

                marg_cov = coverage(y, lo_t, hi_t)
                mean_width = float(np.mean(hi_t - lo_t))
                iscore = interval_score(y, lo_t, hi_t)

                # expression-bin coverage (use the SHARED dataset-level edges)
                hi_mask = y >= hi_q
                lo_mask = y <= lo_q
                hi_cov = coverage(y[hi_mask], lo_t[hi_mask], hi_t[hi_mask])
                lo_cov = coverage(y[lo_mask], lo_t[lo_mask], hi_t[lo_mask])

                row = {
                    "dataset": dataset,
                    "label": label,
                    "method": method,
                    "marginal_coverage_90": marg_cov,
                    "mean_interval_width": mean_width,
                    "interval_score": iscore,
                    "high_expr_coverage": hi_cov,
                    "low_expr_coverage": lo_cov,
                    "n_test": int(len(ti)),
                    "qhat90": r["qhat90"],
                    "cal_unc_error_corr": r["cal_unc_error_corr"],
                    "high_expr_threshold": hi_q,
                    "low_expr_threshold": lo_q,
                    "n_high_expr": int(hi_mask.sum()),
                    "n_low_expr": int(lo_mask.sum()),
                    "wall_s": float(time.time() - tm0),
                }
                rows.append(row)
                per_dataset_test[dataset][method] = {
                    "marg_cov": marg_cov,
                    "width": mean_width,
                    "hi_cov": hi_cov,
                    "lo_cov": lo_cov,
                    "iscore": iscore,
                }
                if method == "A_constant":
                    width_A = mean_width
                log(f"  {method:18s} cov={marg_cov:.4f} "
                    f"width={mean_width:.4f} IS={iscore:.4f} "
                    f"hi_cov={hi_cov:.4f} lo_cov={lo_cov:.4f} "
                    f"qhat={r['qhat90']:.3f} "
                    f"cal_r={r['cal_unc_error_corr']:+.3f} "
                    f"[{time.time()-tm0:.0f}s]", fh)
            except Exception as e:
                log(f"  {method:18s} ERROR: {e}\n{traceback.format_exc()}", fh)

        # attach width_ratio_vs_A in a second pass
        if width_A is not None and width_A > 0:
            for row in rows:
                if row["dataset"] == dataset:
                    row["width_ratio_vs_A"] = (
                        float(row["mean_interval_width"]) / width_A)

        log(f"  dataset wall: {time.time()-t0:.0f}s", fh)

    if not rows:
        log("No rows produced; aborting.", fh)
        fh.close()
        return 1

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "coverage_width_comparison.csv")
    df.to_csv(csv_path, index=False)
    log(f"\nWrote {csv_path}", fh)

    # ---- Summary ----
    def col_mean(df_, m, col):
        sub = df_[(df_["method"] == m) & df_[col].notna()]
        return float(sub[col].mean()) if len(sub) else float("nan")

    summary = {
        "n_datasets": int(df["dataset"].nunique()),
        "datasets": sorted(df["dataset"].unique().tolist()),
        "methods": AUDIT_METHODS,
        "alpha": ALPHA,
        "target_90": TARGET,
        "per_method": {},
    }
    for m in AUDIT_METHODS:
        summary["per_method"][m] = {
            "mean_marginal_coverage": col_mean(df, m, "marginal_coverage_90"),
            "mean_interval_width": col_mean(df, m, "mean_interval_width"),
            "mean_interval_score": col_mean(df, m, "interval_score"),
            "mean_high_expr_coverage": col_mean(df, m, "high_expr_coverage"),
            "mean_low_expr_coverage": col_mean(df, m, "low_expr_coverage"),
            "mean_qhat90": col_mean(df, m, "qhat90"),
            "mean_cal_unc_error_corr": col_mean(df, m, "cal_unc_error_corr"),
        }
    # width ratios B/A and D/A (only where A is finite & > 0)
    ratios_B, ratios_D = [], []
    for ds in df["dataset"].unique():
        sub = df[df["dataset"] == ds].set_index("method")
        if ("A_constant" in sub.index and "B_local_gene" in sub.index
                and pd.notna(sub.loc["A_constant", "mean_interval_width"])
                and sub.loc["A_constant", "mean_interval_width"] > 0):
            ratios_B.append(sub.loc["B_local_gene", "mean_interval_width"]
                            / sub.loc["A_constant", "mean_interval_width"])
        if ("A_constant" in sub.index and "D_residual_spot" in sub.index
                and pd.notna(sub.loc["A_constant", "mean_interval_width"])
                and sub.loc["A_constant", "mean_interval_width"] > 0):
            ratios_D.append(sub.loc["D_residual_spot", "mean_interval_width"]
                            / sub.loc["A_constant", "mean_interval_width"])
    summary["width_ratio_B_over_A"] = {
        "mean": float(np.mean(ratios_B)) if ratios_B else float("nan"),
        "min": float(np.min(ratios_B)) if ratios_B else float("nan"),
        "max": float(np.max(ratios_B)) if ratios_B else float("nan"),
        "per_dataset": list(ratios_B),
    }
    summary["width_ratio_D_over_A"] = {
        "mean": float(np.mean(ratios_D)) if ratios_D else float("nan"),
        "per_dataset": list(ratios_D),
    }

    # High-expr improvement: B minus A, averaged.
    hi_A = summary["per_method"]["A_constant"]["mean_high_expr_coverage"]
    hi_B = summary["per_method"]["B_local_gene"]["mean_high_expr_coverage"]
    hi_D = summary["per_method"]["D_residual_spot"]["mean_high_expr_coverage"]
    summary["high_expr_improvement_B_minus_A"] = hi_B - hi_A
    summary["high_expr_improvement_D_minus_A"] = hi_D - hi_A

    # Does B maintain marginal coverage?  (all datasets within 0.05 of 0.90)
    b_covs = df[df["method"] == "B_local_gene"]["marginal_coverage_90"].dropna()
    summary["B_coverage_within_0.05_of_target_all_datasets"] = bool(
        np.all(np.abs(b_covs - TARGET) <= 0.05)) if len(b_covs) else False
    summary["B_coverage_min"] = float(b_covs.min()) if len(b_covs) else float("nan")
    summary["B_coverage_max"] = float(b_covs.max()) if len(b_covs) else float("nan")

    # VERDICT
    mean_ratio_B = summary["width_ratio_B_over_A"]["mean"]
    b_cov_mean = summary["per_method"]["B_local_gene"]["mean_marginal_coverage"]
    b_cov_ok = bool(abs(b_cov_mean - TARGET) <= 0.05)
    if not b_cov_ok:
        verdict = ("BAD: Method B marginal coverage "
                   f"({b_cov_mean:.3f}) deviates from target {TARGET} by >0.05 "
                   "-> heteroscedastic noise appears to break the conformal "
                   "exchangeability assumption.")
    elif np.isfinite(mean_ratio_B) and mean_ratio_B > 5.0:
        verdict = (f"BAD: Method B intervals are {mean_ratio_B:.1f}x wider than "
                   "Method A -> intervals too wide to be useful despite "
                   "correct coverage.")
    elif np.isfinite(mean_ratio_B) and mean_ratio_B < 2.0:
        verdict = ("GOOD: Method B maintains marginal coverage "
                   f"({b_cov_mean:.3f} ~= {TARGET}) AND width ratio B/A = "
                   f"{mean_ratio_B:.2f} (<2) -> sharper or comparable intervals "
                   "with correct coverage; the r=0.55 improvement IS useful.")
    else:
        verdict = ("MIXED: Method B maintains coverage "
                   f"({b_cov_mean:.3f}) but width ratio B/A = "
                   f"{mean_ratio_B:.2f} -- moderate inflation; acceptable if "
                   "high-expr coverage gain justifies it.")
    summary["verdict"] = verdict

    json_path = os.path.join(OUT_DIR, "summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    log(f"Wrote {json_path}", fh)

    # ---- Figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ds_order = [s[0] for s in tum.SAMPLES
                    if s[0] in per_dataset_test]
        method_colors = {
            "A_constant":      "#999999",
            "B_local_gene":    "#4C72B0",
            "D_residual_spot": "#C44E52",
        }
        fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

        # (a) marginal coverage by method (box over datasets)
        ax = axes[0]
        data_a = [df.loc[df["method"] == m,
                         "marginal_coverage_90"].dropna().values
                  for m in AUDIT_METHODS]
        pos = np.arange(len(AUDIT_METHODS))
        ax.boxplot(data_a, positions=pos, widths=0.5, patch_artist=True,
                   boxprops=dict(facecolor="#9ecae1", alpha=0.6))
        for i, m in enumerate(AUDIT_METHODS):
            vals = df.loc[df["method"] == m,
                          "marginal_coverage_90"].dropna().values
            ax.scatter(np.full_like(vals, i, dtype=float)
                       + np.random.default_rng(0).uniform(-0.08, 0.08, len(vals)),
                       vals, color=method_colors[m], s=30, alpha=0.9, zorder=3)
        ax.axhline(TARGET, color="k", ls="--", lw=1.0,
                   label=f"target {TARGET:.2f}")
        ax.set_xticks(pos)
        ax.set_xticklabels(AUDIT_METHODS, rotation=15)
        ax.set_ylim(0.70, 1.02)
        ax.set_ylabel("Marginal 90% coverage")
        ax.set_title("(a) Marginal coverage by method\n"
                     "(conformal should hold ~0.90 for all)")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)

        # (b) mean interval width by method (box over datasets)
        ax = axes[1]
        data_b = [df.loc[df["method"] == m,
                         "mean_interval_width"].dropna().values
                  for m in AUDIT_METHODS]
        ax.boxplot(data_b, positions=pos, widths=0.5, patch_artist=True,
                   boxprops=dict(facecolor="#a1d99b", alpha=0.6))
        for i, m in enumerate(AUDIT_METHODS):
            vals = df.loc[df["method"] == m,
                          "mean_interval_width"].dropna().values
            ax.scatter(np.full_like(vals, i, dtype=float)
                       + np.random.default_rng(1).uniform(-0.08, 0.08, len(vals)),
                       vals, color=method_colors[m], s=30, alpha=0.9, zorder=3)
        ax.set_xticks(pos)
        ax.set_xticklabels(AUDIT_METHODS, rotation=15)
        ax.set_ylabel("Mean interval width (hi - lo)")
        ax.set_title(f"(b) Interval width by method\n"
                     f"mean B/A = {summary['width_ratio_B_over_A']['mean']:.2f}")
        ax.grid(axis="y", alpha=0.3)

        # (c) high-expr coverage by method (box over datasets) + target line
        ax = axes[2]
        data_c = [df.loc[df["method"] == m,
                         "high_expr_coverage"].dropna().values
                  for m in AUDIT_METHODS]
        ax.boxplot(data_c, positions=pos, widths=0.5, patch_artist=True,
                   boxprops=dict(facecolor="#fdae6b", alpha=0.6))
        for i, m in enumerate(AUDIT_METHODS):
            vals = df.loc[df["method"] == m,
                          "high_expr_coverage"].dropna().values
            ax.scatter(np.full_like(vals, i, dtype=float)
                       + np.random.default_rng(2).uniform(-0.08, 0.08, len(vals)),
                       vals, color=method_colors[m], s=30, alpha=0.9, zorder=3)
        ax.axhline(TARGET, color="k", ls="--", lw=1.0,
                   label=f"target {TARGET:.2f}")
        ax.set_xticks(pos)
        ax.set_xticklabels(AUDIT_METHODS, rotation=15)
        ax.set_ylim(0.60, 1.02)
        ax.set_ylabel("Coverage on top-20% expression bin")
        ax.set_title(f"(c) High-expression coverage by method\n"
                     f"B-A improvement = "
                     f"{summary['high_expr_improvement_B_minus_A']:+.3f}")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)

        fig.suptitle("AUDIT 4: Coverage & interval width under the new "
                     "local_noise model\n"
                     "(does r=0.55 buy better calibration, or just wider "
                     "intervals?)", fontsize=12, y=1.04)
        fig.tight_layout()
        png = os.path.join(OUT_DIR, "figure.png")
        fig.savefig(png, dpi=140, bbox_inches="tight")
        log(f"Wrote {png}", fh)
    except Exception as e:
        log(f"figure error: {e}\n{traceback.format_exc()}", fh)

    # ---- Console summary ----
    log("", fh)
    log("=" * 70, fh)
    log("AUDIT 4 SUMMARY", fh)
    log("=" * 70, fh)
    log(f"datasets: {summary['n_datasets']}  target_90: {TARGET}", fh)
    for m in AUDIT_METHODS:
        pm = summary["per_method"][m]
        log(f"  {m:18s} cov={pm['mean_marginal_coverage']:.4f} "
            f"width={pm['mean_interval_width']:.4f} "
            f"IS={pm['mean_interval_score']:.4f} "
            f"hi_cov={pm['mean_high_expr_coverage']:.4f} "
            f"cal_r={pm['mean_cal_unc_error_corr']:+.4f}", fh)
    log(f"width ratio B/A: mean={summary['width_ratio_B_over_A']['mean']:.3f} "
        f"min={summary['width_ratio_B_over_A']['min']:.3f} "
        f"max={summary['width_ratio_B_over_A']['max']:.3f}", fh)
    log(f"high-expr coverage improvement B-A: "
        f"{summary['high_expr_improvement_B_minus_A']:+.4f}", fh)
    log(f"B coverage range across datasets: "
        f"[{summary['B_coverage_min']:.4f}, {summary['B_coverage_max']:.4f}]", fh)
    log(f"B within 0.05 of target on ALL datasets: "
        f"{summary['B_coverage_within_0.05_of_target_all_datasets']}", fh)
    log("", fh)
    log(f"VERDICT: {verdict}", fh)

    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

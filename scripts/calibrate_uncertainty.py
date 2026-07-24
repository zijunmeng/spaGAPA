#!/usr/bin/env python3
"""
Uncertainty calibration benchmark for spaGAPA's sparse GP.

Demonstrates that the raw sparse-GP posterior std is poorly calibrated
(unc-error corr ~0, 2-sigma coverage ~1.0 = far too conservative) and that
post-hoc **split conformal** calibration produces intervals with rigorous,
target marginal coverage at 80/90/95%.

Protocol
--------
1. Load GSE263789 binned APA matrix + coordinates.
2. Mask 20% of observed entries per gene.
3. Fit ``SparseGPImputerBatch`` on the masked matrix; predict held-out spots.
4. Split the held-out points into (a) a calibration set used to fit the
   conformal quantile and (b) an independent test set on which we report
   coverage / corr / RMSE.  Splitting is essential: reporting coverage on the
   same points used to fit ``q_hat`` would over-state it.
5. Report BEFORE vs AFTER conformal:
     - uncertainty-error Pearson + Spearman correlation
     - interval coverage at 80% / 90% / 95%
     - RMSE on the full test set and after filtering the highest-uncertainty
       spots (a downstream use-case: drop the 20% most-uncertain imputations).
6. Both ``global`` and ``locally_adaptive`` conformal modes are reported.

Usage
-----
python scripts/calibrate_uncertainty.py \
    --apa-matrix  pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200/apa_matrix.csv \
    --coordinates pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200/coordinates.csv \
    --output      pipeline_output/uncertainty_calibration \
    [--mask-fraction 0.2] [--cal-fraction 0.5] [--seed 42]
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
from scipy.stats import pearsonr, spearmanr


def _safe_corr(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if np.std(a) > 0 and np.std(b) > 0:
        r, _ = pearsonr(a, b)
        rho, _ = spearmanr(a, b)
        return float(r), float(rho)
    return float("nan"), float("nan")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apa-matrix", required=True)
    parser.add_argument("--coordinates", required=True)
    parser.add_argument(
        "--output", default="pipeline_output/uncertainty_calibration",
        help="Output directory (JSON + report).")
    parser.add_argument("--mask-fraction", type=float, default=0.2,
                        help="Fraction of observed entries masked per gene.")
    parser.add_argument("--cal-fraction", type=float, default=0.5,
                        help="Fraction of held-out points used for conformal "
                             "calibration (rest = test).")
    parser.add_argument("--filter-frac", type=float, default=0.2,
                        help="Fraction of highest-uncertainty spots to drop "
                             "when reporting RMSE-after-filtering.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-inducing", type=int, default=None,
                        help="Override inducing-point count (default: "
                             "min(500, max(100, n_spots//100))).")
    args = parser.parse_args()

    if not (0.0 < args.mask_fraction < 1.0):
        parser.error("--mask-fraction must be in (0, 1)")
    if not (0.0 < args.cal_fraction < 1.0):
        parser.error("--cal-fraction must be in (0, 1)")

    import pandas as pd
    from spagapa.imputation import SparseGPImputer
    from spagapa.imputation.calibration import (
        ConformalCalibrator,
        evaluate_coverage,
    )
    from scipy.spatial import cKDTree

    t_total = time.time()

    # ------------------------------------------------------------- load data
    t0 = time.time()
    apa = pd.read_csv(args.apa_matrix, index_col=0)
    coords = pd.read_csv(args.coordinates, index_col=0)
    values = apa.values  # (n_genes, n_spots)
    if "x" in coords.columns and "y" in coords.columns:
        xy = coords[["x", "y"]].values
    else:
        xy = coords.iloc[:, :2].values.astype(float)
    # Best-effort spot alignment
    if values.shape[1] != xy.shape[0]:
        common = apa.columns.intersection(coords.index)
        if len(common) > 0:
            apa = apa.loc[:, common]
            coords = coords.loc[common]
            xy = coords[["x", "y"]].values if {"x", "y"} <= set(coords.columns) \
                else coords.iloc[:, :2].values.astype(float)
            values = apa.values
    if values.shape[1] != xy.shape[0]:
        raise ValueError(
            f"Spot count mismatch: apa={values.shape[1]}, coords={xy.shape[0]}"
        )
    n_genes, n_spots = values.shape
    print(f"[load] {n_genes} genes x {n_spots} spots, "
          f"observed frac={(values > 0).mean():.3f} "
          f"({time.time()-t0:.1f}s)")

    # --------------------------------------------------------------- masking
    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    observed = values > 0
    masked_values = values.copy()
    held_out: list[tuple[int, np.ndarray, np.ndarray]] = []  # (gene, idx, truth)
    for g in range(n_genes):
        obs_idx = np.where(observed[g])[0]
        if len(obs_idx) < 5:
            continue
        n_mask = max(1, int(len(obs_idx) * args.mask_fraction))
        mask_idx = rng.choice(obs_idx, size=n_mask, replace=False)
        masked_values[g, mask_idx] = 0.0
        held_out.append((g, mask_idx, values[g, mask_idx]))
    n_held = sum(len(idx) for _, idx, _ in held_out)
    print(f"[mask] held out {n_held} entries across {len(held_out)} genes "
          f"(frac={args.mask_fraction}) ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------- fit sparse GP
    t0 = time.time()
    _kdt = cKDTree(xy)
    _nn, _ = _kdt.query(xy, k=2)
    nn_dist = float(np.median(_nn[:, 1]))
    length_scale = nn_dist * 5
    n_inducing = (args.n_inducing if args.n_inducing is not None
                  else min(500, max(100, n_spots // 100)))
    print(f"[gp] length_scale={length_scale:.1f} (nn={nn_dist:.1f}), "
          f"n_inducing={n_inducing}")
    base = SparseGPImputer(
        n_inducing=n_inducing, length_scale=length_scale, noise_level=0.1,
    )
    train_mask = masked_values > 0
    batch = base.fit_batch(xy, masked_values, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    gp_fit_time = time.time() - t0
    print(f"[gp] fit+impute done ({gp_fit_time:.1f}s)")

    # ---------------------------------- gather held-out pred / std / truth
    truths_all, preds_all, unc_all = [], [], []
    for g, idx, truth in held_out:
        truths_all.append(truth)
        preds_all.append(gp_pred[g, idx])
        unc_all.append(gp_unc[g, idx])
    truths_all = np.concatenate(truths_all)
    preds_all = np.concatenate(preds_all)
    unc_all = np.concatenate(unc_all)
    errors_all = np.abs(truths_all - preds_all)
    print(f"[gather] {len(truths_all)} held-out points")

    # --------------------------------------- split calibration / test (50/50)
    perm = rng.permutation(len(truths_all))
    n_cal = int(len(truths_all) * args.cal_fraction)
    cal_idx = perm[:n_cal]
    test_idx = perm[n_cal:]
    cal_err = errors_all[cal_idx]
    cal_std = unc_all[cal_idx]
    test_truth = truths_all[test_idx]
    test_pred = preds_all[test_idx]
    test_std = unc_all[test_idx]
    test_err = errors_all[test_idx]

    results: dict = {
        "config": {
            "apa_matrix": os.path.abspath(args.apa_matrix),
            "coordinates": os.path.abspath(args.coordinates),
            "n_genes": int(n_genes),
            "n_spots": int(n_spots),
            "mask_fraction": float(args.mask_fraction),
            "cal_fraction": float(args.cal_fraction),
            "n_held_out": int(len(truths_all)),
            "n_calibration": int(n_cal),
            "n_test": int(len(test_idx)),
            "filter_frac": float(args.filter_frac),
            "seed": int(args.seed),
            "gp_length_scale": float(length_scale),
            "gp_n_inducing": int(n_inducing),
            "gp_nn_distance": float(nn_dist),
        },
        "timings": {
            "gp_fit_impute_s": float(gp_fit_time),
        },
    }

    # ============================================================ BEFORE
    # Raw GP intervals: [pred - z*std, pred + z*std] for z = 1.282/1.645/1.96
    z_levels = {0.80: 1.2816, 0.90: 1.6449, 0.95: 1.96}
    before = {
        "unc_error_pearson": _safe_corr(unc_all, errors_all)[0],
        "unc_error_spearman": _safe_corr(unc_all, errors_all)[1],
        "coverage": {},
    }
    for level, z in z_levels.items():
        half = z * test_std
        cov = evaluate_coverage(test_pred - half, test_pred + half, test_truth)
        before["coverage"][f"{int(level*100)}pct"] = cov
    # RMSE before (full test set). RMSE-after-filtering is computed below
    # once the raw-std interval half-width is available for fair comparison.
    before["rmse"] = float(np.sqrt(np.mean(test_err ** 2)))
    print("\n=== BEFORE conformal (raw GP std) ===")
    print(f"  unc-error Pearson : {before['unc_error_pearson']:.4f}")
    print(f"  unc-error Spearman: {before['unc_error_spearman']:.4f}")
    for level in z_levels:
        print(f"  {int(level*100)}% coverage     : "
              f"{before['coverage'][f'{int(level*100)}pct']:.4f} "
              f"(target {level:.2f})")

    # ============================================================ AFTER
    after: dict = {"modes": {}}
    for mode in ("global", "locally_adaptive"):
        mode_res = {"coverage": {}, "q_hat": {}, "calibration_unc_error_corr": None}
        for alpha, level in [(0.20, 80), (0.10, 90), (0.05, 95)]:
            cal = ConformalCalibrator(alpha=alpha, mode=mode)
            cal.fit(cal_err, cal_std if mode == "locally_adaptive" else None)
            lo, hi = cal.predict(test_pred, test_std if mode == "locally_adaptive" else None)
            cov = evaluate_coverage(lo, hi, test_truth)
            mode_res["coverage"][f"{level}pct"] = cov
            mode_res["q_hat"][f"{level}pct"] = cal.interval_.q_hat
            if alpha == 0.10:
                mode_res["calibration_unc_error_corr"] = (
                    cal.interval_.calibration_unc_error_corr
                )
                # RMSE-after-filtering: drop top `filter_frac` by interval
                # half-width (a calibrated proxy for uncertainty). For
                # locally-adaptive this is q_hat*std; for global it is the
                # constant q_hat (filtering is then random -> no improvement,
                # which is itself a useful diagnostic).
                half = (hi - lo) / 2.0
                mode_res["rmse"] = float(np.sqrt(np.mean(test_err ** 2)))
                mode_res["rmse_after_filter_top_unc"] = _rmse_filter(
                    test_err, half, args.filter_frac)
        after["modes"][mode] = mode_res

        print(f"\n=== AFTER conformal ({mode}) ===")
        if mode_res["calibration_unc_error_corr"] is not None:
            c = mode_res["calibration_unc_error_corr"]
            print(f"  calibration-set unc-error corr: "
                  f"{c:.4f}" if not np.isnan(c) else
                  "  calibration-set unc-error corr: NaN")
        for level in (80, 90, 95):
            print(f"  {level}% coverage: "
                  f"{mode_res['coverage'][f'{level}pct']:.4f} "
                  f"(target {level/100:.2f}, "
                  f"q_hat={mode_res['q_hat'][f'{level}pct']:.4f})")
        print(f"  RMSE (all test)            : {mode_res['rmse']:.4f}")
        print(f"  RMSE after dropping top {args.filter_frac:.0%} unc: "
              f"{mode_res['rmse_after_filter_top_unc']:.4f}")

    # Report BEFORE RMSE-after-filter too for fair comparison (filter by raw std)
    raw_half_90 = z_levels[0.90] * test_std
    before["rmse"] = float(np.sqrt(np.mean(test_err ** 2)))
    before["rmse_after_filter_top_unc"] = _rmse_filter(
        test_err, raw_half_90, args.filter_frac)

    results["before"] = before
    results["after"] = after
    results["timings"]["total_s"] = float(time.time() - t_total)

    # ----------------------------------------------------------------- save
    os.makedirs(args.output, exist_ok=True)
    out_path = os.path.join(args.output, "uncertainty_calibration.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    report_path = os.path.join(args.output, "uncertainty_calibration_report.txt")
    with open(report_path, "w") as f:
        f.write(_format_report(results))
    print(f"\n[done] JSON  -> {out_path}")
    print(f"       report -> {report_path}")
    print(f"       total wall time: {results['timings']['total_s']:.1f}s")


def _rmse_filter(err: np.ndarray, half_width: np.ndarray, frac: float) -> float:
    """RMSE after dropping the top `frac` points by `half_width` (uncertainty).

    Returns RMSE of the *remaining* (1-frac) points.  A calibrated uncertainty
    should make the dropped set have higher error, so this RMSE should drop
    relative to the unfiltered RMSE.  If uncertainty is uninformative (e.g.
    global conformal where half-width is constant), this ~ equals the
    unfiltered RMSE.
    """
    err = np.asarray(err, dtype=float)
    half_width = np.asarray(half_width, dtype=float)
    n = len(err)
    n_keep = int(np.ceil(n * (1.0 - frac)))
    order = np.argsort(half_width)  # ascending uncertainty
    keep = order[:n_keep]
    return float(np.sqrt(np.mean(err[keep] ** 2)))


def _format_report(r: dict) -> str:
    cfg = r["config"]
    bef = r["before"]
    aft = r["after"]
    lines = []
    lines.append("=" * 70)
    lines.append("spaGAPA Uncertainty Calibration (split conformal)")
    lines.append("=" * 70)
    lines.append(
        f"Data: {cfg['n_genes']} genes x {cfg['n_spots']} spots "
        f"(observed frac via mask frac={cfg['mask_fraction']})"
    )
    lines.append(
        f"Held-out: {cfg['n_held_out']}  "
        f"calibration={cfg['n_calibration']}  test={cfg['n_test']}"
    )
    lines.append(
        f"GP: length_scale={cfg['gp_length_scale']:.1f}, "
        f"n_inducing={cfg['gp_n_inducing']}"
    )
    lines.append("")
    lines.append("BEFORE conformal (raw GP std):")
    lines.append(f"  unc-error Pearson : {bef['unc_error_pearson']:.4f}")
    lines.append(f"  unc-error Spearman: {bef['unc_error_spearman']:.4f}")
    for level in (80, 90, 95):
        lines.append(
            f"  {level}% coverage: {bef['coverage'][f'{level}pct']:.4f}"
        )
    lines.append(f"  RMSE              : {bef['rmse']:.4f}")
    lines.append(f"  RMSE after dropping top {cfg['filter_frac']:.0%} unc (raw std): "
                 f"{bef['rmse_after_filter_top_unc']:.4f}")
    lines.append("")
    for mode in ("global", "locally_adaptive"):
        m = aft["modes"][mode]
        lines.append(f"AFTER conformal ({mode}):")
        c = m["calibration_unc_error_corr"]
        cs = "NaN" if (c is None or np.isnan(c)) else f"{c:.4f}"
        lines.append(f"  calibration-set unc-error corr: {cs}")
        for level in (80, 90, 95):
            lines.append(
                f"  {level}% coverage: {m['coverage'][f'{level}pct']:.4f} "
                f"(target {level/100:.2f}, q_hat={m['q_hat'][f'{level}pct']:.4f})"
            )
        lines.append(f"  RMSE              : {m['rmse']:.4f}")
        lines.append(
            f"  RMSE after dropping top {cfg['filter_frac']:.0%} unc: "
            f"{m['rmse_after_filter_top_unc']:.4f}"
        )
        lines.append("")
    lines.append(f"Total wall time: {r['timings']['total_s']:.1f}s")
    lines.append(f"GP fit+impute  : {r['timings']['gp_fit_impute_s']:.1f}s")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
AUDIT 2: Leave-One-Dataset-Out Cross-Validation (LOOCV) of noise-method selection.

EXPERT ATTACK (selection bias)
------------------------------
"Methods B-D were selected on the SAME 5 datasets used for evaluation -> the
reported r=0.554 for Method B is optimistic; with a different / held-out
dataset B might not have won."

This script tests that concern directly.  For each of the 5 datasets used in
test_uncertainty_noise_methods.py, we:

  1. Compute the mean Pearson r of Methods A/B/C/D on the OTHER 4 datasets
     (the *training fold*).
  2. Select the method with the highest training-fold mean r.
  3. Evaluate the SELECTED method's r on the HELD-OUT dataset (the *test fold*).

We then report:
  - selection frequency of each method (who would be chosen across folds?)
  - held-out r of the method selected on the training fold (the honest,
    no-peeking estimate)
  - held-out r of EVERY method (so we can see whether B would have won had we
    just always picked it)

VERDICT RULES
-------------
  If Method B wins on held-out data in 4/5 folds  -> selection is ROBUST.
  If Method B wins on 2/5 folds                   -> selection is FRAGILE,
                                                      qualify in the paper.

Important design choices
------------------------
* REUSES the per-(dataset, method) Pearson r already computed by
  scripts/test_uncertainty_noise_methods.py (pipeline_output/
  uncertainty_corr_improvement/per_dataset_corr.csv).  No GP refit needed --
  the r values are deterministic functions of (data, SEED=42, method) and we
  are only asking a QUESTION ABOUT SELECTION, not re-estimating r.  This keeps
  the audit cheap and exactly reproducible.
* Selection tie-break: alphabetical on method name (only triggers if two
  methods have identical mean r -- effectively never).
* We also report an "oracle" row: the held-out r of the best method on the
  held-out fold itself (the upper bound a peeking selector would achieve).

Outputs -> pipeline_output/uncertainty_loocv/
    loocv_results.csv : fold, held_out_dataset, method_selected_on_train,
                        held_out_r_<method> x4, selected_method_held_out_r,
                        oracle_method, oracle_held_out_r
    per_fold_train_means.csv : fold, method, train_fold_mean_r,
                        train_fold_median_r
    summary.json : selection frequency, mean held-out r per method,
                    mean selected-method held-out r, mean oracle held-out r,
                    robustness verdict
    run.log
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT_DIR = os.path.join(REPO, "pipeline_output", "uncertainty_loocv")
SOURCE_CSV = os.path.join(REPO, "pipeline_output", "uncertainty_corr_improvement",
                          "per_dataset_corr.csv")

# Same 5 datasets as test_uncertainty_noise_methods.py.
DATASETS = ["GSE183456", "GSE220442", "GSE169749", "GSE338525", "GSE237183"]
METHODS = ["A_constant", "B_local_gene", "C_spatial_spot", "D_residual_spot"]


def log(msg: str, fh) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = os.path.join(OUT_DIR, "run.log")
    fh = open(log_path, "w")

    if not os.path.exists(SOURCE_CSV):
        log(f"FATAL: source CSV not found: {SOURCE_CSV}", fh)
        log("Run scripts/test_uncertainty_noise_methods.py first.", fh)
        fh.close()
        return 1

    df = pd.read_csv(SOURCE_CSV)
    # Restrict to the 5 LOOCV datasets (defensive -- the source CSV already
    # contains exactly these, but filtering makes the script robust).
    df = df[df["dataset"].isin(DATASETS)].copy()
    # Drop any NaN corr rows (failed runs).
    df = df.dropna(subset=["unc_error_pearson"]).copy()

    present = set(df["dataset"].unique())
    missing = set(DATASETS) - present
    if missing:
        log(f"WARNING: datasets missing from source CSV: {sorted(missing)}", fh)

    log(f"LOOCV over {len(DATASETS)} datasets, methods={METHODS}", fh)
    log(f"Source: {SOURCE_CSV}", fh)
    log("=" * 70, fh)

    loocv_rows = []
    train_means_rows = []
    selection_counts = {m: 0 for m in METHODS}
    held_out_r_per_method = {m: [] for m in METHODS}
    selected_held_out_r = []
    oracle_held_out_r = []

    for held in DATASETS:
        train_set = [d for d in DATASETS if d != held]
        # Training-fold mean r per method (over the 4 held-in datasets).
        train_df = df[df["dataset"].isin(train_set)]
        train_means = {}
        train_medians = {}
        for m in METHODS:
            sub = train_df[train_df["method"] == m]["unc_error_pearson"]
            if len(sub) == 0:
                train_means[m] = float("nan")
                train_medians[m] = float("nan")
            else:
                train_means[m] = float(sub.mean())
                train_medians[m] = float(sub.median())
            train_means_rows.append({
                "fold_held_out": held,
                "method": m,
                "train_fold_mean_r": train_means[m],
                "train_fold_median_r": train_medians[m],
                "n_train_datasets": int(len(sub)),
            })

        # Select best method on training fold (highest mean r).
        finite = {m: v for m, v in train_means.items()
                  if np.isfinite(v)}
        if not finite:
            log(f"  fold held={held}: no finite train means, SKIP", fh)
            continue
        selected = max(finite, key=lambda m: (finite[m], ))  # value then name
        selection_counts[selected] += 1

        # Held-out r per method on the test fold (the held-out dataset).
        held_df = df[df["dataset"] == held]
        held_r = {}
        for m in METHODS:
            sub = held_df[held_df["method"] == m]["unc_error_pearson"]
            if len(sub) == 0:
                held_r[m] = float("nan")
            else:
                held_r[m] = float(sub.iloc[0])
                held_out_r_per_method[m].append(held_r[m])

        sel_held_r = held_r[selected]
        selected_held_out_r.append(sel_held_r)

        # Oracle: the method that WOULD have been best on the held-out fold
        # itself (upper bound -- what a peeking selector would get).
        finite_held = {m: v for m, v in held_r.items() if np.isfinite(v)}
        oracle_m = (max(finite_held, key=lambda m: finite_held[m])
                    if finite_held else None)
        oracle_r = finite_held.get(oracle_m, float("nan")) if oracle_m else float("nan")
        oracle_held_out_r.append(oracle_r)

        row = {
            "fold": len(loocv_rows) + 1,
            "held_out_dataset": held,
            "n_train_datasets": len(train_set),
            "method_selected_on_train": selected,
            "selected_method_train_mean_r": train_means[selected],
        }
        for m in METHODS:
            row[f"held_out_r_{m}"] = held_r[m]
        row["selected_method_held_out_r"] = sel_held_r
        row["oracle_method"] = oracle_m
        row["oracle_held_out_r"] = oracle_r
        loocv_rows.append(row)

        # Pretty per-fold log.
        train_str = " ".join(f"{m}={train_means[m]:+.3f}" for m in METHODS)
        held_str = " ".join(f"{m}={held_r[m]:+.3f}" for m in METHODS)
        log(f"  fold held={held} ({len(train_set)} train datasets)", fh)
        log(f"    train mean r: {train_str}", fh)
        log(f"    SELECTED on train: {selected}  "
            f"(train mean={train_means[selected]:+.4f})", fh)
        log(f"    held-out r:    {held_str}", fh)
        log(f"    selected held-out r = {sel_held_r:+.4f}   "
            f"oracle={oracle_m}({oracle_r:+.4f})", fh)
        log("", fh)

    loocv_df = pd.DataFrame(loocv_rows)
    train_df_out = pd.DataFrame(train_means_rows)

    csv1 = os.path.join(OUT_DIR, "loocv_results.csv")
    csv2 = os.path.join(OUT_DIR, "per_fold_train_means.csv")
    loocv_df.to_csv(csv1, index=False)
    train_df_out.to_csv(csv2, index=False)
    log(f"Wrote {csv1}", fh)
    log(f"Wrote {csv2}", fh)

    # ---- Summary ----
    # B-wins-on-held-out frequency (the key robustness metric).
    b_wins_held = 0
    for r in loocv_rows:
        held_r = {m: r[f"held_out_r_{m}"] for m in METHODS}
        finite_held = {m: v for m, v in held_r.items() if np.isfinite(v)}
        if finite_held:
            best = max(finite_held, key=lambda m: finite_held[m])
            if best == "B_local_gene":
                b_wins_held += 1

    mean_held_per_method = {
        m: (float(np.mean(v)) if v else float("nan"))
        for m, v in held_out_r_per_method.items()
    }
    mean_selected = float(np.mean(selected_held_out_r)) if selected_held_out_r else float("nan")
    mean_oracle = float(np.mean(oracle_held_out_r)) if oracle_held_out_r else float("nan")

    # Always-B strategy: pick B regardless of training fold -> its mean held-out r.
    always_b = mean_held_per_method.get("B_local_gene", float("nan"))

    # Regret of the LOOCV-selected method vs always-B vs oracle.
    regret_vs_oracle = mean_oracle - mean_selected
    regret_always_b_vs_oracle = mean_oracle - always_b

    # Verdict.
    n_folds = len(loocv_rows)
    b_win_frac = b_wins_held / n_folds if n_folds else 0.0
    if b_win_frac >= 0.8:
        verdict = ("ROBUST: Method B (local_noise) wins on held-out data in "
                   f"{b_wins_held}/{n_folds} folds -> selection is NOT an "
                   "artifact of in-sample optimisation.")
    elif b_win_frac >= 0.4:
        verdict = ("PARTIALLY ROBUST: Method B wins on held-out data in "
                   f"{b_wins_held}/{n_folds} folds -> selection is mostly "
                   "stable but should be qualified.")
    else:
        verdict = ("FRAGILE: Method B wins on held-out data in only "
                   f"{b_wins_held}/{n_folds} folds -> the r=0.554 is "
                   "in-sample optimistic; reconsider method choice.")

    summary = {
        "n_folds": n_folds,
        "datasets": DATASETS,
        "methods": METHODS,
        "selection_frequency_on_train_fold": selection_counts,
        "b_wins_on_held_out_count": b_wins_held,
        "b_wins_on_held_out_fraction": b_win_frac,
        "mean_held_out_r_per_method": mean_held_per_method,
        "mean_loocv_selected_held_out_r": mean_selected,
        "mean_always_B_held_out_r": always_b,
        "mean_oracle_held_out_r": mean_oracle,
        "regret_selected_vs_oracle": regret_vs_oracle,
        "regret_always_B_vs_oracle": regret_always_b_vs_oracle,
        "selected_method_per_fold": [
            {"fold": r["fold"], "held_out": r["held_out_dataset"],
             "selected": r["method_selected_on_train"],
             "selected_held_out_r": r["selected_method_held_out_r"],
             "oracle_method": r["oracle_method"],
             "oracle_held_out_r": r["oracle_held_out_r"]}
            for r in loocv_rows
        ],
        "verdict": verdict,
    }
    json_path = os.path.join(OUT_DIR, "summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    log(f"Wrote {json_path}", fh)

    # ---- Console summary ----
    log("", fh)
    log("=" * 70, fh)
    log("LOOCV SUMMARY", fh)
    log("=" * 70, fh)
    log(f"folds: {n_folds}", fh)
    log("selection frequency (method chosen on training fold):", fh)
    for m in METHODS:
        log(f"  {m:18s} {selection_counts[m]}/{n_folds}", fh)
    log("mean held-out r per method (averaged over the 5 folds):", fh)
    for m in METHODS:
        log(f"  {m:18s} {mean_held_per_method[m]:+.4f}", fh)
    log(f"B wins on held-out in {b_wins_held}/{n_folds} folds "
        f"(fraction={b_win_frac:.2f})", fh)
    log(f"mean r of LOOCV-selected method on held-out : {mean_selected:+.4f}", fh)
    log(f"mean r of always-B strategy on held-out    : {always_b:+.4f}", fh)
    log(f"mean r of oracle (peeking) on held-out     : {mean_oracle:+.4f}", fh)
    log(f"regret (selected vs oracle)   : {regret_vs_oracle:+.4f}", fh)
    log(f"regret (always-B vs oracle)   : {regret_always_b_vs_oracle:+.4f}", fh)
    log("", fh)
    log(f"VERDICT: {verdict}", fh)

    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

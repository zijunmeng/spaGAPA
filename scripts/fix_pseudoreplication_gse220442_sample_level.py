#!/usr/bin/env python3
"""Fix pseudoreplication in GSE220442 differential APA: SAMPLE-LEVEL test.

PROBLEM
-------
The existing analysis (`pipeline_output/gse220442_differential_apa_unified/`)
pools all ~29000 spots from 6 samples and treats EACH SPOT as an independent
observation in a t-test (n_control=14147, n_AD=14955). With biological
replicates = 3 control + 3 AD INDIVIDUALS, this is pseudoreplication: the
effective n is 6, not ~29000, and the per-spot t-test massively inflates the
degrees of freedom -> artificially tiny p-values (91 "significant" genes).

CORRECT APPROACH
----------------
The biological replicate unit is the SAMPLE (individual brain section), not the
spot. We therefore:

  1. For each gene in each sample, compute the SAMPLE-LEVEL mean distal-usage
     (aggregate across all that sample's spots). These per-sample means are
     already stored in `results.csv` columns `control_vals` and `ad_vals`
     (semicolon-joined means for the 3 control / 3 AD samples), produced by the
     unified analysis' `sample_level()` helper.
  2. Run a proper two-sample t-test on n=3 control means vs n=3 AD means, per
     gene. (n=3 per arm -> 4 degrees of freedom; appropriately low power.)
  3. BH-FDR correction across all tested genes.
  4. Compare to the spot-level call set: how many spot-sig genes survive at the
     sample level? (overlap = "robust core").
  5. Direction agreement: do all 3 AD samples shift the same way relative to
     the 3 control samples? (min(AD) > max(control) or vice versa.)

OUTPUT
------
pipeline_output/gse220442_sample_level_differential/
  sample_level_results.csv        per-gene: ctrl/ad sample means, sample_t,
                                  sample_p, sample_padj, direction_agree
  spot_vs_sample_comparison.csv   per-gene: spot_padj, sample_padj,
                                  retained_at_sample_level
  summary.json                    n_spot_sig, n_sample_sig, overlap,
                                  direction-agreement rate, top genes

Usage:
  ~/anaconda3/envs/spagapa/bin/python \
      scripts/fix_pseudoreplication_gse220442_sample_level.py
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PKG = Path(__file__).resolve().parents[1]
SPOT_RESULTS = (PKG / "pipeline_output/gse220442_differential_apa_unified"
                / "results.csv")
OUT = PKG / "pipeline_output/gse220442_sample_level_differential"
OUT.mkdir(parents=True, exist_ok=True)

# sample order matches the unified analysis: control = GSM6801751/2/3,
# AD = GSM6801754/5/6
CONTROL_SAMPLES = ["gsm6801751", "gsm6801752", "gsm6801753"]
AD_SAMPLES = ["gsm6801754", "gsm6801755", "gsm6801756"]

PADJ_THRESH = 0.05
DELTA_THRESH = 0.05  # effect-size threshold (mirrors spot-level analysis)


def parse_sample_means(results: pd.DataFrame):
    """Extract per-sample mean distal-usage for every gene.

    results.csv columns:
      control_vals  "0.686;0.673;0.748"   <- ctrl sample means (3)
      ad_vals       "0.676;0.658;0.706"   <- AD   sample means (3)
    Returns a DataFrame indexed by gene with columns:
      ctrl_mean_1, ctrl_mean_2, ctrl_mean_3,
      ad_mean_1,   ad_mean_2,   ad_mean_3
    """
    def _split(col):
        arr = results[col].fillna("").astype(str).str.split(";")
        out = np.full((len(arr), 3), np.nan, dtype=float)
        for i, parts in enumerate(arr):
            vals = []
            for p in parts:
                p = p.strip()
                if p == "" or p.lower() == "nan":
                    vals.append(np.nan)
                else:
                    try:
                        vals.append(float(p))
                    except ValueError:
                        vals.append(np.nan)
            for j in range(min(3, len(vals))):
                out[i, j] = vals[j]
        return out

    ctrl = _split("control_vals")
    ad = _split("ad_vals")
    out = pd.DataFrame({
        "ctrl_mean_1": ctrl[:, 0], "ctrl_mean_2": ctrl[:, 1],
        "ctrl_mean_3": ctrl[:, 2],
        "ad_mean_1": ad[:, 0], "ad_mean_2": ad[:, 1], "ad_mean_3": ad[:, 2],
    }, index=results["gene"].values)
    return out


def main():
    t0 = time.time()
    print("=== GSE220442 sample-level (anti-pseudoreplication) "
          "differential APA ===")
    print(f"  reading spot-level results: {SPOT_RESULTS}")
    spot = pd.read_csv(SPOT_RESULTS)
    print(f"  genes in spot-level table: {len(spot)}")

    # per-sample mean distal-usage per gene
    means = parse_sample_means(spot)
    # a gene is testable at the sample level only if all 6 sample means are
    # finite (NaN -> gene not observed in that sample's spots)
    testable = means.dropna()
    print(f"  genes with all 6 sample means finite (testable): {len(testable)}")

    ctrl_cols = ["ctrl_mean_1", "ctrl_mean_2", "ctrl_mean_3"]
    ad_cols = ["ad_mean_1", "ad_mean_2", "ad_mean_3"]
    ctrl = testable[ctrl_cols].values
    ad = testable[ad_cols].values

    # ---- per-gene sample-level t-test (n=3 vs n=3, Welch's = default of scipy
    #      ttest_ind with equal_var=False; report both Welch and pooled) ----
    # Welch (does not assume equal variances) is the more honest default at n=3.
    # IMPORTANT: many genes have near-constant distal-usage (e.g. all 1.0) ->
    # zero within-group variance -> scipy returns NaN t/p ("catastrophic
    # cancellation"). Such genes carry NO differential signal; assign p=1.0
    # (and a degenerate t=0) so BH-FDR treats them as non-significant rather
    # than propagating NaN.
    t_welch = np.zeros(len(testable))
    p_welch = np.ones(len(testable))
    t_pool = np.zeros(len(testable))
    p_pool = np.ones(len(testable))
    for i in range(len(testable)):
        c = ctrl[i]
        a = ad[i]
        # both groups constant -> no testable difference
        c_var = np.var(c, ddof=1)
        a_var = np.var(a, ddof=1)
        if c_var == 0 and a_var == 0:
            if np.mean(a) == np.mean(c):
                continue  # identical -> p stays 1.0, t stays 0
            else:
                # constant but different means -> infinitely significant.
                # At n=3 this is degenerate; report p=0 but flag.
                p_welch[i] = 0.0
                p_pool[i] = 0.0
                t_welch[i] = np.inf * np.sign(np.mean(a) - np.mean(c))
                t_pool[i] = t_welch[i]
                continue
        with np.errstate(invalid="ignore"):
            tw, pw = stats.ttest_ind(a, c, equal_var=False)
            tp, pp = stats.ttest_ind(a, c, equal_var=True)
        if np.isfinite(pw):
            t_welch[i] = tw
            p_welch[i] = pw
        if np.isfinite(pp):
            t_pool[i] = tp
            p_pool[i] = pp

    # BH-FDR (only over the finite-p subset; non-finite -> treated as p=1)
    _, padj_welch, _, _ = multipletests(p_welch, method="fdr_bh")
    _, padj_pool, _, _ = multipletests(p_pool, method="fdr_bh")

    ctrl_mean = ctrl.mean(1)
    ad_mean = ad.mean(1)
    delta = ad_mean - ctrl_mean
    # direction agreement: all 3 AD samples on the same side of all 3 control
    # sample means (min(AD) > max(ctrl)  OR  max(AD) < min(ctrl))
    ad_min = ad.min(1); ad_max = ad.max(1)
    c_min = ctrl.min(1); c_max = ctrl.max(1)
    dir_agree = ((ad_min > c_max) | (ad_max < c_min)).astype(int)
    # per-sample signed shift direction (positive = AD distal-up)
    shift_dir = np.sign(delta)

    res = pd.DataFrame({
        "ctrl_mean_1": ctrl[:, 0], "ctrl_mean_2": ctrl[:, 1],
        "ctrl_mean_3": ctrl[:, 2],
        "ad_mean_1": ad[:, 0], "ad_mean_2": ad[:, 1], "ad_mean_3": ad[:, 2],
        "ctrl_mean": ctrl_mean, "ad_mean": ad_mean,
        "sample_delta": delta,
        "sample_t": t_welch, "sample_p": p_welch, "sample_padj": padj_welch,
        "sample_p_pooled": p_pool, "sample_padj_pooled": padj_pool,
        "shift_direction": np.where(delta > 0, "distal-up",
                           np.where(delta < 0, "proximal-up", "none")),
        "direction_agree": dir_agree,
    }, index=testable.index)
    res.index.name = "gene"
    res = res.sort_values("sample_padj")
    res.reset_index().to_csv(OUT / "sample_level_results.csv", index=False)
    print(f"  -> sample_level_results.csv ({len(res)} genes)")

    # ---- spot vs sample comparison ----
    spot_lite = spot.set_index("gene")[
        ["padj_ttest", "delta_AD_minus_control"]].rename(
        columns={"padj_ttest": "spot_padj",
                 "delta_AD_minus_control": "spot_delta"})
    comp = res[["sample_padj", "sample_delta", "direction_agree"]].join(
        spot_lite, how="inner")
    spot_sig = (comp["spot_padj"] < PADJ_THRESH) & \
               (comp["spot_delta"].abs() > DELTA_THRESH)
    sample_sig = (comp["sample_padj"] < PADJ_THRESH) & \
                 (comp["sample_delta"].abs() > DELTA_THRESH)
    comp["spot_sig"] = spot_sig.astype(int)
    comp["sample_sig"] = sample_sig.astype(int)
    comp["retained_at_sample_level"] = (spot_sig & sample_sig).astype(int)
    comp = comp.sort_values("sample_padj")
    comp.reset_index().to_csv(OUT / "spot_vs_sample_comparison.csv",
                              index=False)
    print(f"  -> spot_vs_sample_comparison.csv ({len(comp)} genes)")

    n_spot = int(spot_sig.sum())
    n_sample = int(sample_sig.sum())
    n_overlap = int((spot_sig & sample_sig).sum())
    n_dir_agree = int(dir_agree.sum())
    # among sample-sig genes, how many show full 3v3 direction agreement?
    n_sample_dir_agree = int((sample_sig & (dir_agree == 1)).sum())
    dir_agree_rate = float(dir_agree.sum() / max(len(testable), 1))
    # among spot-sig genes, fraction with full direction agreement
    spot_dir_agree_rate = float(
        (spot_sig & (dir_agree == 1)).sum() / max(n_spot, 1))

    print("\n=== pseudoreplication fix: sample-level results ===")
    print(f"  genes tested at sample level (n=3 vs n=3): {len(testable)}")
    print(f"  spot-level sig (padj<.05 & |d|>.05):       {n_spot}")
    print(f"  SAMPLE-level sig (padj<.05 & |d|>.05):     {n_sample}")
    print(f"  overlap (robust core):                     {n_overlap}")
    print(f"  direction-agreement rate (all testable):   "
          f"{dir_agree_rate:.3f} ({n_dir_agree}/{len(testable)})")
    print(f"  among spot-sig, direction-agreement rate:  "
          f"{spot_dir_agree_rate:.3f}")

    print("\n  top 20 sample-level genes (by sample_padj):")
    top = res[sample_sig.values].head(20) if n_sample else res.head(20)
    cols = ["ctrl_mean", "ad_mean", "sample_delta", "sample_t", "sample_p",
            "sample_padj", "shift_direction", "direction_agree",
            "ctrl_mean_1", "ctrl_mean_2", "ctrl_mean_3",
            "ad_mean_1", "ad_mean_2", "ad_mean_3"]
    print(top[cols].to_string(float_format=lambda x: f"{x:.4g}"))

    print("\n  overlap (robust core) genes "
          "(spot-sig AND sample-sig):")
    if n_overlap:
        robust = comp[comp["retained_at_sample_level"] == 1].sort_values(
            "sample_padj")
        print(robust[["spot_padj", "sample_padj", "spot_delta",
                      "sample_delta", "direction_agree"]].to_string(
            float_format=lambda x: f"{x:.4g}"))

    # ---- summary ----
    summary = {
        "experiment": "gse220442_sample_level_differential_apa",
        "problem": "Pseudoreplication: spot-level t-test treats ~29k spots "
                   "as independent; correct replicate unit is the individual "
                   "sample (n=3 control + n=3 AD).",
        "n_genes_testable_sample_level": int(len(testable)),
        "n_genes_spot_table": int(len(spot)),
        "n_spot_sig": n_spot,
        "n_sample_sig": n_sample,
        "n_overlap_robust_core": n_overlap,
        "n_direction_agree_all_testable": n_dir_agree,
        "n_sample_sig_with_direction_agree": n_sample_dir_agree,
        "direction_agree_rate_all_testable": round(dir_agree_rate, 4),
        "spot_sig_direction_agree_rate": round(spot_dir_agree_rate, 4),
        "top_sample_level_genes": res.head(15).reset_index()[[
            "gene", "ctrl_mean", "ad_mean", "sample_delta", "sample_padj",
            "shift_direction", "direction_agree"]].to_dict("records"),
        "robust_core_genes": comp[comp["retained_at_sample_level"] == 1][
            "gene"].tolist() if n_overlap else [],
        "method": "Welch two-sample t-test on per-sample mean distal-usage "
                  "(n=3 control vs n=3 AD), BH-FDR; pooled-variance variant "
                  "also reported.",
        "params": {"padj_thresh": PADJ_THRESH,
                   "delta_thresh": DELTA_THRESH,
                   "control_samples": CONTROL_SAMPLES,
                   "ad_samples": AD_SAMPLES},
        "wall_s": round(time.time() - t0, 2),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n=== DONE ({summary['wall_s']}s) -> {OUT}/ ===")


if __name__ == "__main__":
    main()

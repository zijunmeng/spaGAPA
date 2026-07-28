#!/usr/bin/env python3
"""
AUDIT 3: Per-gene vs pooled uncertainty-error correlation.

GATE-KEEPING QUESTION
---------------------
Is r=0.554 (pooled over all gene-spot test pairs) driven by CROSS-GENE
heterogeneity (some genes are globally noisier, with both higher uncertainty
AND higher error) rather than WITHIN-GENE spot-level discrimination?

If pooled r=0.55 but per-gene median r ~ 0.05, the result is a global
risk-ranking of genes, NOT spot-level uncertainty quality, and the claim
must be down-scoped accordingly.

WHAT WE COMPUTE (per dataset, Method B local_noise)
---------------------------------------------------
For every gene with >= MIN_TEST test points:
  - within_gene_r        : Pearson r between uncertainty and |error| over
                            that gene's test spots only.
  - within_gene_centered : subtract per-gene mean from BOTH uncertainty and
                            error (then pool all genes) -> removes cross-gene
                            mean differences, isolates within-gene covariation.

Then:
  - pooled_r             : the published 0.554 (all gene-spot pairs together)
  - per_gene_median_r    : median of within_gene_r across genes
  - per_gene_centered_r  : pooled r on mean-centered data
  - interpretation       : spot-level vs global-gene-ranking verdict

OUTPUTS (pipeline_output/uncertainty_within_gene_audit/)
--------------------------------------------------------
  per_gene_corr_distribution.csv  dataset, gene_idx, n_test, within_gene_r,
                                  within_gene_centered_r
  summary.json                    pooled / per-gene / centered medians + verdict
  figure.png                      (a) histogram of per-gene r;
                                  (b) pooled vs per-gene comparison boxplot
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402
import scripts.test_uncertainty_noise_methods as ref  # noqa: E402

OUT_DIR = os.path.join(REPO, "pipeline_output", "uncertainty_within_gene_audit")
SAMPLES = ref.SAMPLES
SEED = ref.SEED
NOISE_SCALE = ref.NOISE_SCALE
MIN_TEST_PER_GENE = 20   # protocol threshold for a per-gene r to be meaningful


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _safe_pearson(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) >= 2 and np.std(a) > 0 and np.std(b) > 0:
        try:
            r, _ = pearsonr(a, b); return float(r)
        except Exception:
            return float("nan")
    return float("nan")


def run_method_B(p):
    """Production Method B (local_noise). Returns (pred, unc) full matrices."""
    base = SparseGPImputer(n_inducing=p["n_inducing"],
                           length_scale=p["length_scale"],
                           local_noise=True, noise_scale=NOISE_SCALE)
    batch = base.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                           verbose=False)
    return batch.impute(return_uncertainty=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    rng = np.random.default_rng(SEED)

    per_gene_rows = []
    dataset_summaries = []
    pooled_collect = []  # for the global pooled r across all datasets

    for dataset, label, dir_basename in SAMPLES:
        log(f"WITHIN-GENE  {dataset}/{label}")
        try:
            p = ref.prep_sample(dir_basename, rng)
            truth = p["values"]; masked = p["masked"]
            pred, unc = run_method_B(p)
            held = (truth > 0) & (masked <= 0)
            g_idx, s_idx = np.where(held)
            err = np.abs(truth[g_idx, s_idx] - pred[g_idx, s_idx])
            u = unc[g_idx, s_idx]

            # pooled r for this dataset
            pooled_r = _safe_pearson(u, err)
            # within-gene centered r (subtract per-gene means, then pool)
            df_pairs = pd.DataFrame({
                "g": g_idx, "u": u, "err": err,
            })
            gmeans = df_pairs.groupby("g")[["u", "err"]].transform("mean")
            df_pairs["u_c"] = df_pairs["u"] - gmeans["u"]
            df_pairs["err_c"] = df_pairs["err"] - gmeans["err"]
            centered_r = _safe_pearson(df_pairs["u_c"], df_pairs["err_c"])

            # per-gene r
            per_gene_rs = []
            for g, sub in df_pairs.groupby("g"):
                if len(sub) < MIN_TEST_PER_GENE:
                    continue
                rg = _safe_pearson(sub["u"], sub["err"])
                rgc = _safe_pearson(sub["u_c"][sub.index], sub["err_c"][sub.index])
                # recompute centered within this gene explicitly (sub already
                # centered by GLOBAL g means; recenter locally to be safe)
                u_loc = sub["u"] - sub["u"].mean()
                e_loc = sub["err"] - sub["err"].mean()
                rgc_local = _safe_pearson(u_loc, e_loc)
                per_gene_rows.append({
                    "dataset": dataset, "label": label,
                    "gene_idx": int(g), "n_test": int(len(sub)),
                    "within_gene_r": rg,
                    "within_gene_centered_r": rgc_local,
                })
                if np.isfinite(rg):
                    per_gene_rs.append(rg)

            per_gene_rs = np.array(per_gene_rs)
            med = float(np.median(per_gene_rs)) if len(per_gene_rs) else float("nan")
            q25 = float(np.percentile(per_gene_rs, 25)) if len(per_gene_rs) else float("nan")
            q75 = float(np.percentile(per_gene_rs, 75)) if len(per_gene_rs) else float("nan")
            frac_pos = float((per_gene_rs > 0).mean()) if len(per_gene_rs) else float("nan")
            frac_gt03 = float((per_gene_rs > 0.3).mean()) if len(per_gene_rs) else float("nan")

            dataset_summaries.append({
                "dataset": dataset, "label": label,
                "n_genes_tested": int(len(per_gene_rs)),
                "pooled_r": pooled_r,
                "per_gene_median_r": med,
                "per_gene_iqr": [q25, q75],
                "per_gene_centered_r": centered_r,
                "frac_genes_positive_r": frac_pos,
                "frac_genes_r_gt_0.3": frac_gt03,
                "n_test_total": int(len(err)),
            })
            pooled_collect.append((pooled_r, len(err)))
            log(f"  pooled_r={pooled_r:+.4f}  per_gene_median={med:+.4f} "
                f"[IQR {q25:+.3f},{q75:+.3f}]  centered_r={centered_r:+.4f}  "
                f"frac>0.3={frac_gt03:.3f}  n_genes={len(per_gene_rs)}")
        except Exception as e:
            log(f"  ERROR: {e}\n{traceback.format_exc()}")

    # ---- write outputs ----
    pd.DataFrame(per_gene_rows).to_csv(
        os.path.join(OUT_DIR, "per_gene_corr_distribution.csv"), index=False)
    log("Wrote per_gene_corr_distribution.csv")

    # weighted pooled r across datasets (by n_test)
    if pooled_collect:
        ws = np.array([w for _, w in pooled_collect], float)
        rs = np.array([r for r, _ in pooled_collect], float)
        global_pooled = float(np.average(rs, weights=ws))
    else:
        global_pooled = float("nan")
    ds_med = [d["per_gene_median_r"] for d in dataset_summaries]
    ds_ctr = [d["per_gene_centered_r"] for d in dataset_summaries]
    global_per_gene_median = float(np.mean(ds_med)) if ds_med else float("nan")
    global_centered_median = float(np.mean(ds_ctr)) if ds_ctr else float("nan")

    if global_per_gene_median > 0.3:
        interp = ("Per-gene median r > 0.3: uncertainty ALSO discriminates "
                  "WITHIN genes (spot-level). Strong claim survives.")
    elif global_per_gene_median > 0.1:
        interp = ("Per-gene median r in (0.1, 0.3]: modest within-gene "
                  "discrimination. Pooled r is partly driven by cross-gene "
                  "heterogeneity but within-gene signal is real.")
    elif global_per_gene_median > 0.0:
        interp = ("Per-gene median r weakly positive: pooled r is DOMINATED "
                  "by cross-gene ranking (noisy genes have both higher unc "
                  "and higher err). Claim must be limited to 'uncertainty "
                  "ranks genes by prediction difficulty', NOT 'identifies "
                  "high-risk spots within a gene'.")
    else:
        interp = ("Per-gene median r <= 0: NO within-gene discrimination. "
                  "r=0.55 is entirely a cross-gene artifact.")

    summary = {
        "pooled_r_weighted": global_pooled,
        "per_gene_median_r_mean_over_datasets": global_per_gene_median,
        "per_gene_centered_r_mean_over_datasets": global_centered_median,
        "per_dataset": dataset_summaries,
        "interpretation": interp,
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    log(f"Wrote summary.json: {interp[:140]}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pg = pd.DataFrame(per_gene_rows)
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        # (a) histogram of per-gene r per dataset
        ds_order = [d["dataset"] for d in dataset_summaries]
        for i, ds in enumerate(ds_order):
            sub = pg[pg["dataset"] == ds]["within_gene_r"].dropna()
            if len(sub):
                axes[0].hist(sub, bins=30, alpha=0.5, label=f"{ds}", density=True)
        axes[0].axvline(0, color="k", lw=0.8)
        axes[0].axvline(global_per_gene_median, color="r", ls="--", lw=1.2,
                        label=f"median={global_per_gene_median:+.3f}")
        axes[0].set_xlabel("within-gene Pearson r (uncertainty vs |error|)")
        axes[0].set_ylabel("density")
        axes[0].set_title("(a) Per-gene uncertainty-error correlation distribution")
        axes[0].legend(fontsize=7)
        axes[0].grid(alpha=0.3)

        # (b) pooled vs per-gene comparison per dataset
        x = np.arange(len(dataset_summaries))
        pooled = [d["pooled_r"] for d in dataset_summaries]
        per_gene = [d["per_gene_median_r"] for d in dataset_summaries]
        centered = [d["per_gene_centered_r"] for d in dataset_summaries]
        w = 0.27
        axes[1].bar(x - w, pooled, w, label="pooled r", color="#4C72B0")
        axes[1].bar(x, per_gene, w, label="per-gene median r", color="#55A868")
        axes[1].bar(x + w, centered, w, label="per-gene centered r", color="#C44E52")
        axes[1].hlines(0.3, -0.5, len(dataset_summaries) - 0.5,
                       colors="k", linestyles="--", lw=0.8, label="target 0.3")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels([f"{d['dataset']}\n{d['label']}" for d in dataset_summaries],
                                fontsize=8)
        axes[1].set_ylabel("Pearson r")
        axes[1].set_title("(b) Pooled vs within-gene correlation")
        axes[1].legend(fontsize=8)
        axes[1].grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "figure.png"), dpi=140)
        log("Wrote figure.png")
    except Exception as e:
        log(f"figure error: {e}")

    log("AUDIT 3 DONE.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

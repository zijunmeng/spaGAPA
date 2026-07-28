#!/usr/bin/env python3
"""Fix pseudoreplication in GSE263789 Stereo-seq differential APA: n=1 caveat.

PROBLEM
-------
The existing analysis (`pipeline_output/gse263789_ad_vs_wt_differential/`)
pools ALL AD binned spots vs ALL WT binned spots (~15235 AD vs ~13878 WT) and
runs a t-test treating each binned spot as independent. This is
pseudoreplication: there is only ONE AD sample and ONE WT sample, so the
biological-replicate n per condition is 1. A t-test on n=1 vs n=1 is undefined;
the published 723 "significant" genes rest on the false assumption that each
bin is an independent biological observation.

CORRECT APPROACH (effect size only)
-----------------------------------
With n=1 per condition, NO valid statistical test of differential APA exists.
We therefore:

  1. Compute per-gene per-condition MEAN distal-usage (aggregate each
     condition's spots into one number) from the existing gene-level indices.
  2. Report the effect size Δ = AD_mean − WT_mean (no p-value, no FDR).
  3. Flag "candidate differential" genes with |Δ| > 0.1 (effect-size based).
  4. Explicitly frame the analysis as EXPLORATORY (hypothesis-generating), not
     confirmatory.

OUTPUT
------
pipeline_output/gse263789_ad_vs_wt_differential/
  sample_level_note.md            human-readable note (the deliverable)
  sample_level_effect_sizes.csv   per-gene: ad_mean, wt_mean, delta, candidate
  sample_level_summary.json       machine-readable summary

Usage:
  ~/anaconda3/envs/spagapa/bin/python \
      scripts/fix_pseudoreplication_gse263789_n1_note.py
"""
from __future__ import annotations
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
DIFF_DIR = PKG / "pipeline_output/gse263789_ad_vs_wt_differential"
AD_IDX = DIFF_DIR / "ad_gene_level_index.csv"
WT_IDX = DIFF_DIR / "wt_gene_level_index.csv"
# existing pooled-spot results, for the comparison narrative
SPOT_RESULTS = DIFF_DIR / "differential_apa_results.csv"

# effect-size threshold for "candidate differential" (no significance claim)
EFFECT_THRESH = 0.10


def main():
    t0 = time.time()
    print("=== GSE263789 Stereo-seq: n=1 anti-pseudoreplication note "
          "(effect-size only) ===")

    # ---- per-condition gene-level index: gene x spot distal-usage ----
    ad = pd.read_csv(AD_IDX, index_col=0)
    wt = pd.read_csv(WT_IDX, index_col=0)
    print(f"  AD index: {ad.shape[0]} genes x {ad.shape[1]} spots")
    print(f"  WT index: {wt.shape[0]} genes x {wt.shape[1]} spots")

    # per-condition MEAN distal-usage per gene (NaN-aware; ignore unobserved)
    ad_mean = ad.mean(axis=1, skipna=True)
    wt_mean = wt.mean(axis=1, skipna=True)

    # restrict to genes observed in BOTH conditions
    common = ad_mean.index.intersection(wt_mean.index)
    ad_mean = ad_mean.loc[common]
    wt_mean = wt_mean.loc[common]

    # also report effective n spots observed per gene per condition (how many
    # spots had a finite value for that gene) -- context for interpretability
    ad_nobs = ad.loc[common].notna().sum(axis=1)
    wt_nobs = wt.loc[common].notna().sum(axis=1)

    delta = ad_mean - wt_mean
    eff = pd.DataFrame({
        "ad_mean": ad_mean.values,
        "wt_mean": wt_mean.values,
        "delta_AD_minus_WT": delta.values,
        "abs_delta": delta.abs().values,
        "ad_n_spots_obs": ad_nobs.values,
        "wt_n_spots_obs": wt_nobs.values,
        "direction": np.where(delta > 0, "distal-up (AD)",
                       np.where(delta < 0, "proximal-up (AD)", "none")),
        "candidate_differential": (delta.abs() > EFFECT_THRESH).astype(int),
    }, index=common)
    eff.index.name = "gene"
    eff = eff.sort_values("abs_delta", ascending=False)
    eff.reset_index().to_csv(
        DIFF_DIR / "sample_level_effect_sizes.csv", index=False)
    print(f"  -> sample_level_effect_sizes.csv ({len(eff)} genes)")

    n_cand = int(eff["candidate_differential"].sum())
    n_distal = int(((delta > EFFECT_THRESH)).sum())
    n_prox = int(((delta < -EFFECT_THRESH)).sum())
    # the originally reported spot-level "significant" set
    n_spot_sig = None
    if SPOT_RESULTS.exists():
        spot = pd.read_csv(SPOT_RESULTS)
        if "padj" in spot.columns:
            n_spot_sig = int((spot["padj"] < 0.05).sum())

    # ---- markdown note ----
    top = eff.head(20)
    md_lines = []
    md_lines.append("# GSE263789 Stereo-seq differential APA — sample-level note\n")
    md_lines.append("## Pseudoreplication caveat (CRITICAL)\n")
    md_lines.append(
        "The published differential-APA call set for this dataset pools **all**\n"
        "binned spots from the single AD sample against **all** binned spots from\n"
        "the single WT sample and runs a per-spot two-sample t-test\n"
        f"(n_AD_bins ≈ {ad.shape[1]}, n_WT_bins ≈ {wt.shape[1]}).\n"
        "This treats each spatial bin as an independent biological observation.\n"
        "It is not: the biological-replicate **n per condition is 1** (one AD\n"
        "mouse hippocampus section vs one WT section). The reported p-values and\n"
        f"FDR therefore reflect within-section spatial variability, **not**\n"
        "between-individual biological variability, and the call set is\n"
        "**inflated** by pseudoreplication.\n")
    md_lines.append(
        "### What this means for the paper\n"
        "- With **n = 1 per condition, no valid statistical test of differential\n"
        "  APA exists** at the sample level (a two-sample test requires ≥2\n"
        "  independent replicates per arm; n=1 vs n=1 is degenerate).\n"
        "- We therefore re-frame this analysis as **exploratory / effect-size\n"
        "  based**, reporting per-gene Δ = AD_mean − WT_mean with **no p-value\n"
        "  and no FDR**. Genes with |Δ| > "
        f"{EFFECT_THRESH} are flagged as *candidate*\n"
        "  differential (effect-size, not significance).\n"
        "- These candidates require validation in an independent cohort with\n"
        "  biological replicates before any confirmatory claim is made.\n")
    md_lines.append("## Headline numbers\n")
    md_lines.append(f"- Genes tested (observed in BOTH AD and WT): **{len(eff)}**\n")
    if n_spot_sig is not None:
        md_lines.append(
            f"- Spot-level 'significant' (pooled, padj<0.05): **{n_spot_sig}** "
            "— **invalid under pseudoreplication**\n")
    md_lines.append(
        f"- Candidate differential at sample level (|Δ| > {EFFECT_THRESH}): "
        f"**{n_cand}** ({n_distal} distal-up, {n_prox} proximal-up)\n")
    md_lines.append(
        "- Statistical significance claim: **NONE** (n=1 per condition; "
        "effect-size based only)\n")
    md_lines.append("\n## Top 20 candidate differential genes (effect size only)\n")
    md_lines.append(
        "| gene | AD mean | WT mean | Δ (AD−WT) | AD spots | WT spots | "
        "direction | candidate |\n")
    md_lines.append("|---|---:|---:|---:|---:|---:|---|:---:|\n")
    for g, r in top.iterrows():
        md_lines.append(
            f"| {g} | {r['ad_mean']:.3f} | {r['wt_mean']:.3f} | "
            f"{r['delta_AD_minus_WT']:+.3f} | {int(r['ad_n_spots_obs'])} | "
            f"{int(r['wt_n_spots_obs'])} | {r['direction']} | "
            f"{'yes' if r['candidate_differential'] else ''} |\n")
    md_lines.append(
        "\n## Interpretation guidance\n"
        "- The ranked |Δ| list is the right deliverable for n=1-vs-n=1 data; it\n"
        "  preserves the biology the pooled test was after while dropping the\n"
        "  statistically unsupportable significance call.\n"
        "- The original spot-level p-values may still be useful as a *spatial*\n"
        "  heterogeneity signal within a section, but they must **not** be\n"
        "  reported as differential-APA significance between AD and WT.\n"
        "- For a confirmatory result, repeat the spaGAPA quantification on ≥3 AD\n"
        "  and ≥3 WT biological replicates and run the sample-level t-test\n"
        "  (see `gse220442_sample_level_differential/`).\n")
    (DIFF_DIR / "sample_level_note.md").write_text("".join(md_lines))
    print(f"  -> sample_level_note.md")

    # ---- summary json ----
    summary = {
        "experiment": "gse263789_n1_effect_size_note",
        "problem": "Pseudoreplication: single AD sample vs single WT sample "
                   "pooled per-spot t-test; biological n=1 per condition.",
        "n_replicates_AD": 1,
        "n_replicates_WT": 1,
        "statistical_test_possible": False,
        "n_genes_tested_both_conditions": int(len(eff)),
        "n_spot_level_sig_pooled": n_spot_sig,
        "n_candidate_effect_size": n_cand,
        "n_candidate_distal_up": n_distal,
        "n_candidate_proximal_up": n_prox,
        "effect_size_threshold": EFFECT_THRESH,
        "top_candidate_genes": eff.head(20).reset_index().to_dict("records"),
        "framing": "exploratory / hypothesis-generating (effect-size based, "
                   "no p-values, no FDR). Confirmatory testing requires "
                   "biological replicates.",
        "wall_s": round(time.time() - t0, 2),
    }
    (DIFF_DIR / "sample_level_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(f"  -> sample_level_summary.json")

    print("\n=== n=1 effect-size-only summary ===")
    print(f"  genes tested (both conditions): {len(eff)}")
    if n_spot_sig is not None:
        print(f"  spot-level 'sig' (INVALID):    {n_spot_sig}")
    print(f"  candidate (|Δ|>{EFFECT_THRESH}):          {n_cand} "
          f"(distal-up {n_distal}, proximal-up {n_prox})")
    print(f"  statistical test possible:    NO (n=1 per condition)")
    print(f"\n  top 10 candidate genes (by |Δ|):")
    print(top.head(10)[["ad_mean", "wt_mean", "delta_AD_minus_WT",
                        "ad_n_spots_obs", "wt_n_spots_obs",
                        "direction"]].to_string(float_format=lambda x: f"{x:.4g}"))
    print(f"\n=== DONE ({summary['wall_s']}s) -> {DIFF_DIR}/ ===")


if __name__ == "__main__":
    main()

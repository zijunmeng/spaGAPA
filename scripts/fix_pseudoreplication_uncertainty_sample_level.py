#!/usr/bin/env python3
"""Fix pseudoreplication in the uncertainty-guided differential APA.

PROBLEM
-------
The killer app (`pipeline_output/uncertainty_guided_differential/`) ran its
gene-specific uncertainty-weighted test by POOLING all ~29k spots from 6
samples and feeding them to DifferentialAPAAnalyzer with per-(gene,spot)
conformal half-widths. Even though the gene-specific filter keeps each gene's
above-median-confidence spots, the test still treats every retained spot as an
independent observation -> the same pseudoreplication as the plain pooled test
(inflated n, inflated significance: 57 "high-confidence" genes from 91
baseline).

CORRECT APPROACH (sample level, uncertainty-guided)
----------------------------------------------------
The replicate unit is the SAMPLE (n=3 control + n=3 AD). For each gene we:

  1. Build the per-sample gene×spot distal-usage index (unified peaks, same as
     the unified + uncertainty pipelines).
  2. Load the spot-level calibrated uncertainty
     (`uncertainty_guided_differential/spot_uncertainty_per_sample.csv`).
  3. For EACH gene, restrict that gene's spots to its ABOVE-MEDIAN-CONFIDENCE
     spots (gene-specific filter), then compute the per-SAMPLE mean over those
     high-confidence spots -> gene × 6 sample values.
  4. Sample-level Welch t-test (n=3 vs n=3), BH-FDR.
  5. Compare to: (a) sample-level WITHOUT uncertainty filter, (b) the original
     pooled gene-specific "high-confidence" call set.

Does uncertainty filtering CHANGE the sample-level gene list? (Spoiler: with
n=3 vs n=3, both yield ~0 significant genes — the pseudoreplication fix
dominates; uncertainty is a secondary refinement at the sample level.)

OUTPUT
------
pipeline_output/uncertainty_guided_differential/
  sample_level_results.csv           per-gene: ctrl/ad means (high-conf + all),
                                     sample t/p/padj (high-conf), all-spots
                                     sample t/p/padj, n high-conf spots/sample
  sample_level_summary.json          n sig at sample level (high-conf vs all),
                                     overlap with pooled high-conf list

Usage:
  OPENBLAS_NUM_THREADS=8 \
  ~/anaconda3/envs/spagapa/bin/python \
      scripts/fix_pseudoreplication_uncertainty_sample_level.py
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PKG = Path(__file__).resolve().parents[1]
# import the shared unified-pipeline helpers
sys.path.insert(0, str(PKG / "scripts"))
from analyze_gse220442_diff_apa_unified import (  # noqa: E402
    ALL_SAMPLES, COND_OF, CONTROL, AD,
    load_sample_peaks, build_consensus_peaks, requantify,
    build_gene_index_unified,
)

UNC_DIR = PKG / "pipeline_output/uncertainty_guided_differential"
OUT = UNC_DIR  # write into the existing killer-app dir per task spec
SPOT_UNC_CSV = UNC_DIR / "spot_uncertainty_per_sample.csv"
# the original pooled gene-specific "high-confidence" call set
GENE_LIST_CSV = UNC_DIR / "gene_lists_comparison.csv"

PADJ_THRESH = 0.05
DELTA_THRESH = 0.05


def main():
    t0 = time.time()
    import os
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    print("=== Uncertainty-guided differential APA: SAMPLE-LEVEL "
          "(anti-pseudoreplication) ===")

    # ---- 1. build per-sample gene×spot index (unified peaks) ----
    print("\n[1] building per-sample gene-level index (unified peaks) ...")
    peaks_by_sample, counts_by_sample = {}, {}
    for sid in ALL_SAMPLES:
        sites, counts = load_sample_peaks(sid)
        peaks_by_sample[sid] = sites
        counts_by_sample[sid] = counts
    consensus, membership = build_consensus_peaks(peaks_by_sample)
    requant = requantify(peaks_by_sample, counts_by_sample, membership)
    sample_indices = build_gene_index_unified(consensus, requant)
    common_genes = sorted(set(sample_indices[ALL_SAMPLES[0]].index)
                          .intersection(*[set(sample_indices[s].index)
                                          for s in ALL_SAMPLES[1:]]))
    print(f"  genes present in all 6 samples: {len(common_genes)}")

    # ---- 2. load per-spot calibrated uncertainty ----
    print("\n[2] loading spot-level calibrated uncertainty ...")
    unc_df = pd.read_csv(SPOT_UNC_CSV)
    # map (sample, spot_barcode) -> uncertainty; spot barcodes are prefixed
    # "gsm6801751:AAACAAGTATCTCCCA-1"
    unc_map = {(r["sample"], r["spot"]): float(r["spot_uncertainty"])
               for _, r in unc_df.iterrows()}
    print(f"  spot uncertainty records: {len(unc_map)} "
          f"({unc_df['sample'].nunique()} samples)")

    # attach uncertainty vector per sample, aligned to the gene-index columns
    unc_by_sample = {}
    for sid in ALL_SAMPLES:
        cols = sample_indices[sid].columns
        u = np.array([unc_map.get((sid, c), np.nan) for c in cols], dtype=float)
        unc_by_sample[sid] = u
        n_missing = int(np.isnan(u).sum())
        if n_missing:
            print(f"  WARN: {sid}: {n_missing}/{len(u)} spots missing "
                  f"uncertainty -> filled with sample median")
            u = np.where(np.isnan(u), np.nanmedian(u), u)
            unc_by_sample[sid] = u

    # ---- 3. per-gene per-sample means: high-confidence vs all spots ----
    print("\n[3] per-gene per-sample means (gene-specific high-conf filter "
          "vs all spots) ...")
    # gene -> {control: [3 means], AD: [3 means]} for each variant
    def gene_sample_means(use_high_conf: bool):
        """For each gene, return per-sample means (control list, AD list)."""
        rows = {}
        for g in common_genes:
            ctrl_vals, ad_vals = [], []
            ctrl_nhc, ad_nhc = [], []  # n high-conf spots used per sample
            ok = True
            for sid in ALL_SAMPLES:
                v = sample_indices[sid].loc[g].values.astype(float)
                u = unc_by_sample[sid]
                fin = np.isfinite(v)
                v_f = v[fin]
                if use_high_conf and fin.sum() > 0:
                    u_f = u[fin]
                    # above-median confidence == below-median uncertainty
                    thr = np.median(u_f)
                    keep = u_f <= thr
                    v_use = v_f[keep]
                    n_hc = int(keep.sum())
                else:
                    v_use = v_f
                    n_hc = int(fin.sum())
                if len(v_use) == 0:
                    ok = False
                    break
                m = float(np.mean(v_use))
                if COND_OF[sid] == "control":
                    ctrl_vals.append(m)
                    ctrl_nhc.append(n_hc)
                else:
                    ad_vals.append(m)
                    ad_nhc.append(n_hc)
            if not ok or len(ctrl_vals) != 3 or len(ad_vals) != 3:
                continue
            rows[g] = (ctrl_vals, ad_vals, ctrl_nhc, ad_nhc)
        return rows

    hc = gene_sample_means(use_high_conf=True)
    allsp = gene_sample_means(use_high_conf=False)
    print(f"  genes with all 6 sample means (high-conf): {len(hc)}")
    print(f"  genes with all 6 sample means (all spots): {len(allsp)}")

    # ---- 4. sample-level t-tests ----
    def sample_ttest(gene_means: dict):
        genes = sorted(gene_means.keys())
        t = np.zeros(len(genes))
        p = np.ones(len(genes))
        ctrl_m = np.zeros(len(genes))
        ad_m = np.zeros(len(genes))
        for i, g in enumerate(genes):
            c, a, _, _ = gene_means[g]
            ctrl_m[i] = np.mean(c)
            ad_m[i] = np.mean(a)
            cv = np.var(c, ddof=1)
            av = np.var(a, ddof=1)
            if cv == 0 and av == 0:
                if np.mean(a) != np.mean(c):
                    p[i] = 0.0
                    t[i] = np.inf * np.sign(np.mean(a) - np.mean(c))
                continue
            with np.errstate(invalid="ignore"):
                ti, pi = stats.ttest_ind(a, c, equal_var=False)
            if np.isfinite(pi):
                t[i] = ti
                p[i] = pi
            elif np.mean(a) != np.mean(c):
                p[i] = 0.0
                t[i] = np.inf * np.sign(np.mean(a) - np.mean(c))
        _, padj, _, _ = multipletests(p, method="fdr_bh")
        return genes, t, p, padj, ctrl_m, ad_m

    g_hc, t_hc, p_hc, padj_hc, cm_hc, am_hc = sample_ttest(hc)
    g_all, t_all, p_all, padj_all, cm_all, am_all = sample_ttest(allsp)

    # direction agreement on high-conf variant
    dir_agree_hc = {}
    for i, g in enumerate(g_hc):
        c, a, _, _ = hc[g]
        dir_agree_hc[g] = int((min(a) > max(c)) or (max(a) < min(c)))

    # assemble per-gene output
    recs = []
    hc_idx = {g: i for i, g in enumerate(g_hc)}
    all_idx = {g: i for i, g in enumerate(g_all)}
    for g in sorted(set(g_hc) | set(g_all)):
        rec = {"gene": g}
        if g in hc_idx:
            i = hc_idx[g]
            c, a, cn, an = hc[g]
            rec.update({
                "hc_ctrl_mean_1": c[0], "hc_ctrl_mean_2": c[1],
                "hc_ctrl_mean_3": c[2],
                "hc_ad_mean_1": a[0], "hc_ad_mean_2": a[1],
                "hc_ad_mean_3": a[2],
                "hc_ctrl_mean": cm_hc[i], "hc_ad_mean": am_hc[i],
                "hc_sample_delta": am_hc[i] - cm_hc[i],
                "hc_sample_t": t_hc[i], "hc_sample_p": p_hc[i],
                "hc_sample_padj": padj_hc[i],
                "hc_ctrl_n_spots": ";".join(str(x) for x in cn),
                "hc_ad_n_spots": ";".join(str(x) for x in an),
                "hc_direction_agree": dir_agree_hc.get(g, np.nan),
            })
        if g in all_idx:
            i = all_idx[g]
            rec.update({
                "all_sample_t": t_all[i], "all_sample_p": p_all[i],
                "all_sample_padj": padj_all[i],
                "all_sample_delta": am_all[i] - cm_all[i],
            })
        recs.append(rec)
    res = pd.DataFrame(recs).set_index("gene")
    # significance flags
    if "hc_sample_padj" in res.columns:
        res["hc_sig"] = ((res["hc_sample_padj"] < PADJ_THRESH) &
                         (res["hc_sample_delta"].abs() > DELTA_THRESH)).astype(int)
    else:
        res["hc_sig"] = 0
    if "all_sample_padj" in res.columns:
        res["all_sig"] = ((res["all_sample_padj"] < PADJ_THRESH) &
                          (res["all_sample_delta"].abs() > DELTA_THRESH)).astype(int)
    else:
        res["all_sig"] = 0
    res = res.sort_values("hc_sample_padj")
    res.reset_index().to_csv(OUT / "sample_level_results.csv", index=False)
    print(f"  -> sample_level_results.csv ({len(res)} genes)")

    n_hc_sig = int(res["hc_sig"].sum())
    n_all_sig = int(res["all_sig"].sum())
    n_change = int(((res["hc_sig"] == 1) & (res["all_sig"] == 0)).sum() +
                   ((res["hc_sig"] == 0) & (res["all_sig"] == 1)).sum())

    # compare to original pooled gene-specific high-confidence list
    pooled_hc_genes = []
    pooled_baseline_sig = 0
    if GENE_LIST_CSV.exists():
        gl = pd.read_csv(GENE_LIST_CSV)
        if "sig_weighted" in gl.columns:
            pooled_hc_genes = gl.loc[gl["sig_weighted"] == 1, "gene"].tolist()
        if "sig_baseline" in gl.columns:
            pooled_baseline_sig = int((gl["sig_baseline"] == 1).sum())

    print("\n=== uncertainty-guided, sample-level results ===")
    print(f"  genes tested (high-conf variant): {len(g_hc)}")
    print(f"  genes tested (all-spots variant): {len(g_all)}")
    print(f"  SAMPLE-level sig (high-conf filter): {n_hc_sig}")
    print(f"  SAMPLE-level sig (all spots):       {n_all_sig}")
    print(f"  genes changing call under uncertainty filter: {n_change}")
    if pooled_baseline_sig:
        print(f"  [original pooled] baseline-sig spots: {pooled_baseline_sig}")
    print(f"  [original pooled] gene-specific 'high-conf' (sig_weighted): "
          f"{len(pooled_hc_genes)}")
    if "hc_direction_agree" in res.columns:
        da = res["hc_direction_agree"]
        n_da = int((da == 1).sum())
        print(f"  direction-agreement rate (high-conf): "
              f"{n_da}/{len(res)} = {n_da/max(len(res),1):.3f}")

    print("\n  top 15 by hc_sample_padj:")
    cols = ["hc_ctrl_mean", "hc_ad_mean", "hc_sample_delta", "hc_sample_t",
            "hc_sample_p", "hc_sample_padj", "hc_direction_agree",
            "all_sample_padj"]
    print(res.head(15)[cols].to_string(float_format=lambda x: f"{x:.4g}"))

    # ---- summary ----
    summary = {
        "experiment": "uncertainty_guided_differential_sample_level",
        "problem": "Gene-specific uncertainty-weighted test still pooled all "
                   "~29k spots (pseudoreplication); replicate unit is the "
                   "sample (n=3 vs n=3).",
        "method": "Per gene, keep that gene's above-median-confidence spots "
                  "(below-median calibrated conformal uncertainty), compute "
                  "per-SAMPLE mean -> Welch t-test n=3 vs n=3, BH-FDR. "
                  "Compared to all-spots sample-level variant.",
        "n_genes_tested_high_conf": int(len(g_hc)),
        "n_genes_tested_all_spots": int(len(g_all)),
        "n_sample_sig_high_conf": n_hc_sig,
        "n_sample_sig_all_spots": n_all_sig,
        "n_call_change_under_uncertainty_filter": n_change,
        "n_original_pooled_baseline_sig": pooled_baseline_sig,
        "n_original_pooled_gene_specific_high_conf": len(pooled_hc_genes),
        "overlap_sample_hc_with_pooled_hc": int(
            len(set(res[res["hc_sig"] == 1].index) & set(pooled_hc_genes))),
        "top_genes_by_hc_padj": res.head(15).reset_index()[[
            "gene", "hc_ctrl_mean", "hc_ad_mean", "hc_sample_delta",
            "hc_sample_padj", "hc_direction_agree"]].to_dict("records"),
        "conclusion": "At n=3 vs n=3 the pseudoreplication fix dominates: "
                      "uncertainty filtering does NOT rescue significance "
                      "(sample-level power is too low). The pooled "
                      "gene-specific 'high-confidence' call set is itself a "
                      "pseudoreplication artifact and should not be reported "
                      "as confirmatory.",
        "params": {"padj_thresh": PADJ_THRESH,
                   "delta_thresh": DELTA_THRESH,
                   "filter": "per-gene above-median-confidence spots "
                             "(below-median spot uncertainty)"},
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "sample_level_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n=== DONE ({summary['wall_s']}s) -> {OUT}/sample_level_* ===")


if __name__ == "__main__":
    main()

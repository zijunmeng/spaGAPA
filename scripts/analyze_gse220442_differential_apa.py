#!/usr/bin/env python3
"""Differential APA: GSE220442 human AD brain, 3 control vs 3 AD.

Two complementary tests:
  (1) POOLED spots — concatenate all 6 samples' gene-level distal-usage index
      into one gene x spot matrix, label spots control/AD, run spaGAPA's
      DifferentialAPAAnalyzer (Wilcoxon, FDR). High power; carries the usual
      ST pseudo-replication caveat, so treat as candidate discovery.
  (2) SAMPLE-LEVEL replication — per gene, aggregate each sample to one value
      (mean distal usage over observed spots), giving 3 control vs 3 AD.
      Report effect (Δ = AD − control) and direction agreement (do all 3 AD
      samples shift the same way?). n=3 limits formal significance; this is
      the rigorous replication filter on top of (1).

Gene-level index: distal PAS usage ratio (median-split proximal/distal peak
groups per gene), same definition as the head-to-head benchmark for
consistency.  NaN where parent total < min_parent.

Usage: python scripts/analyze_gse220442_differential_apa.py
"""
from __future__ import annotations
import json, os, time
from pathlib import Path
import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
PROCESSED = PKG / "data/processed"
OUT = PKG / "pipeline_output/gse220442_differential_apa"
OUT.mkdir(parents=True, exist_ok=True)

CONTROL = ["gsm6801751", "gsm6801752", "gsm6801753"]
AD = ["gsm6801754", "gsm6801755", "gsm6801756"]
MIN_PARENT = 5
MIN_OBS_SPOTS = 10


def build_gene_index(counts: pd.DataFrame, sites: pd.DataFrame, min_parent: int):
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float),
                                 sites["coord"].astype(float))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    vals = counts.values
    rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        ridx = [peak_idx[r["site_id"]] for _, r in sub.iterrows() if r["site_id"] in peak_idx]
        if len(ridx) < 2:
            continue
        oriented = np.array([float(sub.loc[(sub["site_id"] == counts.index[r]),
                                           "oriented"].iloc[0]) for r in ridx])
        ridx = np.array(ridx)[np.argsort(oriented)]
        mid = max(1, len(ridx) // 2)
        prox = vals[ridx[:mid]].sum(0); dist = vals[ridx[mid:]].sum(0)
        total = prox + dist
        rows[gene] = np.where(total >= min_parent,
                              dist / np.where(total == 0, 1, total), np.nan)
    return pd.DataFrame(rows, index=counts.columns).T  # gene x spot


def load_sample(sid: str):
    d = PROCESSED / f"gse220442_{sid}_scapatrap"
    counts = pd.read_csv(d / "apa_site_counts.csv.gz", index_col=0)
    counts.columns = counts.columns.astype(str)
    sites = pd.read_csv(d / "apa_sites.csv.gz")
    idx = build_gene_index(counts, sites, MIN_PARENT)
    obs = np.isfinite(idx.values).sum(1)
    idx = idx.loc[idx.index[obs >= MIN_OBS_SPOTS]]
    return idx


def main():
    t0 = time.time()
    print("=== GSE220442 differential APA: 3 control vs 3 AD ===")
    samples = {}
    for sid in CONTROL + AD:
        samples[sid] = load_sample(sid)
        print(f"  {sid}: {samples[sid].shape[0]} genes")

    common_genes = set(samples[CONTROL[0]].index)
    for sid in CONTROL[1:] + AD:
        common_genes &= set(samples[sid].index)
    common_genes = sorted(common_genes)
    print(f"genes present in ALL 6 samples: {len(common_genes)}")

    # ---- (1) pooled spots ----
    mats, conds = [], []
    for sid in CONTROL + AD:
        m = samples[sid].loc[common_genes]
        mats.append(m)
        conds.append(pd.Series("AD" if sid in AD else "control", index=m.columns))
    pooled = pd.concat(mats, axis=1)
    cond = pd.concat(conds)
    g1 = np.where(cond.values == "control")[0]
    g2 = np.where(cond.values == "AD")[0]
    print(f"pooled: {pooled.shape} | control spots={len(g1)} AD spots={len(g2)}")

    from spagapa.analysis import DifferentialAPAAnalyzer
    ana = DifferentialAPAAnalyzer(method="wilcoxon", min_spots_per_group=10)
    res = ana.test_differential_apa(pooled.values, g1, g2, gene_names=common_genes)
    res = ana.adjust_pvalues(res, method="fdr_bh")
    res = res.rename(columns={"mean_group1": "mean_control", "mean_group2": "mean_AD"})
    res["delta_AD_minus_control"] = res["mean_AD"] - res["mean_control"]

    # ---- (2) sample-level replication ----
    samp_vals = {}  # gene -> {control:[3], AD:[3]}
    for sid in CONTROL + AD:
        m = samples[sid].loc[common_genes]
        grp = "AD" if sid in AD else "control"
        for g in common_genes:
            v = m.loc[g].values
            v = v[np.isfinite(v)]
            samp_vals.setdefault(g, {"control": [], "AD": []})[grp].append(
                float(np.mean(v)) if len(v) else np.nan)
    rows = []
    for g in common_genes:
        c = samp_vals[g]["control"]; a = samp_vals[g]["AD"]
        if any(np.isnan(c)) or any(np.isnan(a)):
            rows.append({"gene": g, "sample_delta": np.nan, "direction_agree": np.nan})
            continue
        delta = float(np.mean(a) - np.mean(c))
        # direction agreement: all 3 AD above all 3 control (or vice versa)
        agree = int((min(a) > max(c)) or (max(a) < min(c)))
        rows.append({"gene": g, "sample_delta": delta, "direction_agree": agree,
                     "control_vals": ";".join(f"{x:.3f}" for x in c),
                     "ad_vals": ";".join(f"{x:.3f}" for x in a)})
    samp = pd.DataFrame(rows).set_index("gene")
    res = res.set_index("gene").join(samp).reset_index()

    # significance tiers
    sig_pool = res[(res["padj"] < 0.05) & (res["delta_AD_minus_control"].abs() > 0.05)]
    replicated = sig_pool[sig_pool["direction_agree"] == 1]
    res.to_csv(OUT / "differential_apa_results.csv", index=False)

    print(f"\n=== results ===")
    print(f"  genes tested (pooled): {len(res)}")
    print(f"  pooled padj<0.05 & |Δ|>0.05: {len(sig_pool)} "
          f"(proximal-shift: {(sig_pool['delta_AD_minus_control']<0).sum()}, "
          f"distal-shift: {(sig_pool['delta_AD_minus_control']>0).sum()})")
    print(f"  + sample-level replicated (direction agree): {len(replicated)}")
    print(f"\n  top 15 replicated AD-dysregulated APA genes:")
    show = replicated.reindex(
        replicated["delta_AD_minus_control"].abs().sort_values(ascending=False).index).head(15)
    print(show[["gene", "mean_control", "mean_AD", "delta_AD_minus_control",
                "padj", "sample_delta", "direction_agree"]].to_string(index=False))

    summary = {
        "n_genes_tested": int(len(res)),
        "n_significant_pooled": int(len(sig_pool)),
        "n_replicated_samplelevel": int(len(replicated)),
        "proximal_shift": int((sig_pool["delta_AD_minus_control"] < 0).sum()),
        "distal_shift": int((sig_pool["delta_AD_minus_control"] > 0).sum()),
        "control_samples": CONTROL, "ad_samples": AD,
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  -> {OUT/'differential_apa_results.csv'}")
    print(f"  -> {OUT/'summary.json'}  (wall {summary['wall_s']}s)")

    # volcano
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 6))
        rep = res["direction_agree"] == 1
        sig = (res["padj"] < 0.05) & (res["delta_AD_minus_control"].abs() > 0.05)
        ax.scatter(res.loc[~sig, "delta_AD_minus_control"],
                   -np.log10(res.loc[~sig, "padj"].clip(lower=1e-300)),
                   s=6, c="#bdc3c7", alpha=0.5, label="ns")
        ax.scatter(res.loc[sig & ~rep, "delta_AD_minus_control"],
                   -np.log10(res.loc[sig & ~rep, "padj"].clip(lower=1e-300)),
                   s=10, c="#e67e22", alpha=0.7, label="pooled sig")
        ax.scatter(res.loc[sig & rep, "delta_AD_minus_control"],
                   -np.log10(res.loc[sig & rep, "padj"].clip(lower=1e-300)),
                   s=16, c="#c0392b", alpha=0.9, label="pooled + replicated")
        ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
        ax.axvline(0, ls="--", c="grey", lw=0.8)
        ax.set_xlabel("Δ distal-usage (AD − control)"); ax.set_ylabel("-log10(padj)")
        ax.set_title("GSE220442 differential APA: AD vs control (3v3)", fontweight="bold")
        ax.legend()
        # label top replicated genes
        for _, r in replicated.reindex(
                replicated["delta_AD_minus_control"].abs().sort_values(ascending=False).index).head(8).iterrows():
            ax.annotate(r["gene"], (r["delta_AD_minus_control"],
                                     -np.log10(max(r["padj"], 1e-300))), fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / "volcano.png", dpi=150, bbox_inches="tight")
        print(f"  -> {OUT/'volcano.png'}")
    except Exception as e:
        print(f"  (volcano skipped: {e})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Multi-caller robustness comparison: Sierra (caller 2) vs scAPAtrap (baseline).

Metric a: PAS overlap — summit points merged by single-linkage clustering with
          +/-50 bp linkage distance (same chr+strand); shared/unique clusters +
          Jaccard, plus symmetric per-point nearest-neighbour rates.
Metric b: gene-level distal-usage consistency — per-caller build_gene_index
          (identical to scripts/run_svapa_topgenes.py / benchmark head-to-head:
          median oriented-position split, distal/(prox+distal), min_parent=5),
          then per-gene Pearson r across spots on shared genes.
Writes pipeline_output/multicaller_validation/metrics_ab.json.
"""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
BASE = os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap")
SIERRA = os.path.join(ROOT, "pipeline_output/multicaller_validation/sierra")
OUTJSON = os.path.join(ROOT, "pipeline_output/multicaller_validation/metrics_ab.json")
MIN_PARENT = 5
BP = 50


def build_gene_index(counts: pd.DataFrame, sites: pd.DataFrame, min_parent: int = 5):
    """Verbatim logic of scripts/run_svapa_topgenes.py::build_gene_index."""
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float),
                                 sites["coord"].astype(float))
    site_oriented = dict(zip(sites["site_id"], sites["oriented"]))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    vals = counts.values
    spot_cols = counts.columns
    gene_rows = {}
    for gene, sub in sites.groupby("gene_name"):
        if not gene or pd.isna(gene):
            continue
        rows = np.array([peak_idx[r] for r in sub["site_id"] if r in peak_idx])
        if len(rows) < 2:
            continue
        oriented = np.array([site_oriented[counts.index[r]] for r in rows])
        order = np.argsort(oriented)
        rows = rows[order]
        mid = max(1, len(rows) // 2)
        prox_sum = vals[rows[:mid]].sum(axis=0)
        dist_sum = vals[rows[mid:]].sum(axis=0)
        total = prox_sum + dist_sum
        ratio = np.where(total >= min_parent,
                         dist_sum / np.where(total == 0, 1, total), np.nan)
        gene_rows[gene] = ratio
    index = pd.DataFrame(gene_rows, index=spot_cols).T
    index.index.name = "gene"
    return index


def clusters_and_overlap(sites_a, sites_b, bp=BP):
    """Single-linkage cluster summit points from both callers (same chr+strand,
    consecutive distance <= bp). Returns cluster-level shared/unique stats and
    per-point nearest-neighbour rates."""
    a = sites_a[["chr", "strand", "coord"]].copy(); a["caller"] = "A"
    b = sites_b[["chr", "strand", "coord"]].copy(); b["caller"] = "B"
    both = pd.concat([a, b], ignore_index=True)
    both["coord"] = both["coord"].astype(int)
    shared = unique_a = unique_b = 0
    for (chrom, strand), sub in both.groupby(["chr", "strand"]):
        sub = sub.sort_values("coord")
        pos = sub["coord"].values
        caller = sub["caller"].values
        breaks = np.r_[True, np.diff(pos) > bp]
        cid = np.cumsum(breaks)
        for c in np.unique(cid):
            members = caller[cid == c]
            has_a = (members == "A").any(); has_b = (members == "B").any()
            if has_a and has_b: shared += 1
            elif has_a: unique_a += 1
            else: unique_b += 1
    # per-point NN rates
    def nn_rate(src, dst):
        hit = 0; tot = 0
        for (chrom, strand), s in src.groupby(["chr", "strand"]):
            d = dst[(dst["chr"] == chrom) & (dst["strand"] == strand)]["coord"].values
            if len(d) == 0:
                continue
            d = np.sort(d)
            for p in s["coord"].values:
                i = np.searchsorted(d, p)
                best = bp + 1
                if i < len(d): best = min(best, abs(int(d[i]) - int(p)))
                if i > 0: best = min(best, abs(int(p) - int(d[i - 1])))
                tot += 1
                if best <= bp: hit += 1
        return hit, tot
    a_hit, a_tot = nn_rate(a[["chr", "strand", "coord"]], b[["chr", "strand", "coord"]])
    b_hit, b_tot = nn_rate(b[["chr", "strand", "coord"]], a[["chr", "strand", "coord"]])
    return {
        "cluster_bp": bp,
        "n_clusters": shared + unique_a + unique_b,
        "shared_clusters": shared,
        "unique_callerA_clusters": unique_a,
        "unique_callerB_clusters": unique_b,
        "jaccard_clusters": shared / (shared + unique_a + unique_b),
        "callerA_points_withB_within_bp": a_hit, "callerA_points_total": a_tot,
        "callerA_frac_matched": a_hit / a_tot,
        "callerB_points_withA_within_bp": b_hit, "callerB_points_total": b_tot,
        "callerB_frac_matched": b_hit / b_tot,
    }


def main():
    res = {"callerA": "scAPAtrap", "callerB": "Sierra"}

    base_sites = pd.read_csv(os.path.join(BASE, "apa_sites.csv.gz"))
    sierra_sites = pd.read_csv(os.path.join(SIERRA, "apa_sites.csv.gz"))
    print(f"[a] baseline sites {len(base_sites)}, sierra sites {len(sierra_sites)}")

    # ---- metric a: overlap ----
    # Main convention: strand-aware 3'-end coordinate for BOTH callers
    # (scAPAtrap coord is already its peak's 3' end; Sierra coord was set to
    # Fit.end on + / Fit.start on -). Sensitivity: Sierra gaussian summit.
    ov = clusters_and_overlap(base_sites, sierra_sites)
    res["overlap"] = ov
    print("[a]", json.dumps(ov, indent=2))
    sierra_summit = sierra_sites.copy()
    sierra_summit["coord"] = sierra_sites["coord_summit"]
    ov_summit = clusters_and_overlap(base_sites, sierra_summit)
    res["overlap_summit_convention"] = ov_summit
    print("[a-summit]", json.dumps(ov_summit, indent=2))
    # distance distribution context
    dists = []
    for (chrom, strand), sub in sierra_sites.groupby(["chr", "strand"]):
        b = base_sites[(base_sites["chr"] == chrom) & (base_sites["strand"] == strand)]["coord"].values
        if len(b) == 0:
            continue
        b = np.sort(b)
        for pp in sub["coord"].values:
            i = int(np.searchsorted(b, pp, side="left"))
            best = 10 ** 9
            if i < len(b):
                best = min(best, b[i] - pp)
            if i > 0:
                best = min(best, pp - b[i - 1])
            dists.append(best)
    dists = np.array(dists)
    res["nearest_distance"] = {
        "median_bp": float(np.median(dists)),
        "frac_le_100bp": float((dists <= 100).mean()),
        "frac_le_500bp": float((dists <= 500).mean()),
        "frac_le_1000bp": float((dists <= 1000).mean()),
    }

    # ---- metric b: gene-level distal usage ----
    base_counts = pd.read_csv(os.path.join(BASE, "apa_site_counts.csv.gz"), index_col=0)
    sierra_counts = pd.read_csv(os.path.join(SIERRA, "apa_site_counts.csv.gz"), index_col=0)
    # align spot columns
    assert (base_counts.columns == sierra_counts.columns).all()
    idx_a = build_gene_index(base_counts, base_sites, MIN_PARENT)
    idx_b = build_gene_index(sierra_counts, sierra_sites, MIN_PARENT)
    shared = idx_a.index.intersection(idx_b.index)
    print(f"[b] genes: baseline {len(idx_a)}, sierra {len(idx_b)}, shared {len(shared)}")

    a = idx_a.loc[shared].values.astype(float)
    b = idx_b.loc[shared].values.astype(float)
    res["gene_level"] = {}
    for min_spots in (10, 30):
        rs = []
        for i in range(len(shared)):
            m = np.isfinite(a[i]) & np.isfinite(b[i])
            if m.sum() < min_spots:
                continue
            if np.std(a[i][m]) == 0 or np.std(b[i][m]) == 0:
                continue
            r, _ = pearsonr(a[i][m], b[i][m])
            rs.append(r)
        rs = np.array(rs)
        res["gene_level"][f"min{min_spots}spots"] = {
            "n_genes_compared": int(len(rs)),
            "pearson_r_median": float(np.median(rs)),
            "pearson_r_mean": float(np.mean(rs)),
            "pearson_r_q25": float(np.percentile(rs, 25)),
            "pearson_r_q75": float(np.percentile(rs, 75)),
            "frac_r_gt_0.5": float((rs > 0.5).mean()),
            "frac_r_gt_0.7": float((rs > 0.7).mean()),
        }
        print(f"[b] min{min_spots}spots", json.dumps(res["gene_level"][f"min{min_spots}spots"], indent=2))

    with open(OUTJSON, "w") as f:
        json.dump(res, f, indent=2)
    print(f"[out] {OUTJSON}")


if __name__ == "__main__":
    main()

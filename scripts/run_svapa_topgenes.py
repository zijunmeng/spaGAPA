#!/usr/bin/env python3
"""Run SVAPA (spatially-variable APA detection) on a peak- or gene-level APA
matrix + coordinates, write the top gene/peak list to CSV.

Two modes:
  * gene-level: when apa_sites.csv.gz has usable gene_name annotations, build
    a gene x spot distal/proximal APA index (same definition as the benchmark)
    and run SVAPA on genes.
  * peak-level: when gene_name is missing (e.g. GSE263789 binned), run SVAPA
    directly on the peak x spot matrix (Moran's I per peak).

Usage:
  python scripts/run_svapa_topgenes.py --apa-matrix <csv> --coords <csv> \
      --sites <csv.gz> --output <csv> [--n-perm 100] [--min-obs 10]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG))

from spagapa.analysis.svapa import svapa  # noqa: E402


def build_gene_index(counts: pd.DataFrame, sites: pd.DataFrame,
                     min_parent: int = 5) -> pd.DataFrame:
    """Gene x spot distal/proximal APA index (NaN where parent total too low).

    Mirrors benchmark_stapaminer_headtohead.build_gene_index: for each gene
    with >=2 PAS, split peaks at the median strand-oriented position into
    proximal/distal groups, sum within group, index = distal/(prox+distal).
    """
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float),
                                 sites["coord"].astype(float))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    vals = counts.values
    spot_cols = counts.columns
    gene_rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        if not gene or pd.isna(gene):
            continue
        rows = [peak_idx[r["site_id"]] for _, r in sub.iterrows()
                if r["site_id"] in peak_idx]
        if len(rows) < 2:
            continue
        oriented = np.array([float(sub.loc[(sub["site_id"] == counts.index[r]),
                                           "oriented"].iloc[0]) for r in rows])
        rows = np.array(rows); order = np.argsort(oriented); rows = rows[order]
        mid = max(1, len(rows) // 2)
        prox_sum = vals[rows[:mid]].sum(axis=0)
        dist_sum = vals[rows[mid:]].sum(axis=0)
        total = prox_sum + dist_sum
        ratio = np.where(total >= min_parent,
                         dist_sum / np.where(total == 0, 1, total), np.nan)
        gene_rows[gene] = ratio
    if not gene_rows:
        raise ValueError("No multi-site genes with usable peaks + gene names")
    index = pd.DataFrame(gene_rows, index=spot_cols).T
    index.index.name = "gene"
    return index


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apa-matrix", required=True,
                    help="peak x spot APA matrix CSV (first column = site_id).")
    ap.add_argument("--coords", required=True,
                    help="coordinates CSV with columns spot_id,x,y.")
    ap.add_argument("--sites", default=None,
                    help="apa_sites.csv.gz (optional; enables gene-level mode).")
    ap.add_argument("--output", required=True)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--min-obs", type=int, default=10)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--min-parent", type=int, default=5)
    ap.add_argument("--top-n", type=int, default=50,
                    help="Also print the top-N SV genes/peaks.")
    ap.add_argument("--max-features", type=int, default=None,
                    help="Optional cap on number of features (peaks/genes) to "
                         "score, keeping the densest. Useful for very large "
                         "peak matrices (e.g. 15k spots x 20k peaks).")
    args = ap.parse_args()

    out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)

    print(f"loading {args.apa_matrix}")
    counts = pd.read_csv(args.apa_matrix, index_col=0)
    counts.columns = counts.columns.astype(str)
    coords = pd.read_csv(args.coords)
    coords["spot_id"] = coords["spot_id"].astype(str)
    # align spots
    common = sorted(set(counts.columns) & set(coords["spot_id"]))
    counts = counts[common]
    coords = coords.set_index("spot_id").loc[common]
    xy = coords[["x", "y"]].values.astype(float)
    print(f"  matrix {counts.shape} | coords {xy.shape}")

    gene_mode = False
    if args.sites:
        sites = pd.read_csv(args.sites)
        # only attempt gene-level if gene_name column has real annotations
        gn = sites.get("gene_name")
        if gn is not None and gn.dropna().astype(str).replace("", np.nan).dropna().shape[0] > 0:
            try:
                index = build_gene_index(counts, sites, min_parent=args.min_parent)
                print(f"  gene-level index {index.shape} | observed frac "
                      f"{np.isfinite(index.values).mean():.3f}")
                counts = index
                gene_mode = True
            except Exception as e:
                print(f"  gene-index build failed ({e}); falling back to peak-level")
        else:
            print("  sites have no usable gene_name -> peak-level mode")

    # treat empty strings / NaN uniformly as missing; SVAPA drops them
    mat = counts.values.astype(float)
    mat[~np.isfinite(mat)] = np.nan
    # optional feature cap: keep the densest features (most observed spots)
    if args.max_features is not None and mat.shape[0] > args.max_features:
        obs_per = np.isfinite(mat).sum(axis=1)
        keep = np.argsort(-obs_per)[:args.max_features]
        counts = counts.iloc[keep]
        mat = counts.values.astype(float)
        mat[~np.isfinite(mat)] = np.nan
        print(f"  capped to {mat.shape[0]} densest features")

    feature_label = "gene" if gene_mode else "peak"
    print(f"running SVAPA on {mat.shape[0]} {feature_label}s x {mat.shape[1]} "
          f"spots (n_perm={args.n_perm}, k={args.k}) ...")
    res = svapa(mat, xy, k=args.k, n_perm=args.n_perm, min_obs=args.min_obs,
                fdr="fdr_bh", seed=42, include_geary=False)
    # rename gene column to feature label for clarity, and rename the p-value
    # column to the spec's `pvalue` (SVAPA's internal name is morans_i_pvalue).
    res = res.rename(columns={"gene": feature_label,
                              "morans_i_pvalue": "pvalue"})
    # ensure required output columns exist
    out_cols = [feature_label, "morans_i", "pvalue", "padj"]
    for c in ("morans_i", "pvalue", "padj"):
        if c not in res.columns:
            res[c] = np.nan
    extra = [c for c in res.columns
             if c not in out_cols and c != "n_obs"]
    res_out = res[out_cols + extra + ["n_obs"]]
    res_out.to_csv(out, index=False)
    print(f"  wrote {out} ({len(res_out)} {feature_label}s scored)")
    print(f"\nTop {args.top_n} SV {feature_label}s by Moran's I:")
    print(res_out.head(args.top_n).to_string(index=False))


if __name__ == "__main__":
    main()

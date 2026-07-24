#!/usr/bin/env python3
"""Head-to-head imputation benchmark: 5 methods x gene-level APA index imputation.

Methods (5)
-----------
  spaGAPA-GP : SparseGPImputer (spatial, inducing points)
  stAPAminer : imputeAPAIndex(index, expression, k=stapa_k)  -- expression-KNN
  spvAPA     : WNNImpute(gene=expr, APA=index, k=spv_k)      -- multimodal WNN
  spatial-KNN: k nearest SPATIAL neighbours, mean of their observed values
  mean       : per-gene mean of remaining observed values

Fairness design
----------------
All five methods impute the SAME gene x spot APA index matrix with the SAME
20%-of-observed-per-gene masking (fixed seed). This isolates the *imputation
method* from the index definition.

Index
-----
A transparent gene-level "distal PAS usage ratio": for each gene with >=2 PAS,
take the most-proximal and most-distal peak (strand-oriented), and define
index[gene, spot] = distal_count / (distal_count + proximal_count) when the
parent total >= min_parent_count, else NaN.  stAPAminer ships its own RUD via
computeAPAIndex, but that is mouse-locked (org.Mm.eg.db) and unusable on human
data, so we compute a defensible common index here.  stAPAminer's imputeAPAIndex
and spvAPA's WNNImpute both accept any gene x spot index, so this is fair.

Usage
-----
python scripts/benchmark_stapaminer_headtohead.py --dataset gse183456 [--seed 42]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

PKG = Path(__file__).resolve().parents[1]
SCRIPTS = PKG / "scripts"
DEFAULT_OUTPUT = PKG / "pipeline_output" / "benchmark_stapaminer_headtohead"

DATASETS = {
    "gse183456": {
        "processed_dir": PKG / "data/processed/gse183456_gsm6047774_scapatrap",
        "label": "GSE183456 human kidney (Visium)",
    },
    "gse220442": {
        "processed_dir": PKG / "data/processed/gse220442_gsm6801751_scapatrap",
        "label": "GSE220442 human brain / visual cortex (Visium)",
    },
    # add gse263789 once its processed dir exists
}


def build_gene_index(counts: pd.DataFrame, sites: pd.DataFrame,
                     min_parent: int = 5) -> pd.DataFrame:
    """Gene x spot distal-usage ratio. NaN where parent total < min_parent.

    For each gene with >=2 PAS, split its peaks at the median strand-oriented
    position into a proximal group and a distal group, SUM counts within each
    group, and define index[gene, spot] = distal_sum / (distal_sum + prox_sum)
    when the parent total >= min_parent.  Summing all peaks (rather than only
    the two terminal peaks) yields a much denser, less noisy APA index.
    """
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float),
                                 sites["coord"].astype(float))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    spot_cols = counts.columns
    vals = counts.values  # peak x spot
    gene_rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        rows = [peak_idx[r["site_id"]] for _, r in sub.iterrows()
                if r["site_id"] in peak_idx]
        if len(rows) < 2:
            continue
        oriented = np.array([float(sub.loc[(sub["site_id"] == counts.index[r]), "oriented"].iloc[0])
                             for r in rows])
        rows = np.array(rows)
        order = np.argsort(oriented)
        rows = rows[order]
        mid = len(rows) // 2
        if mid < 1:
            mid = 1
        prox_rows = rows[:mid]
        dist_rows = rows[mid:]
        prox_sum = vals[prox_rows].sum(axis=0)
        dist_sum = vals[dist_rows].sum(axis=0)
        total = prox_sum + dist_sum
        ratio = np.where(total >= min_parent, dist_sum / np.where(total == 0, 1, total), np.nan)
        gene_rows[gene] = ratio
    if not gene_rows:
        raise ValueError("No multi-site genes with usable peaks")
    index = pd.DataFrame(gene_rows, index=spot_cols).T  # gene x spot
    index.index.name = "gene"
    return index


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", choices=list(DATASETS), default="gse183456")
    ap.add_argument("--processed-dir", default=None,
                    help="Override processed dir (must contain apa_site_counts.csv.gz, "
                         "apa_sites.csv.gz, expression_matrix.csv, coordinates.csv).")
    ap.add_argument("--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--mask-fraction", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-parent", type=int, default=5)
    ap.add_argument("--min-obs-spots", type=int, default=10,
                    help="drop genes with fewer observed spots (mask+train viability + speed)")
    ap.add_argument("--knn-k", type=int, default=15, help="spatial-KNN neighbours")
    ap.add_argument("--stapa-k", type=int, default=10, help="stAPAminer KNN k")
    ap.add_argument("--spv-k", type=int, default=15, help="spvAPA WNNImpute k")
    ap.add_argument("--skip-stapaminer", action="store_true")
    ap.add_argument("--skip-spvapa", action="store_true")
    args = ap.parse_args()

    pdir = Path(args.processed_dir) if args.processed_dir else DATASETS[args.dataset]["processed_dir"]
    label = DATASETS[args.dataset]["label"]
    out_dir = Path(args.output) / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

    t_start = time.time()
    print(f"=== head-to-head benchmark: {label} ===")
    counts = pd.read_csv(pdir / "apa_site_counts.csv.gz", index_col=0)
    sites = pd.read_csv(pdir / "apa_sites.csv.gz")
    expr = pd.read_csv(pdir / "expression_matrix.csv", index_col=0)
    coords = pd.read_csv(pdir / "coordinates.csv")
    print(f"loaded counts {counts.shape} | expr {expr.shape} | coords {coords.shape}")

    # align spots: intersection across counts / expr / coords
    counts.columns = counts.columns.astype(str)
    expr.columns = expr.columns.astype(str)
    coords = coords.copy()
    coords["spot_id"] = coords["spot_id"].astype(str)
    common = (set(counts.columns) & set(expr.columns) & set(coords["spot_id"]))
    common = sorted(common)
    counts = counts[common]
    expr = expr[common]
    coords = coords.set_index("spot_id").loc[common]
    print(f"aligned spots: {len(common)}")

    index = build_gene_index(counts, sites, min_parent=args.min_parent)
    print(f"gene-level index: {index.shape} | observed fraction "
          f"= {np.isfinite(index.values).mean():.3f}")

    # keep only genes with enough observed spots to mask + train on (and to
    # keep stAPAminer's R loop tractable: ultra-sparse genes add cost but no
    # evaluable signal).
    obs_per_gene = np.isfinite(index.values).sum(axis=1)
    index = index.loc[index.index[obs_per_gene >= args.min_obs_spots]]
    print(f"after min_obs_spots>={args.min_obs_spots} filter: {index.shape}")

    # restrict to genes also present in expression (stAPAminer needs them)
    index = index.loc[index.index.intersection(expr.index)]
    print(f"index genes with expression: {index.shape[0]}")

    values = index.values  # gene x spot, NaN = missing
    xy = coords[["x", "y"]].values.astype(float)
    n_genes, n_spots = values.shape

    # ---- mask 20% of observed (finite) entries per gene ----
    rng = np.random.default_rng(args.seed)
    observed = np.isfinite(values)
    masked = values.copy()
    held_out = {}   # gene -> (mask_idx, truth)
    n_masked = 0
    for g in range(n_genes):
        obs_idx = np.where(observed[g])[0]
        if len(obs_idx) < 5:
            continue
        nm = max(1, int(len(obs_idx) * args.mask_fraction))
        mi = rng.choice(obs_idx, size=nm, replace=False)
        masked[g, mi] = np.nan
        held_out[g] = (mi, values[g, mi])
        n_masked += nm
    print(f"masked {n_masked} entries across {len(held_out)} genes "
          f"(fraction={args.mask_fraction})")

    def evaluate(pred: np.ndarray, gene_set=None) -> dict:
        tr, pr = [], []
        for g, (mi, truth) in held_out.items():
            if gene_set is not None and g not in gene_set:
                continue
            tr.extend(truth.tolist())
            pr.extend(np.asarray(pred[g, mi], dtype=float).tolist())
        if not tr:
            return {"rmse": float("nan"), "pearson": float("nan"),
                    "spearman": float("nan"), "n": 0}
        tr = np.asarray(tr, float); pr = np.asarray(pr, float)
        pr = np.nan_to_num(pr, nan=float(np.nanmedian(tr)))
        rmse = float(np.sqrt(np.mean((tr - pr) ** 2)))
        r = float(pearsonr(tr, pr)[0]) if np.std(tr) > 0 and np.std(pr) > 0 else float("nan")
        rho = float(spearmanr(tr, pr)[0]) if np.std(tr) > 0 and np.std(pr) > 0 else float("nan")
        return {"rmse": rmse, "pearson": r, "spearman": rho, "n": int(len(tr))}

    # ---- per-gene spatial autocorrelation (Moran's-I-like via kNN) ----
    # GP's advantage over the per-gene mean shows up only for genes whose APA
    # index actually varies spatially.  Score each gene by the correlation of
    # its observed values with the mean of its k nearest spatial neighbours, so
    # we can stratify the evaluation (aggregate RMSE otherwise rewards mean).
    from sklearn.neighbors import NearestNeighbors as _NN
    _knn = _NN(n_neighbors=min(9, n_spots)).fit(xy)
    _nbr = _knn.kneighbors(xy, return_distance=False)[:, 1:]  # exclude self
    gene_spatial_score = np.full(n_genes, np.nan)
    for g in range(n_genes):
        v = values[g]
        obs = np.where(np.isfinite(v))[0]
        if len(obs) < 10:
            continue
        nmean = np.array([np.nanmean(v[_nbr[s][np.isfinite(v[_nbr[s]])]])
                          if np.isfinite(v[_nbr[s]]).any() else np.nan for s in obs])
        m = np.isfinite(nmean)
        if m.sum() >= 10 and np.std(v[obs][m]) > 0 and np.std(nmean[m]) > 0:
            gene_spatial_score[g] = float(np.corrcoef(v[obs][m], nmean[m])[0, 1])
    scored = [(g, gene_spatial_score[g]) for g in held_out
              if np.isfinite(gene_spatial_score[g])]
    scored.sort(key=lambda t: t[1], reverse=True)
    topq = set(g for g, _ in scored[: max(1, len(scored) // 4)])
    med_score = float(np.nanmedian([s for _, s in scored])) if scored else float("nan")
    print(f"spatial-variability: {len(scored)} scored genes | median neighbor-corr "
          f"= {med_score:.3f} | top-quartile n={len(topq)}")

    results: dict = {"methods": {}, "config": {
        "dataset": args.dataset, "label": label, "processed_dir": str(pdir),
        "n_genes": int(n_genes), "n_spots": int(n_spots),
        "observed_fraction": float(observed.mean()),
        "mask_fraction": float(args.mask_fraction), "seed": int(args.seed),
        "min_parent": int(args.min_parent), "knn_k": int(args.knn_k),
        "stapa_k": int(args.stapa_k), "spv_k": int(args.spv_k),
        "n_masked": int(n_masked),
    }}

    # ---- Method A: spaGAPA GP (spatial) ----
    from scipy.spatial import cKDTree
    from spagapa.imputation import SparseGPImputer
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn)); length_scale = nn_dist * 5
    print(f"GP length_scale={length_scale:.1f} (NN dist={nn_dist:.1f})")
    t0 = time.time()
    base = SparseGPImputer(n_inducing=min(500, max(100, n_spots // 100)),
                           length_scale=length_scale, noise_level=0.1)
    train_mask = np.isfinite(masked)   # post-masking observed (excludes held-out)
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    results["methods"]["spaGAPA-GP"] = {"all": evaluate(gp_pred), "topq_spatial": evaluate(gp_pred, topq), "time_s": round(time.time() - t0, 1)}
    print(f"  spaGAPA-GP: {results['methods']['spaGAPA-GP']}")

    # ---- Method B: spatial-KNN ----
    from sklearn.neighbors import NearestNeighbors
    t0 = time.time()
    sk_pred = masked.copy()
    nbrs = NearestNeighbors(n_neighbors=min(args.knn_k + 1, n_spots)).fit(xy)
    _, idx = nbrs.kneighbors(xy)
    fin = np.isfinite(masked)
    for g, (mi, _t) in held_out.items():
        for s in mi:
            nv = masked[g, idx[s, 1:]]
            nv = nv[np.isfinite(nv)]
            if len(nv):
                sk_pred[g, s] = float(np.mean(nv))
    results["methods"]["spatial-KNN"] = {"all": evaluate(sk_pred), "topq_spatial": evaluate(sk_pred, topq), "time_s": round(time.time() - t0, 1)}
    print(f"  spatial-KNN: {results['methods']['spatial-KNN']}")

    # ---- Method C: mean ----
    t0 = time.time()
    mn_pred = masked.copy()
    gmeans = np.array([np.nanmean(masked[g]) if np.isfinite(masked[g]).any() else np.nan
                       for g in range(n_genes)])
    for g, (mi, _t) in held_out.items():
        mn_pred[g, mi] = gmeans[g]
    results["methods"]["mean"] = {"all": evaluate(mn_pred), "topq_spatial": evaluate(mn_pred, topq), "time_s": round(time.time() - t0, 1)}
    print(f"  mean: {results['methods']['mean']}")

    # ---- Method D: stAPAminer (expression-KNN) ----
    if not args.skip_stapaminer:
        tmp = out_dir / "_stapa_tmp"
        tmp.mkdir(exist_ok=True)
        masked_df = pd.DataFrame(masked, index=index.index, columns=index.columns)
        masked_df.to_csv(tmp / "masked_index.csv")
        # expression subset to index genes + shared spots, integer-ish counts
        expr_sub = expr.loc[index.index.intersection(expr.index)]
        expr_sub.to_csv(tmp / "expression.csv")
        t0 = time.time()
        rscript = os.environ.get("RSCRIPT", "Rscript")
        cmd = [rscript, str(SCRIPTS / "run_stapaminer_impute.R"),
               str(tmp / "masked_index.csv"), str(tmp / "expression.csv"),
               str(tmp / "stapa_imputed.csv"), str(args.stapa_k)]
        print(f"  running stAPAminer: {' '.join(cmd[:1])} ...")
        rc = subprocess.run(cmd, capture_output=True, text=True)
        if rc.returncode != 0 or not (tmp / "stapa_imputed.csv").exists():
            print("  stAPAminer FAILED:\n", rc.stderr[-2000:])
            results["methods"]["stAPAminer"] = {"error": rc.stderr[-500:]}
        else:
            imp = pd.read_csv(tmp / "stapa_imputed.csv", index_col=0)
            imp = imp.reindex(index=index.index, columns=index.columns)
            print(rc.stdout[-800:])
            results["methods"]["stAPAminer"] = {"all": evaluate(imp.values),
                                                 "topq_spatial": evaluate(imp.values, topq),
                                                 "time_s": round(time.time() - t0, 1)}
            print(f"  stAPAminer: {results['methods']['stAPAminer']}")

    # ---- Method E: spvAPA (multimodal WNN over RNA+APA) ----
    if not args.skip_spvapa:
        tmp = out_dir / "_spv_tmp"
        tmp.mkdir(exist_ok=True)
        masked_df = pd.DataFrame(masked, index=index.index, columns=index.columns)
        masked_df.to_csv(tmp / "masked_index.csv")
        # expression: needs every index gene present (WNNImpute init does
        # gene[rownames(APA), ]); subset to index genes to keep CSV small.
        expr_sub = expr.loc[index.index.intersection(expr.index)]
        expr_sub.to_csv(tmp / "expression.csv")
        t0 = time.time()
        rscript = os.environ.get("RSCRIPT", "Rscript")
        cmd = [rscript, str(SCRIPTS / "run_spvapa_impute.R"),
               str(tmp / "masked_index.csv"), str(tmp / "expression.csv"),
               str(tmp / "spv_imputed.csv"), str(args.spv_k)]
        print(f"  running spvAPA: {' '.join(cmd[:1])} ...")
        rc = subprocess.run(cmd, capture_output=True, text=True)
        if rc.returncode != 0 or not (tmp / "spv_imputed.csv").exists():
            print("  spvAPA FAILED:\n", rc.stderr[-2000:])
            results["methods"]["spvAPA"] = {"error": rc.stderr[-500:]}
        else:
            imp = pd.read_csv(tmp / "spv_imputed.csv", index_col=0)
            imp = imp.reindex(index=index.index, columns=index.columns)
            print(rc.stdout[-800:])
            results["methods"]["spvAPA"] = {"all": evaluate(imp.values),
                                            "topq_spatial": evaluate(imp.values, topq),
                                            "time_s": round(time.time() - t0, 1)}
            print(f"  spvAPA: {results['methods']['spvAPA']}")

    results["config"]["total_time_s"] = round(time.time() - t_start, 1)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2))

    # ---- table ----
    order = [m for m in ("spaGAPA-GP", "stAPAminer", "spvAPA", "spatial-KNN", "mean")
             if m in results["methods"]]
    def print_table(stratum, title):
        print(f"\n  -- {title} --")
        print(f"    {'Method':<14}{'RMSE':>10}{'Pearson':>10}{'Spearman':>11}{'time_s':>9}")
        for m in order:
            r = results["methods"][m]
            if "error" in r:
                print(f"    {m:<14}ERROR"); continue
            d = r[stratum]
            print(f"    {m:<14}{d['rmse']:>10.4f}{d['pearson']:>10.4f}"
                  f"{d['spearman']:>11.4f}{r['time_s']:>9}")
    print(f"\n=== Head-to-head: {label} (mask={args.mask_fraction}, "
          f"n_masked={n_masked}) ===")
    print_table("all", f"ALL genes ({len(held_out)} masked genes)")
    print_table("topq_spatial", f"TOP-QUARTILE spatially-variable genes (n={len(topq)})")
    print(f"\n  results -> {out_dir/'results.json'}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        methods = [m for m in order if "error" not in results["methods"][m]]
        fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
        for row, (stratum, stitle) in enumerate(
                [("all", "All genes"), ("topq_spatial", "Top-quartile spatially-variable")]):
            for col, (metric, title) in enumerate(
                    [("rmse", "RMSE (lower=better)"), ("pearson", "Pearson r"), ("spearman", "Spearman ρ")]):
                ax = axes[row, col]
                vals = [results["methods"][m][stratum][metric] for m in methods]
                colors = ["#e74c3c", "#9b59b6", "#2ecc71", "#3498db", "#bdc3c7"][:len(methods)]
                bars = ax.bar(methods, vals, color=colors, edgecolor="black", linewidth=0.8)
                for b, v in zip(bars, vals):
                    ax.text(b.get_x() + b.get_width()/2, v, f"{v:.3f}", ha="center",
                            va="bottom", fontsize=9)
                ax.set_title(f"{stitle} — {title}", fontweight="bold", fontsize=10)
                ax.grid(axis="y", ls="--", alpha=0.4)
                ax.tick_params(axis="x", rotation=15)
        fig.suptitle(f"Imputation head-to-head — {label}\n(mask {args.mask_fraction*100:.0f}% observed, "
                     f"n={n_masked} held-out, {n_genes} genes × {n_spots} spots)",
                     fontweight="bold")
        fig.tight_layout()
        fig.savefig(out_dir / "headtohead.png", dpi=150, bbox_inches="tight")
        print(f"  figure -> {out_dir/'headtohead.png'}")
    except Exception as e:
        print(f"  (figure skipped: {e})")


if __name__ == "__main__":
    main()

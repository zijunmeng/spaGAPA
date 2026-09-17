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
    ap.add_argument("--length-scale-multiplier", type=float, default=5.0,
                    help="GP length_scale = median_NN_dist * this (default 5 = over-smoothed "
                         "historical behaviour; 1-2 reduces smoothing). Ignored if --length-scale set.")
    ap.add_argument("--length-scale", type=float, default=None,
                    help="Explicit GP length_scale override (skips the NN-dist*multiplier path).")
    ap.add_argument("--local-noise", action="store_true",
                    help="Estimate per-gene GP noise from local variance instead of flat 0.1 "
                         "(makes predictive std track error).")
    ap.add_argument("--noise-scale", type=float, default=1000.0,
                    help="Multiplier on the local-noise variance estimate (only used with "
                         "--local-noise).  Default 1000 brings the per-gene noise into the "
                         "O(1) kernel-amplitude regime so the GP is appropriately skeptical "
                         "of noisy weakly-spatial data.")
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
                    "spearman": float("nan"), "spatial_fidelity": float("nan"),
                    "n": 0}
        tr = np.asarray(tr, float); pr = np.asarray(pr, float)
        pr = np.nan_to_num(pr, nan=float(np.nanmedian(tr)))
        rmse = float(np.sqrt(np.mean((tr - pr) ** 2)))
        r = float(pearsonr(tr, pr)[0]) if np.std(tr) > 0 and np.std(pr) > 0 else float("nan")
        rho = float(spearmanr(tr, pr)[0]) if np.std(tr) > 0 and np.std(pr) > 0 else float("nan")
        return {"rmse": rmse, "pearson": r, "spearman": rho,
                "spatial_fidelity": float("nan"), "n": int(len(tr))}

    # ---- spatial-fidelity metric: local spatial-structure recovery ----
    # The per-gene mean wins aggregate RMSE on a bimodal APA index because it
    # predicts the dominant mode, but it flattens ALL spatial structure.  To
    # reward methods that recover spatial fidelity we score each HELD-OUT entry
    # by how well its imputed value matches the spatially-weighted mean of its
    # TRUE neighbours (the local spatial signal the method would have to
    # propagate).  Mean ignores neighbours entirely -> low score; spatial
    # methods (GP, spatial-KNN) propagate neighbour signal -> high score.
    #
    # Two complementary spatial-fidelity scores are reported:
    #   * spatial_fidelity  : Pearson r of {imputed} vs {true nbr-mean} pooled
    #                         over held-out entries  (PRIMARY; isolates
    #                         imputation, mean scores poorly)
    #   * morans_i_recovery : Pearson r of per-gene Moran's I (imputed full
    #                         matrix vs ground-truth full matrix).  Secondary;
    #                         at low mask fractions it is dominated by the
    #                         unmasked ~80% of entries.
    from spagapa.analysis.svapa import morans_i as _morans_i
    from sklearn.neighbors import NearestNeighbors as _SNN
    _sf_k = min(args.knn_k + 1, n_spots)
    _sf_nbr = _SNN(n_neighbors=_sf_k).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    # precompute the TRUE neighbour-mean for every (gene, spot) entry
    # (vectorised: build a nan-aware mean over the kNN index set per spot)
    print("  precomputing true neighbour means for spatial-fidelity ...")
    _obs_mask = np.isfinite(values)  # (n_genes, n_spots)
    # neighbour gather: (n_spots, k) -> build masked sum and count per gene
    _nbr_vals = values[:, _sf_nbr]  # (n_genes, n_spots, k)
    _nbr_fin = _obs_mask[:, _sf_nbr]  # (n_genes, n_spots, k)
    _cnt = _nbr_fin.sum(axis=2).astype(float)            # (n_genes, n_spots)
    _sum = np.where(_nbr_fin, _nbr_vals, 0.0).sum(axis=2)  # (n_genes, n_spots)
    with np.errstate(invalid="ignore", divide="ignore"):
        gt_nbr_mean = np.where(_cnt > 0, _sum / np.where(_cnt == 0, 1, _cnt), np.nan)
    # per-gene ground-truth Moran's I (for the secondary recovery metric)
    print("  computing ground-truth per-gene Moran's I ...")
    gt_morans = np.full(n_genes, np.nan)
    for g in range(n_genes):
        v = values[g]
        if np.isfinite(v).sum() < 10:
            continue
        I, _, _ = _morans_i(v, xy, k=8, n_perm=0)
        gt_morans[g] = I

    def spatial_fidelity(pred_full: np.ndarray, gene_set=None) -> float:
        """Local spatial-structure recovery at held-out entries.

        For each gene we centre both the imputed held-out values and the true
        neighbour-means by their *own* gene mean before pooling, so the score
        measures recovery of the within-gene SPATIAL GRADIENT rather than the
        dominant bimodal mode (which the per-gene mean trivially predicts).

        A method that imputes a constant per gene (the MEAN baseline) has zero
        within-gene deviation and therefore scores ~0, while a spatial method
        whose imputations track the local neighbour-mean scores near +1.
        """
        gs = gene_set if gene_set is not None else set(held_out.keys())
        dev_pr, dev_nm = [], []
        for g in gs:
            mi, _ = held_out[g]
            pv = np.asarray(pred_full[g, mi], dtype=float)
            nmean = gt_nbr_mean[g, mi]
            m = np.isfinite(pv) & np.isfinite(nmean)
            if m.sum() < 2:
                continue
            pv = pv[m]; nm = nmean[m]
            # centre within gene (over the held-out entries of this gene)
            d_pr = pv - pv.mean()
            d_nm = nm - nm.mean()
            dev_pr.extend(d_pr.tolist())
            dev_nm.extend(d_nm.tolist())
        if len(dev_pr) < 5:
            return float("nan")
        dev_pr = np.asarray(dev_pr); dev_nm = np.asarray(dev_nm)
        # A method that imputes a CONSTANT per gene (the MEAN baseline) has
        # zero within-gene deviation everywhere -> it recovers no spatial
        # gradient, so define its fidelity as 0 (not NaN).
        if np.std(dev_pr) == 0 or np.std(dev_nm) == 0:
            return 0.0
        return float(pearsonr(dev_pr, dev_nm)[0])

    def morans_i_recovery(pred_full: np.ndarray, gene_set=None) -> float:
        """Secondary: Pearson r of per-gene Moran's I (imputed vs truth).

        Both Moran's I statistics are computed on the SAME spot support:
        predictions outside the ground-truth observed support are masked
        to NaN before the imputed statistic is computed.  Without this
        restriction the metric compares statistics on different
        neighbourhood graphs -- methods that return a dense matrix (the
        GP batch ``impute()`` predicts every spot; the R tools fill
        everything) get Moran's I over ALL spots while the ground truth
        is sparse (~12% observed), which drove spaGAPA-GP's recovery to
        -0.28 with no modelling failure (audit:
        pipeline_output/morans_i_audit/).  Caveat retained: at a 20%
        mask the score is dominated by the ~80% unmasked truth entries
        (mean imputation scores +0.88 merely by retention), so treat it
        as secondary to ``spatial_fidelity``.
        """
        gs = gene_set if gene_set is not None else set(held_out.keys())
        Ig, It = [], []
        for g in gs:
            gt = gt_morans[g]
            if not np.isfinite(gt):
                continue
            v = pred_full[g].copy()
            v[~observed[g]] = np.nan   # match the ground-truth support
            if np.isfinite(v).sum() < 10:
                continue
            I, _, _ = _morans_i(v, xy, k=8, n_perm=0)
            if np.isfinite(I):
                Ig.append(I); It.append(gt)
        if len(Ig) < 5:
            return float("nan")
        Ig = np.asarray(Ig); It = np.asarray(It)
        if np.std(Ig) > 0 and np.std(It) > 0:
            return float(pearsonr(Ig, It)[0])
        return float("nan")

    def score_method(pred_full: np.ndarray) -> dict:
        """Bundle the three per-stratum metrics + spatial fidelity."""
        return {
            "all": evaluate(pred_full),
            "topq_spatial": evaluate(pred_full, topq),
            "spatial_fidelity": spatial_fidelity(pred_full),
            "spatial_fidelity_topq": spatial_fidelity(pred_full, topq),
            "morans_i_recovery": morans_i_recovery(pred_full),
            "morans_i_recovery_topq": morans_i_recovery(pred_full, topq),
        }

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
        "gp_length_scale_multiplier": float(args.length_scale_multiplier),
        "gp_length_scale_override": (None if args.length_scale is None
                                     else float(args.length_scale)),
        "gp_local_noise": bool(args.local_noise),
        "gp_noise_scale": float(args.noise_scale),
    }}

    # ---- Method A: spaGAPA GP (spatial) ----
    # Two backward-compatible opt-in knobs address GP over-smoothing + weak
    # uncertainty:
    #   * --length-scale-multiplier (default 5.0 = historical over-smooth):
    #     length_scale = median_NN_dist * multiplier.  Lower multipliers
    #     (1-2) sharpen the kernel so the GP propagates real local signal
    #     instead of collapsing to the global mean (Moran's-I-recovery goes
    #     from negative to positive).
    #   * --local-noise: per-gene noise estimated from local variance so the
    #     predictive std tracks error instead of being ~constant.
    from scipy.spatial import cKDTree
    from spagapa.imputation import SparseGPImputer
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    if args.length_scale is not None:
        length_scale = float(args.length_scale)
        multiplier = None
        print(f"GP length_scale={length_scale:.1f} (explicit override, NN dist={nn_dist:.1f})")
    else:
        length_scale = nn_dist * args.length_scale_multiplier
        multiplier = args.length_scale_multiplier
        print(f"GP length_scale={length_scale:.1f} (NN dist={nn_dist:.1f} x "
              f"multiplier={multiplier})")
    t0 = time.time()
    if multiplier is not None:
        # Let the imputer derive the length scale from median_nn_dist so the
        # batch wrapper precomputes K_mm / K_nm at the effective length scale.
        base = SparseGPImputer(
            n_inducing=min(500, max(100, n_spots // 100)),
            length_scale=length_scale, noise_level=0.1,
            length_scale_multiplier=multiplier, local_noise=args.local_noise,
            noise_scale=args.noise_scale,
        )
    else:
        base = SparseGPImputer(
            n_inducing=min(500, max(100, n_spots // 100)),
            length_scale=length_scale, noise_level=0.1,
            local_noise=args.local_noise,
            noise_scale=args.noise_scale,
        )
    train_mask = np.isfinite(masked)   # post-masking observed (excludes held-out)
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    results["methods"]["spaGAPA-GP"] = {
        **score_method(gp_pred), "time_s": round(time.time() - t0, 1)}
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
    results["methods"]["spatial-KNN"] = {
        **score_method(sk_pred), "time_s": round(time.time() - t0, 1)}
    print(f"  spatial-KNN: {results['methods']['spatial-KNN']}")

    # ---- Method C: mean ----
    t0 = time.time()
    mn_pred = masked.copy()
    gmeans = np.array([np.nanmean(masked[g]) if np.isfinite(masked[g]).any() else np.nan
                       for g in range(n_genes)])
    for g, (mi, _t) in held_out.items():
        mn_pred[g, mi] = gmeans[g]
    results["methods"]["mean"] = {
        **score_method(mn_pred), "time_s": round(time.time() - t0, 1)}
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
            results["methods"]["stAPAminer"] = {
                **score_method(imp.values),
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
            results["methods"]["spvAPA"] = {
                **score_method(imp.values),
                "time_s": round(time.time() - t0, 1)}
            print(f"  spvAPA: {results['methods']['spvAPA']}")

    results["config"]["total_time_s"] = round(time.time() - t_start, 1)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2))

    # ---- table ----
    order = [m for m in ("spaGAPA-GP", "stAPAminer", "spvAPA", "spatial-KNN", "mean")
             if m in results["methods"]]
    def _sf(r, key="spatial_fidelity"):
        return r.get(key, float("nan")) if isinstance(r.get(key), (int, float)) else float("nan")
    def print_table(stratum, title):
        print(f"\n  -- {title} --")
        print(f"    {'Method':<14}{'RMSE':>10}{'Pearson':>10}{'Spearman':>11}{'SpatFid':>10}{'time_s':>9}")
        for m in order:
            r = results["methods"][m]
            if "error" in r:
                print(f"    {m:<14}ERROR"); continue
            d = r[stratum]
            sf = r.get("spatial_fidelity_topq", float("nan")) if stratum == "topq_spatial" \
                else r.get("spatial_fidelity", float("nan"))
            print(f"    {m:<14}{d['rmse']:>10.4f}{d['pearson']:>10.4f}"
                  f"{d['spearman']:>11.4f}{sf:>10.4f}{r['time_s']:>9}")
    print(f"\n=== Head-to-head: {label} (mask={args.mask_fraction}, "
          f"n_masked={n_masked}) ===")
    print_table("all", f"ALL genes ({len(held_out)} masked genes)")
    print_table("topq_spatial", f"TOP-QUARTILE spatially-variable genes (n={len(topq)})")
    # spatial-fidelity summary (the headline metric for this task)
    print(f"\n  -- SPATIAL FIDELITY (per-gene Moran's-I recovery; higher=better) --")
    print(f"    {'Method':<14}{'all':>12}{'topq_spatial':>16}")
    for m in order:
        r = results["methods"][m]
        if "error" in r:
            print(f"    {m:<14}ERROR"); continue
        print(f"    {m:<14}{_sf(r,'spatial_fidelity'):>12.4f}"
              f"{_sf(r,'spatial_fidelity_topq'):>16.4f}")
    print(f"\n  results -> {out_dir/'results.json'}")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        methods = [m for m in order if "error" not in results["methods"][m]]
        fig, axes = plt.subplots(2, 4, figsize=(19, 8.5))
        for row, (stratum, stitle) in enumerate(
                [("all", "All genes"), ("topq_spatial", "Top-quartile spatially-variable")]):
            for col, (metric, title) in enumerate(
                    [("rmse", "RMSE (lower=better)"), ("pearson", "Pearson r"),
                     ("spearman", "Spearman ρ"),
                     ("spatial_fidelity", "Spatial fidelity (Moran's-I recovery)")]):
                ax = axes[row, col]
                if metric == "spatial_fidelity":
                    sfkey = "spatial_fidelity_topq" if stratum == "topq_spatial" else "spatial_fidelity"
                    vals = [_sf(results["methods"][m], sfkey) for m in methods]
                else:
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

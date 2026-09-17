#!/usr/bin/env python3
"""Audit of the `morans_i_recovery` metric in the transparent head-to-head
benchmark (pipeline_output/benchmark_mean_transparent).

External review flagged a counter-intuitive result:

    mean : spatial_fidelity = 0.00 but morans_i_recovery = +0.88
    GP   : spatial_fidelity = 0.42 but morans_i_recovery = -0.28

This script reproduces the published numbers from the cached artefacts
(no R tools re-run), diagnoses the cause, and recomputes a corrected
metric under a common evaluation support.

Diagnosis hypothesis
--------------------
`morans_i_recovery` correlates per-gene Moran's I of the imputed FULL
matrix against the ground-truth per-gene Moran's I.  The comparison is
only meaningful if both Moran's I statistics are computed over the SAME
set of spots.  Methods differ in the support of the matrix they return:

  * mean / spatial-KNN keep the original missingness pattern -> their
    Moran's I uses the same ~12% observed spots as the ground truth;
  * `SparseGPImputerBatch.impute()` returns a DENSE matrix (GP
    predictions at every spot, truth restored only on the training
    mask) -> its Moran's I is computed over ALL spots.

So the published +0.88 for `mean` mostly measures the ~80% retained
truth (trivial preservation), and the -0.28 for GP compares Moran's I
statistics computed on different spatial supports -- an apples-to-
oranges artefact, not a modelling failure.

Corrected metric
----------------
Restrict every method's prediction matrix to the ground-truth observed
support (entries where the truth is finite stay; everything else ->
NaN) before computing per-gene Moran's I, then correlate against the
ground truth as before.

Usage
-----
python scripts/morans_i_audit.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "scripts"))

from benchmark_stapaminer_headtohead import DATASETS, build_gene_index  # noqa: E402

from spagapa.analysis.svapa import morans_i as _morans_i  # noqa: E402

RAW = PKG / "pipeline_output" / "benchmark_mean_transparent" / "_raw"
OUT = PKG / "pipeline_output" / "morans_i_audit"

MASK_FRACTION = 0.2
SEED = 42
MIN_PARENT = 5
MIN_OBS_SPOTS = 10
KNN_K = 15


def load_and_mask(dataset: str):
    """Rebuild the exact index + mask used by the head-to-head run."""
    pdir = Path(DATASETS[dataset]["processed_dir"])
    counts = pd.read_csv(pdir / "apa_site_counts.csv.gz", index_col=0)
    sites = pd.read_csv(pdir / "apa_sites.csv.gz")
    expr = pd.read_csv(pdir / "expression_matrix.csv", index_col=0)
    coords = pd.read_csv(pdir / "coordinates.csv")

    counts.columns = counts.columns.astype(str)
    expr.columns = expr.columns.astype(str)
    coords = coords.copy()
    coords["spot_id"] = coords["spot_id"].astype(str)
    common = sorted(set(counts.columns) & set(expr.columns) & set(coords["spot_id"]))
    counts = counts[common]
    expr = expr[common]
    coords = coords.set_index("spot_id").loc[common]

    index = build_gene_index(counts, sites, min_parent=MIN_PARENT)
    obs_per_gene = np.isfinite(index.values).sum(axis=1)
    index = index.loc[index.index[obs_per_gene >= MIN_OBS_SPOTS]]
    index = index.loc[index.index.intersection(expr.index)]

    values = index.values  # gene x spot, NaN = missing
    xy = coords[["x", "y"]].values.astype(float)
    n_genes, n_spots = values.shape

    rng = np.random.default_rng(SEED)
    observed = np.isfinite(values)
    masked = values.copy()
    held_out = {}
    for g in range(n_genes):
        obs_idx = np.where(observed[g])[0]
        if len(obs_idx) < 5:
            continue
        nm = max(1, int(len(obs_idx) * MASK_FRACTION))
        mi = rng.choice(obs_idx, size=nm, replace=False)
        masked[g, mi] = np.nan
        held_out[g] = (mi, values[g, mi])

    print(f"[{dataset}] index {values.shape} | observed fraction "
          f"{observed.mean():.3f} | masked {sum(len(m) for m, _ in held_out.values())} "
          f"entries over {len(held_out)} genes")
    return index, values, masked, observed, held_out, xy


def compute_gt_morans(values: np.ndarray, xy: np.ndarray) -> np.ndarray:
    gt_morans = np.full(values.shape[0], np.nan)
    for g in range(values.shape[0]):
        v = values[g]
        if np.isfinite(v).sum() < 10:
            continue
        gt_morans[g] = _morans_i(v, xy, k=8, n_perm=0)[0]
    return gt_morans


def morans_i_recovery(pred_full, observed, held_out, gt_morans, xy,
                      restrict_to_truth_support=False, gene_set=None):
    """Copy of the head-to-head metric (+ optional support restriction)."""
    gs = gene_set if gene_set is not None else set(held_out.keys())
    Ig, It = [], []
    for g in gs:
        gt = gt_morans[g]
        if not np.isfinite(gt):
            continue
        v = pred_full[g]
        if restrict_to_truth_support:
            v = v.copy()
            v[~observed[g]] = np.nan
        if np.isfinite(v).sum() < 10:
            continue
        I = _morans_i(v, xy, k=8, n_perm=0)[0]
        if np.isfinite(I):
            Ig.append(I)
            It.append(gt)
    if len(Ig) < 5:
        return float("nan")
    Ig = np.asarray(Ig)
    It = np.asarray(It)
    if np.std(Ig) > 0 and np.std(It) > 0:
        return float(np.corrcoef(Ig, It)[0, 1])
    return float("nan")


def build_topq(values, xy, held_out):
    """Replicate the top-quartile spatially-variable gene set."""
    from sklearn.neighbors import NearestNeighbors
    n_genes, n_spots = values.shape
    _knn = NearestNeighbors(n_neighbors=min(9, n_spots)).fit(xy)
    _nbr = _knn.kneighbors(xy, return_distance=False)[:, 1:]
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
    return set(g for g, _ in scored[: max(1, len(scored) // 4)])


def build_predictions(dataset: str, index, values, masked, held_out, xy):
    """mean / spatial-KNN / GP recomputed; stAPAminer + spvAPA from cache."""
    n_genes, n_spots = values.shape
    preds = {}

    # ---- mean ----
    mn_pred = masked.copy()
    gmeans = np.array([np.nanmean(masked[g]) if np.isfinite(masked[g]).any() else np.nan
                       for g in range(n_genes)])
    for g, (mi, _t) in held_out.items():
        mn_pred[g, mi] = gmeans[g]
    preds["mean"] = mn_pred

    # ---- spatial-KNN ----
    from sklearn.neighbors import NearestNeighbors
    sk_pred = masked.copy()
    nbrs = NearestNeighbors(n_neighbors=min(KNN_K + 1, n_spots)).fit(xy)
    _, idx = nbrs.kneighbors(xy)
    for g, (mi, _t) in held_out.items():
        for s in mi:
            nv = masked[g, idx[s, 1:]]
            nv = nv[np.isfinite(nv)]
            if len(nv):
                sk_pred[g, s] = float(np.mean(nv))
    preds["spatial-KNN"] = sk_pred

    # ---- spaGAPA-GP (same config as the published run) ----
    from scipy.spatial import cKDTree
    from spagapa.imputation import SparseGPImputer
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    t0 = time.time()
    base = SparseGPImputer(
        n_inducing=min(500, max(100, n_spots // 100)),
        length_scale=nn_dist * 5.0, noise_level=0.1,
        length_scale_multiplier=5.0, local_noise=False, noise_scale=1000.0,
    )
    train_mask = np.isfinite(masked)
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    gp_pred, _ = batch.impute(return_uncertainty=True)
    print(f"[{dataset}] GP reproduced in {time.time() - t0:.1f}s")
    preds["spaGAPA-GP"] = gp_pred

    # ---- stAPAminer / spvAPA from the cached R outputs ----
    for method, tmp in (("stAPAminer", "_stapa_tmp"), ("spvAPA", "_spv_tmp")):
        csv = RAW / dataset / tmp / ("stapa_imputed.csv" if method == "stAPAminer"
                                     else "spv_imputed.csv")
        if csv.exists():
            imp = pd.read_csv(csv, index_col=0)
            preds[method] = imp.reindex(index=index.index, columns=index.columns).values
        else:
            print(f"[{dataset}] WARNING: cached {method} output missing ({csv})")
    return preds


def mask_fraction_sweep(dataset, index, values, xy, gt_morans):
    """Show that mean's high recovery is retention, not recovery.

    As the mask fraction grows (less truth retained in pred_full), the
    mean method's morans_i_recovery collapses toward 0/NaN.
    """
    from sklearn.neighbors import NearestNeighbors  # noqa: F401
    rows = []
    observed = np.isfinite(values)
    for mf in (0.2, 0.5, 0.8, 0.95):
        rng = np.random.default_rng(SEED)
        masked = values.copy()
        held_out = {}
        for g in range(values.shape[0]):
            obs_idx = np.where(observed[g])[0]
            if len(obs_idx) < 5:
                continue
            nm = max(1, int(len(obs_idx) * mf))
            mi = rng.choice(obs_idx, size=nm, replace=False)
            masked[g, mi] = np.nan
            held_out[g] = mi
        mn_pred = masked.copy()
        gmeans = np.array([np.nanmean(masked[g]) if np.isfinite(masked[g]).any()
                           else np.nan for g in range(values.shape[0])])
        for g, mi in held_out.items():
            mn_pred[g, mi] = gmeans[g]
        r_raw = morans_i_recovery(mn_pred, observed, held_out, gt_morans, xy)
        rows.append({"dataset": dataset, "method": "mean", "mask_fraction": mf,
                     "retained_truth_fraction": 1.0 - mf,
                     "morans_i_recovery_raw": r_raw})
        print(f"[{dataset}] mean @mask={mf:.2f}: raw recovery = {r_raw:.4f}")
    return rows


def audit_dataset(dataset: str) -> tuple[list[dict], list[dict]]:
    index, values, masked, observed, held_out, xy = load_and_mask(dataset)
    print(f"[{dataset}] computing ground-truth per-gene Moran's I ...")
    t0 = time.time()
    gt_morans = compute_gt_morans(values, xy)
    print(f"[{dataset}] gt_morans in {time.time() - t0:.1f}s "
          f"({np.isfinite(gt_morans).sum()} genes)")
    topq = build_topq(values, xy, held_out)

    with open(RAW / dataset / "results.json") as fh:
        published = json.load(fh)["methods"]

    preds = build_predictions(dataset, index, values, masked, held_out, xy)

    rows, sweep = [], []
    for method, pred in preds.items():
        n_finite = np.isfinite(pred).sum(axis=1).mean()
        # fraction of truth-observed entries that are exactly the truth
        # (how much of pred_full is just retained input)
        exact = np.mean([np.mean(pred[g][observed[g]] == values[g][observed[g]])
                         for g in list(held_out)[:200]])
        for stratum, gs in (("all", None), ("topq", topq)):
            r_pub = published.get(method, {}).get(
                "morans_i_recovery" if stratum == "all" else "morans_i_recovery_topq",
                float("nan"))
            r_raw = morans_i_recovery(pred, observed, held_out, gt_morans, xy, gs)
            r_fix = morans_i_recovery(pred, observed, held_out, gt_morans, xy,
                                      restrict_to_truth_support=True, gene_set=gs)
            rows.append({
                "dataset": dataset, "method": method, "stratum": stratum,
                "mean_support_n_spots": float(n_finite),
                "truth_exact_fraction": float(exact),
                "morans_i_recovery_published": r_pub,
                "morans_i_recovery_reproduced": r_raw,
                "morans_i_recovery_corrected": r_fix,
            })
            print(f"[{dataset}] {method:<12} {stratum:<5} support={n_finite:7.1f} "
                  f"retained={exact:.3f} pub={r_pub:+.4f} repro={r_raw:+.4f} "
                  f"corrected={r_fix:+.4f}")

    # sanity: oracle (truth) must recover perfectly under both definitions
    r_oracle = morans_i_recovery(values, observed, held_out, gt_morans, xy)
    print(f"[{dataset}] oracle (truth vs itself): {r_oracle:+.4f}")

    sweep = mask_fraction_sweep(dataset, index, values, xy, gt_morans)
    return rows, sweep


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows, all_sweep = [], []
    for dataset in ("gse183456", "gse220442"):
        rows, sweep = audit_dataset(dataset)
        all_rows.extend(rows)
        all_sweep.extend(sweep)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / "morans_i_recovery_audit.csv", index=False)
    pd.DataFrame(all_sweep).to_csv(OUT / "mean_mask_fraction_sweep.csv", index=False)
    print(f"\nwrote {OUT / 'morans_i_recovery_audit.csv'}")
    print(f"wrote {OUT / 'mean_mask_fraction_sweep.csv'}")


if __name__ == "__main__":
    main()

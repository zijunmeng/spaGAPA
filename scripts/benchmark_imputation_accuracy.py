#!/usr/bin/env python3
"""
Imputation accuracy benchmark: mask 20% spots, compare GP vs KNN vs mean.

Masks 20% of the observed entries per gene, imputes them with three methods
(spaGAPA sparse GP, KNN, gene mean), and reports RMSE / Pearson / Spearman
plus GP uncertainty calibration (uncertainty-error correlation, 2-sigma
interval coverage).

Usage
-----
python scripts/benchmark_imputation_accuracy.py \
    --apa-matrix   <path/to/apa_matrix.csv>   \
    --coordinates  <path/to/coordinates.csv>  \
    [--output benchmark_results/imputation_accuracy] \
    [--mask-fraction 0.2] [--seed 42]
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
from scipy.stats import pearsonr, spearmanr


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apa-matrix", required=True,
                        help="CSV (n_genes x n_spots) of APA values.")
    parser.add_argument("--coordinates", required=True,
                        help="CSV with columns ['x','y'], one row per spot.")
    parser.add_argument(
        "--output", default="benchmark_results/imputation_accuracy",
        help="Output directory for results JSON.")
    parser.add_argument("--mask-fraction", type=float, default=0.2,
                        help="Fraction of observed entries to mask per gene.")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if not (0.0 < args.mask_fraction < 1.0):
        parser.error("--mask-fraction must be in (0, 1)")

    import pandas as pd

    t_start = time.time()
    apa = pd.read_csv(args.apa_matrix, index_col=0)
    coords = pd.read_csv(args.coordinates, index_col=0)
    values = apa.values  # (n_genes, n_spots)

    # Resolve coordinate columns robustly: accept either ['x','y'] or the
    # first two numeric columns.
    if "x" in coords.columns and "y" in coords.columns:
        xy = coords[["x", "y"]].values
    else:
        xy = coords.iloc[:, :2].values.astype(float)

    # Align spots: apa columns vs coordinates index.
    if apa.shape[1] != xy.shape[0]:
        # Best-effort alignment by index when shapes disagree.
        common = apa.columns.intersection(coords.index)
        if len(common) > 0:
            apa = apa.loc[:, common]
            coords = coords.loc[common]
            xy = coords[["x", "y"]].values if {"x", "y"} <= set(coords.columns) \
                else coords.iloc[:, :2].values.astype(float)
            values = apa.values
    if values.shape[1] != xy.shape[0]:
        raise ValueError(
            f"Spot count mismatch: apa_matrix has {values.shape[1]} spots, "
            f"coordinates has {xy.shape[0]}"
        )

    print(
        f"Loaded APA matrix: {values.shape[0]} genes x {values.shape[1]} spots, "
        f"observed fraction = {(values > 0).mean():.3f}"
    )

    rng = np.random.default_rng(args.seed)
    observed = values > 0  # mask of observed entries (n_genes, n_spots)

    # For each gene, mask 20% (default) of its observed spots.
    masked_values = values.copy()
    held_out_truth = {}
    n_masked_total = 0
    for g in range(values.shape[0]):
        obs_idx = np.where(observed[g])[0]
        if len(obs_idx) < 5:
            continue
        n_mask = max(1, int(len(obs_idx) * args.mask_fraction))
        mask_idx = rng.choice(obs_idx, size=n_mask, replace=False)
        masked_values[g, mask_idx] = 0.0
        held_out_truth[g] = (mask_idx, values[g, mask_idx])
        n_masked_total += n_mask
    print(f"Masked {n_masked_total} entries across {len(held_out_truth)} genes "
          f"(fraction={args.mask_fraction})")

    results: dict = {}

    # ------------------------------------------------------------------
    # Method 1: GP imputation (spaGAPA sparse GP)
    # ------------------------------------------------------------------
    from spagapa.imputation import SparseGPImputer

    n_spots = values.shape[1]
    # Length scale = (median nearest-neighbor distance) x 5.
    #
    # Two correctness fixes vs. the original formulation:
    #
    # (1) NEAREST-NEIGHBOR DISTANCE. The previous line
    #       nn_dist = np.median(np.linalg.norm(xy[1:] - xy[0]))
    #     did NOT compute the median NN distance -- it computed the median
    #     distance from spot 0 to every other spot, which for a tissue-scale
    #     coordinate layout (values up to ~1e5 here) was ~2.6e6. Multiplying
    #     by 10 gave length_scale ~2.6e7, so the RBF kernel exp(-d^2/2l^2)
    #     evaluated to ~1.0 for every pair of points -> every test point saw
    #     the same inducing-point interpolation -> constant prediction ->
    #     Pearson/Spearman = NaN.
    #     We now use a true k-NN distance (k=2 excludes self) via scipy's
    #     cKDTree, which on this bin-200 grid is exactly 200.
    #
    # (2) TRAINING MASK. We previously passed ``mask=observed`` (the ORIGINAL
    #     pre-masking observed mask) to ``fit_batch``. ``SparseGPImputerBatch.
    #     impute()`` then ran ``predictions[observed] = self.values[observed]``
    #     where ``self.values`` is the MASKED matrix -- so every held-out
    #     (masked) spot still satisfied ``observed[spot] == True`` AND had
    #     ``self.values[spot] == 0``, which overwrote the GP prediction with 0.
    #     Result: every held-out prediction was 0 -> constant -> NaN.
    #     Fix: pass the POST-masking training mask ``masked_values > 0``.
    from scipy.spatial import cKDTree
    _kdt = cKDTree(xy)
    _nn_dists, _ = _kdt.query(xy, k=2)
    nn_dist = float(np.median(_nn_dists[:, 1]))
    length_scale = nn_dist * 5
    print(
        f"  GP length_scale = {length_scale:.1f} "
        f"(median NN dist = {nn_dist:.1f}, multiplier = 5)"
    )
    t0 = time.time()
    base = SparseGPImputer(
        n_inducing=min(500, max(100, n_spots // 100)),
        length_scale=length_scale,
        noise_level=0.1,
    )
    # Post-masking training mask: only spots whose value survives masking
    # are training examples; held-out spots must NOT be overwritten in
    # ``SparseGPImputerBatch.impute()``.
    train_mask = masked_values > 0
    batch = base.fit_batch(xy, masked_values, mask=train_mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=True)
    results["gp_time"] = time.time() - t0
    print(f"  GP imputation: {results['gp_time']:.1f}s")

    # ------------------------------------------------------------------
    # Method 2: KNN imputation (stAPAminer-style spatial neighbors)
    # ------------------------------------------------------------------
    from sklearn.neighbors import NearestNeighbors

    t0 = time.time()
    knn_pred = masked_values.copy()
    nn = NearestNeighbors(n_neighbors=min(16, n_spots)).fit(xy)
    _, indices = nn.kneighbors(xy)
    for g in range(values.shape[0]):
        # Spots we need to impute: masked entries of this gene.
        if g not in held_out_truth:
            continue
        mask_idx, _ = held_out_truth[g]
        for s in mask_idx:
            neighbors = indices[s, 1:]  # exclude self
            neighbor_vals = masked_values[g, neighbors]
            neighbor_obs = neighbor_vals[neighbor_vals > 0]
            if len(neighbor_obs) > 0:
                knn_pred[g, s] = float(np.mean(neighbor_obs))
    results["knn_time"] = time.time() - t0
    print(f"  KNN imputation: {results['knn_time']:.1f}s")

    # ------------------------------------------------------------------
    # Method 3: Mean imputation (per-gene mean of remaining observed)
    # ------------------------------------------------------------------
    t0 = time.time()
    mean_pred = masked_values.copy()
    gene_means = np.array([
        masked_values[g][masked_values[g] > 0].mean()
        if (masked_values[g] > 0).any() else 0.0
        for g in range(values.shape[0])
    ])
    for g, (mask_idx, _truth) in held_out_truth.items():
        mean_pred[g, mask_idx] = gene_means[g]
    results["mean_time"] = time.time() - t0
    print(f"  Mean imputation: {results['mean_time']:.1f}s")

    # ------------------------------------------------------------------
    # Evaluate: gather (truth, pred) pairs over all masked entries.
    # ------------------------------------------------------------------
    def evaluate(pred: np.ndarray) -> dict:
        truths, preds = [], []
        for g, (mask_idx, truth) in held_out_truth.items():
            truths.extend(truth.tolist())
            preds.extend(pred[g, mask_idx].tolist())
        truths = np.asarray(truths, dtype=float)
        preds = np.asarray(preds, dtype=float)
        rmse = float(np.sqrt(np.mean((truths - preds) ** 2)))
        if np.std(truths) > 0 and np.std(preds) > 0:
            r, _ = pearsonr(truths, preds)
            rho, _ = spearmanr(truths, preds)
        else:
            r = rho = float("nan")
        return {"rmse": rmse, "pearson": float(r), "spearman": float(rho)}

    results["gp"] = evaluate(gp_pred)
    results["knn"] = evaluate(knn_pred)
    results["mean"] = evaluate(mean_pred)

    # ------------------------------------------------------------------
    # Uncertainty calibration (GP only)
    # ------------------------------------------------------------------
    gp_truths, gp_preds, gp_uncs = [], [], []
    for g, (mask_idx, truth) in held_out_truth.items():
        gp_truths.extend(truth.tolist())
        gp_preds.extend(gp_pred[g, mask_idx].tolist())
        gp_uncs.extend(gp_unc[g, mask_idx].tolist())
    gp_truths = np.asarray(gp_truths, dtype=float)
    gp_preds = np.asarray(gp_preds, dtype=float)
    gp_uncs = np.asarray(gp_uncs, dtype=float)
    errors = np.abs(gp_truths - gp_preds)
    if np.std(gp_uncs) > 0 and np.std(errors) > 0:
        unc_corr = float(np.corrcoef(gp_uncs, errors)[0, 1])
    else:
        unc_corr = float("nan")
    coverage = float(np.mean(np.abs(gp_truths - gp_preds) <= 2 * gp_uncs))
    results["gp_uncertainty"] = {
        "unc_error_corr": unc_corr,
        "interval_coverage_2sigma": coverage,
    }

    results["config"] = {
        "apa_matrix": os.path.abspath(args.apa_matrix),
        "coordinates": os.path.abspath(args.coordinates),
        "n_genes": int(values.shape[0]),
        "n_spots": int(values.shape[1]),
        "observed_fraction": float(observed.mean()),
        "mask_fraction": float(args.mask_fraction),
        "seed": int(args.seed),
        "n_masked_entries": int(n_masked_total),
        "gp_length_scale": float(length_scale),
        "gp_nn_distance": float(nn_dist),
    }
    results["total_time"] = time.time() - t_start

    os.makedirs(args.output, exist_ok=True)
    out_path = os.path.join(args.output, "imputation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n=== Imputation Accuracy Benchmark ===")
    print(f"  {'Method':<8} {'RMSE':<10} {'Pearson':<10} {'Spearman':<10} {'Time':<10}")
    for method in ("gp", "knn", "mean"):
        r = results[method]
        print(
            f"  {method:<8} {r['rmse']:<10.4f} {r['pearson']:<10.4f} "
            f"{r['spearman']:<10.4f} {results[method + '_time']:<10.1f}s"
        )
    print(f"\n  GP uncertainty-error correlation: "
          f"{results['gp_uncertainty']['unc_error_corr']:.4f}")
    print(f"  GP 2-sigma interval coverage:     "
          f"{results['gp_uncertainty']['interval_coverage_2sigma']:.4f}")
    print(f"\n  Results written to: {out_path}")
    print(f"  Total wall time: {results['total_time']:.1f}s")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Emit parameter_table.csv for the benchmark fairness statement.

One row per method, columns = every tunable parameter (or 'none' for the
parameter-free mean baseline).  Values are the EXACT settings used in the
head-to-head re-run (scripts/benchmark_stapaminer_headtohead.py, seed 42,
mask 0.2).
"""
import csv
from pathlib import Path

OUT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/benchmark_fairness")
OUT.mkdir(parents=True, exist_ok=True)

COLS = [
    "method", "language", "neighbour_basis", "k",
    "n_inducing", "length_scale", "noise", "normalization",
    "other_params", "n_free_params",
]

ROWS = [
    {
        "method": "spaGAPA-GP",
        "language": "Python 3.10 (spagapa env)",
        "neighbour_basis": "sparse GP over spatial coords (RBF + WhiteKernel)",
        "k": "n/a (uses inducing points, not kNN)",
        "n_inducing": "min(500, max(100, n_spots//100))",
        "length_scale": "median_NN_dist x 5 (length_scale_multiplier=5.0)",
        "noise": "0.1 (flat noise_level; local_noise=False)",
        "normalization": "none on APA index (already in [0,1])",
        "other_params": "noise_scale=1000 (only used if local_noise=True; off here)",
        "n_free_params": "3 (n_inducing, length_scale_multiplier, noise_level)",
    },
    {
        "method": "stAPAminer",
        "language": "R 4.4.2 (r442 env)",
        "neighbour_basis": "expression-KNN (imputeAPAIndex, kNN on RNA expression)",
        "k": "10 (stAPAminer README/paper default)",
        "n_inducing": "n/a",
        "length_scale": "n/a",
        "noise": "n/a (no probabilistic model)",
        "normalization": "TenX=TRUE (imputeAPAIndex internal)",
        "other_params": "init=TRUE (zero-init where matching gene has zero expression)",
        "n_free_params": "1 (k)",
    },
    {
        "method": "spvAPA",
        "language": "R 4.4.2 (r442 env)",
        "neighbour_basis": "multimodal WNN (Seurat WNN graph over RNA+APA assays)",
        "k": "15 (spvAPA WNN default; WNNImpute native default is 20)",
        "n_inducing": "n/a",
        "length_scale": "n/a",
        "noise": "n/a (no probabilistic model)",
        "normalization": "SCTransform + PCA (default Seurat RNA workflow inside WNNImpute)",
        "other_params": "init=TRUE; is.weight=FALSE (returns imputed matrix, not weights)",
        "n_free_params": "1 (k)",
    },
    {
        "method": "spatial-KNN",
        "language": "Python 3.10 (spagapa env, sklearn NearestNeighbors)",
        "neighbour_basis": "spatial kNN (on x,y coordinates only)",
        "k": "15",
        "n_inducing": "n/a",
        "length_scale": "n/a",
        "noise": "n/a",
        "normalization": "none",
        "other_params": "impute = mean of k nearest observed spots (nan-aware)",
        "n_free_params": "1 (k)",
    },
    {
        "method": "mean",
        "language": "Python 3.10 (numpy)",
        "neighbour_basis": "none (per-gene constant = mean of observed entries)",
        "k": "n/a",
        "n_inducing": "n/a",
        "length_scale": "n/a",
        "noise": "n/a",
        "normalization": "none",
        "other_params": "prediction[gene] = nanmean(observed_non_heldout[gene])",
        "n_free_params": "0",
    },
]

with open(OUT / "parameter_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    for r in ROWS:
        w.writerow(r)
print(f"wrote {OUT / 'parameter_table.csv'} with {len(ROWS)} methods")

#!/usr/bin/env python
"""Formal real-data benchmark on the prepared stAPAminer MOB dataset.

This benchmark upgrades the earlier smoke test into a multi-seed, multi-mask,
publication-oriented benchmark with saved tables and Matplotlib figures.

Benchmark dimensions
--------------------
1. More genes
2. Multiple random seeds
3. Multiple mask types
   - random
   - spatial_block
   - spatial_block_large
   - ring_sector
   - layer_aware
   - low_coverage
4. Uncertainty calibration for GP
5. Runtime and memory usage
6. Biological consistency using MOB layer labels
7. Saved result tables and figures
"""

from __future__ import annotations

import argparse
import json
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Tuple

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from scipy.stats import pearsonr, spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    mean_absolute_error,
    mean_squared_error,
    normalized_mutual_info_score,
    r2_score,
)
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler

from spagapa.imputation import (
    ExpressionFeatureBuilder,
    ExpressionGPImputer,
    GPImputer,
    GeometryFeatureBuilder,
)


@dataclass
class MethodResult:
    method: str
    gp_kernel: str
    gp_alpha: float
    gp_variant: str
    gp_lambda_expr: float
    gp_product_offset: float
    gp_gate_mode: str
    gp_gate_tau: float
    gp_local_k: int
    gp_layer_gate_mode: str
    expr_n_components: float
    gp_use_theta: bool
    seed: int
    mask_type: str
    rmse: float
    mae: float
    pearson: float
    spearman: float
    r2: float
    n_holdout: int
    runtime_s: float
    peak_rss_mb: float
    layer_ari: float
    layer_nmi: float
    n_domains: int
    uncertainty_coverage_68: float = np.nan
    uncertainty_coverage_95: float = np.nan
    uncertainty_nll: float = np.nan
    uncertainty_error_pearson: float = np.nan
    uncertainty_error_spearman: float = np.nan
    rmse_keep_80pct: float = np.nan
    rmse_keep_50pct: float = np.nan
    rmse_gain_keep_80pct: float = np.nan


METHOD_ORDER = [
    "spagapa_gp",
    "stapaminer_knn_expression",
    "knn_spatial",
    "mean",
    "median",
]

METHOD_COLORS = {
    "spagapa_gp": "#e74c3c",
    "stapaminer_knn_expression": "#8e44ad",
    "knn_spatial": "#3498db",
    "mean": "#95a5a6",
    "median": "#bdc3c7",
}

MASK_ORDER = [
    "random",
    "spatial_block",
    "spatial_block_large",
    "ring_sector",
    "layer_aware",
    "low_coverage",
]

MASK_COLORS = {
    "random": "#2ecc71",
    "spatial_block": "#f39c12",
    "spatial_block_large": "#d35400",
    "ring_sector": "#1abc9c",
    "layer_aware": "#9b59b6",
    "low_coverage": "#34495e",
}

FILTER_RETAIN_FRACTIONS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared stAPAminer MOB directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/stapaminer_mob_formal",
        help="Output directory.",
    )
    parser.add_argument("--n-genes", type=int, default=200)
    parser.add_argument("--min-observed-spots", type=int, default=100)
    parser.add_argument("--mask-fraction", type=float, default=0.2)
    parser.add_argument(
        "--seeds",
        default="42,43,44",
        help="Comma-separated random seeds.",
    )
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--gp-kernel", default="matern")
    parser.add_argument("--gp-alpha", type=float, default=1e-10)
    parser.add_argument(
        "--gp-variant",
        default="spatial",
        choices=[
            "spatial",
            "spatial_radial",
            "expr_additive",
            "expr_product",
            "adaptive_additive",
            "expr_layer_local",
        ],
        help="spaGAPA GP variant to benchmark.",
    )
    parser.add_argument(
        "--gp-lambda-expr",
        type=float,
        default=0.5,
        help="Expression kernel weight for expression-informed GP variants.",
    )
    parser.add_argument(
        "--gp-product-offset",
        type=float,
        default=1.0,
        help="Baseline offset in product kernel: K_space * (offset + lambda * K_expr).",
    )
    parser.add_argument(
        "--gp-gate-mode",
        default="distance",
        choices=["distance"],
        help="Adaptive gating mode for expression-aware GP variants.",
    )
    parser.add_argument(
        "--gp-gate-tau",
        type=float,
        default=None,
        help="Spatial gating length scale for adaptive_additive GP. Defaults to spatial length scale.",
    )
    parser.add_argument(
        "--gp-local-k",
        type=int,
        default=20,
        help="Neighborhood size used to estimate layer-local soft gate scale.",
    )
    parser.add_argument(
        "--gp-layer-gate-mode",
        default="radius",
        choices=["radius", "pseudolayer"],
        help="Coordinate-derived local gate feature for expr_layer_local.",
    )
    parser.add_argument(
        "--expr-n-components",
        type=int,
        default=10,
        help="Number of PCA components for expression embedding.",
    )
    parser.add_argument(
        "--expr-use-hvg",
        action="store_true",
        help="Restrict expression embedding to highly variable genes before PCA.",
    )
    parser.add_argument(
        "--expr-n-top-genes",
        type=int,
        default=2000,
        help="Number of HVGs to keep when --expr-use-hvg is set.",
    )
    parser.add_argument(
        "--gp-use-theta",
        action="store_true",
        help="Include polar angle theta in spatial_radial geometry features.",
    )
    parser.add_argument("--n-domains", type=int, default=5)
    parser.add_argument(
        "--mask-types",
        default="random,spatial_block,spatial_block_large,ring_sector,layer_aware,low_coverage",
        help="Comma-separated mask types.",
    )
    parser.add_argument(
        "--calibration-bins",
        type=int,
        default=8,
        help="Number of uncertainty bins for GP reliability curves.",
    )
    return parser.parse_args()


def load_inputs(data_dir: Path):
    apa = pd.read_csv(data_dir / "apa_matrix.csv", index_col=0)
    expression = pd.read_csv(data_dir / "expression_matrix.csv", index_col=0)
    coords = pd.read_csv(data_dir / "coordinates.csv").set_index("spot_id")
    metadata = pd.read_csv(data_dir / "metadata.csv")
    metadata = metadata.set_index("spot_id").loc[coords.index]
    return apa, expression, coords, metadata


def select_genes(apa: pd.DataFrame, n_genes: int, min_observed_spots: int):
    observed = apa.notna().sum(axis=1)
    has_missing = observed < apa.shape[1]
    eligible = observed[(observed >= min_observed_spots) & has_missing]
    if len(eligible) < n_genes:
        eligible = observed[observed >= min_observed_spots]
    selected = eligible.sort_values(ascending=False).head(n_genes).index
    return apa.loc[selected]


def make_random_mask(values: np.ndarray, mask_fraction: float, seed: int):
    rng = np.random.default_rng(seed)
    observed_idx = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed_idx) * mask_fraction))
    return sample_pairs_to_mask(values.shape, observed_idx, n_mask, rng)


def sample_pairs_to_mask(
    shape: tuple[int, int],
    candidate_pairs: np.ndarray,
    n_mask: int,
    rng: np.random.Generator,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    if len(candidate_pairs) == 0:
        raise ValueError("No candidate pairs available for masking")

    n_select = min(n_mask, len(candidate_pairs))
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
        if weights.sum() == 0:
            weights = None
        else:
            weights = weights / weights.sum()

    chosen_idx = rng.choice(
        len(candidate_pairs),
        size=n_select,
        replace=False,
        p=weights,
    )
    chosen = candidate_pairs[chosen_idx]
    mask = np.zeros(shape, dtype=bool)
    mask[chosen[:, 0], chosen[:, 1]] = True
    return mask


def observed_pairs_from_spots(values: np.ndarray, spot_indices: np.ndarray) -> np.ndarray:
    if len(spot_indices) == 0:
        return np.empty((0, 2), dtype=int)

    observed = np.isfinite(values[:, spot_indices])
    gene_idx, rel_spot_idx = np.where(observed)
    return np.column_stack([gene_idx, spot_indices[rel_spot_idx]]).astype(int)


def collect_nearest_spots(
    values: np.ndarray,
    ordered_spots: np.ndarray,
    target_pairs: int,
) -> np.ndarray:
    chosen_spots: list[int] = []
    n_pairs = 0
    for spot_idx in ordered_spots:
        obs_here = int(np.isfinite(values[:, spot_idx]).sum())
        if obs_here == 0:
            continue
        chosen_spots.append(int(spot_idx))
        n_pairs += obs_here
        if n_pairs >= target_pairs:
            break
    return np.asarray(chosen_spots, dtype=int)


def make_spatial_block_mask(
    values: np.ndarray,
    coords: np.ndarray,
    mask_fraction: float,
    seed: int,
    region_multiplier: float = 1.0,
):
    rng = np.random.default_rng(seed)
    observed_idx = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed_idx) * mask_fraction))

    spot_weights = np.isfinite(values).sum(axis=0).astype(float)
    spot_weights = spot_weights / spot_weights.sum()
    center_spot = rng.choice(values.shape[1], p=spot_weights)
    dists = np.linalg.norm(coords - coords[center_spot], axis=1)
    order = np.argsort(dists).astype(int)
    target_pairs = max(n_mask, int(n_mask * region_multiplier))
    candidate_spots = collect_nearest_spots(values, order, target_pairs)
    candidate_pairs = observed_pairs_from_spots(values, candidate_spots)
    return sample_pairs_to_mask(values.shape, candidate_pairs, n_mask, rng)


def make_low_coverage_mask(values: np.ndarray, mask_fraction: float, seed: int):
    rng = np.random.default_rng(seed)
    observed_idx = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed_idx) * mask_fraction))

    observed_per_gene = np.isfinite(values).sum(axis=1).astype(float)
    gene_weights = np.zeros(values.shape[0], dtype=float)
    valid_genes = observed_per_gene > 0
    gene_weights[valid_genes] = 1.0 / observed_per_gene[valid_genes]

    pair_weights = gene_weights[observed_idx[:, 0]]
    return sample_pairs_to_mask(values.shape, observed_idx, n_mask, rng, weights=pair_weights)


def make_ring_sector_mask(
    values: np.ndarray,
    coords: np.ndarray,
    mask_fraction: float,
    seed: int,
):
    rng = np.random.default_rng(seed)
    observed_idx = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed_idx) * mask_fraction))

    center = coords.mean(axis=0)
    deltas = coords - center
    radii = np.linalg.norm(deltas, axis=1)
    angles = np.arctan2(deltas[:, 1], deltas[:, 0])
    theta0 = rng.uniform(-np.pi, np.pi)
    angle_diff = np.abs(((angles - theta0 + np.pi) % (2 * np.pi)) - np.pi)

    radial_bands = [(0.20, 0.80), (0.10, 0.90), (0.0, 1.0)]
    angle_widths = [np.pi / 2, 2 * np.pi / 3, np.pi]

    candidate_spots = np.array([], dtype=int)
    for q_low, q_high in radial_bands:
        r_low, r_high = np.quantile(radii, [q_low, q_high])
        for width in angle_widths:
            candidate_spots = np.where(
                (radii >= r_low) & (radii <= r_high) & (angle_diff <= width / 2)
            )[0]
            if np.isfinite(values[:, candidate_spots]).sum() >= n_mask:
                break
        if np.isfinite(values[:, candidate_spots]).sum() >= n_mask:
            break

    if len(candidate_spots) == 0:
        return make_spatial_block_mask(values, coords, mask_fraction, seed)

    candidate_pairs = observed_pairs_from_spots(values, candidate_spots.astype(int))
    return sample_pairs_to_mask(values.shape, candidate_pairs, n_mask, rng)


def infer_layer_order(layer_labels: np.ndarray, coords: np.ndarray) -> list[str]:
    center = coords.mean(axis=0)
    radii = np.linalg.norm(coords - center, axis=1)
    layer_df = pd.DataFrame({"layer": layer_labels.astype(str), "radius": radii})
    return (
        layer_df.groupby("layer")["radius"]
        .mean()
        .sort_values()
        .index.astype(str)
        .tolist()
    )


def make_layer_aware_mask(
    values: np.ndarray,
    coords: np.ndarray,
    layer_labels: np.ndarray,
    mask_fraction: float,
    seed: int,
):
    rng = np.random.default_rng(seed)
    observed_idx = np.argwhere(np.isfinite(values))
    n_mask = max(1, int(len(observed_idx) * mask_fraction))

    layer_order = infer_layer_order(layer_labels, coords)
    layer_to_spots = {
        layer: np.where(layer_labels == layer)[0]
        for layer in layer_order
    }
    layer_pair_counts = {
        layer: int(np.isfinite(values[:, spots]).sum())
        for layer, spots in layer_to_spots.items()
    }

    available_layers = [layer for layer in layer_order if layer_pair_counts[layer] > 0]
    if not available_layers:
        return make_random_mask(values, mask_fraction, seed)

    layer_weights = np.array([layer_pair_counts[layer] for layer in available_layers], dtype=float)
    layer_weights = layer_weights / layer_weights.sum()
    target_layer = str(rng.choice(available_layers, p=layer_weights))

    target_idx = layer_order.index(target_layer)
    chosen_layers = [target_layer]
    left = target_idx - 1
    right = target_idx + 1

    def current_pair_count() -> int:
        spots = np.where(np.isin(layer_labels, chosen_layers))[0]
        return int(np.isfinite(values[:, spots]).sum())

    while current_pair_count() < n_mask and (left >= 0 or right < len(layer_order)):
        candidates = []
        if left >= 0:
            candidates.append(layer_order[left])
        if right < len(layer_order):
            candidates.append(layer_order[right])
        if not candidates:
            break
        next_layer = max(candidates, key=lambda layer: layer_pair_counts.get(layer, 0))
        chosen_layers.append(next_layer)
        if left >= 0 and layer_order[left] == next_layer:
            left -= 1
        elif right < len(layer_order) and layer_order[right] == next_layer:
            right += 1

    candidate_spots = np.where(np.isin(layer_labels, chosen_layers))[0]
    candidate_pairs = observed_pairs_from_spots(values, candidate_spots.astype(int))
    return sample_pairs_to_mask(values.shape, candidate_pairs, n_mask, rng)


def build_mask(
    mask_type: str,
    values: np.ndarray,
    coords: np.ndarray,
    layer_labels: np.ndarray,
    mask_fraction: float,
    seed: int,
):
    if mask_type == "random":
        return make_random_mask(values, mask_fraction, seed)
    if mask_type == "spatial_block":
        return make_spatial_block_mask(values, coords, mask_fraction, seed)
    if mask_type == "spatial_block_large":
        return make_spatial_block_mask(values, coords, mask_fraction, seed, region_multiplier=2.0)
    if mask_type == "ring_sector":
        return make_ring_sector_mask(values, coords, mask_fraction, seed)
    if mask_type == "layer_aware":
        return make_layer_aware_mask(values, coords, layer_labels, mask_fraction, seed)
    if mask_type == "low_coverage":
        return make_low_coverage_mask(values, mask_fraction, seed)
    raise ValueError(f"Unknown mask_type: {mask_type}")


def make_train_matrix(values: np.ndarray, holdout_mask: np.ndarray):
    train = values.copy()
    train[holdout_mask] = np.nan
    return train


def impute_mean(train: np.ndarray):
    filled = train.copy()
    gene_means = safe_row_nanmean(filled)
    missing_gene, missing_spot = np.where(~np.isfinite(filled))
    filled[missing_gene, missing_spot] = gene_means[missing_gene]
    return np.clip(filled, 0.0, 1.0)


def impute_median(train: np.ndarray):
    filled = train.copy()
    gene_medians = safe_row_nanmedian(filled)
    missing_gene, missing_spot = np.where(~np.isfinite(filled))
    filled[missing_gene, missing_spot] = gene_medians[missing_gene]
    return np.clip(filled, 0.0, 1.0)


def impute_spatial_knn(train: np.ndarray, coords: np.ndarray, k: int):
    dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    np.fill_diagonal(dist, np.inf)
    neighbors = np.argsort(dist, axis=1)[:, : min(k, dist.shape[0] - 1)]
    filled = train.copy()
    for _ in range(10):
        if np.isfinite(filled).all():
            break
        old_missing = np.isnan(filled).sum()
        for spot_idx in range(filled.shape[1]):
            target = filled[:, spot_idx]
            missing = ~np.isfinite(target)
            if not missing.any():
                continue
            source = filled[:, neighbors[spot_idx]]
            means = safe_row_nanmean(source)
            target[missing] = means[missing]
            filled[:, spot_idx] = target
        if np.isnan(filled).sum() == old_missing:
            break
    gene_means = safe_row_nanmean(filled)
    missing_gene, missing_spot = np.where(~np.isfinite(filled))
    filled[missing_gene, missing_spot] = gene_means[missing_gene]
    return np.clip(filled, 0.0, 1.0)


def impute_stapaminer_expression_knn(
    train: np.ndarray,
    gene_names: List[str],
    spot_names: List[str],
    expression: pd.DataFrame,
    k: int,
):
    expr_all = expression.reindex(columns=spot_names).fillna(0.0)
    scaled = expr_all.T.astype(float)
    scaled = (scaled - scaled.mean(axis=0)) / scaled.std(axis=0).replace(0, np.nan)
    scaled = scaled.fillna(0.0)

    dist = pairwise_distances(scaled.values, metric="euclidean")
    np.fill_diagonal(dist, np.inf)
    neighbors = np.argsort(dist, axis=1)[:, : min(k, dist.shape[0] - 1)]

    filled = train.copy()
    expr_subset = expression.reindex(index=gene_names, columns=spot_names)
    init_zero = (expr_subset.values == 0) & ~np.isfinite(filled)
    filled[init_zero] = 0.0

    for _ in range(10):
        if np.isfinite(filled).all():
            break
        old_missing = np.isnan(filled).sum()
        for spot_idx in range(filled.shape[1]):
            target = filled[:, spot_idx]
            missing = ~np.isfinite(target)
            if not missing.any():
                continue
            source = filled[:, neighbors[spot_idx]]
            means = safe_row_nanmean(source)
            target[missing] = means[missing]
            filled[:, spot_idx] = target
        if np.isnan(filled).sum() == old_missing:
            break

    filled[~np.isfinite(filled)] = 0.0
    return np.clip(filled, 0.0, 1.0)


def impute_spagapa_gp(train: np.ndarray, coords: np.ndarray, kernel: str, alpha: float):
    imputer = GPImputer(kernel_type=kernel, alpha=alpha, n_restarts_optimizer=1)
    batch = imputer.fit_batch(
        coords,
        np.nan_to_num(train, nan=0.0),
        mask=np.isfinite(train),
        n_jobs=1,
        verbose=False,
    )
    imputed, uncertainty = batch.impute(return_uncertainty=True)
    return np.clip(imputed, 0.0, 1.0), uncertainty


def impute_spagapa_gp_variant(
    train: np.ndarray,
    coords: np.ndarray,
    kernel: str,
    alpha: float,
    variant: str,
    expression_embedding: np.ndarray | None = None,
    geometry_features: dict[str, np.ndarray] | None = None,
    lambda_expr: float = 0.5,
    product_offset: float = 1.0,
    gate_mode: str = "distance",
    gate_tau: float | None = None,
    local_k: int = 20,
    layer_gate_mode: str = "radius",
    use_theta: bool = False,
):
    if variant == "spatial":
        return impute_spagapa_gp(train, coords, kernel, alpha)

    if variant == "spatial_radial":
        imputer = ExpressionGPImputer(
            variant="spatial_radial",
            spatial_kernel_type=kernel,
            alpha=alpha,
            use_theta=use_theta,
        )
        batch = imputer.fit_batch(
            coords,
            np.nan_to_num(train, nan=0.0),
            mask=np.isfinite(train),
            geometry_features=geometry_features,
            n_jobs=1,
            verbose=False,
        )
        imputed, uncertainty = batch.impute(return_uncertainty=True)
        return np.clip(imputed, 0.0, 1.0), uncertainty

    if variant == "expr_additive":
        imputer = ExpressionGPImputer(
            variant="additive",
            spatial_kernel_type=kernel,
            expression_kernel_type="rbf",
            alpha=alpha,
            lambda_expr=lambda_expr,
            use_theta=use_theta,
        )
        batch = imputer.fit_batch(
            coords,
            np.nan_to_num(train, nan=0.0),
            mask=np.isfinite(train),
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            n_jobs=1,
            verbose=False,
        )
        imputed, uncertainty = batch.impute(return_uncertainty=True)
        return np.clip(imputed, 0.0, 1.0), uncertainty

    if variant == "expr_product":
        imputer = ExpressionGPImputer(
            variant="product",
            spatial_kernel_type=kernel,
            expression_kernel_type="rbf",
            alpha=alpha,
            lambda_expr=lambda_expr,
            product_offset=product_offset,
            use_theta=use_theta,
        )
        batch = imputer.fit_batch(
            coords,
            np.nan_to_num(train, nan=0.0),
            mask=np.isfinite(train),
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            n_jobs=1,
            verbose=False,
        )
        imputed, uncertainty = batch.impute(return_uncertainty=True)
        return np.clip(imputed, 0.0, 1.0), uncertainty

    if variant == "adaptive_additive":
        imputer = ExpressionGPImputer(
            variant="adaptive_additive",
            spatial_kernel_type=kernel,
            expression_kernel_type="rbf",
            alpha=alpha,
            lambda_expr=lambda_expr,
            gate_mode=gate_mode,
            gate_tau=gate_tau,
            use_theta=use_theta,
        )
        batch = imputer.fit_batch(
            coords,
            np.nan_to_num(train, nan=0.0),
            mask=np.isfinite(train),
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            n_jobs=1,
            verbose=False,
        )
        imputed, uncertainty = batch.impute(return_uncertainty=True)
        return np.clip(imputed, 0.0, 1.0), uncertainty

    if variant == "expr_layer_local":
        imputer = ExpressionGPImputer(
            variant="layer_local",
            spatial_kernel_type=kernel,
            expression_kernel_type="rbf",
            alpha=alpha,
            lambda_expr=lambda_expr,
            product_offset=product_offset,
            gate_tau=gate_tau,
            local_k=local_k,
            layer_gate_mode=layer_gate_mode,
            use_theta=use_theta,
        )
        batch = imputer.fit_batch(
            coords,
            np.nan_to_num(train, nan=0.0),
            mask=np.isfinite(train),
            expression_embedding=expression_embedding,
            geometry_features=geometry_features,
            n_jobs=1,
            verbose=False,
        )
        imputed, uncertainty = batch.impute(return_uncertainty=True)
        return np.clip(imputed, 0.0, 1.0), uncertainty

    raise ValueError(f"Unknown gp variant: {variant}")


def safe_row_nanmean(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    sums = np.where(finite, values, 0.0).sum(axis=1)
    means = np.zeros(values.shape[0], dtype=float)
    valid = counts > 0
    means[valid] = sums[valid] / counts[valid]
    return means


def safe_row_nanmedian(values: np.ndarray) -> np.ndarray:
    medians = np.zeros(values.shape[0], dtype=float)
    for i, row in enumerate(values):
        finite = row[np.isfinite(row)]
        medians[i] = float(np.median(finite)) if finite.size else 0.0
    return medians


def compute_metrics(true: np.ndarray, pred: np.ndarray, holdout: np.ndarray):
    y_true = true[holdout]
    y_pred = pred[holdout]
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    mae = float(mean_absolute_error(y_true, y_pred))
    if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1:
        pearson = float(pearsonr(y_true, y_pred).statistic)
        spearman = float(spearmanr(y_true, y_pred).statistic)
    else:
        pearson = np.nan
        spearman = np.nan
    try:
        r2 = float(r2_score(y_true, y_pred))
    except Exception:
        r2 = np.nan
    return rmse, mae, pearson, spearman, r2


def safe_stat_correlation(x: np.ndarray, y: np.ndarray, kind: str) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3:
        return np.nan
    x_valid = x[valid]
    y_valid = y[valid]
    if len(np.unique(x_valid)) < 2 or len(np.unique(y_valid)) < 2:
        return np.nan
    if kind == "pearson":
        return float(pearsonr(x_valid, y_valid).statistic)
    if kind == "spearman":
        return float(spearmanr(x_valid, y_valid).statistic)
    raise ValueError(f"Unknown correlation kind: {kind}")


def compute_uncertainty_summary(
    true: np.ndarray,
    pred: np.ndarray,
    uncertainty: np.ndarray | None,
    holdout: np.ndarray,
) -> dict[str, float]:
    summary = {
        "uncertainty_coverage_68": np.nan,
        "uncertainty_coverage_95": np.nan,
        "uncertainty_nll": np.nan,
        "uncertainty_error_pearson": np.nan,
        "uncertainty_error_spearman": np.nan,
        "rmse_keep_80pct": np.nan,
        "rmse_keep_50pct": np.nan,
        "rmse_gain_keep_80pct": np.nan,
    }
    if uncertainty is None:
        return summary

    y_true = true[holdout]
    y_pred = pred[holdout]
    y_std = uncertainty[holdout]
    valid = np.isfinite(y_true) & np.isfinite(y_pred) & np.isfinite(y_std)
    if valid.sum() < 3:
        return summary

    y_true = y_true[valid]
    y_pred = y_pred[valid]
    y_std = np.maximum(y_std[valid], 1e-6)
    errors = y_true - y_pred
    abs_errors = np.abs(errors)
    sq_errors = errors ** 2

    summary["uncertainty_coverage_68"] = float(np.mean(abs_errors <= y_std))
    summary["uncertainty_coverage_95"] = float(np.mean(abs_errors <= 1.96 * y_std))
    summary["uncertainty_nll"] = float(
        np.mean(0.5 * np.log(2 * np.pi * y_std ** 2) + sq_errors / (2 * y_std ** 2))
    )
    summary["uncertainty_error_pearson"] = safe_stat_correlation(y_std, abs_errors, "pearson")
    summary["uncertainty_error_spearman"] = safe_stat_correlation(y_std, abs_errors, "spearman")

    order = np.argsort(y_std)
    full_rmse = float(np.sqrt(np.mean(sq_errors)))
    for retain_fraction, key in [(0.8, "rmse_keep_80pct"), (0.5, "rmse_keep_50pct")]:
        keep_n = max(1, int(len(order) * retain_fraction))
        keep_idx = order[:keep_n]
        summary[key] = float(np.sqrt(np.mean(sq_errors[keep_idx])))
    summary["rmse_gain_keep_80pct"] = float(full_rmse - summary["rmse_keep_80pct"])
    return summary


def compute_filtering_curve(
    true: np.ndarray,
    pred: np.ndarray,
    uncertainty: np.ndarray | None,
    holdout: np.ndarray,
    method_name: str,
    seed: int,
    mask_type: str,
) -> list[dict[str, float | int | str]]:
    if uncertainty is None:
        return []

    y_true = true[holdout]
    y_pred = pred[holdout]
    y_std = uncertainty[holdout]
    valid = np.isfinite(y_true) & np.isfinite(y_pred) & np.isfinite(y_std)
    if valid.sum() < 3:
        return []

    y_true = y_true[valid]
    y_pred = y_pred[valid]
    y_std = y_std[valid]
    sq_errors = (y_true - y_pred) ** 2
    order = np.argsort(y_std)

    rows = []
    for retain_fraction in FILTER_RETAIN_FRACTIONS:
        keep_n = max(1, int(len(order) * retain_fraction))
        keep_idx = order[:keep_n]
        rows.append(
            {
                "method": method_name,
                "seed": seed,
                "mask_type": mask_type,
                "retain_fraction": float(retain_fraction),
                "rmse": float(np.sqrt(np.mean(sq_errors[keep_idx]))),
                "mean_uncertainty": float(np.mean(y_std[keep_idx])),
                "n_points": int(keep_n),
            }
        )
    return rows


def compute_reliability_curve(
    true: np.ndarray,
    pred: np.ndarray,
    uncertainty: np.ndarray | None,
    holdout: np.ndarray,
    method_name: str,
    seed: int,
    mask_type: str,
    n_bins: int,
) -> list[dict[str, float | int | str]]:
    if uncertainty is None:
        return []

    y_true = true[holdout]
    y_pred = pred[holdout]
    y_std = uncertainty[holdout]
    valid = np.isfinite(y_true) & np.isfinite(y_pred) & np.isfinite(y_std)
    if valid.sum() < max(3, n_bins):
        return []

    y_true = y_true[valid]
    y_pred = y_pred[valid]
    y_std = y_std[valid]
    sq_errors = (y_true - y_pred) ** 2
    abs_errors = np.abs(y_true - y_pred)
    order = np.argsort(y_std)
    split_idx = np.array_split(order, n_bins)

    rows = []
    for bin_idx, idx in enumerate(split_idx):
        if len(idx) == 0:
            continue
        rows.append(
            {
                "method": method_name,
                "seed": seed,
                "mask_type": mask_type,
                "bin_idx": int(bin_idx),
                "mean_uncertainty": float(np.mean(y_std[idx])),
                "empirical_rmse": float(np.sqrt(np.mean(sq_errors[idx]))),
                "empirical_mae": float(np.mean(abs_errors[idx])),
                "coverage_68": float(np.mean(abs_errors[idx] <= y_std[idx])),
                "coverage_95": float(np.mean(abs_errors[idx] <= 1.96 * y_std[idx])),
                "n_points": int(len(idx)),
            }
        )
    return rows


def evaluate_layer_consistency(
    matrix_gs: np.ndarray,
    layer_codes: np.ndarray,
    n_domains: int,
):
    spot_matrix = matrix_gs.T
    scaler = StandardScaler()
    spot_matrix = scaler.fit_transform(spot_matrix)
    kmeans = KMeans(n_clusters=n_domains, random_state=42, n_init=10)
    pred = kmeans.fit_predict(spot_matrix)
    ari = adjusted_rand_score(layer_codes, pred)
    nmi = normalized_mutual_info_score(layer_codes, pred)
    return float(ari), float(nmi), int(len(np.unique(pred)))


def measure_method(
    method_name: str,
    gp_kernel: str,
    gp_alpha: float,
    gp_variant: str,
    gp_lambda_expr: float,
    gp_product_offset: float,
    gp_gate_mode: str,
    gp_gate_tau: float,
    gp_local_k: int,
    gp_layer_gate_mode: str,
    expr_n_components: float,
    gp_use_theta: bool,
    impute_fn: Callable[[], Tuple[np.ndarray, np.ndarray | None] | np.ndarray],
    true_values: np.ndarray,
    holdout_mask: np.ndarray,
    layer_codes: np.ndarray,
    n_domains: int,
    seed: int,
    mask_type: str,
) -> tuple[MethodResult, np.ndarray, np.ndarray | None]:
    process = psutil.Process(os.getpid())
    peak_rss_mb = process.memory_info().rss / 1024 / 1024
    stop_event = threading.Event()

    def monitor_memory() -> None:
        nonlocal peak_rss_mb
        while not stop_event.is_set():
            rss_mb = process.memory_info().rss / 1024 / 1024
            peak_rss_mb = max(peak_rss_mb, rss_mb)
            time.sleep(0.02)

    monitor = threading.Thread(target=monitor_memory, daemon=True)
    monitor.start()
    t0 = time.perf_counter()
    try:
        out = impute_fn()
    finally:
        runtime = time.perf_counter() - t0
        stop_event.set()
        monitor.join(timeout=1.0)
        peak_rss_mb = max(peak_rss_mb, process.memory_info().rss / 1024 / 1024)

    if isinstance(out, tuple):
        pred, uncertainty = out
    else:
        pred, uncertainty = out, None

    rmse, mae, pearson, spearman, r2 = compute_metrics(true_values, pred, holdout_mask)
    ari, nmi, n_pred_domains = evaluate_layer_consistency(pred, layer_codes, n_domains)

    row = MethodResult(
        method=method_name,
        gp_kernel=gp_kernel,
        gp_alpha=float(gp_alpha),
        gp_variant=gp_variant,
        gp_lambda_expr=float(gp_lambda_expr),
        gp_product_offset=float(gp_product_offset),
        gp_gate_mode=gp_gate_mode,
        gp_gate_tau=float(gp_gate_tau),
        gp_local_k=int(gp_local_k),
        gp_layer_gate_mode=gp_layer_gate_mode,
        expr_n_components=float(expr_n_components),
        gp_use_theta=bool(gp_use_theta),
        seed=seed,
        mask_type=mask_type,
        rmse=rmse,
        mae=mae,
        pearson=pearson,
        spearman=spearman,
        r2=r2,
        n_holdout=int(holdout_mask.sum()),
        runtime_s=float(runtime),
        peak_rss_mb=float(peak_rss_mb),
        layer_ari=ari,
        layer_nmi=nmi,
        n_domains=n_pred_domains,
    )
    return row, pred, uncertainty


def plot_benchmark_summary(results_df: pd.DataFrame, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    metric_panels = [
        ("rmse", "Lower Better"),
        ("pearson", "Higher Better"),
        ("runtime_s", "Lower Better"),
        ("layer_ari", "Higher Better"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.ravel()

    summary = (
        results_df.groupby("method")[["rmse", "pearson", "runtime_s", "layer_ari"]]
        .mean()
        .reindex(METHOD_ORDER)
        .dropna(how="all")
    )

    for ax, (metric, subtitle) in zip(axes, metric_panels):
        vals = summary[metric]
        ax.bar(
            np.arange(len(vals)),
            vals.values,
            color=[METHOD_COLORS[m] for m in vals.index],
            edgecolor="white",
            linewidth=0.6,
        )
        ax.set_xticks(np.arange(len(vals)))
        ax.set_xticklabels(vals.index, rotation=25, ha="right")
        ax.set_title(f"{metric.upper()} ({subtitle})")
        ax.grid(axis="y", alpha=0.25)
        for i, v in enumerate(vals.values):
            ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_summary_panels.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_mask_type_comparison(results_df: pd.DataFrame, output_dir: Path):
    summary = (
        results_df.groupby(["mask_type", "method"])["rmse"]
        .mean()
        .reset_index()
    )
    mask_types = [m for m in MASK_ORDER if m in summary["mask_type"].unique()]
    methods = [m for m in METHOD_ORDER if m in summary["method"].unique()]

    x = np.arange(len(mask_types))
    width = 0.8 / max(len(methods), 1)
    fig, ax = plt.subplots(figsize=(10, 5))

    for i, method in enumerate(methods):
        vals = []
        for mask_type in mask_types:
            row = summary[(summary["mask_type"] == mask_type) & (summary["method"] == method)]
            vals.append(row["rmse"].iloc[0] if len(row) else np.nan)
        ax.bar(
            x + (i - len(methods) / 2 + 0.5) * width,
            vals,
            width * 0.9,
            label=method,
            color=METHOD_COLORS[method],
            edgecolor="white",
            linewidth=0.5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(mask_types)
    ax.set_ylabel("RMSE")
    ax.set_title("RMSE Across Mask Types")
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_mask_types_rmse.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_seed_stability(results_df: pd.DataFrame, output_dir: Path):
    fig, ax = plt.subplots(figsize=(10, 5))
    for method in METHOD_ORDER:
        sub = results_df[results_df["method"] == method]
        if sub.empty:
            continue
        seed_summary = sub.groupby("seed")["rmse"].mean().reset_index()
        ax.plot(
            seed_summary["seed"],
            seed_summary["rmse"],
            marker="o",
            linewidth=2,
            label=method,
            color=METHOD_COLORS[method],
        )
    ax.set_xlabel("Random Seed")
    ax.set_ylabel("Mean RMSE Across Mask Types")
    ax.set_title("Seed Stability")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_seed_stability.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_mask_layouts(
    mask_records: list[dict[str, np.ndarray | str | int]],
    coords: np.ndarray,
    output_dir: Path,
):
    if not mask_records:
        return

    n = len(mask_records)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, record in zip(axes, mask_records):
        mask = record["holdout_mask"]
        spot_fraction = mask.mean(axis=0)
        sc = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=spot_fraction,
            cmap="viridis",
            s=28,
            edgecolor="none",
        )
        ax.set_title(f"{record['mask_type']} (seed {record['seed']})")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes[n:]:
        ax.axis("off")

    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_mask_layouts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_runtime_vs_accuracy(results_df: pd.DataFrame, output_dir: Path):
    summary = (
        results_df.groupby("method")[["rmse", "runtime_s"]]
        .mean()
        .reindex(METHOD_ORDER)
        .dropna(how="all")
    )
    fig, ax = plt.subplots(figsize=(7, 5))
    for method, row in summary.iterrows():
        ax.scatter(
            row["runtime_s"],
            row["rmse"],
            s=90,
            color=METHOD_COLORS[method],
            label=method,
            edgecolor="white",
            linewidth=0.7,
        )
        ax.text(row["runtime_s"], row["rmse"], f" {method}", va="center", fontsize=9)
    ax.set_xlabel("Runtime (s)")
    ax.set_ylabel("RMSE")
    ax.set_title("Runtime vs Accuracy")
    if (summary["runtime_s"] > 0).all():
        ax.set_xscale("log")
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_runtime_vs_accuracy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_memory_usage(results_df: pd.DataFrame, output_dir: Path):
    summary = (
        results_df.groupby("method")["peak_rss_mb"]
        .mean()
        .reindex(METHOD_ORDER)
        .dropna()
    )
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(
        np.arange(len(summary)),
        summary.values,
        color=[METHOD_COLORS[m] for m in summary.index],
        edgecolor="white",
        linewidth=0.6,
    )
    ax.set_xticks(np.arange(len(summary)))
    ax.set_xticklabels(summary.index, rotation=25, ha="right")
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_title("Memory Usage by Method")
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i, value in enumerate(summary.values):
        ax.text(i, value, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_memory_usage.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_biological_consistency(results_df: pd.DataFrame, output_dir: Path):
    summary = (
        results_df.groupby("method")[["layer_ari", "layer_nmi"]]
        .mean()
        .reindex(METHOD_ORDER)
        .dropna(how="all")
    )
    x = np.arange(len(summary))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width / 2, summary["layer_ari"], width, label="Layer ARI", color="#16a085")
    ax.bar(x + width / 2, summary["layer_nmi"], width, label="Layer NMI", color="#2980b9")
    ax.set_xticks(x)
    ax.set_xticklabels(summary.index, rotation=25, ha="right")
    ax.set_ylabel("Score")
    ax.set_title("Biological Consistency with MOB Layer Labels")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_biological_consistency.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gp_uncertainty_panels(
    results_df: pd.DataFrame,
    filter_curve_df: pd.DataFrame,
    reliability_df: pd.DataFrame,
    output_dir: Path,
):
    gp_results = results_df[results_df["method"] == "spagapa_gp"].copy()
    if gp_results.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes = axes.ravel()
    mask_types = [m for m in MASK_ORDER if m in gp_results["mask_type"].unique()]

    coverage = (
        gp_results.groupby("mask_type")[["uncertainty_coverage_68", "uncertainty_coverage_95"]]
        .mean()
        .reindex(mask_types)
    )
    x = np.arange(len(coverage))
    width = 0.35
    axes[0].bar(x - width / 2, coverage["uncertainty_coverage_68"], width, label="Empirical 68%", color="#16a085")
    axes[0].bar(x + width / 2, coverage["uncertainty_coverage_95"], width, label="Empirical 95%", color="#2980b9")
    axes[0].axhline(0.6827, color="#16a085", linestyle="--", linewidth=1)
    axes[0].axhline(0.95, color="#2980b9", linestyle="--", linewidth=1)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(mask_types, rotation=25, ha="right")
    axes[0].set_ylabel("Coverage")
    axes[0].set_title("GP Interval Coverage")
    axes[0].legend(frameon=False)
    axes[0].grid(axis="y", alpha=0.25)

    corr = (
        gp_results.groupby("mask_type")[["uncertainty_error_pearson", "uncertainty_error_spearman"]]
        .mean()
        .reindex(mask_types)
    )
    axes[1].bar(x - width / 2, corr["uncertainty_error_pearson"], width, label="Pearson", color="#c0392b")
    axes[1].bar(x + width / 2, corr["uncertainty_error_spearman"], width, label="Spearman", color="#8e44ad")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(mask_types, rotation=25, ha="right")
    axes[1].set_ylabel("Correlation")
    axes[1].set_title("Uncertainty vs Absolute Error")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", alpha=0.25)

    if not filter_curve_df.empty:
        filter_summary = (
            filter_curve_df[filter_curve_df["method"] == "spagapa_gp"]
            .groupby(["mask_type", "retain_fraction"])["rmse"]
            .mean()
            .reset_index()
        )
        for mask_type in mask_types:
            sub = filter_summary[filter_summary["mask_type"] == mask_type]
            if sub.empty:
                continue
            axes[2].plot(
                sub["retain_fraction"],
                sub["rmse"],
                marker="o",
                linewidth=2,
                label=mask_type,
                color=MASK_COLORS.get(mask_type, "#7f8c8d"),
            )
        axes[2].set_xlabel("Retained Lowest-Uncertainty Fraction")
        axes[2].set_ylabel("RMSE")
        axes[2].set_title("GP Accuracy After High-Uncertainty Filtering")
        axes[2].legend(frameon=False, ncol=2)
        axes[2].grid(alpha=0.25)

    if not reliability_df.empty:
        rel_summary = (
            reliability_df[reliability_df["method"] == "spagapa_gp"]
            .groupby("bin_idx")[["mean_uncertainty", "empirical_rmse"]]
            .mean()
            .reset_index()
        )
        axes[3].plot(
            rel_summary["mean_uncertainty"],
            rel_summary["empirical_rmse"],
            marker="o",
            linewidth=2,
            color="#2c3e50",
        )
        max_val = float(
            max(
                rel_summary["mean_uncertainty"].max(),
                rel_summary["empirical_rmse"].max(),
            )
        )
        axes[3].plot([0, max_val], [0, max_val], linestyle="--", color="#7f8c8d", linewidth=1)
        axes[3].set_xlabel("Mean Predicted Std")
        axes[3].set_ylabel("Empirical RMSE")
        axes[3].set_title("GP Reliability Curve")
        axes[3].grid(alpha=0.25)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_dir / "benchmark_gp_uncertainty_panels.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    mask_types = [x.strip() for x in args.mask_types.split(",") if x.strip()]

    apa, expression, coords_df, metadata = load_inputs(data_dir)
    apa = select_genes(apa, args.n_genes, args.min_observed_spots)
    apa = apa.loc[:, coords_df.index]

    values = apa.values.astype(float)
    coords = coords_df[["x", "y"]].values.astype(float)
    gene_names = apa.index.astype(str).tolist()
    spot_names = apa.columns.astype(str).tolist()
    layer_labels = metadata.loc[spot_names, "layer"].astype(str).values
    layer_codes = pd.Categorical(metadata.loc[spot_names, "layer"]).codes
    expression_embedding = None
    if args.gp_variant in {"expr_additive", "expr_product", "adaptive_additive", "expr_layer_local"}:
        expression_aligned = expression.reindex(columns=spot_names).fillna(0.0)
        expression_embedding = ExpressionFeatureBuilder(
            n_components=args.expr_n_components,
            use_hvg=args.expr_use_hvg,
            n_top_genes=args.expr_n_top_genes,
            orientation="genes_by_spots",
        ).fit_transform(expression_aligned)

    geometry_features = None
    if args.gp_variant == "spatial_radial":
        geometry_features = GeometryFeatureBuilder(use_theta=args.gp_use_theta).build(
            coords,
            layer_labels=layer_labels,
        )
    elif args.gp_variant == "expr_layer_local":
        geometry_features = GeometryFeatureBuilder(use_theta=args.gp_use_theta).build(coords)

    rows: List[MethodResult] = []
    sample_outputs: Dict[str, np.ndarray] = {}
    filter_curve_rows: list[dict[str, float | int | str]] = []
    reliability_rows: list[dict[str, float | int | str]] = []
    mask_records: list[dict[str, np.ndarray | str | int]] = []
    example_seed = seeds[0] if seeds else 42

    for seed in seeds:
        for mask_type in mask_types:
            holdout = build_mask(mask_type, values, coords, layer_labels, args.mask_fraction, seed)
            train = make_train_matrix(values, holdout)
            if seed == example_seed:
                mask_records.append(
                    {
                        "mask_type": mask_type,
                        "seed": seed,
                        "holdout_mask": holdout.copy(),
                    }
                )

            method_specs: Dict[str, Callable[[], np.ndarray | Tuple[np.ndarray, np.ndarray | None]]] = {
                "mean": lambda train=train: impute_mean(train),
                "median": lambda train=train: impute_median(train),
                "knn_spatial": lambda train=train: impute_spatial_knn(train, coords, args.knn_k),
                "stapaminer_knn_expression": lambda train=train: impute_stapaminer_expression_knn(
                    train, gene_names, spot_names, expression, args.knn_k
                ),
                "spagapa_gp": lambda train=train: impute_spagapa_gp_variant(
                    train,
                    coords,
                    args.gp_kernel,
                    args.gp_alpha,
                    args.gp_variant,
                    expression_embedding=expression_embedding,
                    geometry_features=geometry_features,
                    lambda_expr=args.gp_lambda_expr,
                    product_offset=args.gp_product_offset,
                    gate_mode=args.gp_gate_mode,
                    gate_tau=args.gp_gate_tau,
                    local_k=args.gp_local_k,
                    layer_gate_mode=args.gp_layer_gate_mode,
                    use_theta=args.gp_use_theta,
                ),
            }

            for method_name, fn in method_specs.items():
                result, pred, uncertainty = measure_method(
                    method_name=method_name,
                    gp_kernel=args.gp_kernel if method_name == "spagapa_gp" else "na",
                    gp_alpha=args.gp_alpha if method_name == "spagapa_gp" else np.nan,
                    gp_variant=args.gp_variant if method_name == "spagapa_gp" else "na",
                    gp_lambda_expr=args.gp_lambda_expr if method_name == "spagapa_gp" else np.nan,
                    gp_product_offset=args.gp_product_offset if method_name == "spagapa_gp" else np.nan,
                    gp_gate_mode=args.gp_gate_mode if method_name == "spagapa_gp" else "na",
                    gp_gate_tau=args.gp_gate_tau if method_name == "spagapa_gp" and args.gp_gate_tau is not None else np.nan,
                    gp_local_k=args.gp_local_k if method_name == "spagapa_gp" else -1,
                    gp_layer_gate_mode=args.gp_layer_gate_mode if method_name == "spagapa_gp" else "na",
                    expr_n_components=args.expr_n_components if method_name == "spagapa_gp" else np.nan,
                    gp_use_theta=args.gp_use_theta if method_name == "spagapa_gp" else False,
                    impute_fn=fn,
                    true_values=values,
                    holdout_mask=holdout,
                    layer_codes=layer_codes,
                    n_domains=args.n_domains,
                    seed=seed,
                    mask_type=mask_type,
                )

                uncertainty_summary = compute_uncertainty_summary(values, pred, uncertainty, holdout)
                for key, value in uncertainty_summary.items():
                    setattr(result, key, value)
                filter_curve_rows.extend(
                    compute_filtering_curve(values, pred, uncertainty, holdout, method_name, seed, mask_type)
                )
                reliability_rows.extend(
                    compute_reliability_curve(
                        values,
                        pred,
                        uncertainty,
                        holdout,
                        method_name,
                        seed,
                        mask_type,
                        n_bins=args.calibration_bins,
                    )
                )
                rows.append(result)

                sample_key = f"{mask_type}_seed{seed}_{method_name}"
                if sample_key not in sample_outputs:
                    sample_outputs[sample_key] = pred
                    np.save(output_dir / f"{sample_key}_imputed.npy", pred)
                    if uncertainty is not None:
                        np.save(output_dir / f"{sample_key}_uncertainty.npy", uncertainty)

    results_df = pd.DataFrame([r.__dict__ for r in rows])
    results_df.to_csv(output_dir / "benchmark_results_long.csv", index=False)
    filter_curve_df = pd.DataFrame(filter_curve_rows)
    reliability_df = pd.DataFrame(reliability_rows)
    if not filter_curve_df.empty:
        filter_curve_df.to_csv(output_dir / "uncertainty_filtering_curve.csv", index=False)
    if not reliability_df.empty:
        reliability_df.to_csv(output_dir / "uncertainty_reliability_curve.csv", index=False)

    summary_df = (
        results_df.groupby(["mask_type", "method"])[
            [
                "rmse",
                "mae",
                "pearson",
                "spearman",
                "r2",
                "runtime_s",
                "peak_rss_mb",
                "layer_ari",
                "layer_nmi",
                "uncertainty_coverage_68",
                "uncertainty_coverage_95",
                "uncertainty_nll",
                "uncertainty_error_pearson",
                "uncertainty_error_spearman",
                "rmse_keep_80pct",
                "rmse_keep_50pct",
                "rmse_gain_keep_80pct",
            ]
        ]
        .agg(["mean", "std"])
    )
    summary_df.to_csv(output_dir / "benchmark_results_summary.csv")

    overall_df = (
        results_df.groupby("method")[
            [
                "rmse",
                "mae",
                "pearson",
                "spearman",
                "r2",
                "runtime_s",
                "peak_rss_mb",
                "layer_ari",
                "layer_nmi",
                "uncertainty_coverage_68",
                "uncertainty_coverage_95",
                "uncertainty_nll",
                "uncertainty_error_pearson",
                "uncertainty_error_spearman",
                "rmse_keep_80pct",
                "rmse_keep_50pct",
                "rmse_gain_keep_80pct",
            ]
        ]
        .mean()
        .reindex(METHOD_ORDER)
    )
    overall_df.to_csv(output_dir / "benchmark_results_overall.csv")

    plot_benchmark_summary(results_df, figure_dir)
    plot_mask_type_comparison(results_df, figure_dir)
    plot_seed_stability(results_df, figure_dir)
    plot_runtime_vs_accuracy(results_df, figure_dir)
    plot_memory_usage(results_df, figure_dir)
    plot_biological_consistency(results_df, figure_dir)
    plot_mask_layouts(mask_records, coords, figure_dir)
    plot_gp_uncertainty_panels(results_df, filter_curve_df, reliability_df, figure_dir)

    summary = {
        "dataset": "stAPAminer_MOB",
        "n_genes": int(values.shape[0]),
        "n_spots": int(values.shape[1]),
        "raw_missing_rate": float(np.isnan(values).mean()),
        "mask_fraction": float(args.mask_fraction),
        "seeds": seeds,
        "mask_types": mask_types,
        "methods": METHOD_ORDER,
        "gp_kernel": str(args.gp_kernel),
        "gp_alpha": float(args.gp_alpha),
        "gp_variant": str(args.gp_variant),
        "gp_lambda_expr": float(args.gp_lambda_expr),
        "gp_product_offset": float(args.gp_product_offset),
        "gp_gate_mode": str(args.gp_gate_mode),
        "gp_gate_tau": None if args.gp_gate_tau is None else float(args.gp_gate_tau),
        "gp_local_k": int(args.gp_local_k),
        "gp_layer_gate_mode": str(args.gp_layer_gate_mode),
        "expr_n_components": int(args.expr_n_components),
        "expr_use_hvg": bool(args.expr_use_hvg),
        "expr_n_top_genes": int(args.expr_n_top_genes),
        "gp_use_theta": bool(args.gp_use_theta),
        "calibration_bins": int(args.calibration_bins),
        "best_method_by_mean_rmse": str(overall_df["rmse"].idxmin()),
        "overall_metrics": overall_df.reset_index().rename(columns={"index": "method"}).to_dict(orient="records"),
        "figure_dir": str(figure_dir),
    }
    (output_dir / "benchmark_summary.json").write_text(json.dumps(summary, indent=2))

    print(overall_df.to_string())
    print(f"\nSaved figures to: {figure_dir}")


if __name__ == "__main__":
    main()

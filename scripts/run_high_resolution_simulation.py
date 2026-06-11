#!/usr/bin/env python
"""High-resolution pseudo-bin benchmark for spaGAPA-BioML.

This script stress-tests high-resolution readiness without downloading a large
external dataset. It expands a prepared real spatial APA dataset into
pseudo-high-resolution sub-bins, increases missingness/noise, and compares
CPU-friendly methods:

* raw gene-mean fill
* expression KNN imputation
* sparse GP imputation
* sparse GP + BioML multi-view graph recovery

The benchmark is intended as a scalability and robustness test. It should not
be presented as a substitute for true high-resolution spatial APA data.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.neighbors import NearestNeighbors

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from spagapa.bioml import BioMLDomainDetector, GraphRegularizedAPAFactorizer, MultiViewGraphBuilder
from spagapa.imputation import ExpressionFeatureBuilder, SparseGPImputer

from run_external_bioml_validation import (
    evaluate_layer_metrics,
    fill_missing_by_gene_mean,
    load_dataset,
    resolve_project_path,
    select_genes,
    track_runtime_memory,
)


METHOD_ORDER = ["raw", "expression_knn", "sparse_gp", "sparse_bioml"]
METHOD_COLORS = {
    "raw": "#7f8c8d",
    "expression_knn": "#9b59b6",
    "sparse_gp": "#e74c3c",
    "sparse_bioml": "#d35400",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared real dataset directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/highres_simulation_v1",
        help="Output directory.",
    )
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--n-genes", type=int, default=40)
    parser.add_argument("--min-observed-spots", type=int, default=100)
    parser.add_argument(
        "--max-parent-spots",
        type=int,
        default=None,
        help="Optional random subset of parent spots for quick tests.",
    )
    parser.add_argument(
        "--subbins-per-spot",
        default="4",
        help="Comma-separated pseudo-bins per parent spot, e.g. 2,4,8.",
    )
    parser.add_argument("--capture-rate", type=float, default=0.45)
    parser.add_argument("--dropout-rate", type=float, default=0.25)
    parser.add_argument("--measurement-noise", type=float, default=0.08)
    parser.add_argument("--micro-noise", type=float, default=0.04)
    parser.add_argument("--jitter-fraction", type=float, default=0.18)
    parser.add_argument("--expression-noise", type=float, default=0.03)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--methods",
        default="raw,expression_knn,sparse_gp,sparse_bioml",
        help="Comma-separated methods.",
    )
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--expr-n-components", type=int, default=10)
    parser.add_argument("--n-domains", type=int, default=None)
    parser.add_argument("--sparse-n-inducing", type=int, default=60)
    parser.add_argument(
        "--sparse-length-scale",
        default="auto",
        help="Sparse GP RBF length scale, or 'auto' to estimate from high-res coordinates.",
    )
    parser.add_argument(
        "--sparse-length-scale-multiplier",
        type=float,
        default=1.0,
        help="Multiplier applied when --sparse-length-scale=auto.",
    )
    parser.add_argument("--sparse-noise-level", type=float, default=0.08)
    parser.add_argument("--bioml-rank", type=int, default=8)
    parser.add_argument("--bioml-lambda-graph", type=float, default=0.5)
    parser.add_argument("--bioml-lambda-l2", type=float, default=1e-2)
    parser.add_argument("--bioml-max-iter", type=int, default=20)
    parser.add_argument("--bioml-n-neighbors", type=int, default=15)
    parser.add_argument("--bioml-blend", type=float, default=0.1)
    parser.add_argument(
        "--bioml-domain-method",
        default="spectral",
        choices=["spectral", "kmeans"],
    )
    parser.add_argument("--bioml-spatial-weight", type=float, default=0.4)
    parser.add_argument("--bioml-expression-weight", type=float, default=0.4)
    parser.add_argument("--bioml-apa-weight", type=float, default=0.2)
    return parser.parse_args()


def parse_int_list(text: str) -> list[int]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        value = int(item)
        if value <= 0:
            raise ValueError("subbins-per-spot values must be positive")
        values.append(value)
    if not values:
        raise ValueError("At least one subbins-per-spot value is required")
    return values


def sample_parent_spots(
    spot_names: list[str],
    metadata: pd.DataFrame,
    layer_column: str,
    max_parent_spots: int | None,
    seed: int,
) -> list[str]:
    if max_parent_spots is None or max_parent_spots >= len(spot_names):
        return spot_names

    rng = np.random.default_rng(seed)
    if layer_column not in metadata.columns:
        chosen = rng.choice(spot_names, size=max_parent_spots, replace=False)
        return sorted(chosen.tolist(), key=spot_names.index)

    selected = []
    metadata_subset = metadata.loc[spot_names]
    layer_sizes = metadata_subset[layer_column].astype(str).value_counts()
    for layer, size in layer_sizes.items():
        layer_spots = metadata_subset[metadata_subset[layer_column].astype(str) == layer].index.tolist()
        n_layer = max(1, int(round(max_parent_spots * size / len(spot_names))))
        n_layer = min(n_layer, len(layer_spots))
        selected.extend(rng.choice(layer_spots, size=n_layer, replace=False).tolist())

    if len(selected) > max_parent_spots:
        selected = rng.choice(selected, size=max_parent_spots, replace=False).tolist()
    elif len(selected) < max_parent_spots:
        remaining = [spot for spot in spot_names if spot not in set(selected)]
        add = rng.choice(remaining, size=min(max_parent_spots - len(selected), len(remaining)), replace=False)
        selected.extend(add.tolist())

    selected_set = set(selected)
    return [spot for spot in spot_names if spot in selected_set]


def median_nearest_distance(coords: np.ndarray) -> float:
    if len(coords) < 2:
        return 1.0
    nn = NearestNeighbors(n_neighbors=2)
    nn.fit(coords)
    distances, _ = nn.kneighbors(coords)
    positive = distances[:, 1]
    positive = positive[positive > 0]
    return float(np.median(positive)) if positive.size else 1.0


def build_parent_expression_embedding(
    expression: pd.DataFrame | None,
    spot_names: list[str],
    coords: np.ndarray,
    n_components: int,
) -> np.ndarray:
    if expression is None:
        return coords.astype(float)
    expr = expression.reindex(columns=spot_names).fillna(0.0)
    builder = ExpressionFeatureBuilder(
        n_components=n_components,
        orientation="genes_by_spots",
    )
    return builder.fit_transform(expr)


def simulate_pseudo_bins(
    apa: pd.DataFrame,
    expression: pd.DataFrame | None,
    coords: pd.DataFrame,
    metadata: pd.DataFrame,
    subbins_per_spot: int,
    args: argparse.Namespace,
) -> dict[str, Any]:
    rng = np.random.default_rng(args.seed + subbins_per_spot * 1009)
    spot_names = apa.columns.astype(str).tolist()
    parent_coords = coords.loc[spot_names][["x", "y"]].values.astype(float)
    n_genes, n_parent = apa.shape
    n_bins = n_parent * subbins_per_spot

    parent_index = np.repeat(np.arange(n_parent), subbins_per_spot)
    parent_spots = np.repeat(np.asarray(spot_names, dtype=object), subbins_per_spot)
    subbin_ids = [
        f"{spot}_bin{bin_idx}"
        for spot in spot_names
        for bin_idx in range(subbins_per_spot)
    ]

    jitter_scale = median_nearest_distance(parent_coords) * args.jitter_fraction
    highres_coords = np.repeat(parent_coords, subbins_per_spot, axis=0)
    highres_coords = highres_coords + rng.normal(scale=jitter_scale, size=highres_coords.shape)

    parent_truth = fill_missing_by_gene_mean(apa.values.astype(float))
    truth = np.repeat(parent_truth, subbins_per_spot, axis=1)
    truth = truth + rng.normal(scale=args.micro_noise, size=truth.shape)
    truth = np.clip(truth, 0.0, 1.0)

    observed_mask = rng.random(size=truth.shape) < args.capture_rate
    observed_mask &= rng.random(size=truth.shape) >= args.dropout_rate
    observed_values = truth + rng.normal(scale=args.measurement_noise, size=truth.shape)
    observed_values = np.clip(observed_values, 0.0, 1.0)
    observed = np.where(observed_mask, observed_values, np.nan)

    parent_embedding = build_parent_expression_embedding(
        expression,
        spot_names,
        parent_coords,
        args.expr_n_components,
    )
    expr_scale = np.nanstd(parent_embedding, axis=0, keepdims=True)
    expr_scale = np.where(expr_scale > 0, expr_scale, 1.0)
    expression_embedding = np.repeat(parent_embedding, subbins_per_spot, axis=0)
    expression_embedding = expression_embedding + rng.normal(
        scale=args.expression_noise * expr_scale,
        size=expression_embedding.shape,
    )

    layer_labels = None
    layer_codes = None
    if args.layer_column in metadata.columns:
        parent_layers = metadata.loc[spot_names, args.layer_column].astype(str).values
        layer_labels = np.repeat(parent_layers, subbins_per_spot)
        unique_layers = {layer: idx for idx, layer in enumerate(sorted(np.unique(layer_labels)))}
        layer_codes = np.asarray([unique_layers[layer] for layer in layer_labels], dtype=int)

    return {
        "truth": truth,
        "observed": observed,
        "observed_mask": np.isfinite(observed),
        "coords": highres_coords,
        "expression_embedding": expression_embedding,
        "parent_index": parent_index,
        "parent_spots": parent_spots,
        "subbin_ids": subbin_ids,
        "layer_labels": layer_labels,
        "layer_codes": layer_codes,
        "n_parent": n_parent,
        "n_bins": n_bins,
    }


def row_nanmean(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    sums = np.where(finite, values, 0.0).sum(axis=1)
    means = np.divide(sums, np.maximum(counts, 1))
    means[counts == 0] = 0.0
    return means


def expression_knn_impute(values: np.ndarray, expression_embedding: np.ndarray, k: int) -> np.ndarray:
    n_spots = values.shape[1]
    n_neighbors = min(max(1, k), n_spots - 1)
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1)
    nn.fit(expression_embedding)
    _, indices = nn.kneighbors(expression_embedding)
    neighbors = indices[:, 1:]

    filled = values.copy()
    for _ in range(10):
        if np.isfinite(filled).all():
            break
        old_missing = int((~np.isfinite(filled)).sum())
        for spot_idx in range(n_spots):
            missing = ~np.isfinite(filled[:, spot_idx])
            if not missing.any():
                continue
            source = filled[:, neighbors[spot_idx]]
            means = row_nanmean(source)
            filled[missing, spot_idx] = means[missing]
        if int((~np.isfinite(filled)).sum()) == old_missing:
            break

    gene_means = row_nanmean(filled)
    missing_gene, missing_spot = np.where(~np.isfinite(filled))
    filled[missing_gene, missing_spot] = gene_means[missing_gene]
    return np.clip(filled, 0.0, 1.0)


def sparse_gp_impute(
    observed: np.ndarray,
    coords: np.ndarray,
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray | None]:
    length_scale = resolve_sparse_length_scale(coords, args)
    imputer = SparseGPImputer(
        n_inducing=args.sparse_n_inducing,
        inducing_method="kmeans",
        length_scale=length_scale,
        noise_level=args.sparse_noise_level,
        kernel_type="rbf",
    )
    mask = np.isfinite(observed)
    batch = imputer.fit_batch(
        coords,
        np.nan_to_num(observed, nan=0.0),
        mask=mask,
        verbose=False,
    )
    imputed, uncertainty = batch.impute(return_uncertainty=True)
    return np.clip(imputed, 0.0, 1.0), uncertainty


def resolve_sparse_length_scale(coords: np.ndarray, args: argparse.Namespace) -> float:
    """Resolve numeric sparse GP length scale from CLI settings."""
    value = str(args.sparse_length_scale).strip().lower()
    if value not in {"auto", "none"}:
        length_scale = float(args.sparse_length_scale)
        if length_scale <= 0:
            raise ValueError("--sparse-length-scale must be positive or 'auto'")
        return length_scale

    if len(coords) < 2:
        return 1.0
    k = min(max(2, int(args.bioml_n_neighbors)), len(coords))
    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(coords)
    distances, _ = nn.kneighbors(coords)
    positive = distances[:, -1]
    positive = positive[np.isfinite(positive) & (positive > 0)]
    base = float(np.median(positive)) if positive.size else median_nearest_distance(coords)
    return max(1e-6, base * float(args.sparse_length_scale_multiplier))


def confidence_from_uncertainty(uncertainty: np.ndarray | None, mask: np.ndarray) -> np.ndarray | None:
    if uncertainty is None:
        return None
    finite = uncertainty[np.isfinite(uncertainty)]
    fallback = float(np.median(finite)) if finite.size else 1.0
    unc = np.nan_to_num(uncertainty, nan=fallback, posinf=fallback, neginf=fallback)
    confidence = 1.0 / (unc + 1e-6)
    scale_values = confidence[np.isfinite(confidence) & (confidence > 0)]
    scale = float(np.median(scale_values)) if scale_values.size else 1.0
    confidence = np.clip(confidence / max(scale, 1e-6), 0.0, 10.0)
    confidence[~mask] = 0.0
    return confidence


def sparse_bioml_refine(
    observed: np.ndarray,
    coords: np.ndarray,
    expression_embedding: np.ndarray,
    sparse_gp: np.ndarray,
    uncertainty: np.ndarray | None,
    n_domains: int,
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(observed)
    confidence = confidence_from_uncertainty(uncertainty, mask)
    graph = MultiViewGraphBuilder(
        n_neighbors=args.bioml_n_neighbors,
        spatial_weight=args.bioml_spatial_weight,
        expression_weight=args.bioml_expression_weight,
        apa_weight=args.bioml_apa_weight,
    ).build(
        coords,
        expression_embedding=expression_embedding,
        apa_matrix=sparse_gp,
        uncertainty=uncertainty,
    )

    factorizer = GraphRegularizedAPAFactorizer(
        rank=args.bioml_rank,
        lambda_graph=args.bioml_lambda_graph,
        lambda_l2=args.bioml_lambda_l2,
        max_iter=args.bioml_max_iter,
        random_state=42,
        preserve_observed=True,
    )
    bioml_imputed = factorizer.fit_transform(
        observed,
        graph_laplacian=graph.laplacian(),
        mask=mask,
        confidence=confidence,
    )
    blend = float(np.clip(args.bioml_blend, 0.0, 1.0))
    refined = (1.0 - blend) * sparse_gp + blend * bioml_imputed
    refined[mask] = observed[mask]
    refined = np.clip(refined, 0.0, 1.0)

    if args.bioml_domain_method == "spectral":
        domains = BioMLDomainDetector(
            method="spectral",
            n_domains=n_domains,
            random_state=42,
        ).fit_predict(graph=graph.fused)
    else:
        domains = BioMLDomainDetector(
            method="kmeans",
            n_domains=n_domains,
            random_state=42,
        ).fit_predict(spot_factors=factorizer.spot_factors_)

    return refined, domains


def aggregate_to_parent(matrix: np.ndarray, parent_index: np.ndarray, n_parent: int) -> np.ndarray:
    out = np.zeros((matrix.shape[0], n_parent), dtype=float)
    counts = np.zeros(n_parent, dtype=float)
    for parent in range(n_parent):
        idx = parent_index == parent
        counts[parent] = max(1, int(idx.sum()))
        out[:, parent] = np.nanmean(matrix[:, idx], axis=1)
    return np.nan_to_num(out, nan=0.0)


def compute_numeric_metrics(
    truth: np.ndarray,
    pred: np.ndarray,
    observed_mask: np.ndarray,
    parent_index: np.ndarray,
    n_parent: int,
    uncertainty: np.ndarray | None,
) -> dict[str, float]:
    holdout = ~observed_mask
    y_true = truth[holdout]
    y_pred = pred[holdout]
    all_true = truth.ravel()
    all_pred = pred.ravel()
    truth_parent = aggregate_to_parent(truth, parent_index, n_parent)
    pred_parent = aggregate_to_parent(pred, parent_index, n_parent)

    metrics = {
        "rmse_holdout": float(np.sqrt(mean_squared_error(y_true, y_pred))) if y_true.size else np.nan,
        "mae_holdout": float(mean_absolute_error(y_true, y_pred)) if y_true.size else np.nan,
        "rmse_all": float(np.sqrt(mean_squared_error(all_true, all_pred))),
        "mae_all": float(mean_absolute_error(all_true, all_pred)),
        "parent_rmse": float(np.sqrt(mean_squared_error(truth_parent.ravel(), pred_parent.ravel()))),
        "parent_mae": float(mean_absolute_error(truth_parent.ravel(), pred_parent.ravel())),
        "uncertainty_error_spearman": np.nan,
        "mean_uncertainty": np.nan,
    }
    if uncertainty is not None and y_true.size:
        abs_error = np.abs(y_true - y_pred)
        unc = uncertainty[holdout]
        valid = np.isfinite(unc) & np.isfinite(abs_error)
        if valid.sum() > 2 and len(np.unique(unc[valid])) > 1:
            metrics["uncertainty_error_spearman"] = float(
                spearmanr(unc[valid], abs_error[valid]).statistic
            )
        metrics["mean_uncertainty"] = float(np.nanmean(unc))
    return metrics


def run_methods(sim: dict[str, Any], args: argparse.Namespace, n_domains: int) -> dict[str, dict[str, Any]]:
    selected_methods = [method.strip() for method in args.methods.split(",") if method.strip()]
    observed = sim["observed"]
    outputs: dict[str, dict[str, Any]] = {}

    if "raw" in selected_methods:
        outputs["raw"] = {
            "matrix": fill_missing_by_gene_mean(observed),
            "uncertainty": None,
            "domains": None,
            "runtime_s": 0.0,
            "peak_rss_mb": 0.0,
        }

    if "expression_knn" in selected_methods:
        matrix, runtime_s, peak_rss_mb = track_runtime_memory(
            lambda: expression_knn_impute(observed, sim["expression_embedding"], args.knn_k)
        )
        outputs["expression_knn"] = {
            "matrix": matrix,
            "uncertainty": None,
            "domains": None,
            "runtime_s": runtime_s,
            "peak_rss_mb": peak_rss_mb,
        }

    sparse_result = None
    sparse_runtime_s = 0.0
    sparse_peak_rss_mb = 0.0
    if "sparse_gp" in selected_methods or "sparse_bioml" in selected_methods:
        sparse_result, sparse_runtime_s, sparse_peak_rss_mb = track_runtime_memory(
            lambda: sparse_gp_impute(observed, sim["coords"], args)
        )
        sparse_matrix, sparse_uncertainty = sparse_result
        if "sparse_gp" in selected_methods:
            outputs["sparse_gp"] = {
                "matrix": sparse_matrix,
                "uncertainty": sparse_uncertainty,
                "domains": None,
                "runtime_s": sparse_runtime_s,
                "peak_rss_mb": sparse_peak_rss_mb,
            }

    if "sparse_bioml" in selected_methods:
        if sparse_result is None:
            sparse_result, sparse_runtime_s, sparse_peak_rss_mb = track_runtime_memory(
                lambda: sparse_gp_impute(observed, sim["coords"], args)
            )
        sparse_matrix, sparse_uncertainty = sparse_result
        (matrix, domains), runtime_s, peak_rss_mb = track_runtime_memory(
            lambda: sparse_bioml_refine(
                observed,
                sim["coords"],
                sim["expression_embedding"],
                sparse_matrix,
                sparse_uncertainty,
                n_domains,
                args,
            )
        )
        outputs["sparse_bioml"] = {
            "matrix": matrix,
            "uncertainty": sparse_uncertainty,
            "domains": domains,
            "runtime_s": sparse_runtime_s + runtime_s,
            "peak_rss_mb": max(sparse_peak_rss_mb, peak_rss_mb),
        }

    return outputs


def evaluate_outputs(
    outputs: dict[str, dict[str, Any]],
    sim: dict[str, Any],
    scenario_name: str,
    dataset_name: str,
    n_domains: int,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rows = []
    domain_maps = {}
    for method in [m for m in METHOD_ORDER if m in outputs]:
        out = outputs[method]
        layer_metrics, pred_domains, _ = evaluate_layer_metrics(
            out["matrix"],
            sim["layer_codes"],
            n_domains,
            domain_labels=out.get("domains"),
        )
        domain_maps[method] = pred_domains
        numeric = compute_numeric_metrics(
            sim["truth"],
            out["matrix"],
            sim["observed_mask"],
            sim["parent_index"],
            sim["n_parent"],
            out.get("uncertainty"),
        )
        rows.append(
            {
                "dataset": dataset_name,
                "scenario": scenario_name,
                "method": method,
                "n_genes": int(sim["truth"].shape[0]),
                "n_parent_spots": int(sim["n_parent"]),
                "n_bins": int(sim["n_bins"]),
                "observed_fraction": float(sim["observed_mask"].mean()),
                "runtime_s": float(out["runtime_s"]),
                "peak_rss_mb": float(out["peak_rss_mb"]),
                **numeric,
                **layer_metrics,
            }
        )
    return pd.DataFrame(rows), domain_maps


def plot_method_comparison(results: pd.DataFrame, figure_dir: Path) -> None:
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "runtime_s"]
    available = [m for m in metrics if m in results.columns and results[m].notna().any()]
    if not available:
        return
    fig, axes = plt.subplots(1, len(available), figsize=(4.2 * len(available), 4.0))
    if len(available) == 1:
        axes = [axes]
    methods = [method for method in METHOD_ORDER if method in results["method"].values]
    x = np.arange(len(methods))
    for ax, metric in zip(axes, available):
        values = [
            float(results.loc[results["method"] == method, metric].mean())
            for method in methods
        ]
        ax.bar(x, values, color=[METHOD_COLORS.get(method, "#333333") for method in methods])
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_title(metric)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_method_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_scaling(results: pd.DataFrame, figure_dir: Path) -> None:
    if results.empty:
        return
    metrics = ["runtime_s", "peak_rss_mb", "rmse_holdout", "layer_ari"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 4.0))
    methods = [method for method in METHOD_ORDER if method in results["method"].values]
    for ax, metric in zip(axes, metrics):
        if metric not in results.columns:
            ax.axis("off")
            continue
        for method in methods:
            sub = results[results["method"] == method].sort_values("n_bins")
            if sub.empty:
                continue
            ax.plot(sub["n_bins"], sub[metric], marker="o", label=method)
        ax.set_xlabel("Pseudo-bins")
        ax.set_title(metric)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_scaling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_domain_maps(
    sim: dict[str, Any],
    domain_maps: dict[str, np.ndarray],
    figure_dir: Path,
    scenario_name: str,
) -> None:
    panels: list[tuple[str, np.ndarray]] = []
    if sim["layer_labels"] is not None:
        panels.append(("true_layer", sim["layer_labels"]))
    panels.extend((method, labels) for method, labels in domain_maps.items())
    if not panels:
        return

    n_cols = min(3, len(panels))
    n_rows = int(np.ceil(len(panels) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 4.0 * n_rows))
    axes = np.asarray(axes).reshape(-1)
    xy = sim["coords"]
    for ax, (name, labels) in zip(axes, panels):
        labels = np.asarray(labels).astype(str)
        uniques = {label: idx for idx, label in enumerate(sorted(np.unique(labels)))}
        codes = np.asarray([uniques[label] for label in labels])
        ax.scatter(xy[:, 0], xy[:, 1], c=codes, cmap="tab10", s=12, edgecolor="none")
        ax.set_title(name)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(panels) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figure_dir / f"highres_domain_maps_{scenario_name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data_dir = resolve_project_path(args.data_dir)
    output_dir = resolve_project_path(args.output_dir)
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    loaded = load_dataset(data_dir)
    apa = select_genes(loaded["apa"], args.n_genes, args.min_observed_spots)
    parent_spots = sample_parent_spots(
        apa.columns.astype(str).tolist(),
        loaded["metadata"],
        args.layer_column,
        args.max_parent_spots,
        args.seed,
    )
    apa = apa.loc[:, parent_spots]
    coords = loaded["coords"].loc[parent_spots]
    metadata = loaded["metadata"].loc[parent_spots]
    expression = loaded["expression"].reindex(columns=parent_spots) if loaded["expression"] is not None else None

    dataset_name = args.dataset_name
    if dataset_name is None:
        if "dataset" in metadata.columns and metadata["dataset"].notna().any():
            dataset_name = str(metadata["dataset"].dropna().astype(str).iloc[0])
        else:
            dataset_name = data_dir.name

    if args.n_domains is not None:
        n_domains = args.n_domains
    elif args.layer_column in metadata.columns:
        n_domains = int(metadata[args.layer_column].astype(str).nunique())
    else:
        n_domains = min(5, max(2, len(parent_spots) // 20))

    all_results = []
    last_sim = None
    last_maps = None
    for subbins in parse_int_list(args.subbins_per_spot):
        scenario_name = f"subbins{subbins}_capture{args.capture_rate:.2f}_dropout{args.dropout_rate:.2f}"
        print(f"\nRunning high-resolution scenario: {scenario_name}")
        sim = simulate_pseudo_bins(apa, expression, coords, metadata, subbins, args)
        outputs = run_methods(sim, args, n_domains)
        results, domain_maps = evaluate_outputs(outputs, sim, scenario_name, dataset_name, n_domains)
        results.to_csv(output_dir / f"{scenario_name}_results.csv", index=False)
        all_results.append(results)
        last_sim = sim
        last_maps = domain_maps

    results_df = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
    results_df.to_csv(output_dir / "highres_results_summary.csv", index=False)
    plot_method_comparison(results_df, figure_dir)
    plot_scaling(results_df, figure_dir)
    if last_sim is not None and last_maps is not None:
        plot_domain_maps(last_sim, last_maps, figure_dir, "last")

    decision = {
        "output_dir": str(output_dir),
        "data_dir": str(data_dir),
        "dataset": dataset_name,
        "n_parent_spots": int(len(parent_spots)),
        "n_genes": int(apa.shape[0]),
        "subbins_per_spot": parse_int_list(args.subbins_per_spot),
        "summary_table": str(output_dir / "highres_results_summary.csv"),
        "figures": sorted(str(path) for path in figure_dir.glob("*.png")),
    }
    if not results_df.empty and "layer_ari" in results_df.columns:
        best_ari = results_df.sort_values("layer_ari", ascending=False).iloc[0]
        best_rmse = results_df.sort_values("rmse_holdout", ascending=True).iloc[0]
        decision["best_layer_ari"] = best_ari.to_dict()
        decision["best_holdout_rmse"] = best_rmse.to_dict()
    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
    print(json.dumps(decision, indent=2, default=str))


if __name__ == "__main__":
    main()

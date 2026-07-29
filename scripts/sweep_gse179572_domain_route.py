#!/usr/bin/env python
"""Sweep lightweight domain routes on GSE179572 without rerunning GP imputation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD

from spagapa.bioml import BioMLDomainDetector, MultiViewGraphBuilder
from spagapa.imputation.feature_builders import ExpressionFeatureBuilder


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--benchmark-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--imputed-file", default="random_seed42_spagapa_fast_gp_imputed.npy")
    parser.add_argument("--layer-column", default="marker_weak_label")
    parser.add_argument("--n-domains", type=int, default=8)
    parser.add_argument("--expr-n-components", type=int, default=10)
    parser.add_argument("--neighbors", default="8,15,25,40")
    parser.add_argument(
        "--mode",
        default="quick",
        choices=["quick", "full"],
        help="quick avoids repeated spectral clustering and is suitable for iteration.",
    )
    return parser.parse_args()


def read_inputs(data_dir: Path, benchmark_dir: Path, imputed_file: str, expr_n_components: int):
    coords_df = pd.read_csv(data_dir / "coordinates.csv").set_index("spot_id")
    metadata = pd.read_csv(data_dir / "metadata.csv").set_index("spot_id").loc[coords_df.index]
    expression = pd.read_csv(data_dir / "expression_matrix.csv", index_col=0)
    expression = expression.reindex(columns=coords_df.index).fillna(0.0)
    expression_embedding = ExpressionFeatureBuilder(
        n_components=expr_n_components,
        use_hvg=True,
        n_top_genes=2000,
        orientation="genes_by_spots",
    ).fit_transform(expression)
    coords = coords_df[["x", "y"]].values.astype(float)
    gp_imputed = np.load(benchmark_dir / imputed_file)
    return coords, expression_embedding, gp_imputed, metadata


def evaluate(labels: np.ndarray, layer_labels: np.ndarray) -> dict[str, float]:
    labels = np.asarray(labels)
    layer_labels = np.asarray(layer_labels).astype(str)
    layer_codes = pd.Categorical(layer_labels).codes
    out = {
        "layer_ari": float(adjusted_rand_score(layer_codes, labels)),
        "layer_nmi": float(normalized_mutual_info_score(layer_codes, labels)),
        "n_domains_observed": int(len(np.unique(labels))),
    }
    keep = layer_labels != "ambiguous"
    if keep.sum() >= 2 and len(np.unique(layer_labels[keep])) > 1:
        keep_codes = pd.Categorical(layer_labels[keep]).codes
        out["nonambig_layer_ari"] = float(adjusted_rand_score(keep_codes, labels[keep]))
        out["nonambig_layer_nmi"] = float(normalized_mutual_info_score(keep_codes, labels[keep]))
        out["nonambig_n_spots"] = int(keep.sum())
    else:
        out["nonambig_layer_ari"] = np.nan
        out["nonambig_layer_nmi"] = np.nan
        out["nonambig_n_spots"] = int(keep.sum())
    return out


def kmeans_labels(features: np.ndarray, n_domains: int, random_state: int = 42) -> np.ndarray:
    X = StandardScaler().fit_transform(np.asarray(features, dtype=float))
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return KMeans(n_clusters=n_domains, random_state=random_state, n_init=20).fit_predict(X)


def graph_svd_kmeans_labels(graph, n_domains: int, n_components: int = 16) -> np.ndarray:
    n_components = max(n_domains, min(n_components, graph.shape[0] - 1))
    embedding = TruncatedSVD(n_components=n_components, random_state=42).fit_transform(graph)
    return kmeans_labels(embedding, n_domains)


def spatial_knn_spectral(coords: np.ndarray, n_neighbors: int, n_domains: int) -> np.ndarray:
    graph = kneighbors_graph(
        StandardScaler().fit_transform(coords),
        n_neighbors=n_neighbors,
        mode="distance",
        include_self=False,
    )
    distances = graph.data
    sigma = float(np.median(distances[distances > 0])) if np.any(distances > 0) else 1.0
    graph.data = np.exp(-0.5 * (graph.data / max(sigma, 1e-6)) ** 2)
    graph = graph.maximum(graph.T).tocsr()
    return BioMLDomainDetector(
        method="spectral",
        n_domains=n_domains,
        random_state=42,
    ).fit_predict(graph=graph)


def plot(rows: pd.DataFrame, output_dir: Path) -> None:
    top = rows.sort_values(["layer_ari", "nonambig_layer_ari"], ascending=False).head(20)
    labels = top["route"] + "\n" + top["params"]
    x = np.arange(len(top))
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    axes[0].bar(x, top["layer_ari"], color="#2c7fb8")
    axes[0].set_title("All Spots ARI")
    axes[1].bar(x, top["nonambig_layer_ari"], color="#f03b20")
    axes[1].set_title("Non-ambiguous Spots ARI")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=7)
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "domain_route_sweep_top20.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    benchmark_dir = Path(args.benchmark_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    coords, expr, gp_imputed, metadata = read_inputs(
        data_dir,
        benchmark_dir,
        args.imputed_file,
        args.expr_n_components,
    )
    layer_labels = metadata[args.layer_column].astype(str).values
    neighbors = [int(x) for x in args.neighbors.split(",") if x.strip()]

    rows: list[dict[str, object]] = []

    baseline_features = {
        "coords_kmeans": coords,
        "expr_kmeans": expr,
        "gp_apa_kmeans": gp_imputed.T,
        "coords_expr_kmeans": np.column_stack([coords, expr]),
        "coords_expr_gpapa_kmeans": np.column_stack([coords, expr, gp_imputed.T]),
    }
    for name, features in baseline_features.items():
        labels = kmeans_labels(features, args.n_domains)
        row = {"route": name, "params": "kmeans"}
        row.update(evaluate(labels, layer_labels))
        rows.append(row)

    if args.mode == "full":
        for k in neighbors:
            labels = spatial_knn_spectral(coords, k, args.n_domains)
            row = {"route": "spatial_graph", "params": f"k={k},spectral"}
            row.update(evaluate(labels, layer_labels))
            rows.append(row)

    weight_grid = [
        (1.0, 0.0, 0.0),
        (0.8, 0.2, 0.0),
        (0.7, 0.3, 0.0),
        (0.6, 0.4, 0.0),
        (0.5, 0.5, 0.0),
        (0.4, 0.6, 0.0),
        (0.3, 0.7, 0.0),
        (0.2, 0.8, 0.0),
        (0.1, 0.9, 0.0),
        (0.1, 0.7, 0.2),
        (0.2, 0.6, 0.2),
        (0.3, 0.5, 0.2),
        (0.5, 0.3, 0.2),
        (0.7, 0.1, 0.2),
        (0.8, 0.0, 0.2),
    ]
    if args.mode == "quick":
        weight_grid = [
            (1.0, 0.0, 0.0),
            (0.8, 0.2, 0.0),
            (0.6, 0.4, 0.0),
            (0.4, 0.6, 0.0),
            (0.2, 0.8, 0.0),
            (0.1, 0.7, 0.2),
            (0.5, 0.3, 0.2),
        ]
        neighbors = [k for k in neighbors if k in {8, 15, 25}]

    for k in neighbors:
        for sw, ew, aw in weight_grid:
            graph = MultiViewGraphBuilder(
                n_neighbors=k,
                spatial_weight=sw,
                expression_weight=ew,
                apa_weight=aw,
            ).build(
                coords,
                expression_embedding=expr,
                apa_matrix=gp_imputed,
                uncertainty=None,
            )
            if args.mode == "full":
                labels = BioMLDomainDetector(
                    method="spectral",
                    n_domains=args.n_domains,
                    random_state=42,
                ).fit_predict(graph=graph.fused)
                route = "multiview_spectral"
            else:
                labels = graph_svd_kmeans_labels(graph.fused, args.n_domains)
                route = "multiview_graph_svd_kmeans"
            row = {
                "route": route,
                "params": f"k={k},s={sw:g},e={ew:g},a={aw:g}",
                "spatial_weight": sw,
                "expression_weight": ew,
                "apa_weight": aw,
                "n_neighbors": k,
            }
            row.update(evaluate(labels, layer_labels))
            rows.append(row)

    results = pd.DataFrame(rows).sort_values(["layer_ari", "nonambig_layer_ari"], ascending=False)
    results.to_csv(output_dir / "domain_route_sweep.csv", index=False)
    results.head(25).to_csv(output_dir / "domain_route_sweep_top25.csv", index=False)
    plot(results, output_dir)
    print(results.head(15).to_string(index=False))


if __name__ == "__main__":
    main()

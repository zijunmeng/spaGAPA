#!/usr/bin/env python
"""External real-data validation for spaGAPA-BioML.

The runner is designed around the experiments used by stAPAminer and metaAPA:

* stAPAminer-style biological consistency:
  layer/domain recovery, internal clustering metrics, within-layer APA
  correlation, SVAPA/LSAPA/DEAPA counts, and method-wise spatial maps.
* metaAPA-style upstream readiness:
  a site-level validation checklist for multi-caller poly(A) site integration
  and sequence/long-read support. This runner does not force metaAPA into a
  downstream imputation comparison, because its main contribution is upstream
  poly(A) site consensus.

Required dataset files
----------------------
apa_matrix.csv          genes x spots APA usage matrix, with NaN for missing
coordinates.csv         columns spot_id/barcode, x, y
metadata.csv            spot-level annotations; layer column is optional

Optional files
--------------
expression_matrix.csv   genes x spots expression matrix for BioML/KNN
stapaminer_rud_imputed.csv
apa_sites.csv.gz        site annotations for readiness checks
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import Any

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    pairwise_distances,
    silhouette_score,
)
from sklearn.preprocessing import LabelEncoder, StandardScaler

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from spagapa import APADataset, SpaGAPA
from spagapa.analysis import DifferentialAPAAnalyzer, SpatialPatternAnalyzer, find_domain_markers
from spagapa.imputation import ExpressionFeatureBuilder

from run_stapaminer_mob_benchmark import impute_stapaminer_expression_knn


METHOD_ORDER = [
    "raw",
    "stapaminer_original_imputed",
    "stapaminer_knn_expression",
    "spagapa_gp",
    "spagapa_bioml",
]

METHOD_COLORS = {
    "raw": "#7f8c8d",
    "stapaminer_original_imputed": "#6a51a3",
    "stapaminer_knn_expression": "#9b59b6",
    "spagapa_gp": "#e74c3c",
    "spagapa_bioml": "#d35400",
}


def resolve_project_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        action="append",
        default=None,
        help=(
            "Prepared dataset directory. Can be repeated. Default: "
            "spaGAPA/data/processed/stapaminer_mob"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/external_bioml_validation_v1",
        help="Output directory.",
    )
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--n-genes", type=int, default=120)
    parser.add_argument("--min-observed-spots", type=int, default=80)
    parser.add_argument(
        "--methods",
        default="raw,stapaminer_original_imputed,stapaminer_knn_expression,spagapa_gp,spagapa_bioml",
        help="Comma-separated methods to run.",
    )
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--gp-kernel", default="matern", choices=["rbf", "matern", "auto"])
    parser.add_argument("--gp-alpha", type=float, default=1e-3)
    parser.add_argument("--gp-n-restarts", type=int, default=1)
    parser.add_argument("--n-domains", type=int, default=None)
    parser.add_argument("--expr-n-components", type=int, default=10)
    parser.add_argument("--svapa-neighbors", type=int, default=6)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--logfc-threshold", type=float, default=0.5)
    parser.add_argument("--skip-downstream", action="store_true")
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
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output tables for the selected output directory.",
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def load_dataset(data_dir: Path) -> dict[str, Any]:
    apa = read_table(data_dir / "apa_matrix.csv")
    coords = pd.read_csv(data_dir / "coordinates.csv")
    if "spot_id" in coords.columns:
        coords = coords.set_index("spot_id")
    elif "barcode" in coords.columns:
        coords = coords.set_index("barcode")
    else:
        coords = coords.set_index(coords.columns[0])
    coords.index = coords.index.astype(str)
    coords = coords[["x", "y"]]

    metadata_path = data_dir / "metadata.csv"
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path)
        if "spot_id" in metadata.columns:
            metadata = metadata.set_index("spot_id")
        elif "barcode" in metadata.columns:
            metadata = metadata.set_index("barcode")
        else:
            metadata = metadata.set_index(metadata.columns[0])
        metadata.index = metadata.index.astype(str)
    else:
        metadata = pd.DataFrame(index=coords.index)

    common_spots = [spot for spot in apa.columns.astype(str) if spot in coords.index]
    if len(common_spots) < 3:
        raise ValueError(f"{data_dir} has fewer than 3 aligned spots")

    apa.columns = apa.columns.astype(str)
    apa = apa.loc[:, common_spots]
    coords = coords.loc[common_spots]
    metadata = metadata.reindex(common_spots)

    expression = None
    expr_path = data_dir / "expression_matrix.csv"
    if expr_path.exists():
        expression = read_table(expr_path)
        expression.columns = expression.columns.astype(str)
        expression = expression.reindex(columns=common_spots)

    stapa_imputed = None
    stapa_path = data_dir / "stapaminer_rud_imputed.csv"
    if stapa_path.exists():
        stapa_imputed = read_table(stapa_path)
        stapa_imputed.columns = stapa_imputed.columns.astype(str)

    return {
        "apa": apa,
        "coords": coords,
        "metadata": metadata,
        "expression": expression,
        "stapaminer_imputed": stapa_imputed,
    }


def select_genes(apa: pd.DataFrame, n_genes: int, min_observed_spots: int) -> pd.DataFrame:
    observed = apa.notna().sum(axis=1)
    eligible = observed[observed >= min_observed_spots]
    if eligible.empty:
        eligible = observed[observed >= max(3, min_observed_spots // 2)]
    selected = eligible.sort_values(ascending=False).head(n_genes).index
    return apa.loc[selected]


def fill_missing_by_gene_mean(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=float).copy()
    finite = np.isfinite(values)
    counts = finite.sum(axis=1)
    sums = np.where(finite, values, 0.0).sum(axis=1)
    means = np.divide(sums, np.maximum(counts, 1), where=np.maximum(counts, 1) > 0)
    means[counts == 0] = 0.0
    missing_gene, missing_spot = np.where(~finite)
    values[missing_gene, missing_spot] = means[missing_gene]
    return np.clip(values, 0.0, 1.0)


def build_expression_embedding(
    expression: pd.DataFrame | None,
    spot_names: list[str],
    n_components: int,
) -> np.ndarray | None:
    if expression is None:
        return None
    expr = expression.reindex(columns=spot_names).fillna(0.0)
    builder = ExpressionFeatureBuilder(
        n_components=n_components,
        orientation="genes_by_spots",
    )
    return builder.fit_transform(expr)


def encode_layers(metadata: pd.DataFrame, layer_column: str) -> tuple[np.ndarray | None, np.ndarray | None]:
    if layer_column not in metadata.columns:
        return None, None
    labels = metadata[layer_column].astype(str).fillna("NA").values
    if len(np.unique(labels)) < 2:
        return labels, None
    codes = LabelEncoder().fit_transform(labels)
    return labels, codes


def safe_metric(fn, *args, default: float = np.nan) -> float:
    try:
        return float(fn(*args))
    except Exception:
        return default


def pairwise_cluster_jaccard(true_codes: np.ndarray, pred: np.ndarray) -> float:
    if len(true_codes) < 2:
        return np.nan
    true_same = true_codes[:, None] == true_codes[None, :]
    pred_same = pred[:, None] == pred[None, :]
    tri = np.triu_indices(len(true_codes), k=1)
    true_flat = true_same[tri]
    pred_flat = pred_same[tri]
    union = np.logical_or(true_flat, pred_flat).sum()
    if union == 0:
        return np.nan
    return float(np.logical_and(true_flat, pred_flat).sum() / union)


def purity_score(true_codes: np.ndarray, pred: np.ndarray) -> float:
    total = len(true_codes)
    if total == 0:
        return np.nan
    score = 0
    for cluster in np.unique(pred):
        idx = pred == cluster
        if idx.sum() == 0:
            continue
        counts = np.bincount(true_codes[idx])
        score += int(counts.max())
    return float(score / total)


def dunn_index(features: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan
    distances = squareform(pdist(features, metric="euclidean"))
    max_intra = 0.0
    min_inter = np.inf
    for label in unique:
        idx = np.where(labels == label)[0]
        if len(idx) > 1:
            max_intra = max(max_intra, float(distances[np.ix_(idx, idx)].max()))
    for i, label_i in enumerate(unique):
        idx_i = np.where(labels == label_i)[0]
        for label_j in unique[i + 1 :]:
            idx_j = np.where(labels == label_j)[0]
            min_inter = min(min_inter, float(distances[np.ix_(idx_i, idx_j)].min()))
    if max_intra <= 0 or not np.isfinite(min_inter):
        return np.nan
    return float(min_inter / max_intra)


def cluster_from_matrix(matrix_gs: np.ndarray, n_domains: int) -> tuple[np.ndarray, np.ndarray]:
    spot_features = fill_missing_by_gene_mean(matrix_gs).T
    spot_features = StandardScaler().fit_transform(spot_features)
    spot_features = np.nan_to_num(spot_features, nan=0.0)
    n_clusters = min(max(2, n_domains), spot_features.shape[0] - 1)
    labels = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(spot_features)
    return labels, spot_features


def evaluate_layer_metrics(
    matrix_gs: np.ndarray,
    true_codes: np.ndarray | None,
    n_domains: int,
    domain_labels: np.ndarray | None = None,
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    if domain_labels is None:
        pred, features = cluster_from_matrix(matrix_gs, n_domains)
        domain_source = "kmeans_on_matrix"
    else:
        pred = np.asarray(domain_labels)
        features = fill_missing_by_gene_mean(matrix_gs).T
        features = StandardScaler().fit_transform(features)
        features = np.nan_to_num(features, nan=0.0)
        domain_source = "model_domain"

    metrics: dict[str, float] = {
        "domain_source": domain_source,
        "n_domains_pred": int(len(np.unique(pred))),
        "layer_ari": np.nan,
        "layer_nmi": np.nan,
        "layer_purity": np.nan,
        "layer_pairwise_jaccard": np.nan,
        "dbi": np.nan,
        "calinski_harabasz": np.nan,
        "silhouette": np.nan,
        "dunn": np.nan,
    }
    if true_codes is not None:
        metrics["layer_ari"] = float(adjusted_rand_score(true_codes, pred))
        metrics["layer_nmi"] = float(normalized_mutual_info_score(true_codes, pred))
        metrics["layer_purity"] = purity_score(true_codes, pred)
        metrics["layer_pairwise_jaccard"] = pairwise_cluster_jaccard(true_codes, pred)

    if len(np.unique(pred)) > 1 and len(np.unique(pred)) < len(pred):
        metrics["dbi"] = safe_metric(davies_bouldin_score, features, pred)
        metrics["calinski_harabasz"] = safe_metric(calinski_harabasz_score, features, pred)
        metrics["silhouette"] = safe_metric(silhouette_score, features, pred)
        metrics["dunn"] = dunn_index(features, pred)
    return metrics, pred, features


def within_layer_correlations(matrix_gs: np.ndarray, layer_labels: np.ndarray | None) -> list[dict[str, Any]]:
    if layer_labels is None:
        return []
    filled = fill_missing_by_gene_mean(matrix_gs)
    rows = []
    for layer in sorted(np.unique(layer_labels)):
        idx = np.where(layer_labels == layer)[0]
        if len(idx) < 2:
            continue
        spot_profiles = filled[:, idx].T
        corr = np.corrcoef(spot_profiles)
        tri = np.triu_indices(len(idx), k=1)
        values = corr[tri]
        values = values[np.isfinite(values)]
        rows.append(
            {
                "layer": str(layer),
                "n_spots": int(len(idx)),
                "median_pearson": float(np.median(values)) if values.size else np.nan,
                "mean_pearson": float(np.mean(values)) if values.size else np.nan,
            }
        )
    return rows


def track_runtime_memory(fn):
    process = psutil.Process(os.getpid())
    peak_rss_mb = process.memory_info().rss / 1024 / 1024
    stop_event = threading.Event()

    def monitor() -> None:
        nonlocal peak_rss_mb
        while not stop_event.is_set():
            peak_rss_mb = max(peak_rss_mb, process.memory_info().rss / 1024 / 1024)
            time.sleep(0.02)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    t0 = time.perf_counter()
    try:
        out = fn()
    finally:
        runtime = time.perf_counter() - t0
        stop_event.set()
        thread.join(timeout=0.2)
    return out, float(runtime), float(peak_rss_mb)


def run_pipeline_method(
    method: str,
    values: np.ndarray,
    coords: np.ndarray,
    gene_names: list[str],
    spot_names: list[str],
    expression: pd.DataFrame | None,
    args: argparse.Namespace,
    n_domains: int,
):
    dataset = APADataset.from_counts(
        values,
        coords,
        gene_names=gene_names,
        spot_names=spot_names,
    )
    pipeline = SpaGAPA(
        n_neighbors=args.svapa_neighbors,
        kernel_type=args.gp_kernel,
        gp_alpha=args.gp_alpha,
        gp_n_restarts_optimizer=args.gp_n_restarts,
        input_type="apa_index",
        min_spots=3,
        use_bioml=(method == "spagapa_bioml"),
        bioml_rank=args.bioml_rank,
        bioml_lambda_graph=args.bioml_lambda_graph,
        bioml_lambda_l2=args.bioml_lambda_l2,
        bioml_max_iter=args.bioml_max_iter,
        bioml_n_neighbors=args.bioml_n_neighbors,
        bioml_blend=args.bioml_blend,
        bioml_domain_method=args.bioml_domain_method,
        bioml_spatial_weight=args.bioml_spatial_weight,
        bioml_expression_weight=args.bioml_expression_weight,
        bioml_apa_weight=args.bioml_apa_weight,
        expression_n_components=args.expr_n_components,
        verbose=False,
    )
    results = pipeline.run(
        dataset=dataset,
        expression_matrix=expression,
        expression_orientation="genes_by_spots",
        impute=True,
        quantify=True,
        identify_domains=True,
        use_bioml=(method == "spagapa_bioml"),
        differential_analysis=False,
        detect_svapa=False,
        n_domains=n_domains,
        use_uncertainty_weights=True,
    )
    matrix = results["domains"]["imputed_values"] if method == "spagapa_bioml" else results["imputed_values"]
    uncertainty = results.get("uncertainty")
    domains = results["domains"]["labels"]
    bioml_metadata = None
    if method == "spagapa_bioml":
        bioml_metadata = results["domains"].get("metadata")
    return np.clip(matrix, 0.0, 1.0), uncertainty, domains, bioml_metadata


def run_methods(
    selected: pd.DataFrame,
    coords: pd.DataFrame,
    metadata: pd.DataFrame,
    expression: pd.DataFrame | None,
    stapa_original: pd.DataFrame | None,
    args: argparse.Namespace,
    n_domains: int,
) -> dict[str, dict[str, Any]]:
    values = selected.values.astype(float)
    gene_names = selected.index.astype(str).tolist()
    spot_names = selected.columns.astype(str).tolist()
    coords_arr = coords.loc[spot_names][["x", "y"]].values.astype(float)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    outputs: dict[str, dict[str, Any]] = {}
    if "raw" in methods:
        outputs["raw"] = {
            "matrix": fill_missing_by_gene_mean(values),
            "uncertainty": None,
            "domains": None,
            "runtime_s": 0.0,
            "peak_rss_mb": float(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024),
        }

    if "stapaminer_original_imputed" in methods and stapa_original is not None:
        aligned = stapa_original.reindex(index=gene_names, columns=spot_names)
        if aligned.notna().sum().sum() > 0:
            outputs["stapaminer_original_imputed"] = {
                "matrix": fill_missing_by_gene_mean(aligned.values.astype(float)),
                "uncertainty": None,
                "domains": None,
                "runtime_s": 0.0,
                "peak_rss_mb": float(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024),
            }

    if "stapaminer_knn_expression" in methods:
        if expression is None:
            print("Skipping stapaminer_knn_expression: expression_matrix.csv is missing")
        else:
            def run_knn():
                return impute_stapaminer_expression_knn(
                    values.copy(),
                    gene_names,
                    spot_names,
                    expression,
                    args.knn_k,
                )

            matrix, runtime_s, peak_rss_mb = track_runtime_memory(run_knn)
            outputs["stapaminer_knn_expression"] = {
                "matrix": matrix,
                "uncertainty": None,
                "domains": None,
                "runtime_s": runtime_s,
                "peak_rss_mb": peak_rss_mb,
            }

    for method in ["spagapa_gp", "spagapa_bioml"]:
        if method not in methods:
            continue
        if method == "spagapa_bioml" and expression is None:
            print("Skipping spagapa_bioml: expression_matrix.csv is missing")
            continue

        def run_spagapa():
            return run_pipeline_method(
                method,
                values.copy(),
                coords_arr,
                gene_names,
                spot_names,
                expression,
                args,
                n_domains,
            )

        (matrix, uncertainty, domains, bioml_metadata), runtime_s, peak_rss_mb = track_runtime_memory(run_spagapa)
        outputs[method] = {
            "matrix": matrix,
            "uncertainty": uncertainty,
            "domains": domains,
            "runtime_s": runtime_s,
            "peak_rss_mb": peak_rss_mb,
            "bioml_metadata": bioml_metadata,
        }

    return outputs


def run_downstream(
    method: str,
    matrix_gs: np.ndarray,
    uncertainty: np.ndarray | None,
    coords: np.ndarray,
    gene_names: list[str],
    true_layer_labels: np.ndarray | None,
    true_layer_codes: np.ndarray | None,
    pred_domains: np.ndarray,
    out_dir: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    if args.skip_downstream:
        return {
            "n_svapa": np.nan,
            "n_lsapa": np.nan,
            "n_deapa": np.nan,
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    filled = fill_missing_by_gene_mean(matrix_gs)

    spatial = SpatialPatternAnalyzer(alpha=args.fdr, random_state=42)
    svapa = spatial.test_spatial_autocorrelation(
        filled,
        coords,
        gene_names=gene_names,
        n_neighbors=args.svapa_neighbors,
    )
    if not svapa.empty:
        from statsmodels.stats.multitest import multipletests

        svapa = svapa.copy()
        _, svapa["padj"], _, _ = multipletests(svapa["pvalue"].values, method="fdr_bh")
        svapa["significant"] = (svapa["padj"] < args.fdr) & (svapa["morans_i"] > 0)
        svapa = svapa.sort_values(["padj", "morans_i"], ascending=[True, False])
    svapa.to_csv(out_dir / f"{method}_svapa_genes.csv", index=False)
    n_svapa = int(svapa["significant"].sum()) if "significant" in svapa.columns else 0

    markers = find_domain_markers(
        filled,
        pred_domains,
        gene_names=gene_names,
        method="wilcoxon",
        padj_threshold=args.fdr,
        logfc_threshold=args.logfc_threshold,
        uncertainty=uncertainty,
    )
    marker_rows = []
    for domain, df in markers.items():
        tmp = df.copy()
        tmp.insert(0, "domain", domain)
        marker_rows.append(tmp)
    markers_df = pd.concat(marker_rows, ignore_index=True) if marker_rows else pd.DataFrame()
    markers_df.to_csv(out_dir / f"{method}_lsapa_domain_markers.csv", index=False)
    n_lsapa = int(markers_df["gene"].nunique()) if "gene" in markers_df.columns else 0

    deapa_rows = []
    n_deapa = np.nan
    if true_layer_codes is not None:
        analyzer = DifferentialAPAAnalyzer(method="wilcoxon", min_spots_per_group=3)
        pairwise = analyzer.test_all_pairwise(
            filled,
            true_layer_labels if true_layer_labels is not None else true_layer_codes,
            gene_names=gene_names,
            uncertainty=uncertainty,
        )
        sig_genes = set()
        for comp, df in pairwise.items():
            adj = analyzer.adjust_pvalues(df)
            filt = analyzer.filter_results(
                adj,
                padj_threshold=args.fdr,
                logfc_threshold=args.logfc_threshold,
            )
            if not filt.empty:
                filt = filt.copy()
                filt.insert(0, "comparison", comp)
                deapa_rows.append(filt)
                sig_genes.update(filt["gene"].astype(str).tolist())
        n_deapa = len(sig_genes)
    deapa_df = pd.concat(deapa_rows, ignore_index=True) if deapa_rows else pd.DataFrame()
    deapa_df.to_csv(out_dir / f"{method}_deapa_layer_pairwise.csv", index=False)

    return {
        "n_svapa": n_svapa,
        "n_lsapa": n_lsapa,
        "n_deapa": int(n_deapa) if np.isfinite(n_deapa) else np.nan,
    }


def build_site_level_checklist(data_dir: Path) -> dict[str, Any]:
    files = {path.name for path in data_dir.iterdir() if path.is_file()}
    checklist = {
        "data_dir": str(data_dir),
        "apa_sites_available": "apa_sites.csv" in files or "apa_sites.csv.gz" in files,
        "site_counts_available": "apa_site_counts.csv" in files or "apa_site_counts.csv.gz" in files,
        "multi_caller_site_tables_available": False,
        "sequence_context_available": False,
        "long_read_support_available": False,
        "reference_pas_annotation_available": False,
        "metaapa_ready": False,
        "notes": [
            "metaAPA evaluates upstream poly(A) site consensus rather than downstream APA imputation.",
            "To reproduce metaAPA-style validation, add Sierra/polyApipe/SCAPE site tables, genome FASTA/GTF, and optional long-read support.",
        ],
    }
    caller_markers = ["sierra", "polyapipe", "scape", "metaapa"]
    checklist["multi_caller_site_tables_available"] = any(
        any(marker in name.lower() for marker in caller_markers) for name in files
    )
    checklist["metaapa_ready"] = bool(
        checklist["apa_sites_available"]
        and checklist["site_counts_available"]
        and checklist["multi_caller_site_tables_available"]
        and checklist["sequence_context_available"]
    )
    return checklist


def plot_layer_metrics(summary: pd.DataFrame, figure_dir: Path) -> None:
    metrics = ["layer_ari", "layer_nmi", "layer_purity", "layer_pairwise_jaccard"]
    available = [m for m in metrics if m in summary.columns and summary[m].notna().any()]
    if not available:
        return
    fig, axes = plt.subplots(1, len(available), figsize=(4.0 * len(available), 4.0))
    if len(available) == 1:
        axes = [axes]
    methods = [m for m in METHOD_ORDER if m in summary["method"].values]
    x = np.arange(len(methods))
    for ax, metric in zip(axes, available):
        values = [
            float(summary.loc[summary["method"] == method, metric].mean())
            for method in methods
        ]
        ax.bar(x, values, color=[METHOD_COLORS.get(m, "#333333") for m in methods])
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_title(metric)
        ax.set_ylim(0, max(1.0, np.nanmax(values) * 1.1 if values else 1.0))
    fig.tight_layout()
    fig.savefig(figure_dir / "layer_separation_metrics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_internal_metrics(summary: pd.DataFrame, figure_dir: Path) -> None:
    metrics = ["dbi", "calinski_harabasz", "silhouette", "dunn"]
    available = [m for m in metrics if m in summary.columns and summary[m].notna().any()]
    if not available:
        return
    fig, axes = plt.subplots(1, len(available), figsize=(4.0 * len(available), 4.0))
    if len(available) == 1:
        axes = [axes]
    methods = [m for m in METHOD_ORDER if m in summary["method"].values]
    x = np.arange(len(methods))
    for ax, metric in zip(axes, available):
        values = [
            float(summary.loc[summary["method"] == method, metric].mean())
            for method in methods
        ]
        ax.bar(x, values, color=[METHOD_COLORS.get(m, "#333333") for m in methods])
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
        ax.set_title(metric)
    fig.tight_layout()
    fig.savefig(figure_dir / "internal_clustering_metrics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_downstream_counts(summary: pd.DataFrame, figure_dir: Path) -> None:
    metrics = ["n_svapa", "n_lsapa", "n_deapa"]
    available = [m for m in metrics if m in summary.columns and summary[m].notna().any()]
    if not available:
        return
    methods = [m for m in METHOD_ORDER if m in summary["method"].values]
    x = np.arange(len(methods))
    width = 0.8 / len(available)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, metric in enumerate(available):
        values = [
            float(summary.loc[summary["method"] == method, metric].mean())
            for method in methods
        ]
        ax.bar(x + (i - (len(available) - 1) / 2) * width, values, width=width, label=metric)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_ylabel("Gene count")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "spatial_apa_gene_counts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_within_layer_correlation(corr_df: pd.DataFrame, figure_dir: Path) -> None:
    if corr_df.empty:
        return
    methods = [m for m in METHOD_ORDER if m in corr_df["method"].values]
    layers = sorted(corr_df["layer"].astype(str).unique())
    fig, ax = plt.subplots(figsize=(max(7, len(layers) * 1.2), 4.5))
    width = 0.8 / max(len(methods), 1)
    x = np.arange(len(layers))
    for i, method in enumerate(methods):
        vals = []
        for layer in layers:
            sub = corr_df[(corr_df["method"] == method) & (corr_df["layer"].astype(str) == layer)]
            vals.append(float(sub["median_pearson"].mean()) if not sub.empty else np.nan)
        ax.bar(x + (i - (len(methods) - 1) / 2) * width, vals, width=width, label=method)
    ax.set_xticks(x)
    ax.set_xticklabels(layers)
    ax.set_ylabel("Median within-layer Pearson r")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "within_layer_correlation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_domain_maps(
    coords: pd.DataFrame,
    layer_labels: np.ndarray | None,
    domain_maps: dict[str, np.ndarray],
    figure_dir: Path,
) -> None:
    panels = []
    if layer_labels is not None:
        panels.append(("true_layer", layer_labels))
    panels.extend((method, labels) for method, labels in domain_maps.items())
    if not panels:
        return
    n_cols = min(3, len(panels))
    n_rows = int(np.ceil(len(panels) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 4.0 * n_rows))
    axes = np.asarray(axes).reshape(-1)
    xy = coords[["x", "y"]].values
    for ax, (name, labels) in zip(axes, panels):
        codes = LabelEncoder().fit_transform(np.asarray(labels).astype(str))
        scatter = ax.scatter(xy[:, 0], xy[:, 1], c=codes, cmap="tab10", s=26, edgecolor="none")
        ax.set_title(name)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(panels) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figure_dir / "method_domain_maps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_experiment_coverage(checklists: list[dict[str, Any]], figure_dir: Path) -> None:
    rows = []
    for item in checklists:
        rows.append(
            {
                "dataset": Path(item["data_dir"]).name,
                "APA sites": item["apa_sites_available"],
                "Site counts": item["site_counts_available"],
                "Multi-caller": item["multi_caller_site_tables_available"],
                "Sequence": item["sequence_context_available"],
                "Long-read": item["long_read_support_available"],
                "metaAPA-ready": item["metaapa_ready"],
            }
        )
    if not rows:
        return
    df = pd.DataFrame(rows).set_index("dataset")
    bool_df = df.astype(bool).astype(int)
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.7 * len(bool_df) + 1.5)))
    ax.imshow(bool_df.values, cmap="Greens", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(bool_df.shape[1]))
    ax.set_xticklabels(bool_df.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(bool_df.shape[0]))
    ax.set_yticklabels(bool_df.index)
    for i in range(bool_df.shape[0]):
        for j in range(bool_df.shape[1]):
            ax.text(j, i, "Y" if bool_df.iat[i, j] else "N", ha="center", va="center", fontsize=9)
    ax.set_title("metaAPA-style site-level readiness")
    fig.tight_layout()
    fig.savefig(figure_dir / "ref_package_experiment_coverage.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_replicate_overlap(method_gene_sets: dict[str, dict[str, set[str]]]) -> pd.DataFrame:
    rows = []
    datasets = sorted(method_gene_sets)
    for i, d1 in enumerate(datasets):
        for d2 in datasets[i + 1 :]:
            for method in sorted(set(method_gene_sets[d1]) | set(method_gene_sets[d2])):
                s1 = method_gene_sets[d1].get(method, set())
                s2 = method_gene_sets[d2].get(method, set())
                union = s1 | s2
                rows.append(
                    {
                        "dataset_1": d1,
                        "dataset_2": d2,
                        "method": method,
                        "shared_spatial_genes": len(s1 & s2),
                        "union_spatial_genes": len(union),
                        "jaccard": len(s1 & s2) / len(union) if union else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def run_one_dataset(
    data_dir: Path,
    output_dir: Path,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, set[str]]], dict[str, Any]]:
    loaded = load_dataset(data_dir)
    apa = loaded["apa"]
    metadata = loaded["metadata"]
    coords = loaded["coords"]
    expression = loaded["expression"]
    stapa_original = loaded["stapaminer_imputed"]

    selected = select_genes(apa, args.n_genes, args.min_observed_spots)
    spot_names = selected.columns.astype(str).tolist()
    coords = coords.loc[spot_names]
    metadata = metadata.loc[spot_names]
    layer_labels, layer_codes = encode_layers(metadata, args.layer_column)
    if args.n_domains is not None:
        n_domains = args.n_domains
    elif layer_codes is not None:
        n_domains = len(np.unique(layer_codes))
    else:
        n_domains = min(5, max(2, selected.shape[1] // 20))

    dataset_name = args.dataset_name or metadata.get("dataset", pd.Series([data_dir.name])).dropna().astype(str).iloc[0]
    dataset_out = output_dir / str(dataset_name)
    downstream_dir = dataset_out / "downstream"
    dataset_out.mkdir(parents=True, exist_ok=True)

    outputs = run_methods(
        selected,
        coords,
        metadata,
        expression,
        stapa_original,
        args,
        n_domains,
    )

    summary_rows = []
    corr_rows = []
    domain_maps = {}
    method_gene_sets: dict[str, set[str]] = {}
    coords_arr = coords[["x", "y"]].values.astype(float)
    gene_names = selected.index.astype(str).tolist()

    for method in [m for m in METHOD_ORDER if m in outputs]:
        out = outputs[method]
        matrix = out["matrix"]
        layer_metrics, pred_domains, _ = evaluate_layer_metrics(
            matrix,
            layer_codes,
            n_domains,
            domain_labels=out.get("domains"),
        )
        domain_maps[method] = pred_domains

        downstream = run_downstream(
            method,
            matrix,
            out.get("uncertainty"),
            coords_arr,
            gene_names,
            layer_labels,
            layer_codes,
            pred_domains,
            downstream_dir,
            args,
        )

        svapa_file = downstream_dir / f"{method}_svapa_genes.csv"
        if svapa_file.exists():
            svapa = pd.read_csv(svapa_file)
            if "significant" in svapa.columns:
                method_gene_sets[method] = set(svapa.loc[svapa["significant"], "gene"].astype(str))
            else:
                method_gene_sets[method] = set()

        corr = within_layer_correlations(matrix, layer_labels)
        for row in corr:
            row.update({"dataset": dataset_name, "method": method})
            corr_rows.append(row)

        row = {
            "dataset": dataset_name,
            "data_dir": str(data_dir),
            "method": method,
            "n_genes": int(selected.shape[0]),
            "n_spots": int(selected.shape[1]),
            "n_observed": int(np.isfinite(selected.values).sum()),
            "observed_fraction": float(np.isfinite(selected.values).mean()),
            "runtime_s": out["runtime_s"],
            "peak_rss_mb": out["peak_rss_mb"],
            "n_domains_target": int(n_domains),
            **layer_metrics,
            **downstream,
        }
        if out.get("bioml_metadata") is not None:
            row["bioml_effective_spatial_weight"] = out["bioml_metadata"]["graph_weights_effective"].get("spatial", np.nan)
            row["bioml_effective_expression_weight"] = out["bioml_metadata"]["graph_weights_effective"].get("expression", np.nan)
            row["bioml_effective_apa_weight"] = out["bioml_metadata"]["graph_weights_effective"].get("apa", np.nan)
            row["bioml_reconstruction_error"] = out["bioml_metadata"].get("factorization_reconstruction_error", np.nan)
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    corr_df = pd.DataFrame(corr_rows)
    summary.to_csv(dataset_out / "external_validation_summary.csv", index=False)
    corr_df.to_csv(dataset_out / "within_layer_correlation.csv", index=False)

    site_checklist = build_site_level_checklist(data_dir)
    (dataset_out / "site_level_validation_checklist.json").write_text(json.dumps(site_checklist, indent=2))
    plot_domain_maps(coords, layer_labels, domain_maps, dataset_out)

    return summary, corr_df, {str(dataset_name): method_gene_sets}, site_checklist


def main() -> None:
    args = parse_args()
    data_dirs = args.dataset_dir or ["spaGAPA/data/processed/stapaminer_mob"]
    data_dirs = [resolve_project_path(path) for path in data_dirs]
    output_dir = resolve_project_path(args.output_dir)
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    all_summary = []
    all_corr = []
    all_gene_sets: dict[str, dict[str, set[str]]] = {}
    site_checklists = []

    for data_dir in data_dirs:
        print(f"\nRunning external validation on {data_dir}")
        summary, corr_df, gene_sets, checklist = run_one_dataset(data_dir, output_dir, args)
        all_summary.append(summary)
        all_corr.append(corr_df)
        all_gene_sets.update(gene_sets)
        site_checklists.append(checklist)

    summary_df = pd.concat(all_summary, ignore_index=True) if all_summary else pd.DataFrame()
    corr_df = pd.concat(all_corr, ignore_index=True) if all_corr else pd.DataFrame()
    summary_df.to_csv(output_dir / "external_validation_summary.csv", index=False)
    corr_df.to_csv(output_dir / "within_layer_correlation.csv", index=False)

    overlap_df = summarize_replicate_overlap(all_gene_sets)
    overlap_df.to_csv(output_dir / "replicate_spatial_gene_overlap.csv", index=False)

    (output_dir / "site_level_validation_checklists.json").write_text(json.dumps(site_checklists, indent=2))
    plot_layer_metrics(summary_df, figure_dir)
    plot_internal_metrics(summary_df, figure_dir)
    plot_downstream_counts(summary_df, figure_dir)
    plot_within_layer_correlation(corr_df, figure_dir)
    plot_experiment_coverage(site_checklists, figure_dir)

    decision = {
        "output_dir": str(output_dir),
        "datasets": [str(path) for path in data_dirs],
        "summary_table": str(output_dir / "external_validation_summary.csv"),
        "figures": sorted(str(path) for path in figure_dir.glob("*.png")),
        "reference_package_alignment": {
            "stAPAminer": [
                "layer/domain ARI, NMI, purity, pairwise Jaccard",
                "internal clustering metrics: DBI, CH, silhouette, Dunn",
                "within-layer APA Pearson correlation",
                "SVAPA, LSAPA/domain-marker, and DEAPA/layer-pairwise gene tables",
            ],
            "metaAPA": [
                "site-level readiness checklist for multi-caller consensus and sequence validation",
                "not used as a downstream imputation baseline because metaAPA is upstream site integration",
            ],
        },
    }
    if not summary_df.empty:
        best_layer = summary_df.sort_values("layer_ari", ascending=False).head(1)
        if not best_layer.empty:
            decision["best_layer_ari"] = best_layer.iloc[0].to_dict()
    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))

    print("\nExternal BioML validation summary:")
    print(json.dumps(decision, indent=2, default=str))


if __name__ == "__main__":
    main()

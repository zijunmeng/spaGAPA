#!/usr/bin/env python
"""Pipeline-level smoke benchmark for high-resolution spaGAPA presets.

This runner checks that the package pipeline is no longer weaker than the
standalone high-resolution benchmark path. It uses the same pseudo-bin
simulation and compares:

* ``runner_highres_bioml`` from ``run_high_resolution_simulation.py``
* ``pipeline_highres_accuracy`` from ``SpaGAPA(analysis_preset='highres_accuracy')``
* ``pipeline_highres_fast`` from ``SpaGAPA(analysis_preset='highres_fast')``

The goal is not a full formal benchmark; it is a quick reproducibility guard
for the high-resolution mainline integration.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from spagapa import APADataset, SpaGAPA

from run_external_bioml_validation import (
    evaluate_layer_metrics,
    load_dataset,
    resolve_project_path,
    select_genes,
    track_runtime_memory,
)
from run_high_resolution_simulation import (
    compute_numeric_metrics,
    parse_int_list,
    run_methods,
    sample_parent_spots,
    simulate_pseudo_bins,
)


METHOD_COLORS = {
    "runner_highres_bioml": "#1f9d8a",
    "pipeline_highres_accuracy": "#2c7fb8",
    "pipeline_highres_fast": "#f39c12",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="spaGAPA/data/processed/stapaminer_mob")
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/pipeline_highres_smoke_v1",
    )
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--n-genes", type=int, default=24)
    parser.add_argument("--min-observed-spots", type=int, default=80)
    parser.add_argument("--max-parent-spots", type=int, default=80)
    parser.add_argument("--subbins-per-spot", default="4")
    parser.add_argument("--capture-rate", type=float, default=0.45)
    parser.add_argument("--dropout-rate", type=float, default=0.25)
    parser.add_argument("--measurement-noise", type=float, default=0.08)
    parser.add_argument("--micro-noise", type=float, default=0.04)
    parser.add_argument("--jitter-fraction", type=float, default=0.18)
    parser.add_argument("--expression-noise", type=float, default=0.03)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument(
        "--methods",
        default="highres_bioml",
        help="Runner methods to execute; smoke requires highres_bioml.",
    )
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--expr-n-components", type=int, default=10)
    parser.add_argument("--n-domains", type=int, default=None)
    parser.add_argument("--sparse-n-inducing", type=int, default=60)
    parser.add_argument("--sparse-length-scale", default="auto")
    parser.add_argument("--sparse-length-scale-multiplier", type=float, default=1.0)
    parser.add_argument("--sparse-noise-level", type=float, default=0.08)
    parser.add_argument("--bioml-rank", type=int, default=8)
    parser.add_argument("--bioml-lambda-graph", type=float, default=0.5)
    parser.add_argument("--bioml-lambda-l2", type=float, default=1e-2)
    parser.add_argument("--bioml-max-iter", type=int, default=20)
    parser.add_argument("--bioml-n-neighbors", type=int, default=15)
    parser.add_argument("--bioml-blend", type=float, default=0.1)
    parser.add_argument("--bioml-domain-method", default="spectral", choices=["spectral", "kmeans"])
    parser.add_argument("--bioml-spatial-weight", type=float, default=0.4)
    parser.add_argument("--bioml-expression-weight", type=float, default=0.4)
    parser.add_argument("--bioml-apa-weight", type=float, default=0.2)
    parser.add_argument("--highres-bioml-gp-blend", type=float, default=0.3)
    parser.add_argument("--highres-bioml-spatial-weight", type=float, default=0.2)
    parser.add_argument("--highres-bioml-expression-weight", type=float, default=0.6)
    parser.add_argument("--highres-bioml-apa-weight", type=float, default=0.2)
    parser.add_argument("--highres-bioml-neighbor-mode", default="adaptive", choices=["fixed", "adaptive"])
    parser.add_argument("--highres-bioml-adaptive-neighbor-scale", type=float, default=10.0)
    parser.add_argument("--highres-bioml-parent-weight", type=float, default=0.0)
    parser.add_argument("--highres-bioml-parent-neighbors", type=int, default=8)
    parser.add_argument(
        "--highres-bioml-apa-source",
        default="expression_knn",
        choices=["expression_knn", "raw", "sparse_gp", "none"],
    )
    return parser.parse_args()


def make_dataset(sim: Dict[str, Any], gene_names: list[str]) -> APADataset:
    spot_names = [str(x) for x in sim["subbin_ids"]]
    dataset = APADataset.from_counts(
        sim["observed"],
        sim["coords"],
        gene_names=gene_names,
        spot_names=spot_names,
    )
    dataset.adata.obs["parent_index"] = sim["parent_index"].astype(int)
    dataset.adata.obs["parent_spot"] = pd.Categorical(sim["parent_spots"].astype(str))
    if sim.get("layer_labels") is not None:
        dataset.adata.obs["true_layer"] = pd.Categorical(sim["layer_labels"].astype(str))
    return dataset


def run_pipeline_method(
    sim: Dict[str, Any],
    gene_names: list[str],
    preset: str,
    n_domains: int,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    dataset = make_dataset(sim, gene_names)
    pipeline = SpaGAPA(
        analysis_preset=preset,
        n_neighbors=min(6, max(1, sim["coords"].shape[0] - 1)),
        n_inducing=args.sparse_n_inducing,
        sparse_gp_inducing_method="kmeans",
        sparse_gp_length_scale=args.sparse_length_scale,
        sparse_gp_length_scale_multiplier=args.sparse_length_scale_multiplier,
        sparse_gp_noise_level=args.sparse_noise_level,
        bioml_rank=args.bioml_rank,
        bioml_lambda_graph=args.bioml_lambda_graph,
        bioml_lambda_l2=args.bioml_lambda_l2,
        bioml_max_iter=args.bioml_max_iter,
        bioml_n_neighbors=args.bioml_n_neighbors,
        highres_bioml_gp_blend=args.highres_bioml_gp_blend,
        highres_bioml_apa_source=args.highres_bioml_apa_source,
        highres_bioml_expression_knn_k=args.knn_k,
        highres_bioml_neighbor_mode=args.highres_bioml_neighbor_mode,
        highres_bioml_adaptive_neighbor_scale=args.highres_bioml_adaptive_neighbor_scale,
        highres_bioml_parent_weight=args.highres_bioml_parent_weight,
        highres_bioml_parent_neighbors=args.highres_bioml_parent_neighbors,
        verbose=False,
    )
    started = time.perf_counter()
    results = pipeline.run(
        dataset=dataset,
        expression_embedding=sim["expression_embedding"],
        impute=True,
        quantify=True,
        identify_domains=True,
        differential_analysis=False,
        detect_svapa=False,
        n_domains=n_domains,
    )
    runtime_s = time.perf_counter() - started
    domains = results["domains"]
    return {
        "matrix": domains["imputed_values"],
        "uncertainty": results.get("uncertainty"),
        "domains": domains["labels"],
        "runtime_s": runtime_s,
        "peak_rss_mb": np.nan,
        "metadata": domains.get("metadata", {}),
        "analysis_preset": results.get("analysis_preset", {}),
    }


def evaluate_method(
    method: str,
    output: Dict[str, Any],
    sim: Dict[str, Any],
    scenario_name: str,
    dataset_name: str,
    n_domains: int,
) -> Dict[str, Any]:
    layer_metrics, pred_domains, _ = evaluate_layer_metrics(
        output["matrix"],
        sim["layer_codes"],
        n_domains,
        domain_labels=output.get("domains"),
    )
    numeric = compute_numeric_metrics(
        sim["truth"],
        output["matrix"],
        sim["observed_mask"],
        sim["parent_index"],
        sim["n_parent"],
        output.get("uncertainty"),
    )
    return {
        "dataset": dataset_name,
        "scenario": scenario_name,
        "method": method,
        "n_genes": int(sim["truth"].shape[0]),
        "n_parent_spots": int(sim["n_parent"]),
        "n_bins": int(sim["n_bins"]),
        "observed_fraction": float(sim["observed_mask"].mean()),
        "runtime_s": float(output["runtime_s"]),
        "peak_rss_mb": float(output["peak_rss_mb"]) if np.isfinite(output["peak_rss_mb"]) else np.nan,
        **numeric,
        **layer_metrics,
    }, pred_domains


def plot_smoke_summary(results: pd.DataFrame, figure_dir: Path) -> None:
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "layer_nmi", "runtime_s"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 4.0))
    methods = [m for m in METHOD_COLORS if m in results["method"].values]
    x = np.arange(len(methods))
    for ax, metric in zip(axes, metrics):
        values = [
            float(results.loc[results["method"] == method, metric].mean())
            for method in methods
        ]
        ax.bar(x, values, color=[METHOD_COLORS[method] for method in methods])
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(figure_dir / "pipeline_highres_smoke_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_domain_maps(
    sim: Dict[str, Any],
    domain_maps: Dict[str, np.ndarray],
    figure_dir: Path,
) -> None:
    panels = []
    if sim["layer_labels"] is not None:
        panels.append(("true_layer", sim["layer_labels"]))
    panels.extend(domain_maps.items())
    if not panels:
        return
    n_cols = min(4, len(panels))
    n_rows = int(np.ceil(len(panels) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.0 * n_cols, 4.0 * n_rows))
    axes = np.asarray(axes).reshape(-1)
    xy = sim["coords"]
    for ax, (name, labels) in zip(axes, panels):
        labels = np.asarray(labels).astype(str)
        code_map = {label: idx for idx, label in enumerate(sorted(np.unique(labels)))}
        codes = np.asarray([code_map[label] for label in labels])
        ax.scatter(xy[:, 0], xy[:, 1], c=codes, cmap="tab10", s=12, edgecolor="none")
        ax.set_title(name)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(panels) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figure_dir / "pipeline_highres_smoke_domains.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_decision(results: pd.DataFrame, output_dir: Path) -> Dict[str, Any]:
    decision: Dict[str, Any] = {
        "output_dir": str(output_dir),
        "summary_table": str(output_dir / "pipeline_highres_smoke_summary.csv"),
        "figures": sorted(str(path) for path in (output_dir / "figures").glob("*.png")),
        "comparisons": {},
    }
    by_method = results.groupby("method").mean(numeric_only=True)
    runner = by_method.loc["runner_highres_bioml"] if "runner_highres_bioml" in by_method.index else None
    for method in ["pipeline_highres_accuracy", "pipeline_highres_fast"]:
        if runner is None or method not in by_method.index:
            continue
        row = by_method.loc[method]
        decision["comparisons"][method] = {
            "rmse_delta_vs_runner": float(row["rmse_holdout"] - runner["rmse_holdout"]),
            "parent_rmse_delta_vs_runner": float(row["parent_rmse"] - runner["parent_rmse"]),
            "layer_ari_delta_vs_runner": float(row["layer_ari"] - runner["layer_ari"]),
            "layer_nmi_delta_vs_runner": float(row["layer_nmi"] - runner["layer_nmi"]),
            "runtime_ratio_vs_runner": float(row["runtime_s"] / max(runner["runtime_s"], 1e-9)),
        }
    return decision


def main() -> None:
    args = parse_args()
    if "highres_bioml" not in {m.strip() for m in args.methods.split(",")}:
        raise ValueError("--methods must include highres_bioml for pipeline smoke comparison")
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
        dataset_name = str(metadata["dataset"].dropna().astype(str).iloc[0]) if "dataset" in metadata.columns else data_dir.name
    if args.n_domains is not None:
        n_domains = args.n_domains
    elif args.layer_column in metadata.columns:
        n_domains = int(metadata[args.layer_column].astype(str).nunique())
    else:
        n_domains = min(5, max(2, len(parent_spots) // 20))

    rows = []
    domain_maps: Dict[str, np.ndarray] = {}
    metadata_out: Dict[str, Any] = {}
    last_sim = None
    for subbins in parse_int_list(args.subbins_per_spot):
        scenario_name = f"subbins{subbins}_capture{args.capture_rate:.2f}_dropout{args.dropout_rate:.2f}"
        print(f"\nRunning pipeline highres smoke scenario: {scenario_name}")
        sim = simulate_pseudo_bins(apa, expression, coords, metadata, subbins, args)
        gene_names = apa.index.astype(str).tolist()

        runner_outputs, _, _ = track_runtime_memory(
            lambda: run_methods(sim, args, n_domains)
        )
        runner = runner_outputs["highres_bioml"]
        method_outputs = {
            "runner_highres_bioml": runner,
            "pipeline_highres_accuracy": run_pipeline_method(
                sim, gene_names, "highres_accuracy", n_domains, args
            ),
            "pipeline_highres_fast": run_pipeline_method(
                sim, gene_names, "highres_fast", n_domains, args
            ),
        }
        for method, output in method_outputs.items():
            row, domains = evaluate_method(method, output, sim, scenario_name, dataset_name, n_domains)
            rows.append(row)
            domain_maps[method] = domains
            metadata_out[f"{scenario_name}:{method}"] = {
                "metadata": output.get("metadata", {}),
                "analysis_preset": output.get("analysis_preset", {}),
            }
        last_sim = sim

    results = pd.DataFrame(rows)
    results.to_csv(output_dir / "pipeline_highres_smoke_summary.csv", index=False)
    (output_dir / "pipeline_highres_smoke_metadata.json").write_text(
        json.dumps(metadata_out, indent=2, default=str)
    )
    plot_smoke_summary(results, figure_dir)
    if last_sim is not None:
        plot_domain_maps(last_sim, domain_maps, figure_dir)

    decision = build_decision(results, output_dir)
    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
    print(json.dumps(decision, indent=2, default=str))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Run formal BioML benchmark expansion on the prepared MOB dataset.

This suite wraps ``run_stapaminer_mob_benchmark.py`` to make BioML formal
validation reproducible:

1. A larger balanced BioML formal benchmark.
2. A graph-weight sweep over spatial / expression / APA views.
3. Summary tables, decision JSON, and Matplotlib figures.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Union

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_SCRIPT = PACKAGE_ROOT / "scripts" / "run_stapaminer_mob_benchmark.py"


WEIGHT_SWEEPS = [
    {
        "name": "balanced_s040_e040_a020",
        "spatial_weight": 0.4,
        "expression_weight": 0.4,
        "apa_weight": 0.2,
    },
    {
        "name": "spatial_heavy_s060_e030_a010",
        "spatial_weight": 0.6,
        "expression_weight": 0.3,
        "apa_weight": 0.1,
    },
    {
        "name": "expression_heavy_s020_e060_a020",
        "spatial_weight": 0.2,
        "expression_weight": 0.6,
        "apa_weight": 0.2,
    },
    {
        "name": "apa_heavy_s030_e030_a040",
        "spatial_weight": 0.3,
        "expression_weight": 0.3,
        "apa_weight": 0.4,
    },
    {
        "name": "no_apa_s050_e050_a000",
        "spatial_weight": 0.5,
        "expression_weight": 0.5,
        "apa_weight": 0.0,
    },
]


def resolve_project_path(path: Union[str, Path]) -> Path:
    """Resolve relative paths against the project root used by benchmark scripts."""
    path = Path(path)
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/bioml_formal_suite_v1",
        help="Output directory.",
    )
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared MOB dataset directory.",
    )
    parser.add_argument("--formal-n-genes", type=int, default=60)
    parser.add_argument("--formal-min-observed-spots", type=int, default=110)
    parser.add_argument("--formal-mask-fraction", type=float, default=0.15)
    parser.add_argument("--formal-seeds", default="42,43")
    parser.add_argument(
        "--formal-mask-types",
        default="random,spatial_block_large,ring_sector,layer_aware,low_coverage",
    )
    parser.add_argument("--sweep-n-genes", type=int, default=40)
    parser.add_argument("--sweep-min-observed-spots", type=int, default=110)
    parser.add_argument("--sweep-mask-fraction", type=float, default=0.15)
    parser.add_argument("--sweep-seeds", default="42")
    parser.add_argument(
        "--sweep-mask-types",
        default="random,spatial_block_large,ring_sector,layer_aware,low_coverage",
    )
    parser.add_argument("--gp-kernel", default="matern")
    parser.add_argument("--gp-alpha", type=float, default=1e-3)
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--n-domains", type=int, default=5)
    parser.add_argument("--expr-n-components", type=int, default=10)
    parser.add_argument("--calibration-bins", type=int, default=5)
    parser.add_argument("--bioml-rank", type=int, default=8)
    parser.add_argument("--bioml-lambda-graph", type=float, default=0.5)
    parser.add_argument("--bioml-lambda-l2", type=float, default=1e-2)
    parser.add_argument("--bioml-max-iter", type=int, default=20)
    parser.add_argument("--bioml-n-neighbors", type=int, default=15)
    parser.add_argument("--bioml-blend", type=float, default=0.1)
    parser.add_argument(
        "--bioml-domain-method",
        default="spectral",
        choices=["kmeans", "spectral"],
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run benchmark groups even if summary files already exist.",
    )
    parser.add_argument("--skip-formal", action="store_true")
    parser.add_argument("--skip-sweep", action="store_true")
    return parser.parse_args()


def build_command(
    output_dir: Path,
    *,
    data_dir: str,
    n_genes: int,
    min_observed_spots: int,
    mask_fraction: float,
    seeds: str,
    mask_types: str,
    gp_kernel: str,
    gp_alpha: float,
    knn_k: int,
    n_domains: int,
    expr_n_components: int,
    calibration_bins: int,
    bioml_rank: int,
    bioml_lambda_graph: float,
    bioml_lambda_l2: float,
    bioml_max_iter: int,
    bioml_n_neighbors: int,
    bioml_blend: float,
    bioml_domain_method: str,
    spatial_weight: float,
    expression_weight: float,
    apa_weight: float,
) -> list[str]:
    return [
        sys.executable,
        str(BENCHMARK_SCRIPT),
        "--data-dir",
        str(resolve_project_path(data_dir)),
        "--output-dir",
        str(output_dir.resolve()),
        "--n-genes",
        str(n_genes),
        "--min-observed-spots",
        str(min_observed_spots),
        "--mask-fraction",
        str(mask_fraction),
        "--seeds",
        str(seeds),
        "--mask-types",
        str(mask_types),
        "--knn-k",
        str(knn_k),
        "--gp-kernel",
        str(gp_kernel),
        "--gp-alpha",
        str(gp_alpha),
        "--gp-variant",
        "spatial",
        "--expr-n-components",
        str(expr_n_components),
        "--calibration-bins",
        str(calibration_bins),
        "--n-domains",
        str(n_domains),
        "--include-bioml",
        "--bioml-rank",
        str(bioml_rank),
        "--bioml-lambda-graph",
        str(bioml_lambda_graph),
        "--bioml-lambda-l2",
        str(bioml_lambda_l2),
        "--bioml-max-iter",
        str(bioml_max_iter),
        "--bioml-n-neighbors",
        str(bioml_n_neighbors),
        "--bioml-blend",
        str(bioml_blend),
        "--bioml-domain-method",
        str(bioml_domain_method),
        "--bioml-spatial-weight",
        str(spatial_weight),
        "--bioml-expression-weight",
        str(expression_weight),
        "--bioml-apa-weight",
        str(apa_weight),
    ]


def run_one(cmd: list[str]) -> None:
    subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)


def run_exists(run_dir: Path) -> bool:
    return (run_dir / "benchmark_results_overall.csv").exists() and (
        run_dir / "benchmark_results_long.csv"
    ).exists()


def load_run_summary(run_dir: Path, run_name: str, extra: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    overall = pd.read_csv(run_dir / "benchmark_results_overall.csv", index_col=0)
    long_df = pd.read_csv(run_dir / "benchmark_results_long.csv")

    overall_rows: list[dict[str, Any]] = []
    for method, row in overall.iterrows():
        record = row.to_dict()
        record.update(extra)
        record["run_name"] = run_name
        record["method"] = method
        overall_rows.append(record)

    mask_summary = (
        long_df.groupby(["mask_type", "method"])[
            [
                "rmse",
                "layer_ari",
                "layer_nmi",
                "uncertainty_error_spearman",
                "uncertainty_coverage_68",
                "runtime_s",
            ]
        ]
        .mean()
        .reset_index()
    )
    mask_rows = []
    for _, row in mask_summary.iterrows():
        record = row.to_dict()
        record.update(extra)
        record["run_name"] = run_name
        mask_rows.append(record)
    return overall_rows, mask_rows


def summarize_suite(suite_dir: Path, runs: list[tuple[str, Path, dict[str, Any]]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    overall_rows: list[dict[str, Any]] = []
    mask_rows: list[dict[str, Any]] = []
    for run_name, run_dir, extra in runs:
        if not run_exists(run_dir):
            continue
        run_overall, run_masks = load_run_summary(run_dir, run_name, extra)
        overall_rows.extend(run_overall)
        mask_rows.extend(run_masks)

    overall_df = pd.DataFrame(overall_rows)
    mask_df = pd.DataFrame(mask_rows)
    if not overall_df.empty:
        overall_df.to_csv(suite_dir / "bioml_formal_overall_summary.csv", index=False)
    if not mask_df.empty:
        mask_df.to_csv(suite_dir / "bioml_formal_mask_summary.csv", index=False)
    return overall_df, mask_df


def plot_suite(overall_df: pd.DataFrame, sweep_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not overall_df.empty:
        formal = overall_df[overall_df["run_type"] == "formal"]
        if not formal.empty:
            methods = ["spagapa_gp", "spagapa_bioml", "stapaminer_knn_expression", "knn_spatial", "mean"]
            formal = formal[formal["method"].isin(methods)].copy()
            formal["method"] = pd.Categorical(formal["method"], categories=methods, ordered=True)
            formal = formal.sort_values("method")
            fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
            x = np.arange(len(formal))
            axes[0].bar(x, formal["rmse"], color="#2c7fb8")
            axes[0].set_title("Formal RMSE")
            axes[0].set_ylabel("RMSE")
            axes[1].bar(x, formal["layer_ari"], color="#41ab5d")
            axes[1].set_title("Formal Layer ARI")
            axes[1].set_ylabel("ARI")
            axes[2].bar(x, formal["uncertainty_error_spearman"], color="#f16913")
            axes[2].set_title("Formal Uncertainty vs Error")
            axes[2].set_ylabel("Spearman")
            for ax in axes:
                ax.set_xticks(x)
                ax.set_xticklabels(formal["method"], rotation=25, ha="right")
                ax.grid(axis="y", alpha=0.25)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
            plt.tight_layout()
            fig.savefig(output_dir / "bioml_formal_method_comparison.png", dpi=300, bbox_inches="tight")
            plt.close(fig)

    if not sweep_df.empty:
        bioml = sweep_df[sweep_df["method"] == "spagapa_bioml"].copy()
        if not bioml.empty:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            x = np.arange(len(bioml))
            labels = bioml["run_name"].str.replace("_", "\n")
            axes[0].bar(x, bioml["rmse"], color="#3182bd")
            axes[0].set_title("Weight Sweep: RMSE")
            axes[0].set_ylabel("RMSE")
            axes[1].bar(x, bioml["layer_ari"], color="#31a354")
            axes[1].set_title("Weight Sweep: Layer ARI")
            axes[1].set_ylabel("ARI")
            for ax in axes:
                ax.set_xticks(x)
                ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=8)
                ax.grid(axis="y", alpha=0.25)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
            plt.tight_layout()
            fig.savefig(output_dir / "bioml_graph_weight_sweep.png", dpi=300, bbox_inches="tight")
            plt.close(fig)


def build_decision_summary(overall_df: pd.DataFrame, suite_dir: Path) -> dict[str, Any]:
    decision: dict[str, Any] = {
        "suite_date": pd.Timestamp.now().isoformat(),
        "formal_best_by_rmse": None,
        "formal_best_by_layer_ari": None,
        "sweep_best_bioml_by_rmse": None,
        "sweep_best_bioml_by_layer_ari": None,
        "recommendation": "",
    }
    notes: list[str] = []
    if overall_df.empty:
        decision["recommendation"] = "No completed runs found."
        (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2))
        return decision

    formal = overall_df[overall_df["run_type"] == "formal"]
    if not formal.empty:
        best_rmse = formal.sort_values("rmse", ascending=True).iloc[0]
        best_ari = formal.sort_values("layer_ari", ascending=False).iloc[0]
        decision["formal_best_by_rmse"] = best_rmse.to_dict()
        decision["formal_best_by_layer_ari"] = best_ari.to_dict()
        notes.append(f"Formal best RMSE: {best_rmse['method']}")
        notes.append(f"Formal best layer ARI: {best_ari['method']}")

    sweep_bioml = overall_df[
        (overall_df["run_type"] == "weight_sweep") & (overall_df["method"] == "spagapa_bioml")
    ]
    if not sweep_bioml.empty:
        best_sweep_rmse = sweep_bioml.sort_values("rmse", ascending=True).iloc[0]
        best_sweep_ari = sweep_bioml.sort_values("layer_ari", ascending=False).iloc[0]
        decision["sweep_best_bioml_by_rmse"] = best_sweep_rmse.to_dict()
        decision["sweep_best_bioml_by_layer_ari"] = best_sweep_ari.to_dict()
        notes.append(f"Sweep best BioML RMSE: {best_sweep_rmse['run_name']}")
        notes.append(f"Sweep best BioML layer ARI: {best_sweep_ari['run_name']}")

    decision["recommendation"] = " ".join(notes)
    (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2))
    return decision


def main() -> None:
    args = parse_args()
    suite_dir = Path(args.output_dir)
    if not suite_dir.is_absolute():
        suite_dir = resolve_project_path(suite_dir)
    formal_dir = suite_dir / "formal"
    sweep_dir = suite_dir / "weight_sweeps"
    figure_dir = suite_dir / "figures"
    suite_dir.mkdir(parents=True, exist_ok=True)
    formal_dir.mkdir(exist_ok=True)
    sweep_dir.mkdir(exist_ok=True)

    completed_runs: list[tuple[str, Path, dict[str, Any]]] = []
    balanced = WEIGHT_SWEEPS[0]

    if not args.skip_formal:
        run_name = f"formal_{balanced['name']}"
        run_dir = formal_dir / run_name
        cmd = build_command(
            run_dir,
            data_dir=args.data_dir,
            n_genes=args.formal_n_genes,
            min_observed_spots=args.formal_min_observed_spots,
            mask_fraction=args.formal_mask_fraction,
            seeds=args.formal_seeds,
            mask_types=args.formal_mask_types,
            gp_kernel=args.gp_kernel,
            gp_alpha=args.gp_alpha,
            knn_k=args.knn_k,
            n_domains=args.n_domains,
            expr_n_components=args.expr_n_components,
            calibration_bins=args.calibration_bins,
            bioml_rank=args.bioml_rank,
            bioml_lambda_graph=args.bioml_lambda_graph,
            bioml_lambda_l2=args.bioml_lambda_l2,
            bioml_max_iter=args.bioml_max_iter,
            bioml_n_neighbors=args.bioml_n_neighbors,
            bioml_blend=args.bioml_blend,
            bioml_domain_method=args.bioml_domain_method,
            spatial_weight=balanced["spatial_weight"],
            expression_weight=balanced["expression_weight"],
            apa_weight=balanced["apa_weight"],
        )
        print(f"\n[formal] Running {run_name}")
        if run_exists(run_dir) and not args.force:
            print(f"[formal] Skipping completed run {run_name}. Use --force to re-run.")
        else:
            run_one(cmd)
        completed_runs.append((run_name, run_dir, {"run_type": "formal", **balanced}))

    if not args.skip_sweep:
        for spec in WEIGHT_SWEEPS:
            run_name = spec["name"]
            run_dir = sweep_dir / run_name
            cmd = build_command(
                run_dir,
                data_dir=args.data_dir,
                n_genes=args.sweep_n_genes,
                min_observed_spots=args.sweep_min_observed_spots,
                mask_fraction=args.sweep_mask_fraction,
                seeds=args.sweep_seeds,
                mask_types=args.sweep_mask_types,
                gp_kernel=args.gp_kernel,
                gp_alpha=args.gp_alpha,
                knn_k=args.knn_k,
                n_domains=args.n_domains,
                expr_n_components=args.expr_n_components,
                calibration_bins=args.calibration_bins,
                bioml_rank=args.bioml_rank,
                bioml_lambda_graph=args.bioml_lambda_graph,
                bioml_lambda_l2=args.bioml_lambda_l2,
                bioml_max_iter=args.bioml_max_iter,
                bioml_n_neighbors=args.bioml_n_neighbors,
                bioml_blend=args.bioml_blend,
                bioml_domain_method=args.bioml_domain_method,
                spatial_weight=spec["spatial_weight"],
                expression_weight=spec["expression_weight"],
                apa_weight=spec["apa_weight"],
            )
            print(f"\n[sweep] Running {run_name}")
            if run_exists(run_dir) and not args.force:
                print(f"[sweep] Skipping completed run {run_name}. Use --force to re-run.")
            else:
                run_one(cmd)
            completed_runs.append((run_name, run_dir, {"run_type": "weight_sweep", **spec}))

    if args.skip_formal:
        formal_runs = [
            (f"formal_{balanced['name']}", formal_dir / f"formal_{balanced['name']}", {"run_type": "formal", **balanced})
        ]
        completed_runs.extend(formal_runs)
    if args.skip_sweep:
        completed_runs.extend(
            [(spec["name"], sweep_dir / spec["name"], {"run_type": "weight_sweep", **spec}) for spec in WEIGHT_SWEEPS]
        )

    overall_df, mask_df = summarize_suite(suite_dir, completed_runs)
    sweep_df = overall_df[overall_df["run_type"] == "weight_sweep"] if not overall_df.empty else pd.DataFrame()
    plot_suite(overall_df, sweep_df, figure_dir)
    decision = build_decision_summary(overall_df, suite_dir)

    print("\nBioML formal suite decision summary:")
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()

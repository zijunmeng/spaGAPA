#!/usr/bin/env python
"""Run a unified benchmark suite for spaGAPA GP variants on the MOB dataset."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
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


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_SCRIPT = PACKAGE_ROOT / "scripts" / "run_stapaminer_mob_benchmark.py"


FORMAL_VARIANTS = [
    {
        "name": "spatial",
        "dir_name": "formal_v2_spatial",
        "gp_variant": "spatial",
        "gp_use_theta": False,
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "spatial_radial",
        "dir_name": "formal_v2_spatial_radial",
        "gp_variant": "spatial_radial",
        "gp_use_theta": True,
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive",
        "dir_name": "formal_v2_expr_additive",
        "gp_variant": "expr_additive",
        "gp_use_theta": True,
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
]


SWEEP_SPECS = [
    {
        "name": "expr_additive_lambda_0.2_pca10",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_0.5_pca10",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_1.0_pca10",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 1.0,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_2.0_pca10",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 2.0,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_0.5_pca5",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 5,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_0.5_pca15",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 15,
        "gp_use_theta": True,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
    },
    {
        "name": "expr_additive_lambda_0.5_pca10_hvg",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "gp_use_theta": True,
        "expr_use_hvg": True,
        "expr_n_top_genes": 1000,
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/expression_gp_suite",
        help="Output directory for the suite.",
    )
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared MOB dataset directory.",
    )
    parser.add_argument("--n-genes", type=int, default=100)
    parser.add_argument("--min-observed-spots", type=int, default=100)
    parser.add_argument("--mask-fraction", type=float, default=0.2)
    parser.add_argument("--seeds", default="42,43")
    parser.add_argument(
        "--mask-types",
        default="random,spatial_block_large,ring_sector,layer_aware,low_coverage",
    )
    parser.add_argument("--gp-kernel", default="matern")
    parser.add_argument("--gp-alpha", type=float, default=1e-3)
    parser.add_argument("--calibration-bins", type=int, default=8)
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--n-domains", type=int, default=5)
    parser.add_argument(
        "--skip-formal",
        action="store_true",
        help="Skip formal three-variant benchmark runs.",
    )
    parser.add_argument(
        "--skip-sweep",
        action="store_true",
        help="Skip expr_additive parameter sweeps.",
    )
    return parser.parse_args()


def build_command(
    output_dir: Path,
    base_args: argparse.Namespace,
    run_spec: dict[str, Any],
) -> list[str]:
    cmd = [
        sys.executable,
        str(BENCHMARK_SCRIPT),
        "--data-dir",
        str(Path(base_args.data_dir).resolve()),
        "--output-dir",
        str(output_dir.resolve()),
        "--n-genes",
        str(base_args.n_genes),
        "--min-observed-spots",
        str(base_args.min_observed_spots),
        "--mask-fraction",
        str(base_args.mask_fraction),
        "--seeds",
        str(base_args.seeds),
        "--mask-types",
        str(base_args.mask_types),
        "--knn-k",
        str(base_args.knn_k),
        "--gp-kernel",
        str(base_args.gp_kernel),
        "--gp-alpha",
        str(base_args.gp_alpha),
        "--calibration-bins",
        str(base_args.calibration_bins),
        "--n-domains",
        str(base_args.n_domains),
        "--gp-variant",
        str(run_spec["gp_variant"]),
        "--gp-lambda-expr",
        str(run_spec["gp_lambda_expr"]),
        "--gp-product-offset",
        str(run_spec.get("gp_product_offset", 1.0)),
        "--expr-n-components",
        str(run_spec["expr_n_components"]),
        "--expr-n-top-genes",
        str(run_spec["expr_n_top_genes"]),
    ]
    if run_spec.get("gp_use_theta", False):
        cmd.append("--gp-use-theta")
    if run_spec.get("expr_use_hvg", False):
        cmd.append("--expr-use-hvg")
    return cmd


def run_one(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def load_run_summary(run_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    overall = pd.read_csv(run_dir / "benchmark_results_overall.csv", index_col=0)
    summary = json.loads((run_dir / "benchmark_summary.json").read_text())
    return overall, summary


def run_exists(run_dir: Path) -> bool:
    return (run_dir / "benchmark_results_overall.csv").exists() and (run_dir / "benchmark_summary.json").exists()


def summarize_formal_runs(suite_dir: Path, formal_runs: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for run in formal_runs:
        run_dir = suite_dir / "formal" / run["dir_name"]
        if not run_exists(run_dir):
            continue
        overall, summary = load_run_summary(run_dir)
        gp_row = overall.loc["spagapa_gp"].to_dict()
        gp_row["run_name"] = run["name"]
        gp_row["gp_variant"] = summary.get("gp_variant")
        gp_row["gp_lambda_expr"] = summary.get("gp_lambda_expr")
        gp_row["gp_product_offset"] = summary.get("gp_product_offset")
        gp_row["expr_n_components"] = summary.get("expr_n_components")
        gp_row["expr_use_hvg"] = summary.get("expr_use_hvg")
        gp_row["gp_use_theta"] = summary.get("gp_use_theta")
        rows.append(gp_row)
    formal_df = pd.DataFrame(rows)
    if formal_df.empty:
        return formal_df
    formal_df = formal_df[
        [
            "run_name",
            "gp_variant",
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
            "gp_lambda_expr",
            "gp_product_offset",
            "expr_n_components",
            "expr_use_hvg",
            "gp_use_theta",
        ]
    ]
    formal_df.to_csv(suite_dir / "formal_variant_comparison.csv", index=False)
    return formal_df


def summarize_sweep_runs(suite_dir: Path, sweep_specs: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for run in sweep_specs:
        run_dir = suite_dir / "sweeps" / run["name"]
        if not run_exists(run_dir):
            continue
        overall, summary = load_run_summary(run_dir)
        gp_row = overall.loc["spagapa_gp"].to_dict()
        gp_row["run_name"] = run["name"]
        gp_row["gp_variant"] = summary.get("gp_variant")
        gp_row["gp_lambda_expr"] = summary.get("gp_lambda_expr")
        gp_row["gp_product_offset"] = summary.get("gp_product_offset")
        gp_row["expr_n_components"] = summary.get("expr_n_components")
        gp_row["expr_use_hvg"] = summary.get("expr_use_hvg")
        gp_row["gp_use_theta"] = summary.get("gp_use_theta")
        rows.append(gp_row)
    sweep_df = pd.DataFrame(rows)
    if sweep_df.empty:
        return sweep_df
    sweep_df = sweep_df.sort_values(["rmse", "layer_ari"], ascending=[True, False])
    sweep_df.to_csv(suite_dir / "expr_additive_sweep_summary.csv", index=False)
    return sweep_df


def plot_suite(formal_df: pd.DataFrame, sweep_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    if not formal_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        x = np.arange(len(formal_df))
        axes[0].bar(x, formal_df["rmse"], color=["#c0392b", "#2980b9", "#8e44ad"][: len(formal_df)])
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(formal_df["run_name"], rotation=20, ha="right")
        axes[0].set_ylabel("RMSE")
        axes[0].set_title("Formal v2 GP Variant RMSE")

        axes[1].bar(x, formal_df["layer_ari"], color=["#16a085", "#27ae60", "#8e44ad"][: len(formal_df)])
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(formal_df["run_name"], rotation=20, ha="right")
        axes[1].set_ylabel("Layer ARI")
        axes[1].set_title("Formal v2 Biological Consistency")

        axes[2].bar(x, formal_df["uncertainty_error_spearman"], color=["#2c3e50", "#34495e", "#7f8c8d"][: len(formal_df)])
        axes[2].set_xticks(x)
        axes[2].set_xticklabels(formal_df["run_name"], rotation=20, ha="right")
        axes[2].set_ylabel("Spearman")
        axes[2].set_title("Uncertainty vs Error")

        for ax in axes:
            ax.grid(axis="y", alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(output_dir / "formal_variant_comparison.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    if not sweep_df.empty:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        scatter = axes[0].scatter(
            sweep_df["gp_lambda_expr"],
            sweep_df["rmse"],
            c=sweep_df["expr_n_components"],
            cmap="viridis",
            s=85,
            edgecolor="white",
            linewidth=0.6,
        )
        axes[0].set_xlabel("lambda_expr")
        axes[0].set_ylabel("RMSE")
        axes[0].set_title("Expr-additive Sweep: RMSE")
        cbar = fig.colorbar(scatter, ax=axes[0])
        cbar.set_label("Expression PCA dims")

        axes[1].scatter(
            sweep_df["gp_lambda_expr"],
            sweep_df["layer_ari"],
            c=sweep_df["expr_n_components"],
            cmap="plasma",
            s=85,
            edgecolor="white",
            linewidth=0.6,
        )
        axes[1].set_xlabel("lambda_expr")
        axes[1].set_ylabel("Layer ARI")
        axes[1].set_title("Expr-additive Sweep: Layer ARI")

        for ax in axes:
            ax.grid(alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(output_dir / "expr_additive_sweep.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def build_decision_summary(
    suite_dir: Path,
    formal_df: pd.DataFrame,
    sweep_df: pd.DataFrame,
    args: argparse.Namespace,
) -> dict[str, Any]:
    decision: dict[str, Any] = {
        "suite_date": pd.Timestamp.now().isoformat(),
        "dataset": "stAPAminer_MOB",
        "n_genes": int(args.n_genes),
        "seeds": [int(x.strip()) for x in str(args.seeds).split(",") if x.strip()],
        "mask_types": [x.strip() for x in str(args.mask_types).split(",") if x.strip()],
        "gp_kernel": str(args.gp_kernel),
        "gp_alpha": float(args.gp_alpha),
        "formal_best_variant_by_rmse": None,
        "formal_best_variant_by_layer_ari": None,
        "best_sweep_run_by_rmse": None,
        "best_sweep_run_by_layer_ari": None,
        "recommendation": "",
    }

    if not formal_df.empty:
        rmse_best = formal_df.sort_values("rmse", ascending=True).iloc[0]
        ari_best = formal_df.sort_values("layer_ari", ascending=False).iloc[0]
        decision["formal_best_variant_by_rmse"] = rmse_best.to_dict()
        decision["formal_best_variant_by_layer_ari"] = ari_best.to_dict()

    if not sweep_df.empty:
        decision["best_sweep_run_by_rmse"] = sweep_df.sort_values("rmse", ascending=True).iloc[0].to_dict()
        decision["best_sweep_run_by_layer_ari"] = sweep_df.sort_values("layer_ari", ascending=False).iloc[0].to_dict()

    recommendation = []
    if not formal_df.empty:
        rmse_best_variant = str(decision["formal_best_variant_by_rmse"]["gp_variant"])
        ari_best_variant = str(decision["formal_best_variant_by_layer_ari"]["gp_variant"])
        recommendation.append(f"Best RMSE variant: {rmse_best_variant}")
        recommendation.append(f"Best layer ARI variant: {ari_best_variant}")
        if rmse_best_variant == "spatial" and ari_best_variant != "expr_additive":
            recommendation.append("Current evidence does not justify switching the mainline model to expression-informed GP yet.")
        elif ari_best_variant == "expr_additive":
            recommendation.append("Expression-informed GP shows biological promise, but should be judged against its RMSE/runtime cost.")
    if not sweep_df.empty:
        best_sweep = decision["best_sweep_run_by_rmse"]
        recommendation.append(
            f"Best expr_additive sweep by RMSE: {best_sweep['run_name']} (lambda={best_sweep['gp_lambda_expr']}, pca={best_sweep['expr_n_components']})"
        )
    decision["recommendation"] = " ".join(recommendation)

    (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2))
    return decision


def main() -> None:
    args = parse_args()
    suite_dir = Path(args.output_dir)
    suite_dir.mkdir(parents=True, exist_ok=True)

    formal_dir = suite_dir / "formal"
    sweep_dir = suite_dir / "sweeps"
    formal_dir.mkdir(exist_ok=True)
    sweep_dir.mkdir(exist_ok=True)

    if not args.skip_formal:
        for spec in FORMAL_VARIANTS:
            run_dir = formal_dir / spec["dir_name"]
            cmd = build_command(run_dir, args, spec)
            print(f"\n[formal] Running {spec['name']}")
            run_one(cmd, cwd=PROJECT_ROOT)

    if not args.skip_sweep:
        for spec in SWEEP_SPECS:
            run_dir = sweep_dir / spec["name"]
            cmd = build_command(run_dir, args, spec)
            print(f"\n[sweep] Running {spec['name']}")
            run_one(cmd, cwd=PROJECT_ROOT)

    formal_df = summarize_formal_runs(suite_dir, FORMAL_VARIANTS)
    sweep_df = summarize_sweep_runs(suite_dir, SWEEP_SPECS)
    plot_suite(formal_df, sweep_df, suite_dir / "figures")
    decision = build_decision_summary(suite_dir, formal_df, sweep_df, args)

    print("\nSuite decision summary:")
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()

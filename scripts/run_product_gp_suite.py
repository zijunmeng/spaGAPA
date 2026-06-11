#!/usr/bin/env python
"""Run a focused product-kernel GP benchmark suite on the MOB dataset."""

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


PILOT_VARIANTS = [
    {
        "name": "spatial",
        "dir_name": "pilot_spatial",
        "gp_variant": "spatial",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": False,
    },
    {
        "name": "expr_additive",
        "dir_name": "pilot_expr_additive",
        "gp_variant": "expr_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product",
        "dir_name": "pilot_expr_product",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_additive",
        "dir_name": "pilot_adaptive_additive",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_layer_local",
        "dir_name": "pilot_expr_layer_local",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 20,
        "gp_layer_gate_mode": "radius",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
]


PRODUCT_SWEEPS = [
    {
        "name": "expr_product_lambda_0.2_offset_0.5",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_0.2_offset_1.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_0.2_offset_2.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 2.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_0.5_offset_0.5",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_0.5_offset_1.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_0.5_offset_2.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 2.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_1.0_offset_0.5",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 1.0,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_1.0_offset_1.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 1.0,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "expr_product_lambda_1.0_offset_2.0",
        "gp_variant": "expr_product",
        "gp_lambda_expr": 1.0,
        "gp_product_offset": 2.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
]


ADAPTIVE_SWEEPS = [
    {
        "name": "adaptive_lambda_0.2_tau_auto",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_lambda_0.5_tau_auto",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_lambda_1.0_tau_auto",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 1.0,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_lambda_0.5_tau_0.5",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": 0.5,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_lambda_0.5_tau_1.0",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": 1.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "adaptive_lambda_0.5_tau_2.0",
        "gp_variant": "adaptive_additive",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 1.0,
        "gp_gate_mode": "distance",
        "gp_gate_tau": 2.0,
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
]


LAYER_LOCAL_SWEEPS = [
    {
        "name": "layer_local_lambda_0.2_k_10_radius",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 10,
        "gp_layer_gate_mode": "radius",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "layer_local_lambda_0.5_k_10_radius",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 10,
        "gp_layer_gate_mode": "radius",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "layer_local_lambda_0.2_k_20_radius",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 20,
        "gp_layer_gate_mode": "radius",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "layer_local_lambda_0.5_k_20_radius",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 20,
        "gp_layer_gate_mode": "radius",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "layer_local_lambda_0.2_k_20_pseudolayer",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.2,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 20,
        "gp_layer_gate_mode": "pseudolayer",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
    {
        "name": "layer_local_lambda_0.5_k_20_pseudolayer",
        "gp_variant": "expr_layer_local",
        "gp_lambda_expr": 0.5,
        "gp_product_offset": 0.5,
        "gp_gate_mode": "distance",
        "gp_gate_tau": None,
        "gp_local_k": 20,
        "gp_layer_gate_mode": "pseudolayer",
        "expr_n_components": 10,
        "expr_use_hvg": False,
        "expr_n_top_genes": 2000,
        "gp_use_theta": True,
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/product_gp_suite",
        help="Output directory.",
    )
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared MOB dataset directory.",
    )
    parser.add_argument("--pilot-n-genes", type=int, default=30)
    parser.add_argument("--pilot-min-observed-spots", type=int, default=110)
    parser.add_argument("--pilot-mask-fraction", type=float, default=0.15)
    parser.add_argument("--pilot-seeds", default="42")
    parser.add_argument(
        "--pilot-mask-types",
        default="random,ring_sector,layer_aware",
    )
    parser.add_argument("--sweep-n-genes", type=int, default=30)
    parser.add_argument("--sweep-min-observed-spots", type=int, default=110)
    parser.add_argument("--sweep-mask-fraction", type=float, default=0.15)
    parser.add_argument("--sweep-seeds", default="42")
    parser.add_argument(
        "--sweep-mask-types",
        default="random,ring_sector,layer_aware",
    )
    parser.add_argument("--gp-kernel", default="matern")
    parser.add_argument("--gp-alpha", type=float, default=1e-3)
    parser.add_argument("--calibration-bins", type=int, default=5)
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--n-domains", type=int, default=5)
    parser.add_argument("--skip-pilot", action="store_true")
    parser.add_argument("--skip-sweep", action="store_true")
    parser.add_argument(
        "--sweep-family",
        default="product",
        choices=["product", "adaptive", "layer_local"],
        help="Which sweep family to run.",
    )
    return parser.parse_args()


def build_command(
    output_dir: Path,
    run_spec: dict[str, Any],
    *,
    data_dir: str,
    n_genes: int,
    min_observed_spots: int,
    mask_fraction: float,
    seeds: str,
    mask_types: str,
    gp_kernel: str,
    gp_alpha: float,
    calibration_bins: int,
    knn_k: int,
    n_domains: int,
) -> list[str]:
    cmd = [
        sys.executable,
        str(BENCHMARK_SCRIPT),
        "--data-dir",
        str(Path(data_dir).resolve()),
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
        "--calibration-bins",
        str(calibration_bins),
        "--n-domains",
        str(n_domains),
        "--gp-variant",
        str(run_spec["gp_variant"]),
        "--gp-lambda-expr",
        str(run_spec["gp_lambda_expr"]),
        "--gp-product-offset",
        str(run_spec["gp_product_offset"]),
        "--gp-gate-mode",
        str(run_spec.get("gp_gate_mode", "distance")),
        "--gp-local-k",
        str(run_spec.get("gp_local_k", 20)),
        "--gp-layer-gate-mode",
        str(run_spec.get("gp_layer_gate_mode", "radius")),
        "--expr-n-components",
        str(run_spec["expr_n_components"]),
        "--expr-n-top-genes",
        str(run_spec["expr_n_top_genes"]),
    ]
    if run_spec.get("gp_gate_tau") is not None:
        cmd.extend(["--gp-gate-tau", str(run_spec["gp_gate_tau"])])
    if run_spec.get("gp_use_theta", False):
        cmd.append("--gp-use-theta")
    if run_spec.get("expr_use_hvg", False):
        cmd.append("--expr-use-hvg")
    return cmd


def run_one(cmd: list[str]) -> None:
    subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)


def run_exists(run_dir: Path) -> bool:
    return (run_dir / "benchmark_results_overall.csv").exists() and (run_dir / "benchmark_summary.json").exists()


def load_summary(run_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    overall = pd.read_csv(run_dir / "benchmark_results_overall.csv", index_col=0)
    summary = json.loads((run_dir / "benchmark_summary.json").read_text())
    return overall, summary


def summarize_runs(
    base_dir: Path,
    specs: list[dict[str, Any]],
    summary_name: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for spec in specs:
        run_dir = base_dir / spec["name"]
        if not run_exists(run_dir):
            run_dir = base_dir / spec.get("dir_name", spec["name"])
        if not run_exists(run_dir):
            continue
        overall, summary = load_summary(run_dir)
        gp_row = overall.loc["spagapa_gp"].to_dict()
        gp_row["run_name"] = spec["name"]
        gp_row["gp_variant"] = summary.get("gp_variant")
        gp_row["gp_lambda_expr"] = summary.get("gp_lambda_expr")
        gp_row["gp_product_offset"] = summary.get("gp_product_offset")
        gp_row["gp_gate_mode"] = summary.get("gp_gate_mode")
        gp_row["gp_gate_tau"] = summary.get("gp_gate_tau")
        gp_row["gp_local_k"] = summary.get("gp_local_k")
        gp_row["gp_layer_gate_mode"] = summary.get("gp_layer_gate_mode")
        gp_row["expr_n_components"] = summary.get("expr_n_components")
        gp_row["expr_use_hvg"] = summary.get("expr_use_hvg")
        gp_row["gp_use_theta"] = summary.get("gp_use_theta")
        rows.append(gp_row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df.to_csv(base_dir.parent / summary_name, index=False)
    return df


def plot_results(pilot_df: pd.DataFrame, sweep_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    if not pilot_df.empty:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        x = np.arange(len(pilot_df))
        rmse_colors = ["#c0392b", "#8e44ad", "#2980b9", "#16a085", "#d35400"]
        ari_colors = ["#16a085", "#27ae60", "#9b59b6", "#2c3e50", "#c0392b"]
        unc_colors = ["#2c3e50", "#34495e", "#7f8c8d", "#95a5a6", "#bdc3c7"]
        axes[0].bar(x, pilot_df["rmse"], color=rmse_colors[: len(pilot_df)])
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(pilot_df["run_name"], rotation=20, ha="right")
        axes[0].set_ylabel("RMSE")
        axes[0].set_title("GP Variant Pilot: RMSE")

        axes[1].bar(x, pilot_df["layer_ari"], color=ari_colors[: len(pilot_df)])
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(pilot_df["run_name"], rotation=20, ha="right")
        axes[1].set_ylabel("Layer ARI")
        axes[1].set_title("GP Variant Pilot: Layer ARI")

        axes[2].bar(x, pilot_df["uncertainty_error_spearman"], color=unc_colors[: len(pilot_df)])
        axes[2].set_xticks(x)
        axes[2].set_xticklabels(pilot_df["run_name"], rotation=20, ha="right")
        axes[2].set_ylabel("Spearman")
        axes[2].set_title("GP Variant Pilot: Uncertainty vs Error")

        for ax in axes:
            ax.grid(axis="y", alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(output_dir / "gp_variant_pilot_comparison.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    if not sweep_df.empty:
        is_product = (sweep_df["gp_variant"] == "expr_product").all()
        is_adaptive = (sweep_df["gp_variant"] == "adaptive_additive").all()
        if is_product:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            scatter = axes[0].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["rmse"],
                c=sweep_df["gp_product_offset"],
                cmap="viridis",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[0].set_xlabel("lambda_expr")
            axes[0].set_ylabel("RMSE")
            axes[0].set_title("Product Sweep: RMSE")
            cbar = fig.colorbar(scatter, ax=axes[0])
            cbar.set_label("product_offset")

            scatter2 = axes[1].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["layer_ari"],
                c=sweep_df["gp_product_offset"],
                cmap="plasma",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[1].set_xlabel("lambda_expr")
            axes[1].set_ylabel("Layer ARI")
            axes[1].set_title("Product Sweep: Layer ARI")
            cbar2 = fig.colorbar(scatter2, ax=axes[1])
            cbar2.set_label("product_offset")
            figure_name = "product_sweep.png"
        elif is_adaptive:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            tau_series = sweep_df["gp_gate_tau"].copy()
            if tau_series.isna().any():
                fallback_tau = float(sweep_df["gp_lambda_expr"].median())
                tau_series = tau_series.fillna(fallback_tau)
            scatter = axes[0].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["rmse"],
                c=tau_series,
                cmap="viridis",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[0].set_xlabel("lambda_expr")
            axes[0].set_ylabel("RMSE")
            axes[0].set_title("Adaptive Sweep: RMSE")
            cbar = fig.colorbar(scatter, ax=axes[0])
            cbar.set_label("gate_tau")

            scatter2 = axes[1].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["layer_ari"],
                c=tau_series,
                cmap="plasma",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[1].set_xlabel("lambda_expr")
            axes[1].set_ylabel("Layer ARI")
            axes[1].set_title("Adaptive Sweep: Layer ARI")
            cbar2 = fig.colorbar(scatter2, ax=axes[1])
            cbar2.set_label("gate_tau")
            figure_name = "adaptive_sweep.png"
        else:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
            scatter = axes[0].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["rmse"],
                c=sweep_df["gp_local_k"],
                cmap="viridis",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[0].set_xlabel("lambda_expr")
            axes[0].set_ylabel("RMSE")
            axes[0].set_title("Layer-Local Sweep: RMSE")
            cbar = fig.colorbar(scatter, ax=axes[0])
            cbar.set_label("local_k")

            scatter2 = axes[1].scatter(
                sweep_df["gp_lambda_expr"],
                sweep_df["layer_ari"],
                c=sweep_df["gp_local_k"],
                cmap="plasma",
                s=90,
                edgecolor="white",
                linewidth=0.6,
            )
            axes[1].set_xlabel("lambda_expr")
            axes[1].set_ylabel("Layer ARI")
            axes[1].set_title("Layer-Local Sweep: Layer ARI")
            cbar2 = fig.colorbar(scatter2, ax=axes[1])
            cbar2.set_label("local_k")
            figure_name = "layer_local_sweep.png"

        for ax in axes:
            ax.grid(alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(output_dir / figure_name, dpi=300, bbox_inches="tight")
        plt.close(fig)


def build_decision_summary(
    pilot_df: pd.DataFrame,
    sweep_df: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    decision = {
        "suite_date": pd.Timestamp.now().isoformat(),
        "pilot_best_by_rmse": None,
        "pilot_best_by_layer_ari": None,
        "best_sweep_by_rmse": None,
        "best_sweep_by_layer_ari": None,
        "recommendation": "",
    }
    rec: list[str] = []

    if not pilot_df.empty:
        best_rmse = pilot_df.sort_values("rmse", ascending=True).iloc[0]
        best_ari = pilot_df.sort_values("layer_ari", ascending=False).iloc[0]
        decision["pilot_best_by_rmse"] = best_rmse.to_dict()
        decision["pilot_best_by_layer_ari"] = best_ari.to_dict()
        rec.append(f"Pilot best RMSE: {best_rmse['run_name']}")
        rec.append(f"Pilot best layer ARI: {best_ari['run_name']}")

    if not sweep_df.empty:
        best_sweep_rmse = sweep_df.sort_values("rmse", ascending=True).iloc[0]
        best_sweep_ari = sweep_df.sort_values("layer_ari", ascending=False).iloc[0]
        decision["best_sweep_by_rmse"] = best_sweep_rmse.to_dict()
        decision["best_sweep_by_layer_ari"] = best_sweep_ari.to_dict()
        if best_sweep_rmse["gp_variant"] == "expr_product":
            rec.append(
                f"Best product sweep by RMSE: {best_sweep_rmse['run_name']} (lambda={best_sweep_rmse['gp_lambda_expr']}, offset={best_sweep_rmse['gp_product_offset']})"
            )
            rec.append(
                f"Best product sweep by layer ARI: {best_sweep_ari['run_name']} (lambda={best_sweep_ari['gp_lambda_expr']}, offset={best_sweep_ari['gp_product_offset']})"
            )
        elif best_sweep_rmse["gp_variant"] == "adaptive_additive":
            rec.append(
                f"Best adaptive sweep by RMSE: {best_sweep_rmse['run_name']} (lambda={best_sweep_rmse['gp_lambda_expr']}, tau={best_sweep_rmse['gp_gate_tau']})"
            )
            rec.append(
                f"Best adaptive sweep by layer ARI: {best_sweep_ari['run_name']} (lambda={best_sweep_ari['gp_lambda_expr']}, tau={best_sweep_ari['gp_gate_tau']})"
            )
        else:
            rec.append(
                f"Best layer-local sweep by RMSE: {best_sweep_rmse['run_name']} (lambda={best_sweep_rmse['gp_lambda_expr']}, k={best_sweep_rmse['gp_local_k']}, gate={best_sweep_rmse['gp_layer_gate_mode']})"
            )
            rec.append(
                f"Best layer-local sweep by layer ARI: {best_sweep_ari['run_name']} (lambda={best_sweep_ari['gp_lambda_expr']}, k={best_sweep_ari['gp_local_k']}, gate={best_sweep_ari['gp_layer_gate_mode']})"
            )

    decision["recommendation"] = " ".join(rec)
    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2))
    return decision


def main() -> None:
    args = parse_args()
    suite_dir = Path(args.output_dir)
    suite_dir.mkdir(parents=True, exist_ok=True)

    pilot_dir = suite_dir / "pilot"
    sweep_dir = suite_dir / "sweeps"
    pilot_dir.mkdir(exist_ok=True)
    sweep_dir.mkdir(exist_ok=True)

    if not args.skip_pilot:
        for spec in PILOT_VARIANTS:
            run_dir = pilot_dir / spec["dir_name"]
            cmd = build_command(
                run_dir,
                spec,
                data_dir=args.data_dir,
                n_genes=args.pilot_n_genes,
                min_observed_spots=args.pilot_min_observed_spots,
                mask_fraction=args.pilot_mask_fraction,
                seeds=args.pilot_seeds,
                mask_types=args.pilot_mask_types,
                gp_kernel=args.gp_kernel,
                gp_alpha=args.gp_alpha,
                calibration_bins=args.calibration_bins,
                knn_k=args.knn_k,
                n_domains=args.n_domains,
            )
            print(f"\n[pilot] Running {spec['name']}")
            run_one(cmd)

    sweep_specs_by_family = {
        "product": PRODUCT_SWEEPS,
        "adaptive": ADAPTIVE_SWEEPS,
        "layer_local": LAYER_LOCAL_SWEEPS,
    }
    sweep_specs = sweep_specs_by_family[args.sweep_family]

    if not args.skip_sweep:
        for spec in sweep_specs:
            run_dir = sweep_dir / spec["name"]
            cmd = build_command(
                run_dir,
                spec,
                data_dir=args.data_dir,
                n_genes=args.sweep_n_genes,
                min_observed_spots=args.sweep_min_observed_spots,
                mask_fraction=args.sweep_mask_fraction,
                seeds=args.sweep_seeds,
                mask_types=args.sweep_mask_types,
                gp_kernel=args.gp_kernel,
                gp_alpha=args.gp_alpha,
                calibration_bins=args.calibration_bins,
                knn_k=args.knn_k,
                n_domains=args.n_domains,
            )
            print(f"\n[sweep] Running {spec['name']}")
            run_one(cmd)

    pilot_df = summarize_runs(pilot_dir, PILOT_VARIANTS, "pilot_comparison.csv")
    sweep_summary_name = {
        "product": "product_sweep_summary.csv",
        "adaptive": "adaptive_sweep_summary.csv",
        "layer_local": "layer_local_sweep_summary.csv",
    }[args.sweep_family]
    sweep_df = summarize_runs(sweep_dir, sweep_specs, sweep_summary_name)
    plot_results(pilot_df, sweep_df, suite_dir / "figures")
    decision = build_decision_summary(pilot_df, sweep_df, suite_dir)

    print("\nGP variant suite decision summary:")
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()

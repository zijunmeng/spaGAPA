#!/usr/bin/env python
"""Validate current vs candidate high-resolution spaGAPA defaults.

This runner uses prepared real spatial APA data and the pipeline-level
pseudo-high-resolution smoke benchmark to compare two user-facing defaults:

* current_default: the conservative highres setting currently used by spaGAPA.
* candidate_default: the best-supported setting from the highres formal suite.

The goal is to decide whether candidate defaults are strong enough to promote
into the main package defaults before changing user-facing behavior.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_pipeline_highres_suite import (
    SUMMARY_METRICS,
    build_smoke_command,
    load_run,
    parse_float_list,
    resolve_project_path,
    run_exists,
    run_one,
)


DEFAULT_CONFIGS = [
    {
        "name": "current_default",
        "spatial_weight": 0.2,
        "expression_weight": 0.6,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.3,
    },
    {
        "name": "candidate_default",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.1,
    },
]


CONFIG_COLORS = {
    "current_default": "#7f8c8d",
    "candidate_default": "#2c7fb8",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/highres_default_validation_v1",
    )
    parser.add_argument("--data-dir", default="spaGAPA/data/processed/stapaminer_mob")
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--n-genes", type=int, default=24)
    parser.add_argument("--min-observed-spots", type=int, default=80)
    parser.add_argument("--max-parent-spots", type=int, default=80)
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--subbins-per-spot", default="2,4,8")
    parser.add_argument("--capture-rate", type=float, default=0.45)
    parser.add_argument("--dropout-rate", type=float, default=0.25)
    parser.add_argument(
        "--dropout-rates",
        default=None,
        help="Optional comma-separated dropout stress levels. Defaults to --dropout-rate only.",
    )
    parser.add_argument("--measurement-noise", type=float, default=0.08)
    parser.add_argument("--micro-noise", type=float, default=0.04)
    parser.add_argument("--jitter-fraction", type=float, default=0.18)
    parser.add_argument("--expression-noise", type=float, default=0.03)
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
    parser.add_argument("--highres-bioml-neighbor-mode", default="adaptive", choices=["fixed", "adaptive"])
    parser.add_argument("--highres-bioml-adaptive-neighbor-scale", type=float, default=10.0)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def dropout_values(args: argparse.Namespace) -> list[float]:
    if args.dropout_rates is None:
        return [float(args.dropout_rate)]
    return parse_float_list(args.dropout_rates)


def run_default_validation(args: argparse.Namespace, output_dir: Path) -> list[tuple[Path, dict[str, Any]]]:
    completed: list[tuple[Path, dict[str, Any]]] = []
    for config in DEFAULT_CONFIGS:
        for dropout in dropout_values(args):
            dropout_label = f"dropout_{dropout:.2f}".replace(".", "p")
            run_name = f"{config['name']}_{dropout_label}"
            run_dir = output_dir / "runs" / run_name
            cmd = build_smoke_command(
                run_dir,
                args,
                seeds=args.seeds,
                subbins=args.subbins_per_spot,
                dropout_rate=dropout,
                graph_spec=config,
            )
            print(f"\n[{config['name']}] Running {dropout_label}")
            if run_exists(run_dir) and not args.force:
                print(f"[{config['name']}] Skipping completed run. Use --force to re-run.")
            else:
                run_one(cmd)
            completed.append(
                (
                    run_dir,
                    {
                        "default_config": config["name"],
                        "suite_dropout_rate": float(dropout),
                        "highres_bioml_spatial_weight": config["spatial_weight"],
                        "highres_bioml_expression_weight": config["expression_weight"],
                        "highres_bioml_apa_weight": config["apa_weight"],
                        "highres_bioml_apa_source": config["apa_source"],
                        "highres_bioml_gp_blend": config["gp_blend"],
                    },
                )
            )
    return completed


def summarize(completed: list[tuple[Path, dict[str, Any]]], output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for run_dir, extra in completed:
        rows.extend(load_run(run_dir, extra))

    long_df = pd.DataFrame(rows)
    if long_df.empty:
        return long_df, pd.DataFrame()

    long_df.to_csv(output_dir / "highres_default_validation_results_long.csv", index=False)
    group_cols = [
        "default_config",
        "suite_dropout_rate",
        "subbins_per_spot",
        "method",
        "highres_bioml_spatial_weight",
        "highres_bioml_expression_weight",
        "highres_bioml_apa_weight",
        "highres_bioml_gp_blend",
    ]
    present = [col for col in group_cols if col in long_df.columns]
    metrics = [col for col in SUMMARY_METRICS if col in long_df.columns]
    summary = long_df.groupby(present, dropna=False)[metrics].agg(["mean", "std", "count"]).reset_index()
    summary.columns = [
        "_".join([str(part) for part in col if str(part)])
        if isinstance(col, tuple)
        else str(col)
        for col in summary.columns
    ]
    summary.to_csv(output_dir / "highres_default_validation_summary.csv", index=False)
    return long_df, summary


def plot_metric_bars(long_df: pd.DataFrame, figure_dir: Path) -> None:
    primary = long_df[long_df["method"] == "pipeline_highres_accuracy"].copy()
    if primary.empty:
        return
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "layer_nmi", "runtime_s"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 4.0))
    configs = [cfg["name"] for cfg in DEFAULT_CONFIGS if cfg["name"] in primary["default_config"].values]
    x = np.arange(len(configs))
    for ax, metric in zip(axes, metrics):
        means = []
        stds = []
        for config in configs:
            vals = primary.loc[primary["default_config"] == config, metric].dropna()
            means.append(float(vals.mean()) if not vals.empty else np.nan)
            stds.append(float(vals.std(ddof=0)) if len(vals) > 1 else 0.0)
        ax.bar(x, means, color=[CONFIG_COLORS.get(config, "#333333") for config in configs])
        ax.errorbar(x, means, yerr=stds, fmt="none", ecolor="#333333", capsize=3)
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(configs, rotation=35, ha="right")
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_default_validation_metrics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_scaling(long_df: pd.DataFrame, figure_dir: Path) -> None:
    primary = long_df[long_df["method"] == "pipeline_highres_accuracy"].copy()
    if primary.empty:
        return
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "layer_nmi"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 4.0))
    for ax, metric in zip(axes, metrics):
        for config in [cfg["name"] for cfg in DEFAULT_CONFIGS]:
            sub = primary[primary["default_config"] == config]
            grouped = sub.groupby("subbins_per_spot")[metric].agg(["mean", "std"]).reset_index()
            if grouped.empty:
                continue
            ax.errorbar(
                grouped["subbins_per_spot"],
                grouped["mean"],
                yerr=grouped["std"].fillna(0.0),
                marker="o",
                capsize=3,
                label=config,
                color=CONFIG_COLORS.get(config, "#333333"),
            )
        ax.set_xlabel("Pseudo-bins per parent spot")
        ax.set_title(metric)
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_default_validation_scaling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_decision(long_df: pd.DataFrame, output_dir: Path) -> dict[str, Any]:
    decision: dict[str, Any] = {
        "output_dir": str(output_dir),
        "summary_table": str(output_dir / "highres_default_validation_results_long.csv"),
        "overall_summary": str(output_dir / "highres_default_validation_summary.csv"),
        "figures": sorted(str(path) for path in (output_dir / "figures").glob("*.png")),
        "primary_method": "pipeline_highres_accuracy",
        "comparisons": {},
        "recommendation": "",
    }
    primary = long_df[long_df["method"] == "pipeline_highres_accuracy"].copy()
    if primary.empty:
        decision["recommendation"] = "No pipeline_highres_accuracy results found."
        (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
        return decision

    means = primary.groupby("default_config")[SUMMARY_METRICS].mean(numeric_only=True)
    if {"current_default", "candidate_default"}.issubset(means.index):
        current = means.loc["current_default"]
        candidate = means.loc["candidate_default"]
        comparison = {
            "rmse_delta_candidate_minus_current": float(candidate["rmse_holdout"] - current["rmse_holdout"]),
            "parent_rmse_delta_candidate_minus_current": float(candidate["parent_rmse"] - current["parent_rmse"]),
            "layer_ari_delta_candidate_minus_current": float(candidate["layer_ari"] - current["layer_ari"]),
            "layer_nmi_delta_candidate_minus_current": float(candidate["layer_nmi"] - current["layer_nmi"]),
            "runtime_ratio_candidate_vs_current": float(candidate["runtime_s"] / max(current["runtime_s"], 1e-9)),
            "current_default": current.to_dict(),
            "candidate_default": candidate.to_dict(),
        }
        decision["comparisons"]["candidate_vs_current"] = comparison
        if (
            comparison["rmse_delta_candidate_minus_current"] <= 0.002
            and comparison["layer_ari_delta_candidate_minus_current"] >= 0.0
            and comparison["layer_nmi_delta_candidate_minus_current"] >= 0.0
        ):
            decision["recommendation"] = (
                "Promote candidate highres defaults: gp_blend=0.1, "
                "spatial=0.1, expression=0.7, apa=0.2."
            )
        else:
            decision["recommendation"] = (
                "Keep current highres defaults for now; candidate did not meet the "
                "RMSE and biological-consistency gate."
            )

    by_subbin = {}
    for subbins, sub in primary.groupby("subbins_per_spot"):
        means_sub = sub.groupby("default_config")[SUMMARY_METRICS].mean(numeric_only=True)
        if {"current_default", "candidate_default"}.issubset(means_sub.index):
            current = means_sub.loc["current_default"]
            candidate = means_sub.loc["candidate_default"]
            by_subbin[str(int(subbins))] = {
                "rmse_delta_candidate_minus_current": float(candidate["rmse_holdout"] - current["rmse_holdout"]),
                "layer_ari_delta_candidate_minus_current": float(candidate["layer_ari"] - current["layer_ari"]),
                "layer_nmi_delta_candidate_minus_current": float(candidate["layer_nmi"] - current["layer_nmi"]),
            }
    decision["comparisons_by_subbins"] = by_subbin

    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
    return decision


def main() -> None:
    args = parse_args()
    output_dir = resolve_project_path(args.output_dir)
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    completed = run_default_validation(args, output_dir)
    long_df, summary = summarize(completed, output_dir)
    if not long_df.empty:
        plot_metric_bars(long_df, figure_dir)
        plot_scaling(long_df, figure_dir)
    decision = build_decision(long_df, output_dir)

    print("\nHigh-resolution default validation decision summary:")
    print(json.dumps(decision, indent=2, default=str))
    print(f"\nOutputs: {output_dir}")
    if not summary.empty:
        print(f"Summary: {output_dir / 'highres_default_validation_summary.csv'}")


if __name__ == "__main__":
    main()

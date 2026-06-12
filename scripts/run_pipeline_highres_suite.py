#!/usr/bin/env python
"""Pipeline-level high-resolution spaGAPA suite.

This suite promotes the pipeline smoke runner into a compact formal benchmark:

1. Multi-seed pipeline takeover validation.
2. 2x / 4x / 8x scaling.
3. Dropout stress.
4. High-resolution graph-weight sweep.
5. GP-blend / fast-mode tradeoff.

All method outputs come through ``SpaGAPA(analysis_preset='highres_accuracy')``
or ``SpaGAPA(analysis_preset='highres_fast')`` inside
``run_pipeline_highres_smoke.py``.
"""

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
SMOKE_SCRIPT = PACKAGE_ROOT / "scripts" / "run_pipeline_highres_smoke.py"


GRAPH_SWEEPS = [
    {
        "name": "default_s010_e070_a020_exprknn",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.1,
    },
    {
        "name": "legacy_default_s020_e060_a020_exprknn",
        "spatial_weight": 0.2,
        "expression_weight": 0.6,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.3,
    },
    {
        "name": "spatial_heavy_s040_e040_a020_exprknn",
        "spatial_weight": 0.4,
        "expression_weight": 0.4,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.3,
    },
    {
        "name": "expr_heavy_s010_e070_a020_exprknn",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.3,
    },
    {
        "name": "no_apa_s040_e060_a000",
        "spatial_weight": 0.4,
        "expression_weight": 0.6,
        "apa_weight": 0.0,
        "apa_source": "none",
        "gp_blend": 0.3,
    },
    {
        "name": "raw_apa_s020_e060_a020",
        "spatial_weight": 0.2,
        "expression_weight": 0.6,
        "apa_weight": 0.2,
        "apa_source": "raw",
        "gp_blend": 0.3,
    },
]


GP_BLEND_SWEEPS = [
    {
        "name": "gp_blend_0.0_fast_like",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.0,
    },
    {
        "name": "gp_blend_0.1",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.1,
    },
    {
        "name": "gp_blend_0.3",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.3,
    },
    {
        "name": "gp_blend_0.5",
        "spatial_weight": 0.1,
        "expression_weight": 0.7,
        "apa_weight": 0.2,
        "apa_source": "expression_knn",
        "gp_blend": 0.5,
    },
]


SUMMARY_METRICS = [
    "rmse_holdout",
    "parent_rmse",
    "layer_ari",
    "layer_nmi",
    "runtime_s",
]


def resolve_project_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    if path.parts and path.parts[0] == PACKAGE_ROOT.name:
        return (PROJECT_ROOT / path).resolve()
    cwd_path = (Path.cwd() / path).resolve()
    if cwd_path.exists():
        return cwd_path
    package_path = (PACKAGE_ROOT / path).resolve()
    if package_path.exists():
        return package_path
    if Path.cwd().resolve() == PACKAGE_ROOT:
        return package_path
    return (PROJECT_ROOT / path).resolve()


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("Expected at least one integer value")
    return values


def parse_float_list(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("Expected at least one float value")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="spaGAPA/benchmark_results/real/pipeline_highres_suite_v1")
    parser.add_argument("--data-dir", default="spaGAPA/data/processed/stapaminer_mob")
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--n-genes", type=int, default=32)
    parser.add_argument("--min-observed-spots", type=int, default=80)
    parser.add_argument("--max-parent-spots", type=int, default=80)
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--formal-subbins", default="4")
    parser.add_argument("--scaling-subbins", default="2,4,8")
    parser.add_argument("--capture-rate", type=float, default=0.45)
    parser.add_argument("--dropout-rate", type=float, default=0.25)
    parser.add_argument("--dropout-rates", default="0.10,0.25,0.40")
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
    parser.add_argument("--skip-formal", action="store_true")
    parser.add_argument("--skip-scaling", action="store_true")
    parser.add_argument("--skip-dropout", action="store_true")
    parser.add_argument("--skip-graph-sweep", action="store_true")
    parser.add_argument("--skip-blend-sweep", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def build_smoke_command(
    run_dir: Path,
    args: argparse.Namespace,
    *,
    seeds: str,
    subbins: str,
    dropout_rate: float,
    graph_spec: dict[str, Any],
) -> list[str]:
    cmd = [
        sys.executable,
        str(SMOKE_SCRIPT),
        "--data-dir",
        str(resolve_project_path(args.data_dir)),
        "--output-dir",
        str(run_dir.resolve()),
        "--n-genes",
        str(args.n_genes),
        "--min-observed-spots",
        str(args.min_observed_spots),
        "--max-parent-spots",
        str(args.max_parent_spots),
        "--subbins-per-spot",
        str(subbins),
        "--seeds",
        str(seeds),
        "--capture-rate",
        str(args.capture_rate),
        "--dropout-rate",
        str(dropout_rate),
        "--measurement-noise",
        str(args.measurement_noise),
        "--micro-noise",
        str(args.micro_noise),
        "--jitter-fraction",
        str(args.jitter_fraction),
        "--expression-noise",
        str(args.expression_noise),
        "--knn-k",
        str(args.knn_k),
        "--expr-n-components",
        str(args.expr_n_components),
        "--sparse-n-inducing",
        str(args.sparse_n_inducing),
        "--sparse-length-scale",
        str(args.sparse_length_scale),
        "--sparse-length-scale-multiplier",
        str(args.sparse_length_scale_multiplier),
        "--sparse-noise-level",
        str(args.sparse_noise_level),
        "--bioml-rank",
        str(args.bioml_rank),
        "--bioml-lambda-graph",
        str(args.bioml_lambda_graph),
        "--bioml-lambda-l2",
        str(args.bioml_lambda_l2),
        "--bioml-max-iter",
        str(args.bioml_max_iter),
        "--bioml-n-neighbors",
        str(args.bioml_n_neighbors),
        "--highres-bioml-neighbor-mode",
        str(args.highres_bioml_neighbor_mode),
        "--highres-bioml-adaptive-neighbor-scale",
        str(args.highres_bioml_adaptive_neighbor_scale),
        "--highres-bioml-gp-blend",
        str(graph_spec["gp_blend"]),
        "--highres-bioml-spatial-weight",
        str(graph_spec["spatial_weight"]),
        "--highres-bioml-expression-weight",
        str(graph_spec["expression_weight"]),
        "--highres-bioml-apa-weight",
        str(graph_spec["apa_weight"]),
        "--highres-bioml-apa-source",
        str(graph_spec["apa_source"]),
    ]
    if args.dataset_name is not None:
        cmd.extend(["--dataset-name", str(args.dataset_name)])
    if args.n_domains is not None:
        cmd.extend(["--n-domains", str(args.n_domains)])
    return cmd


def run_exists(run_dir: Path) -> bool:
    return (run_dir / "pipeline_highres_smoke_summary.csv").exists()


def run_one(cmd: list[str]) -> None:
    result = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)


def iter_run_specs(args: argparse.Namespace) -> list[tuple[str, str, dict[str, Any], str, str, float]]:
    default = GRAPH_SWEEPS[0]
    first_seed = str(parse_int_list(args.seeds)[0])
    specs: list[tuple[str, str, dict[str, Any], str, str, float]] = []
    if not args.skip_formal:
        specs.append(("formal", "formal_all_seeds", default, args.seeds, args.formal_subbins, args.dropout_rate))
    if not args.skip_scaling:
        specs.append(("scaling", "scaling_all_seeds", default, args.seeds, args.scaling_subbins, args.dropout_rate))
    if not args.skip_dropout:
        for dropout in parse_float_list(args.dropout_rates):
            name = f"dropout_{dropout:.2f}".replace(".", "p")
            specs.append(("dropout", name, default, first_seed, args.formal_subbins, dropout))
    if not args.skip_graph_sweep:
        for graph_spec in GRAPH_SWEEPS:
            specs.append(("graph_sweep", graph_spec["name"], graph_spec, first_seed, args.formal_subbins, args.dropout_rate))
    if not args.skip_blend_sweep:
        for graph_spec in GP_BLEND_SWEEPS:
            specs.append(("gp_blend_sweep", graph_spec["name"], graph_spec, first_seed, args.formal_subbins, args.dropout_rate))
    return specs


def load_run(run_dir: Path, extra: dict[str, Any]) -> list[dict[str, Any]]:
    path = run_dir / "pipeline_highres_smoke_summary.csv"
    if not path.exists():
        return []
    df = pd.read_csv(path)
    records = []
    for _, row in df.iterrows():
        record = row.to_dict()
        record.update(extra)
        records.append(record)
    return records


def summarize_suite(suite_dir: Path, completed: list[tuple[Path, dict[str, Any]]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for run_dir, extra in completed:
        rows.extend(load_run(run_dir, extra))
    long_df = pd.DataFrame(rows)
    if long_df.empty:
        return long_df, pd.DataFrame()

    long_df.to_csv(suite_dir / "pipeline_highres_suite_results_long.csv", index=False)
    group_cols = [
        "run_type",
        "config_name",
        "method",
        "highres_bioml_spatial_weight",
        "highres_bioml_expression_weight",
        "highres_bioml_apa_weight",
        "highres_bioml_apa_source",
        "highres_bioml_gp_blend",
    ]
    present = [col for col in group_cols if col in long_df.columns]
    metrics = [col for col in SUMMARY_METRICS if col in long_df.columns]
    overall = long_df.groupby(present, dropna=False)[metrics].agg(["mean", "std", "count"]).reset_index()
    overall.columns = [
        "_".join([str(part) for part in col if str(part)])
        if isinstance(col, tuple)
        else str(col)
        for col in overall.columns
    ]
    overall.to_csv(suite_dir / "pipeline_highres_suite_overall_summary.csv", index=False)
    return long_df, overall


def method_mean(df: pd.DataFrame, method: str) -> pd.Series | None:
    sub = df[df["method"] == method]
    if sub.empty:
        return None
    return sub.select_dtypes(include=[np.number]).mean()


def best_record(df: pd.DataFrame, metric: str, ascending: bool) -> dict[str, Any] | None:
    if metric not in df.columns:
        return None
    valid = df[df[metric].notna()]
    if valid.empty:
        return None
    return valid.sort_values(metric, ascending=ascending).iloc[0].to_dict()


def plot_method_bars(formal: pd.DataFrame, figure_dir: Path) -> None:
    if formal.empty:
        return
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "layer_nmi", "runtime_s"]
    methods = ["runner_highres_bioml", "pipeline_highres_accuracy", "pipeline_highres_fast"]
    colors = ["#1f9d8a", "#2c7fb8", "#f39c12"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 4.2))
    for ax, metric in zip(axes, metrics):
        vals = []
        errs = []
        labels = []
        bar_colors = []
        for method, color in zip(methods, colors):
            sub = formal[formal["method"] == method]
            if sub.empty:
                continue
            vals.append(float(sub[metric].mean()))
            errs.append(float(sub[metric].std(ddof=0)) if len(sub) > 1 else 0.0)
            labels.append(method)
            bar_colors.append(color)
        x = np.arange(len(vals))
        ax.bar(x, vals, color=bar_colors)
        ax.errorbar(x, vals, yerr=errs, fmt="none", ecolor="#333333", capsize=3)
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "pipeline_highres_suite_formal.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_scaling(long_df: pd.DataFrame, figure_dir: Path) -> None:
    scaling = long_df[long_df["run_type"] == "scaling"].copy()
    if scaling.empty:
        return
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "runtime_s"]
    methods = ["runner_highres_bioml", "pipeline_highres_accuracy", "pipeline_highres_fast"]
    colors = {
        "runner_highres_bioml": "#1f9d8a",
        "pipeline_highres_accuracy": "#2c7fb8",
        "pipeline_highres_fast": "#f39c12",
    }
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.2 * len(metrics), 4.0))
    for ax, metric in zip(axes, metrics):
        for method in methods:
            sub = scaling[scaling["method"] == method]
            grouped = sub.groupby("subbins_per_spot")[metric].agg(["mean", "std"]).reset_index()
            if grouped.empty:
                continue
            ax.errorbar(
                grouped["subbins_per_spot"],
                grouped["mean"],
                yerr=grouped["std"].fillna(0.0),
                marker="o",
                capsize=3,
                label=method,
                color=colors[method],
            )
        ax.set_xlabel("Pseudo-bins per parent spot")
        ax.set_title(metric)
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "pipeline_highres_suite_scaling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sweep(long_df: pd.DataFrame, run_type: str, figure_name: str, figure_dir: Path) -> None:
    sweep = long_df[
        (long_df["run_type"] == run_type)
        & (long_df["method"].isin(["pipeline_highres_accuracy", "pipeline_highres_fast"]))
    ].copy()
    if sweep.empty:
        return
    summary = (
        sweep.groupby(["config_name", "method"])[["rmse_holdout", "layer_ari", "runtime_s"]]
        .mean()
        .reset_index()
    )
    labels = summary["config_name"] + "\n" + summary["method"].str.replace("pipeline_highres_", "")
    x = np.arange(len(summary))
    fig, axes = plt.subplots(1, 3, figsize=(max(12, 0.7 * len(summary)), 4.2))
    for ax, metric, color in zip(axes, ["rmse_holdout", "layer_ari", "runtime_s"], ["#3182bd", "#31a354", "#f16913"]):
        ax.bar(x, summary[metric], color=color)
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
        ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / figure_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_suite(long_df: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plot_method_bars(long_df[long_df["run_type"] == "formal"], figure_dir)
    plot_scaling(long_df, figure_dir)
    plot_sweep(long_df, "dropout", "pipeline_highres_suite_dropout.png", figure_dir)
    plot_sweep(long_df, "graph_sweep", "pipeline_highres_suite_graph_sweep.png", figure_dir)
    plot_sweep(long_df, "gp_blend_sweep", "pipeline_highres_suite_gp_blend.png", figure_dir)


def build_decision_summary(long_df: pd.DataFrame, suite_dir: Path) -> dict[str, Any]:
    decision: dict[str, Any] = {
        "suite_date": pd.Timestamp.now().isoformat(),
        "output_dir": str(suite_dir),
        "summary_table": str(suite_dir / "pipeline_highres_suite_results_long.csv"),
        "overall_summary": str(suite_dir / "pipeline_highres_suite_overall_summary.csv"),
        "figures": sorted(str(path) for path in (suite_dir / "figures").glob("*.png")),
        "comparisons": {},
        "best_graph_sweep_by_layer_ari": None,
        "best_graph_sweep_by_rmse": None,
        "best_gp_blend_by_rmse": None,
        "recommendation": "",
    }
    if long_df.empty:
        decision["recommendation"] = "No completed pipeline high-resolution suite runs found."
        (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
        return decision

    formal = long_df[long_df["run_type"] == "formal"]
    runner = method_mean(formal, "runner_highres_bioml")
    for method in ["pipeline_highres_accuracy", "pipeline_highres_fast"]:
        row = method_mean(formal, method)
        if row is None or runner is None:
            continue
        decision["comparisons"][method] = {
            "rmse_delta_vs_runner": float(row["rmse_holdout"] - runner["rmse_holdout"]),
            "parent_rmse_delta_vs_runner": float(row["parent_rmse"] - runner["parent_rmse"]),
            "layer_ari_delta_vs_runner": float(row["layer_ari"] - runner["layer_ari"]),
            "layer_nmi_delta_vs_runner": float(row["layer_nmi"] - runner["layer_nmi"]),
            "runtime_ratio_vs_runner": float(row["runtime_s"] / max(runner["runtime_s"], 1e-9)),
        }

    for run_type, key in [
        ("graph_sweep", "best_graph_sweep"),
        ("gp_blend_sweep", "best_gp_blend"),
    ]:
        sub = long_df[(long_df["run_type"] == run_type) & (long_df["method"] == "pipeline_highres_accuracy")]
        if sub.empty:
            continue
        grouped = sub.groupby("config_name", as_index=False)[SUMMARY_METRICS].mean(numeric_only=True)
        decision[f"{key}_by_rmse"] = best_record(grouped, "rmse_holdout", True)
        if run_type == "graph_sweep":
            decision[f"{key}_by_layer_ari"] = best_record(grouped, "layer_ari", False)

    notes = []
    comp = decision["comparisons"].get("pipeline_highres_accuracy")
    if comp:
        notes.append(
            "pipeline_highres_accuracy vs runner: "
            f"delta_RMSE={comp['rmse_delta_vs_runner']:.4g}, "
            f"delta_ARI={comp['layer_ari_delta_vs_runner']:.4g}."
        )
    fast = decision["comparisons"].get("pipeline_highres_fast")
    if fast:
        notes.append(
            "pipeline_highres_fast tradeoff: "
            f"delta_RMSE={fast['rmse_delta_vs_runner']:.4g}, "
            f"delta_ARI={fast['layer_ari_delta_vs_runner']:.4g}, "
            f"runtime_ratio={fast['runtime_ratio_vs_runner']:.3g}."
        )
    best_graph = decision.get("best_graph_sweep_by_layer_ari")
    if best_graph:
        notes.append(f"Best graph sweep by layer ARI: {best_graph['config_name']}.")
    best_blend = decision.get("best_gp_blend_by_rmse")
    if best_blend:
        notes.append(f"Best GP blend by RMSE: {best_blend['config_name']}.")
    decision["recommendation"] = " ".join(notes)
    (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
    return decision


def main() -> None:
    args = parse_args()
    suite_dir = resolve_project_path(args.output_dir)
    figure_dir = suite_dir / "figures"
    suite_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(exist_ok=True)

    completed: list[tuple[Path, dict[str, Any]]] = []
    for run_type, run_name, graph_spec, seeds, subbins, dropout_rate in iter_run_specs(args):
        run_dir = suite_dir / run_type / run_name
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        cmd = build_smoke_command(
            run_dir,
            args,
            seeds=seeds,
            subbins=subbins,
            dropout_rate=dropout_rate,
            graph_spec=graph_spec,
        )
        print(f"\n[{run_type}] Running {run_name}")
        if run_exists(run_dir) and not args.force:
            print(f"[{run_type}] Skipping completed run {run_name}. Use --force to re-run.")
        else:
            run_one(cmd)
        extra = {
            "run_type": run_type,
            "run_name": run_name,
            "config_name": graph_spec["name"],
            "suite_seeds": seeds,
            "suite_subbins": subbins,
            "suite_dropout_rate": dropout_rate,
            "highres_bioml_spatial_weight": graph_spec["spatial_weight"],
            "highres_bioml_expression_weight": graph_spec["expression_weight"],
            "highres_bioml_apa_weight": graph_spec["apa_weight"],
            "highres_bioml_apa_source": graph_spec["apa_source"],
            "highres_bioml_gp_blend": graph_spec["gp_blend"],
        }
        completed.append((run_dir, extra))

    long_df, overall_df = summarize_suite(suite_dir, completed)
    if not long_df.empty:
        plot_suite(long_df, figure_dir)
    decision = build_decision_summary(long_df, suite_dir)
    print("\nPipeline high-resolution suite decision summary:")
    print(json.dumps(decision, indent=2, default=str))
    print(f"\nSuite outputs: {suite_dir}")
    if not overall_df.empty:
        print(f"Overall summary: {suite_dir / 'pipeline_highres_suite_overall_summary.csv'}")
    if not long_df.empty:
        print(f"Long summary: {suite_dir / 'pipeline_highres_suite_results_long.csv'}")


if __name__ == "__main__":
    main()

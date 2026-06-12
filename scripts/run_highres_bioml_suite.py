#!/usr/bin/env python
"""Run a formal high-resolution BioML benchmark suite.

This suite wraps ``run_high_resolution_simulation.py`` and focuses on the
decoupled ``highres_bioml`` route:

1. Formal multi-seed comparison against raw, expression KNN, sparse GP, and
   direct sparse BioML.
2. Pseudo-bin scaling and dropout stress tests.
3. Graph-weight and GP-blend sweeps for accuracy, biological consistency, and
   runtime trade-offs.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
HIGHRES_SCRIPT = PACKAGE_ROOT / "scripts" / "run_high_resolution_simulation.py"


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
    {
        "name": "sparsegp_apa_s020_e060_a020",
        "spatial_weight": 0.2,
        "expression_weight": 0.6,
        "apa_weight": 0.2,
        "apa_source": "sparse_gp",
        "gp_blend": 0.3,
    },
]


GP_BLEND_SWEEPS = [
    {
        "name": "gp_blend_0.0",
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
    "layer_purity",
    "uncertainty_error_spearman",
    "runtime_s",
    "peak_rss_mb",
]


def resolve_project_path(path: str | Path) -> Path:
    """Resolve paths relative to the project root."""
    path = Path(path)
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_int_list(text: str) -> list[int]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    if not values:
        raise ValueError("Expected at least one integer value")
    return values


def parse_float_list(text: str) -> list[float]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    if not values:
        raise ValueError("Expected at least one float value")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/benchmark_results/real/highres_bioml_suite_v1",
        help="Suite output directory.",
    )
    parser.add_argument(
        "--data-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Prepared real dataset directory.",
    )
    parser.add_argument("--n-genes", type=int, default=40)
    parser.add_argument("--min-observed-spots", type=int, default=100)
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
    parser.add_argument("--bioml-blend", type=float, default=0.1)
    parser.add_argument(
        "--highres-bioml-neighbor-mode",
        default="adaptive",
        choices=["fixed", "adaptive"],
    )
    parser.add_argument("--highres-bioml-adaptive-neighbor-scale", type=float, default=10.0)
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument(
        "--formal-methods",
        default="raw,expression_knn,sparse_gp,sparse_bioml,highres_bioml",
        help="Methods for formal comparison.",
    )
    parser.add_argument(
        "--sweep-methods",
        default="raw,expression_knn,highres_bioml",
        help="Methods for graph/blend/dropout/scaling sweeps.",
    )
    parser.add_argument("--skip-formal", action="store_true")
    parser.add_argument("--skip-scaling", action="store_true")
    parser.add_argument("--skip-dropout", action="store_true")
    parser.add_argument("--skip-graph-sweep", action="store_true")
    parser.add_argument("--skip-blend-sweep", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run completed run directories.",
    )
    return parser.parse_args()


def build_command(
    run_dir: Path,
    args: argparse.Namespace,
    *,
    seed: int,
    subbins: str,
    dropout_rate: float,
    methods: str,
    graph_spec: dict[str, Any],
) -> list[str]:
    cmd = [
        sys.executable,
        str(HIGHRES_SCRIPT),
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
        "--seed",
        str(seed),
        "--methods",
        methods,
        "--layer-column",
        str(args.layer_column),
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
        "--bioml-blend",
        str(args.bioml_blend),
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
    if args.n_domains is not None:
        cmd.extend(["--n-domains", str(args.n_domains)])
    return cmd


def run_exists(run_dir: Path) -> bool:
    return (run_dir / "highres_results_summary.csv").exists()


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


def iter_formal_runs(args: argparse.Namespace) -> Iterable[tuple[str, str, dict[str, Any], int, str, float, str]]:
    default = GRAPH_SWEEPS[0]
    for seed in parse_int_list(args.seeds):
        yield (
            "formal",
            f"formal_seed{seed}",
            default,
            seed,
            args.formal_subbins,
            args.dropout_rate,
            args.formal_methods,
        )


def iter_scaling_runs(args: argparse.Namespace) -> Iterable[tuple[str, str, dict[str, Any], int, str, float, str]]:
    default = GRAPH_SWEEPS[0]
    seed = parse_int_list(args.seeds)[0]
    for subbins in parse_int_list(args.scaling_subbins):
        yield (
            "scaling",
            f"scaling_subbins{subbins}_seed{seed}",
            default,
            seed,
            str(subbins),
            args.dropout_rate,
            args.sweep_methods,
        )


def iter_dropout_runs(args: argparse.Namespace) -> Iterable[tuple[str, str, dict[str, Any], int, str, float, str]]:
    default = GRAPH_SWEEPS[0]
    seed = parse_int_list(args.seeds)[0]
    for dropout in parse_float_list(args.dropout_rates):
        name = f"dropout_{dropout:.2f}".replace(".", "p")
        yield (
            "dropout",
            f"{name}_seed{seed}",
            default,
            seed,
            args.formal_subbins,
            dropout,
            args.sweep_methods,
        )


def iter_graph_sweep_runs(args: argparse.Namespace) -> Iterable[tuple[str, str, dict[str, Any], int, str, float, str]]:
    seed = parse_int_list(args.seeds)[0]
    for spec in GRAPH_SWEEPS:
        yield (
            "graph_sweep",
            f"{spec['name']}_seed{seed}",
            spec,
            seed,
            args.formal_subbins,
            args.dropout_rate,
            args.sweep_methods,
        )


def iter_blend_sweep_runs(args: argparse.Namespace) -> Iterable[tuple[str, str, dict[str, Any], int, str, float, str]]:
    seed = parse_int_list(args.seeds)[0]
    for spec in GP_BLEND_SWEEPS:
        yield (
            "gp_blend_sweep",
            f"{spec['name']}_seed{seed}",
            spec,
            seed,
            args.formal_subbins,
            args.dropout_rate,
            args.sweep_methods,
        )


def load_run(run_dir: Path, extra: dict[str, Any]) -> list[dict[str, Any]]:
    if not run_exists(run_dir):
        return []
    df = pd.read_csv(run_dir / "highres_results_summary.csv")
    rows = []
    for _, row in df.iterrows():
        record = row.to_dict()
        record.update(extra)
        rows.append(record)
    return rows


def summarize_runs(suite_dir: Path, completed: list[tuple[Path, dict[str, Any]]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for run_dir, extra in completed:
        rows.extend(load_run(run_dir, extra))

    long_df = pd.DataFrame(rows)
    if long_df.empty:
        return long_df, pd.DataFrame()

    long_df.to_csv(suite_dir / "highres_bioml_suite_results_long.csv", index=False)

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
    present_group_cols = [col for col in group_cols if col in long_df.columns]
    metric_cols = [col for col in SUMMARY_METRICS if col in long_df.columns]
    overall = (
        long_df.groupby(present_group_cols, dropna=False)[metric_cols]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    overall.columns = [
        "_".join([part for part in col if part]) if isinstance(col, tuple) else col
        for col in overall.columns
    ]
    overall.to_csv(suite_dir / "highres_bioml_suite_overall_summary.csv", index=False)
    return long_df, overall


def metric_mean(df: pd.DataFrame, metric: str) -> float:
    return float(df[metric].mean()) if metric in df.columns and df[metric].notna().any() else np.nan


def barplot_methods(formal: pd.DataFrame, figure_dir: Path) -> None:
    if formal.empty:
        return
    methods = ["raw", "expression_knn", "sparse_gp", "sparse_bioml", "highres_bioml"]
    metrics = ["rmse_holdout", "parent_rmse", "layer_ari", "layer_nmi", "runtime_s"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 4.2))
    for ax, metric in zip(axes, metrics):
        vals = []
        labels = []
        for method in methods:
            sub = formal[formal["method"] == method]
            if sub.empty:
                continue
            vals.append(metric_mean(sub, metric))
            labels.append(method)
        ax.bar(np.arange(len(vals)), vals, color="#1f9d8a")
        ax.set_title(metric)
        ax.set_xticks(np.arange(len(vals)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.grid(axis="y", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_formal_method_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_sweep(df: pd.DataFrame, run_type: str, figure_name: str, figure_dir: Path) -> None:
    sweep = df[(df["run_type"] == run_type) & (df["method"] == "highres_bioml")].copy()
    if sweep.empty:
        return
    summary = (
        sweep.groupby("config_name")[["rmse_holdout", "layer_ari", "runtime_s"]]
        .mean()
        .reset_index()
        .sort_values("layer_ari", ascending=False)
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    x = np.arange(len(summary))
    labels = summary["config_name"].str.replace("_", "\n")
    for ax, metric, color in zip(axes, ["rmse_holdout", "layer_ari", "runtime_s"], ["#3182bd", "#31a354", "#f16913"]):
        ax.bar(x, summary[metric], color=color)
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.grid(axis="y", alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_dir / figure_name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_scaling(df: pd.DataFrame, figure_dir: Path) -> None:
    scaling = df[(df["run_type"] == "scaling") & (df["method"] == "highres_bioml")].copy()
    if scaling.empty:
        return
    scaling = scaling.sort_values("n_bins")
    metrics = ["rmse_holdout", "layer_ari", "runtime_s", "peak_rss_mb"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.0 * len(metrics), 4.0))
    for ax, metric in zip(axes, metrics):
        ax.plot(scaling["n_bins"], scaling[metric], marker="o", color="#1f9d8a")
        ax.set_xlabel("Pseudo-bins")
        ax.set_title(metric)
        ax.grid(alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(figure_dir / "highres_bioml_scaling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_suite(long_df: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    formal = long_df[long_df["run_type"] == "formal"]
    barplot_methods(formal, figure_dir)
    plot_sweep(long_df, "dropout", "highres_dropout_stress.png", figure_dir)
    plot_sweep(long_df, "graph_sweep", "highres_graph_weight_sweep.png", figure_dir)
    plot_sweep(long_df, "gp_blend_sweep", "highres_gp_blend_sweep.png", figure_dir)
    plot_scaling(long_df, figure_dir)


def best_record(df: pd.DataFrame, metric: str, ascending: bool) -> dict[str, Any] | None:
    valid = df[df[metric].notna()] if metric in df.columns else pd.DataFrame()
    if valid.empty:
        return None
    return valid.sort_values(metric, ascending=ascending).iloc[0].to_dict()


def method_mean(formal: pd.DataFrame, method: str) -> pd.Series | None:
    sub = formal[formal["method"] == method]
    if sub.empty:
        return None
    numeric = sub.select_dtypes(include=[np.number])
    return numeric.mean()


def build_decision_summary(long_df: pd.DataFrame, suite_dir: Path) -> dict[str, Any]:
    decision: dict[str, Any] = {
        "suite_date": pd.Timestamp.now().isoformat(),
        "formal_best_by_rmse": None,
        "formal_best_by_layer_ari": None,
        "highres_bioml_vs_expression_knn": None,
        "best_graph_sweep_by_layer_ari": None,
        "best_graph_sweep_by_rmse": None,
        "best_gp_blend_by_layer_ari": None,
        "best_gp_blend_by_rmse": None,
        "recommendation": "",
    }
    if long_df.empty:
        decision["recommendation"] = "No completed high-resolution suite runs found."
        (suite_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))
        return decision

    formal = long_df[long_df["run_type"] == "formal"]
    if not formal.empty:
        formal_method = formal.groupby("method", as_index=False)[SUMMARY_METRICS].mean(numeric_only=True)
        decision["formal_best_by_rmse"] = best_record(formal_method, "rmse_holdout", True)
        decision["formal_best_by_layer_ari"] = best_record(formal_method, "layer_ari", False)
        highres = method_mean(formal, "highres_bioml")
        expr = method_mean(formal, "expression_knn")
        raw = method_mean(formal, "raw")
        if highres is not None and expr is not None:
            decision["highres_bioml_vs_expression_knn"] = {
                "rmse_delta": float(highres["rmse_holdout"] - expr["rmse_holdout"]),
                "layer_ari_delta": float(highres["layer_ari"] - expr["layer_ari"]),
                "layer_nmi_delta": float(highres["layer_nmi"] - expr["layer_nmi"]),
                "runtime_ratio": float(highres["runtime_s"] / max(expr["runtime_s"], 1e-9)),
            }
        if highres is not None and raw is not None:
            decision["highres_bioml_vs_raw"] = {
                "rmse_delta": float(highres["rmse_holdout"] - raw["rmse_holdout"]),
                "parent_rmse_delta": float(highres["parent_rmse"] - raw["parent_rmse"]),
                "layer_ari_delta": float(highres["layer_ari"] - raw["layer_ari"]),
            }

    for run_type, key_prefix in [
        ("graph_sweep", "best_graph_sweep"),
        ("gp_blend_sweep", "best_gp_blend"),
    ]:
        sub = long_df[(long_df["run_type"] == run_type) & (long_df["method"] == "highres_bioml")]
        if sub.empty:
            continue
        grouped = sub.groupby("config_name", as_index=False)[SUMMARY_METRICS].mean(numeric_only=True)
        decision[f"{key_prefix}_by_layer_ari"] = best_record(grouped, "layer_ari", False)
        decision[f"{key_prefix}_by_rmse"] = best_record(grouped, "rmse_holdout", True)

    notes = []
    comp = decision.get("highres_bioml_vs_expression_knn")
    raw_comp = decision.get("highres_bioml_vs_raw")
    if comp:
        notes.append(
            "highres_bioml vs expression_knn: "
            f"delta_ARI={comp['layer_ari_delta']:.4f}, "
            f"delta_RMSE={comp['rmse_delta']:.4f}, "
            f"runtime_ratio={comp['runtime_ratio']:.1f}x."
        )
    if raw_comp:
        notes.append(
            "highres_bioml vs raw: "
            f"delta_RMSE={raw_comp['rmse_delta']:.4f}, "
            f"delta_ARI={raw_comp['layer_ari_delta']:.4f}."
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

    run_specs: list[tuple[str, str, dict[str, Any], int, str, float, str]] = []
    if not args.skip_formal:
        run_specs.extend(iter_formal_runs(args))
    if not args.skip_scaling:
        run_specs.extend(iter_scaling_runs(args))
    if not args.skip_dropout:
        run_specs.extend(iter_dropout_runs(args))
    if not args.skip_graph_sweep:
        run_specs.extend(iter_graph_sweep_runs(args))
    if not args.skip_blend_sweep:
        run_specs.extend(iter_blend_sweep_runs(args))

    completed: list[tuple[Path, dict[str, Any]]] = []
    for run_type, run_name, graph_spec, seed, subbins, dropout_rate, methods in run_specs:
        run_dir = suite_dir / run_type / run_name
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        cmd = build_command(
            run_dir,
            args,
            seed=seed,
            subbins=subbins,
            dropout_rate=dropout_rate,
            methods=methods,
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
            "seed": seed,
            "suite_subbins": str(subbins),
            "suite_dropout_rate": dropout_rate,
            "suite_methods": methods,
            "highres_bioml_spatial_weight": graph_spec["spatial_weight"],
            "highres_bioml_expression_weight": graph_spec["expression_weight"],
            "highres_bioml_apa_weight": graph_spec["apa_weight"],
            "highres_bioml_apa_source": graph_spec["apa_source"],
            "highres_bioml_gp_blend": graph_spec["gp_blend"],
        }
        completed.append((run_dir, extra))

    long_df, overall_df = summarize_runs(suite_dir, completed)
    if not long_df.empty:
        plot_suite(long_df, figure_dir)
    decision = build_decision_summary(long_df, suite_dir)

    print("\nHigh-resolution BioML suite decision summary:")
    print(json.dumps(decision, indent=2, default=str))
    print(f"\nSuite outputs: {suite_dir}")
    if not overall_df.empty:
        print(f"Overall summary: {suite_dir / 'highres_bioml_suite_overall_summary.csv'}")
    if not long_df.empty:
        print(f"Long summary: {suite_dir / 'highres_bioml_suite_results_long.csv'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""BIB-oriented real-data benchmark suite orchestrator.

The first responsibility of this script is readiness: discover prepared real
datasets and report which BIB benchmark tracks can run. Optional flags execute
the current real-data validation runners on ready datasets.
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

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from spagapa.benchmark import (  # noqa: E402
    check_prepared_dataset,
    discover_prepared_datasets,
    statuses_to_dataframe,
    summarize_bib_readiness,
)


EXTERNAL_SCRIPT = PACKAGE_ROOT / "scripts" / "run_external_bioml_validation.py"
HIGHRES_DEFAULT_SCRIPT = PACKAGE_ROOT / "scripts" / "run_highres_default_validation.py"


def resolve_path(path: str | Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    cwd_path = (Path.cwd() / path).resolve()
    if cwd_path.exists():
        return cwd_path
    package_path = (PACKAGE_ROOT / path).resolve()
    if package_path.exists():
        return package_path
    if path.parts and path.parts[0] == PACKAGE_ROOT.name:
        return (PROJECT_ROOT / path).resolve()
    return package_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", default="data/processed")
    parser.add_argument("--dataset-dir", action="append", default=None)
    parser.add_argument("--output-dir", default="benchmark_results/real/real_data_suite_v1")
    parser.add_argument("--layer-column", default="layer")
    parser.add_argument("--run-external", action="store_true")
    parser.add_argument("--run-highres-default", action="store_true")
    parser.add_argument("--force", action="store_true")

    parser.add_argument("--n-genes", type=int, default=24)
    parser.add_argument("--min-observed-spots", type=int, default=80)
    parser.add_argument("--max-parent-spots", type=int, default=80)
    parser.add_argument("--seeds", default="42,43,44")
    parser.add_argument("--subbins-per-spot", default="2,4,8")
    parser.add_argument("--dropout-rate", type=float, default=0.25)
    parser.add_argument("--dropout-rates", default="0.10,0.25,0.40")
    parser.add_argument("--sparse-n-inducing", type=int, default=60)
    parser.add_argument("--bioml-n-neighbors", type=int, default=15)
    parser.add_argument("--bioml-rank", type=int, default=8)
    parser.add_argument("--bioml-max-iter", type=int, default=20)
    parser.add_argument(
        "--external-methods",
        default="raw,stapaminer_original_imputed,stapaminer_knn_expression,spagapa_gp,spagapa_bioml",
    )
    parser.add_argument("--skip-downstream", action="store_true")
    return parser.parse_args()


def collect_statuses(args: argparse.Namespace):
    if args.dataset_dir:
        return [
            check_prepared_dataset(resolve_path(path), layer_column=args.layer_column)
            for path in args.dataset_dir
        ]
    return discover_prepared_datasets(resolve_path(args.processed_root), layer_column=args.layer_column)


def run_command(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)


def maybe_run_external(statuses, args: argparse.Namespace, output_dir: Path) -> list[dict[str, Any]]:
    rows = []
    for status in statuses:
        if not status.external_validation_ready:
            continue
        dataset_out = output_dir / "external_validation" / status.dataset
        summary_file = dataset_out / "external_validation_summary.csv"
        if summary_file.exists() and not args.force:
            rows.append({"dataset": status.dataset, "track": "external_validation", "status": "skipped_existing"})
            continue
        cmd = [
            sys.executable,
            str(EXTERNAL_SCRIPT),
            "--dataset-dir",
            status.data_dir,
            "--output-dir",
            str(dataset_out),
            "--dataset-name",
            status.dataset,
            "--n-genes",
            str(args.n_genes),
            "--min-observed-spots",
            str(args.min_observed_spots),
            "--methods",
            args.external_methods,
        ]
        if args.skip_downstream:
            cmd.append("--skip-downstream")
        run_command(cmd, PACKAGE_ROOT)
        rows.append({"dataset": status.dataset, "track": "external_validation", "status": "completed"})
    return rows


def maybe_run_highres_default(statuses, args: argparse.Namespace, output_dir: Path) -> list[dict[str, Any]]:
    rows = []
    for status in statuses:
        if not status.highres_validation_ready:
            continue
        dataset_out = output_dir / "highres_default_validation" / status.dataset
        summary_file = dataset_out / "highres_default_validation_results_long.csv"
        if summary_file.exists() and not args.force:
            rows.append({"dataset": status.dataset, "track": "highres_default_validation", "status": "skipped_existing"})
            continue
        cmd = [
            sys.executable,
            str(HIGHRES_DEFAULT_SCRIPT),
            "--data-dir",
            status.data_dir,
            "--output-dir",
            str(dataset_out),
            "--dataset-name",
            status.dataset,
            "--n-genes",
            str(args.n_genes),
            "--min-observed-spots",
            str(args.min_observed_spots),
            "--max-parent-spots",
            str(args.max_parent_spots),
            "--seeds",
            args.seeds,
            "--subbins-per-spot",
            args.subbins_per_spot,
            "--dropout-rate",
            str(args.dropout_rate),
            "--dropout-rates",
            args.dropout_rates,
            "--sparse-n-inducing",
            str(args.sparse_n_inducing),
            "--bioml-n-neighbors",
            str(args.bioml_n_neighbors),
            "--bioml-rank",
            str(args.bioml_rank),
            "--bioml-max-iter",
            str(args.bioml_max_iter),
        ]
        if args.force:
            cmd.append("--force")
        run_command(cmd, PACKAGE_ROOT)
        rows.append({"dataset": status.dataset, "track": "highres_default_validation", "status": "completed"})
    return rows


def main() -> None:
    args = parse_args()
    output_dir = resolve_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    statuses = collect_statuses(args)
    readiness = statuses_to_dataframe(statuses)
    readiness.to_csv(output_dir / "real_data_readiness.csv", index=False)

    bib_summary = summarize_bib_readiness(statuses)
    (output_dir / "real_data_readiness_summary.json").write_text(
        json.dumps(bib_summary, indent=2, default=str)
    )

    run_rows: list[dict[str, Any]] = []
    if args.run_external:
        run_rows.extend(maybe_run_external(statuses, args, output_dir))
    if args.run_highres_default:
        run_rows.extend(maybe_run_highres_default(statuses, args, output_dir))
    run_manifest = pd.DataFrame(run_rows)
    run_manifest.to_csv(output_dir / "real_data_suite_run_manifest.csv", index=False)

    decision = {
        "output_dir": str(output_dir),
        "readiness_table": str(output_dir / "real_data_readiness.csv"),
        "readiness_summary": str(output_dir / "real_data_readiness_summary.json"),
        "run_manifest": str(output_dir / "real_data_suite_run_manifest.csv"),
        "bib_readiness": bib_summary,
        "ready_datasets": {
            "external_validation": [
                status.dataset for status in statuses if status.external_validation_ready
            ],
            "highres_validation": [
                status.dataset for status in statuses if status.highres_validation_ready
            ],
        },
    }
    (output_dir / "decision_summary.json").write_text(json.dumps(decision, indent=2, default=str))

    print("\nReal-data suite decision summary:")
    print(json.dumps(decision, indent=2, default=str))


if __name__ == "__main__":
    main()

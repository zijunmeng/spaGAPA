#!/usr/bin/env python3
"""Phase 3 Task 1 (B) rerun + merge driver.

Reruns the runtime/scalability benchmark at 42k and 100k for all four methods
using the FIXED harness (OPENBLAS_NUM_THREADS=8; highres_fast skips the
factorizer; highres_accuracy marks SKIPPED above 15k), then merges the new
high-scale rows with the existing 1k/5k/15k rows into a single
runtime_table.csv, regenerates scaling.png / mem_scaling.png, and writes
summary.json.

Run in the background; writes everything to disk so the parent agent can
verify via the filesystem.  Does NOT git commit.

Usage:
  python scripts/benchmark_runtime_rerun_and_merge.py [--wall-cap 1200]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

# ---- env FIRST (mirror the harness; OPENBLAS_NUM_THREADS=8 is load-bearing) -
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import pandas as pd

PKG = Path(__file__).resolve().parents[1]
OUT = PKG / "pipeline_output" / "benchmark_runtime"
HARNESS = PKG / "scripts" / "benchmark_runtime_scalability.py"
PY = os.path.expanduser("~/anaconda3/envs/spagapa/bin/python")

METHODS = [
    "spaGAPA-highres_fast",
    "spaGAPA-highres_accuracy",
    "stAPAminer",
    "spvAPA",
]


def run_harness_highscale(wall_cap: int) -> Path:
    """Run the harness for 42k+100k, all four methods, into a side output dir.

    Returns the path to the side runtime_table.csv.  We run into a side dir so
    the existing 1k/5k/15k rows + plots in OUT are untouched until the merge
    step (partial-failure isolation).
    """
    side = OUT / "_highscale_rerun"
    side.mkdir(parents=True, exist_ok=True)
    # reuse the already-generated 42k synthetic data + generate 100k by pointing
    # synthetic_data under the side dir via symlinks for 42k so we don't regen.
    side_data = side / "synthetic_data"
    side_data.mkdir(parents=True, exist_ok=True)
    src_42k = OUT / "synthetic_data" / "42000"
    dst_42k = side_data / "42000"
    if src_42k.exists() and not dst_42k.exists():
        dst_42k.symlink_to(src_42k)

    cmd = [
        PY, str(HARNESS),
        "--scales", "42000", "100000",
        "--methods", *METHODS,
        "--gene-cap", "500",
        "--wall-cap", str(wall_cap),
        "--output", str(side),
    ]
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "8"
    env["OMP_NUM_THREADS"] = "8"
    env["MKL_NUM_THREADS"] = "8"
    env["TMPDIR"] = "/s3/mengzijun/tmp"
    print(f"[rerun] launching harness: {' '.join(cmd)}", flush=True)
    t0 = time.time()
    # stream output to a log file AND inherit so Monitor/background can see it
    log_path = side / "_harness.log"
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.time() - t0
    print(f"[rerun] harness exited rc={proc.returncode} in {elapsed:.0f}s", flush=True)
    return side / "runtime_table.csv"


def merge_tables(high_csv: Path) -> pd.DataFrame:
    """Keep existing 1k/5k/15k rows; replace 42k/100k with the rerun rows."""
    main_csv = OUT / "runtime_table.csv"
    keep_scales = {1000, 5000, 15000}
    if main_csv.exists():
        main = pd.read_csv(main_csv)
        keep = main[main["n_spots"].isin(keep_scales)].copy()
    else:
        keep = pd.DataFrame(columns=[
            "method", "n_spots", "n_genes_used", "status",
            "time_s", "peak_mem_mb", "reason"])
    new_rows = []
    if high_csv.exists():
        new_rows.append(pd.read_csv(high_csv))
    # also keep any 42k/100k rows the rerun produced (handles partial reruns)
    table = pd.concat([keep] + new_rows, ignore_index=True)
    # dedupe: prefer the rerun (last occurrence) for a given (method, n_spots)
    table = table.drop_duplicates(subset=["method", "n_spots"], keep="last")
    # sort: method order first, then n_spots ascending
    method_order = {m: i for i, m in enumerate(METHODS)}
    table["_mo"] = table["method"].map(method_order)
    table = table.sort_values(["_mo", "n_spots"]).drop(columns=["_mo"])
    # reorder columns (drop the helper 'rc' col if present)
    cols = ["method", "n_spots", "n_genes_used", "status",
            "time_s", "peak_mem_mb", "reason"]
    table = table[[c for c in cols if c in table.columns]]
    table.to_csv(main_csv, index=False)
    print(f"[merge] wrote {main_csv} ({len(table)} rows)", flush=True)
    return table


def regen_plots(table: pd.DataFrame) -> None:
    """Regenerate scaling.png + mem_scaling.png from the merged table."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("brs", HARNESS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.make_plots(table, OUT)
    print(f"[plots] regenerated scaling.png + mem_scaling.png in {OUT}",
          flush=True)


def write_summary(table: pd.DataFrame) -> None:
    """Per-method: max scale COMPLETED + time/mem at 100k or failure reason."""
    summary = {}
    for m in METHODS:
        sub = table[table["method"] == m].sort_values("n_spots")
        completed = sub[sub["status"] == "COMPLETED"]
        max_done = (int(completed["n_spots"].max())
                    if len(completed) else None)
        entry = {"max_scale_completed": max_done}
        # 100k row (or largest scale attempted)
        for n in (100000, 42000):
            row = sub[sub["n_spots"] == n]
            if len(row):
                r = row.iloc[0]
                entry[f"at_{n}"] = {
                    "status": r["status"],
                    "time_s": float(r["time_s"]) if r["time_s"] else None,
                    "peak_mem_mb": (float(r["peak_mem_mb"])
                                    if r["peak_mem_mb"] else None),
                    "reason": (str(r["reason"]) if isinstance(r["reason"], str)
                               and r["reason"] else None),
                }
                break
        summary[m] = entry
    out = OUT / "summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"[summary] wrote {out}", flush=True)
    print(json.dumps(summary, indent=2), flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wall-cap", type=int, default=1200)
    ap.add_argument("--skip-rerun", action="store_true",
                    help="only merge existing tables + regen plots/summary")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    if not args.skip_rerun:
        high_csv = run_harness_highscale(args.wall_cap)
    else:
        high_csv = OUT / "_highscale_rerun" / "runtime_table.csv"

    table = merge_tables(high_csv)
    try:
        regen_plots(table)
    except Exception as exc:
        print(f"[plots] SKIPPED: {exc!r}", flush=True)
    write_summary(table)
    print("[rerun+merge] DONE", flush=True)


if __name__ == "__main__":
    main()

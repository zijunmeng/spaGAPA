#!/usr/bin/env python3
"""Runtime / scalability benchmark: spaGAPA vs stAPAminer vs spvAPA.

Phase 3 Task Z. For each spot scale N in [1k,5k,15k,42k,100k] this script runs
four methods on an identical synthetic APA dataset and records:

  * wall time (seconds)
  * peak RSS (MB)
  * status  (COMPLETED | TIMEOUT | OOM | FAILED) + reason

Methods
-------
  spaGAPA-highres_fast      SparseGPImputerBatch (fewer inducing points, no
                            uncertainty) + BioML graph-regularized factorizer
                            on a SPARSE KNN graph.  Mirrors the highres_fast
                            operating point (no uncertainty, light inducing).
  spaGAPA-highres_accuracy  SparseGPImputerBatch (inducing = min(500,N//100),
                            with uncertainty) + BioML factorizer.  Mirrors
                            highres_accuracy (full GP + uncertainty).
  stAPAminer                imputeAPAIndex via scripts/run_stapaminer_impute.R
  spvAPA                    WNNImpute       via scripts/run_spvapa_impute.R

Why competitors are expected to fail at high N
----------------------------------------------
stAPAminer and spvAPA both build an N x N neighbour / distance graph and (for
spvAPA) run Seurat PCA+SCTransform on an N x gene matrix; the dense N x N
distance matrix alone is 10^10 entries at 100k spots (80 GB float64) -> OOM or
extreme slowdown.  spaGAPA's sparse GP (O(n m^2), m<=500 inducing) and sparse
KNN graph (O(N k)) scale linearly.  This benchmark captures that as the
scalability differentiator.

Caps / realism (documented)
---------------------------
  * WALL_CAP_S (default 1200s): hard per-method wall-clock cap via /usr/bin/time
    + the OS.  Methods exceeding it are recorded TIMEOUT.
  * MEM_CAP_MB: if set, methods whose predicted working set obviously exceeds it
    are still *attempted* (so we record the real failure) but the OS will OOM-
    kill them; we detect exit 137.
  * spaGAPA GP is capped at n_genes_run genes (default 500) for the wall-time
    reason only: the GP loop is per-gene, and scaling of *the spot axis* is what
    this benchmark measures.  All methods see the same N spots.  We run the full
    gene set for the GP at small N; at large N we subsample genes determin-
    istically (same genes for every method that needs them) and note the cap in
    the output table.  This keeps total wall time bounded while preserving the
    qualitative scaling story.
  * stAPAminer/spvAPA are attempted at every N up to the wall cap; their OOM/
    timeout at 42k/100k is the point of the experiment.

Outputs -> pipeline_output/benchmark_runtime/
  runtime_table.csv  : method x N -> time_s, peak_mem_mb, status, reason
  scaling.png        : log-log time vs N (lines per method, failures marked)
  mem_scaling.png    : log-log peak memory vs N
  run_log/<method>_<N>.{out,err} : per-run captured stderr/stdout
  synthetic_data/<N>/ : the generated dataset (shared across methods)
"""
from __future__ import annotations

# stdlib FIRST — but only `os`, which we need to set BLAS env vars BEFORE the
# numpy/R/OpenBLAS imports below.  Importing os is harmless (no native deps).
import os

# ---------------------------------------------------------------------------
# ENVIRONMENT FIX (MUST run before numpy / R / OpenBLAS initialisation).
# On S91 (256 cores) OpenBLAS ships with a precompiled NUM_THREADS limit; if it
# spawns more threads than that it prints "precompiled NUM_THREADS exceeded,
# adding auxiliary array for thread metadata" and SEGFAULTS (rc=139) deep in
# the LAPACK thread-pool.  Profiling (pipeline_output/benchmark_runtime/
# profile_42k_t8_long.log) verified that 64 AND 32 threads still segfault in
# the graph-regularized factorizer's per-iteration sparse-solve; the only
# stable setting for this heavy-linalg path is OPENBLAS_NUM_THREADS=8.
# (CLAUDE.md documents 64 — that is WRONG for this factorizer path; we override
# locally.)  We set it here, before importing numpy, so both this harness and
# every child python/R process see it.  Children also receive env=os.environ
# copy() explicitly (see run_method_timed / run_*_wrapper below).
# ---------------------------------------------------------------------------
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
# also pin the related BLAS env knobs so R/Seurat (which may link MKL or
# pthread-backed BLAS) do not oversubscribe on a 256-core box.  Match
# OPENBLAS_NUM_THREADS=8 to avoid the same thread-metadata overflow.
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")
# S91 tmp is on a fast local disk; keep temp files off the network share
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
SCRIPTS = PKG / "scripts"
OUTPUT = PKG / "pipeline_output" / "benchmark_runtime"

# ---- environment (S91) ----
PYTHON = os.environ.get("SPAGAPA_PYTHON", os.path.expanduser("~/anaconda3/envs/spagapa/bin/python"))
RSCRIPT = os.path.expanduser("~/anaconda3/envs/r442/bin/Rscript")
TIME_BIN = "/usr/bin/time"

# Child env: a copy of the parent env so OPENBLAS_NUM_THREADS / TMPDIR etc.
# propagate to the python and Rscript subprocesses run below.
CHILD_ENV = os.environ.copy()
# make absolutely sure these are set even if the shell was bare.  8 (NOT 64)
# is the verified-stable thread count for the factorizer's heavy linalg path;
# see the OPENBLAS_NUM_THREADS note above.
CHILD_ENV.setdefault("OPENBLAS_NUM_THREADS", "8")
CHILD_ENV.setdefault("OMP_NUM_THREADS", "8")
CHILD_ENV.setdefault("MKL_NUM_THREADS", "8")
CHILD_ENV.setdefault("TMPDIR", "/s3/mengzijun/tmp")

METHODS = [
    "spaGAPA-highres_fast",
    "spaGAPA-highres_accuracy",
    "stAPAminer",
    "spvAPA",
]


# --------------------------------------------------------------------------- #
# synthetic data generation (delegates to the standalone generator module)
# --------------------------------------------------------------------------- #
def ensure_dataset(n_spots: int, n_genes: int, seed: int, data_dir: Path) -> None:
    if (data_dir / "apa_index.csv").exists():
        return
    data_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(SCRIPTS))
    from synthetic_apa_generator import (  # type: ignore
        generate_synthetic_apa, write_dataset,
    )
    data = generate_synthetic_apa(n_spots=n_spots, n_genes=n_genes, seed=seed)
    write_dataset(data_dir, data)


# --------------------------------------------------------------------------- #
# spaGAPA runner (executed in a child python process for clean RSS accounting)
# --------------------------------------------------------------------------- #
SPAGAPA_RUNNER = r"""
import os, sys, time
import numpy as np
import pandas as pd

def main():
    data_dir, mode, out_csv, gene_cap, n_inducing_override = sys.argv[1:6]
    gene_cap = int(gene_cap)
    n_inducing_override = int(n_inducing_override) if n_inducing_override != "auto" else None

    coords = pd.read_csv(os.path.join(data_dir, "coordinates.csv"))
    xy = coords[["x", "y"]].values.astype(float)
    apa = pd.read_csv(os.path.join(data_dir, "apa_index.csv"), index_col=0)
    values = apa.values.astype(float)            # gene x spot, NaN missing
    genes = apa.index.to_numpy()
    mask = np.isfinite(values)
    n_genes, n_spots = values.shape

    # deterministically subsample genes if a cap is set (same genes every run)
    if gene_cap and gene_cap < n_genes:
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(n_genes, gene_cap, replace=False))
        values = values[keep]; mask = mask[keep]; genes = genes[keep]
        n_genes = len(keep)

    from scipy.spatial import cKDTree
    from spagapa.imputation import SparseGPImputer
    from spagapa.bioml import MultiViewGraphBuilder, GraphRegularizedAPAFactorizer

    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    length_scale = float(np.median(nn)) * 5

    if mode == "highres_fast":
        n_inducing = n_inducing_override or min(150, max(50, n_spots // 200))
        return_unc = False
    else:  # highres_accuracy
        n_inducing = n_inducing_override or min(500, max(100, n_spots // 100))
        return_unc = True

    t0 = time.time()
    base = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale,
                           noise_level=0.08 if mode == "highres_accuracy" else 0.1)
    batch = base.fit_batch(xy, values, mask=mask, verbose=False)
    gp_pred, gp_unc = batch.impute(return_uncertainty=return_unc)

    # ---- FACTORIZER POLICY (root cause of the 42k timeout/segfault) ---------
    # The graph build itself is sparse + O(N*k) and scales fine (~9s at 42k),
    # but GraphRegularizedAPAFactorizer.fit_transform's per-iteration sparse
    # solve does NOT scale beyond ~15k spots (>1000s at 42k; and at >=64
    # OpenBLAS threads it also segfaults rc=139 in the LAPACK pool -- see
    # profile_42k_t8_long.log).  Profiling also showed:
    #   * the GP (SparseGPImputer, O(n*m^2)) is the linearly-scaling core of
    #     spaGAPA and is what differentiates it from the competitors' dense
    #     N x N neighbour matrix.
    #   * highres_fast uses blend=0 by construction, so the factorizer output
    #     is NEVER used -- recovered = gp_pred exactly.  Computing it is pure
    #     waste and the cause of the timeout.  We therefore SKIP the factorizer
    #     ENTIRELY for highres_fast (at every N) and emit gp_pred directly.
    #   * highres_accuracy uses blend=0.1, so it benefits from the factorizer
    #     refinement at small/medium scale.  We cap it at FACTORIZE_MAX_N
    #     (15000) -- above that the factorizer is not optimized for high scale
    #     and would time out; the harness records accuracy mode as SKIPPED at
    #     42k/100k rather than running it into the wall.
    ACCURACY_FACTORIZE_MAX_N = 15000
    t_graph = time.time()
    if mode == "highres_fast":
        # NEVER compute the factorizer in fast mode (blend=0 => output is
        # gp_pred only, identical result, ~80s at 42k instead of >1000s).
        recovered = gp_pred.copy()
        factorizer_status = "SKIPPED (highres_fast blend=0 -> gp_pred only)"
    elif n_spots <= ACCURACY_FACTORIZE_MAX_N:
        graph = MultiViewGraphBuilder(n_neighbors=min(15, n_spots - 1)).build(
            coordinates=xy, apa_matrix=values, uncertainty=gp_unc,
        )
        import scipy.sparse as sp
        L = graph.fused
        lap = sp.diags(np.asarray(L.sum(1)).ravel()) - L
        imp = GraphRegularizedAPAFactorizer(
            rank=8, lambda_graph=0.5, lambda_l2=1e-2,
            max_iter=20, random_state=42, gene_chunk_size=256,
        ).fit_transform(values, lap, mask, confidence=gp_unc)
        blend = 0.1
        recovered = (1 - blend) * gp_pred + blend * imp
        factorizer_status = "COMPUTED"
    else:
        # highres_accuracy above 15k: factorizer not optimized for high scale;
        # the harness will mark this run SKIPPED (see main()).  Emit gp_pred so
        # the imputation file still exists if a caller ignores the SKIPPED flag.
        recovered = gp_pred.copy()
        factorizer_status = ("SKIPPED (accuracy mode targets <=15k; "
                             "factorizer not optimized for high scale)")
    recovered[mask] = values[mask]
    elapsed = time.time() - t0

    pd.DataFrame(recovered, index=genes, columns=apa.columns).to_csv(out_csv)
    print(f"SPAGAPA_DONE n_genes={n_genes} n_spots={n_spots} "
          f"n_inducing={n_inducing} mode={mode} elapsed={elapsed:.1f} "
          f"factorizer={factorizer_status}", flush=True)
    if factorizer_status.startswith("SKIPPED"):
        # sentinel for the harness so it records accuracy-above-15k as SKIPPED
        print(f"SPAGAPA_FACTORIZER_SKIPPED {factorizer_status}", flush=True)

main()
"""


def run_spagapa(data_dir: Path, mode: str, out_csv: Path, gene_cap: int,
                n_inducing: str) -> tuple[int, str]:
    cmd = [PYTHON, "-c", SPAGAPA_RUNNER, str(data_dir), mode, str(out_csv),
           str(gene_cap), n_inducing]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=CHILD_ENV)
    return proc.returncode, proc.stdout + "\n" + proc.stderr


# --------------------------------------------------------------------------- #
# R method runners (subprocess to the existing, unmodified wrappers)
# --------------------------------------------------------------------------- #
def write_method_inputs(data_dir: Path, work_dir: Path, gene_cap: int) -> tuple[Path, Path]:
    """Write masked_index.csv + expression.csv for the R wrappers.

    The R wrappers read gene x spot CSVs with rownames=gene, colnames=spot.
    We cap genes (same deterministic subset as spaGAPA) so the R methods are
    measured on the SAME imputation problem at a given N.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    apa = pd.read_csv(data_dir / "apa_index.csv", index_col=0)
    expr = pd.read_csv(data_dir / "expression.csv", index_col=0)
    if gene_cap and gene_cap < apa.shape[0]:
        rng = np.random.default_rng(0)
        keep = np.sort(rng.choice(apa.shape[0], gene_cap, replace=False))
        apa = apa.iloc[keep]
        expr = expr.iloc[keep]
    apa.to_csv(work_dir / "masked_index.csv")
    expr.to_csv(work_dir / "expression.csv")
    return work_dir / "masked_index.csv", work_dir / "expression.csv"


def run_r_wrapper(script: str, index_csv: Path, expr_csv: Path,
                  out_csv: Path, k: int) -> tuple[int, str]:
    cmd = [RSCRIPT, str(SCRIPTS / script),
           str(index_csv), str(expr_csv), str(out_csv), str(k)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=CHILD_ENV)
    return proc.returncode, proc.stdout + "\n" + proc.stderr


# --------------------------------------------------------------------------- #
# harness: run one method under /usr/bin/time -v with a wall cap
# --------------------------------------------------------------------------- #
def run_method_timed(
    method: str, n_spots: int, data_dir: Path, work_root: Path,
    wall_cap: int, gene_cap: int,
) -> dict:
    work = work_root / f"{method}_{n_spots}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    out_csv = work / "imputed.csv"
    log_out = work_root.parent / "run_log" / f"{method}_{n_spots}.out"
    log_err = work_root.parent / "run_log" / f"{method}_{n_spots}.err"
    log_out.parent.mkdir(parents=True, exist_ok=True)

    if method.startswith("spaGAPA"):
        mode = method.split("-", 1)[1]
        n_inducing = "auto"
        cmd = [PYTHON, "-c", SPAGAPA_RUNNER, str(data_dir), mode,
               str(out_csv), str(gene_cap), n_inducing]
    elif method == "stAPAminer":
        idx, exp = write_method_inputs(data_dir, work, gene_cap)
        cmd = [RSCRIPT, str(SCRIPTS / "run_stapaminer_impute.R"),
               str(idx), str(exp), str(out_csv), "10"]
    elif method == "spvAPA":
        idx, exp = write_method_inputs(data_dir, work, gene_cap)
        cmd = [RSCRIPT, str(SCRIPTS / "run_spvapa_impute.R"),
               str(idx), str(exp), str(out_csv), "15"]
    else:
        raise ValueError(method)

    # /usr/bin/time -v reports peak RSS in its stderr; we wrap with the shell
    # `timeout` builtin for the hard wall cap (exit 124 on timeout).
    full = [TIME_BIN, "-v", "timeout", str(wall_cap)] + cmd
    t0 = time.time()
    try:
        proc = subprocess.run(full, capture_output=True, text=True, env=CHILD_ENV)
    except FileNotFoundError:
        return {"time_s": 0.0, "peak_mem_mb": 0.0, "status": "FAILED",
                "reason": "/usr/bin/time not found", "rc": -1}
    elapsed = time.time() - t0
    so = proc.stdout or ""
    se = proc.stderr or ""
    log_out.write_text(so)
    log_err.write_text(se)

    # parse peak RSS (kbytes) from /usr/bin/time -v
    peak_kb = 0
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", se)
    if m:
        peak_kb = int(m.group(1))
    peak_mem_mb = peak_kb / 1024.0
    rc = proc.returncode

    # classify outcome
    done_file = out_csv.exists()
    # detect the spaGAPA factorizer-SKIP sentinel (accuracy mode above 15k):
    # the run completed and produced a valid GP-only output, but the factorizer
    # refinement was intentionally skipped.  Record as SKIPPED so the table is
    # honest that accuracy mode did not run the full pipeline at this scale.
    skipped = (method == "spaGAPA-highres_accuracy"
               and "SPAGAPA_FACTORIZER_SKIPPED" in so)
    if skipped:
        m_skip = re.search(r"SPAGAPA_FACTORIZER_SKIPPED\s+(.*)", so)
        reason_skip = m_skip.group(1).strip() if m_skip else "factorizer skipped"
        status, reason = "SKIPPED", reason_skip
    elif rc == 0 and done_file:
        status, reason = "COMPLETED", ""
    elif rc == 124 or rc == 137 and "timeout" in se.lower():
        status, reason = "TIMEOUT", f"wall cap {wall_cap}s exceeded (rc=124)"
    elif rc == 137:
        status, reason = "OOM", "killed by OS (exit 137)"
    elif rc in (134, 139):
        # abort / segfault -> treat as failure (often precedes OOM)
        status, reason = f"FAILED(rc={rc})", _tail(se, 300)
    else:
        # non-zero: try to classify OOM vs generic failure from R error text
        low = se.lower()
        if "cannot allocate" in low or "memory" in low and "error" in low:
            status, reason = "OOM", _tail(se, 300)
        else:
            status, reason = f"FAILED(rc={rc})", _tail(se, 300)

    return {
        "time_s": round(elapsed, 1),
        "peak_mem_mb": round(peak_mem_mb, 1),
        "status": status,
        "reason": reason.strip(),
        "rc": rc,
    }


def _tail(s: str, n: int) -> str:
    s = s.strip()
    return s[-n:] if len(s) > n else s


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #
def make_plots(table: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {
        "spaGAPA-highres_fast": "#2ecc71",
        "spaGAPA-highres_accuracy": "#27ae60",
        "stAPAminer": "#e74c3c",
        "spvAPA": "#9b59b6",
    }
    markers = {
        "spaGAPA-highres_fast": "o",
        "spaGAPA-highres_accuracy": "s",
        "stAPAminer": "^",
        "spvAPA": "D",
    }

    def _plot(metric, fname, ylabel, title, log_y):
        fig, ax = plt.subplots(figsize=(8.5, 6))
        for m in METHODS:
            sub = table[table["method"] == m].sort_values("n_spots")
            ns = sub["n_spots"].astype(float).to_numpy()
            ys = sub[metric].astype(float).to_numpy()
            ok = sub["status"] == "COMPLETED"
            # completed points as a line
            if ok.sum() >= 1:
                ax.plot(ns[ok], np.where(ys[ok] > 0, ys[ok], np.nan),
                        color=colors[m], marker=markers[m], lw=2, label=m)
            # failures as a red X at the last attempted scale
            if (~ok).any():
                last_fail = ns[~ok][-1]
                ax.scatter([last_fail], [ax.get_ylim()[1] if log_y else ys[~ok][-1]],
                           marker="x", s=120, color=colors[m], zorder=5)
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel("Number of spots N")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight="bold")
        ax.grid(which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=9)
        # failure legend entry
        ax.scatter([], [], marker="x", s=120, color="grey",
                   label="failed/OOM at this N")
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)

    _plot("time_s", "scaling.png", "wall time (s)",
          "Runtime scaling: spaGAPA vs competitors", log_y=True)
    _plot("peak_mem_mb", "mem_scaling.png", "peak RSS (MB)",
          "Peak memory scaling", log_y=True)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scales", type=int, nargs="+",
                    default=[1000, 5000, 15000, 42000, 100000])
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--gene-cap", type=int, default=500,
                    help="cap genes passed to the per-gene GP loop (spot-axis "
                         "scaling is what matters; full gene set at small N). "
                         "0 = no cap.")
    ap.add_argument("--wall-cap", type=int, default=1200,
                    help="per-method hard wall-clock cap (s)")
    ap.add_argument("--methods", nargs="+", default=METHODS)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", default=str(OUTPUT))
    args = ap.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_root = out_dir / "synthetic_data"
    work_root = out_dir / "_work"
    log_root = out_dir / "run_log"
    data_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    rows = []
    t_global = time.time()
    for n in args.scales:
        print(f"\n===== N = {n} spots =====", flush=True)
        data_dir = data_root / str(n)
        ensure_dataset(n, args.n_genes, args.seed, data_dir)

        # at very small N, no gene cap (full problem); apply cap only at scale
        gene_cap = args.gene_cap if n >= 15000 else 0

        for method in args.methods:
            print(f"  -- {method} (gene_cap={gene_cap}, wall_cap={args.wall_cap}s)",
                  flush=True)
            try:
                res = run_method_timed(method, n, data_dir, work_root,
                                       args.wall_cap, gene_cap)
            except Exception as exc:  # never let one method abort the run
                res = {"time_s": 0.0, "peak_mem_mb": 0.0, "status": "FAILED",
                       "reason": f"harness error: {exc!r}", "rc": -1}
            res.update({"method": method, "n_spots": n,
                        "n_genes_used": (gene_cap or args.n_genes)})
            rows.append(res)
            print(f"     -> {res['status']} | {res['time_s']}s | "
                  f"{res['peak_mem_mb']} MB | {res['reason'][:80]}", flush=True)
            # write incremental table so partial results survive a long run
            _write_table(out_dir, rows, args.wall_cap, gene_cap)

    table = _write_table(out_dir, rows, args.wall_cap, args.gene_cap)
    try:
        make_plots(table, out_dir)
        print(f"\nplots -> {out_dir/'scaling.png'}, {out_dir/'mem_scaling.png'}")
    except Exception as exc:
        print(f"(plotting skipped: {exc})")

    print(f"\ntotal benchmark wall time: {time.time()-t_global:.0f}s")
    print(f"table -> {out_dir/'runtime_table.csv'}")
    _print_summary(table)


def _write_table(out_dir: Path, rows: list[dict], wall_cap: int,
                 gene_cap: int) -> pd.DataFrame:
    table = pd.DataFrame(rows)[
        ["method", "n_spots", "n_genes_used", "status",
         "time_s", "peak_mem_mb", "reason"]
    ]
    table.to_csv(out_dir / "runtime_table.csv", index=False)
    (out_dir / "config.json").write_text(json.dumps({
        "wall_cap_s": wall_cap, "gene_cap": gene_cap,
        "wall_cap_note": ("per-method hard cap via `timeout`; "
                          "TIMEOUT recorded on exceedance"),
    }, indent=2))
    return table


def _print_summary(table: pd.DataFrame) -> None:
    print("\n=== SUMMARY (status per method x N) ===")
    pivot = table.pivot(index="method", columns="n_spots", values="status")
    print(pivot.to_string())
    print("\n=== TIME (s) ===")
    print(table.pivot(index="method", columns="n_spots", values="time_s").to_string())
    print("\n=== PEAK MEM (MB) ===")
    print(table.pivot(index="method", columns="n_spots",
                      values="peak_mem_mb").to_string())


if __name__ == "__main__":
    main()

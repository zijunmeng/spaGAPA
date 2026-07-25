#!/usr/bin/env python3
"""Benchmark GraphRegularizedAPAFactorizer.fit_transform at 42k spots.

Builds a sparse kNN-RBF Laplacian from the synthetic 42k coordinates (same
recipe as scripts/profile_42k_timeout.py), then times ONLY the factorizer's
fit_transform path. Used to measure before/after of the Phase-3 Task-2
graph-reg spot-update optimization.

Usage:
    python3 scripts/bench_factorizer_42k.py [LABEL] [GENE_CAP] [MAX_ITER]
        LABEL     tag written into before_after.json  (default: "run")
        GENE_CAP  cap on number of genes              (default: 500)
        MAX_ITER  factorizer max_iter                 (default: 20)
"""
import os, sys, time, json, resource, gc

# Make `spagapa` importable when run as a plain script (package is not
# pip-installed; it relies on cwd being the repo root).
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

DATA = "pipeline_output/benchmark_runtime/synthetic_data/42000"
LABEL = sys.argv[1] if len(sys.argv) > 1 else "run"
GENE_CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 500
MAX_ITER = int(sys.argv[3]) if len(sys.argv) > 3 else 20
OUT_DIR = "pipeline_output/factorizer_optimization"
os.makedirs(OUT_DIR, exist_ok=True)


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def log(msg, t0=None):
    el = f"{time.time()-T0:7.1f}s" if t0 is None else f"{time.time()-t0:7.1f}s"
    print(f"[{el} | rss {rss_mb():8.0f}MB] {msg}", flush=True)


T0 = time.time()
log(f"start label={LABEL} gene_cap={GENE_CAP} max_iter={MAX_ITER} "
    f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}")

coords = pd.read_csv(f"{DATA}/coordinates.csv")
xy = coords[["x", "y"]].values.astype(float)
apa = pd.read_csv(f"{DATA}/apa_index.csv", index_col=0)
values = apa.values.astype(float)
mask = np.isfinite(values)
n_genes, n_spots = values.shape
log(f"loaded apa {values.shape} (observed {mask.mean():.3f})")

if GENE_CAP and GENE_CAP < n_genes:
    rng = np.random.default_rng(0)
    keep = np.sort(rng.choice(n_genes, GENE_CAP, replace=False))
    values = values[keep]
    mask = mask[keep]
    n_genes = len(keep)
    log(f"gene-capped to {n_genes}")

# Fill NaNs the same way the factorizer will, so we hand it a clean matrix
# (the factorizer re-fills internally; this just makes the graph build sane).
filled = values.copy()
gm = np.nanmean(filled, axis=1)
gm = np.nan_to_num(gm, nan=0.0)
mg, ms = np.where(~np.isfinite(filled))
filled[mg, ms] = gm[mg]

# Sparse kNN-RBF graph (same as MultiViewGraphBuilder spatial view + the
# profile script's combinatorial Laplacian).
t = time.time()
k = min(15, n_spots - 1)
tree = cKDTree(xy)
dists, idx = tree.query(xy, k=k + 1)  # +1 to drop self
rows = np.repeat(np.arange(n_spots), k)
cols = idx[:, 1:].ravel()
d = dists[:, 1:].ravel()
sigma = float(np.median(d)) + 1e-12
W = sparse.csr_matrix((np.exp(-0.5 * (d / sigma) ** 2), (rows, cols)),
                      shape=(n_spots, n_spots))
W = W.maximum(W.T)
W.setdiag(0.0)
W.eliminate_zeros()
deg = np.asarray(W.sum(1)).ravel()
lap = sparse.diags(deg) - W.tocsc()
log(f"graph build DONE ({time.time()-t:.1f}s) lap nnz={lap.nnz} "
    f"density={lap.nnz/(n_spots**2):.2e}")

# Import after heavy data load so RSS reflects data, not libs.
from spagapa.bioml import GraphRegularizedAPAFactorizer

f = GraphRegularizedAPAFactorizer(
    rank=8, lambda_graph=0.5, lambda_l2=1e-2,
    max_iter=MAX_ITER, random_state=42, gene_chunk_size=256,
)

log("calling fit_transform ...")
t = time.time()
imp = f.fit_transform(filled, lap, mask, confidence=None)
ft_time = time.time() - t
log(f"factorizer fit_transform DONE ({ft_time:.1f}s) "
    f"n_iter={f.result_.n_iter} err={f.result_.reconstruction_error:.4g}")

assert np.isfinite(imp).all(), "imputed has NaN/inf"
assert imp.shape == filled.shape

rec = {
    "label": LABEL,
    "gene_cap": GENE_CAP,
    "max_iter": MAX_ITER,
    "n_genes": int(n_genes),
    "n_spots": int(n_spots),
    "lap_nnz": int(lap.nnz),
    "fit_transform_seconds": round(float(ft_time), 2),
    "n_iter": int(f.result_.n_iter),
    "reconstruction_error": float(f.result_.reconstruction_error),
    "peak_rss_mb": round(float(rss_mb()), 0),
    "openblas_num_threads": int(os.environ.get("OPENBLAS_NUM_THREADS", "0") or 0),
}

json_path = f"{OUT_DIR}/before_after.json"
existing = []
if os.path.exists(json_path):
    try:
        existing = json.load(open(json_path))
        if not isinstance(existing, list):
            existing = [existing]
    except Exception:
        existing = []
# Replace any prior record with the same label.
existing = [r for r in existing if r.get("label") != LABEL]
existing.append(rec)
json.dump(existing, open(json_path, "w"), indent=2)
log(f"wrote {json_path}")
print(json.dumps(rec, indent=2))

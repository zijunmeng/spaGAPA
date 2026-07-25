#!/usr/bin/env python3
"""Microbench: smoother per-rank-column loop vs single multi-RHS solve,
on the REAL 42k kNN Laplacian (structured, factors cleanly)."""
import os, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
_REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.sparse.linalg import factorized

coords = pd.read_csv(f"{_REPO}/pipeline_output/benchmark_runtime/synthetic_data/42000/coordinates.csv")
xy = coords[["x", "y"]].values.astype(float)
n = xy.shape[0]
k = 15
tree = cKDTree(xy)
dists, idx = tree.query(xy, k=k + 1)
rows = np.repeat(np.arange(n), k)
cols = idx[:, 1:].ravel()
d = dists[:, 1:].ravel()
sigma = float(np.median(d)) + 1e-12
W = sparse.csr_matrix((np.exp(-0.5 * (d / sigma) ** 2), (rows, cols)), shape=(n, n))
W = W.maximum(W.T); W.setdiag(0.0); W.eliminate_zeros()
deg = np.asarray(W.sum(1)).ravel()
lap = sparse.diags(deg) - W.tocsc()
system = sparse.eye(n, format="csc") + 0.5 * lap.tocsc()
print(f"n={n} lap nnz={lap.nnz}", flush=True)

t = time.time(); f = factorized(system); t_fac = time.time() - t
print(f"LU factorization: {t_fac:.2f}s", flush=True)

rng = np.random.default_rng(0)
rank = 8
z_raw = rng.random((n, rank))

# Warm up
_ = f(z_raw[:, 0])

# Current: per-column loop + column_stack
Nrep = 20
t = time.time()
for _ in range(Nrep):
    z1 = np.column_stack([f(z_raw[:, d]) for d in range(rank)])
t_loop = (time.time() - t) / Nrep
print(f"per-column loop ({rank} solves) + column_stack: {t_loop*1000:.2f}ms/iter", flush=True)

# Optimized: single multi-RHS solve
t = time.time()
for _ in range(Nrep):
    z2 = f(z_raw)
t_multi = (time.time() - t) / Nrep
print(f"single multi-RHS solve: {t_multi*1000:.2f}ms/iter", flush=True)
print(f"equivalent: {np.allclose(z1, z2)}  per-iter speedup: {t_loop/t_multi:.1f}x")
print(f"over 20 iters: loop={t_loop*20*1000:.0f}ms  multi={t_multi*20*1000:.0f}ms  saves {(t_loop-t_multi)*20*1000:.0f}ms")

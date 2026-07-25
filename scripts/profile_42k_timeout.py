#!/usr/bin/env python3
"""Profile spaGAPA-highres_fast at 42k spots to find the timeout/memory hotspot.

Reuses the runtime benchmark's synthetic 42k data. Times + RSS-tracks each step:
load -> GP fit_batch -> GP impute -> graph build -> factorizer fit_transform.
"""
import os, sys, time, resource, gc
os.environ.setdefault("OPENBLAS_NUM_THREADS", "64")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")
import numpy as np
import pandas as pd

DATA = sys.argv[1] if len(sys.argv) > 1 else "pipeline_output/benchmark_runtime/synthetic_data/42000"
GENE_CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 500

def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

def stamp(label, t0):
    print(f"[{time.time()-T0:7.1f}s | rss {rss_mb():8.0f}MB] {label}", flush=True)

T0 = time.time()
stamp(f"start (gene_cap={GENE_CAP}, data={DATA})", 0)

coords = pd.read_csv(f"{DATA}/coordinates.csv")
xy = coords[["x", "y"]].values.astype(float)
apa = pd.read_csv(f"{DATA}/apa_index.csv", index_col=0)
values = apa.values.astype(float)
mask = np.isfinite(values)
n_genes, n_spots = values.shape
stamp(f"loaded apa {values.shape} (observed {mask.mean():.3f})", 0)

if GENE_CAP and GENE_CAP < n_genes:
    rng = np.random.default_rng(0)
    keep = np.sort(rng.choice(n_genes, GENE_CAP, replace=False))
    values = values[keep]; mask = mask[keep]; n_genes = len(keep)
    stamp(f"gene-capped to {n_genes}", 0)

from scipy.spatial import cKDTree
from spagapa.imputation import SparseGPImputer
from spagapa.bioml import MultiViewGraphBuilder, GraphRegularizedAPAFactorizer

nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
length_scale = float(np.median(nn)) * 5
n_inducing = min(150, max(50, n_spots // 200))
stamp(f"GP setup: n_inducing={n_inducing} length_scale={length_scale:.1f} n_spots={n_spots}", 0)

t = time.time()
base = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale, noise_level=0.1)
batch = base.fit_batch(xy, values, mask=mask, verbose=False)
stamp(f"GP fit_batch DONE ({time.time()-t:.1f}s)", 0)

t = time.time()
gp_pred = batch.impute(return_uncertainty=False)
stamp(f"GP impute DONE ({time.time()-t:.1f})", 0)
del batch; gc.collect()

t = time.time()
graph = MultiViewGraphBuilder(n_neighbors=min(15, n_spots - 1)).build(
    coordinates=xy, apa_matrix=values, uncertainty=None)
import scipy.sparse as sp
L = graph.fused
stamp(f"graph build DONE ({time.time()-t:.1f}s) fused nnz={L.nnz} density={L.nnz/(n_spots**2):.2e}", 0)

t = time.time()
lap = sp.diags(np.asarray(L.sum(1)).ravel()) - L
stamp(f"laplacian DONE ({time.time()-t:.1f}s) lap nnz={lap.nnz}", 0)

t = time.time()
imp = GraphRegularizedAPAFactorizer(
    rank=8, lambda_graph=0.5, lambda_l2=1e-2,
    max_iter=20, random_state=42, gene_chunk_size=256,
).fit_transform(values, lap, mask, confidence=None)
stamp(f"factorizer fit_transform DONE ({time.time()-t:.1f}s)", 0)

stamp(f"ALL DONE — total {time.time()-T0:.1f}s peak rss {rss_mb():.0f}MB", 0)

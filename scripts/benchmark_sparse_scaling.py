#!/usr/bin/env python3
"""Scaling benchmark: run factorizer at 3 spot scales, record runtime/memory."""
import time, json, os, resource
import numpy as np
from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer

def benchmark_factorizer(n_genes, n_spots, rank=8, max_iter=10):
    rng = np.random.default_rng(42)
    apa = rng.random((n_genes, n_spots))
    apa[apa < 0.9] = np.nan  # 90% sparse
    f = GraphRegularizedAPAFactorizer(rank=rank, max_iter=max_iter, random_state=42)
    start = time.time()
    f.fit(apa)
    elapsed = time.time() - start
    peak_mem_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # Linux: KB -> GB
    return {
        "n_genes": n_genes, "n_spots": n_spots, "rank": rank,
        "runtime_s": round(elapsed, 2),
        "peak_mem_gb": round(peak_mem_gb, 2),
        "completed": True,
        "recon_error": round(f.result_.reconstruction_error, 6),
    }

if __name__ == "__main__":
    results = []
    for n_spots in [12344, 42438, 100000]:
        print(f"Benchmarking n_spots={n_spots}...", flush=True)
        try:
            r = benchmark_factorizer(n_genes=8659, n_spots=n_spots)
            results.append(r)
            print(f"  OK: {r['runtime_s']}s, {r['peak_mem_gb']}GB, error={r['recon_error']}")
        except Exception as e:
            results.append({"n_spots": n_spots, "completed": False, "error": str(e)})
            print(f"  FAIL: {e}")
    out_dir = "benchmark_results/sparse_upgrade_scaling"
    os.makedirs(out_dir, exist_ok=True)
    with open(f"{out_dir}/scaling_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {out_dir}/scaling_results.json")

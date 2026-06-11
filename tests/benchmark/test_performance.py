"""
Performance benchmarks: Task 3.5 — 1000 genes × 2000 spots.

Run with:
    pytest tests/benchmark/test_performance.py -v --no-cov -s
"""

import time
import numpy as np
import pytest
from spagapa.benchmark import simulate_spatial_apa
from spagapa.imputation import GPImputer, SparseGPImputer
from spagapa.analysis import (
    identify_spatial_domains,
    test_differential_apa,
    identify_svapa_genes,
)
from spagapa.benchmark.evaluator import knn_imputation


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_data(n_spots, n_genes, dropout=0.4, seed=0):
    data = simulate_spatial_apa(
        n_spots=n_spots, n_genes=n_genes,
        dropout_rate=dropout, random_state=seed
    )
    apa = data['observed_apa'].T   # (n_genes, n_spots)
    coords = data['coordinates']
    return apa, coords


def _impute_all(apa_matrix, coords, imputer_cls, **kwargs):
    """Impute all genes using given imputer class."""
    n_genes = apa_matrix.shape[0]
    imputer = imputer_cls(**kwargs)
    imputed = np.zeros_like(apa_matrix)
    for i in range(n_genes):
        result, _ = imputer.impute(coords, apa_matrix[i])
        imputed[i] = result
    return imputed


# ── imputation performance ────────────────────────────────────────────────────

class TestImputationPerformance:

    @pytest.mark.parametrize("n_spots,n_genes", [
        (500,  100),
        (1000, 200),
        (2000, 100),
    ])
    def test_knn_imputation_time(self, n_spots, n_genes):
        apa, coords = _make_data(n_spots, n_genes)
        t0 = time.time()
        result = knn_imputation(apa.T, coords, k=10)
        elapsed = time.time() - t0
        assert not np.any(np.isnan(result))
        print(f"\nKNN {n_spots}×{n_genes}: {elapsed:.2f}s")
        assert elapsed < 120, f"KNN too slow: {elapsed:.1f}s"

    @pytest.mark.parametrize("n_spots,n_genes", [
        (200, 5),
        (500, 5),
    ])
    def test_gp_imputation_time(self, n_spots, n_genes):
        apa, coords = _make_data(n_spots, n_genes)
        imputer = GPImputer(kernel_type='matern')
        t0 = time.time()
        for i in range(n_genes):
            imputer.impute(coords, apa[i])
        elapsed = time.time() - t0
        per_gene = elapsed / n_genes
        print(f"\nGP {n_spots}×{n_genes}: {elapsed:.2f}s ({per_gene:.2f}s/gene)")
        assert per_gene < 30, f"GP too slow: {per_gene:.1f}s/gene"

    @pytest.mark.parametrize("n_spots,n_genes", [
        (500,  5),
        (1000, 5),
    ])
    def test_sparse_gp_imputation_time(self, n_spots, n_genes):
        apa, coords = _make_data(n_spots, n_genes)
        imputer = SparseGPImputer(n_inducing=30, inducing_method='kmeans')
        t0 = time.time()
        for i in range(n_genes):
            imputer.impute(coords, apa[i])
        elapsed = time.time() - t0
        per_gene = elapsed / n_genes
        print(f"\nSparseGP {n_spots}×{n_genes}: {elapsed:.2f}s ({per_gene:.2f}s/gene)")
        assert per_gene < 30, f"SparseGP too slow: {per_gene:.1f}s/gene"


# ── analysis performance ──────────────────────────────────────────────────────

class TestAnalysisPerformance:

    def test_domain_identification_1000x2000(self):
        """Domain identification on 1000 genes × 2000 spots < 60s."""
        # Use imputed (no NaN) data for domain ID
        data = simulate_spatial_apa(n_spots=2000, n_genes=1000,
                                    dropout_rate=0.0, random_state=0)
        apa = data['true_apa'].T   # no NaN
        coords = data['coordinates']

        t0 = time.time()
        labels, stats = identify_spatial_domains(
            apa, coords, method='kmeans', n_clusters=5, refine=False
        )
        elapsed = time.time() - t0
        print(f"\nDomain ID 2000×1000: {elapsed:.2f}s")
        assert elapsed < 60
        assert len(labels) == 2000

    def test_differential_apa_1000x2000(self):
        """Differential APA on 1000 genes × 2000 spots < 120s."""
        data = simulate_spatial_apa(n_spots=2000, n_genes=1000,
                                    dropout_rate=0.0, random_state=0)
        apa = data['true_apa'].T
        g1 = np.arange(1000)
        g2 = np.arange(1000, 2000)

        t0 = time.time()
        results = test_differential_apa(apa, g1, g2)
        elapsed = time.time() - t0
        print(f"\nDiff APA 2000×1000: {elapsed:.2f}s")
        assert elapsed < 120
        assert len(results) == 1000

    def test_svapa_detection_500x500(self):
        """SVAPA detection on 500 genes × 500 spots < 60s."""
        data = simulate_spatial_apa(n_spots=500, n_genes=500,
                                    dropout_rate=0.0, random_state=0)
        apa = data['true_apa'].T
        coords = data['coordinates']

        t0 = time.time()
        svapa, results = identify_svapa_genes(apa, coords)
        elapsed = time.time() - t0
        print(f"\nSVAPA 500×500: {elapsed:.2f}s, found {len(svapa)} genes")
        assert elapsed < 60


# ── memory usage ──────────────────────────────────────────────────────────────

class TestMemoryUsage:

    def test_knn_memory_reasonable(self):
        """KNN imputation on 2000×500 should not crash."""
        apa, coords = _make_data(2000, 500)
        result = knn_imputation(apa.T, coords, k=10)
        assert result.shape == (2000, 500)

    def test_domain_id_memory_reasonable(self):
        """Domain ID on 2000×1000 should not crash."""
        data = simulate_spatial_apa(n_spots=2000, n_genes=1000,
                                    dropout_rate=0.0, random_state=0)
        apa = data['true_apa'].T
        coords = data['coordinates']
        labels, _ = identify_spatial_domains(
            apa, coords, method='kmeans', n_clusters=5, refine=False
        )
        assert len(labels) == 2000

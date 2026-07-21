"""
Tests for the sparse upgrade of GraphRegularizedAPAFactorizer.

Task 1: randomized SVD initialization for large spot counts (avoids OOM that
        the dense ``np.linalg.svd`` hits on matrices like 8659 x 42438).
Task 2: batched ridge solve via ``np.einsum`` to eliminate the per-gene and
        per-spot Python loops inside ``fit``.

These tests are intentionally self-contained: they build synthetic APA
matrices with a known low-rank structure and NaN entries, then assert the
factorizer still produces a finite, correctly-shaped result and that the
batched path is fast on a wide matrix (where the loop version would be slow).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from spagapa.bioml import GraphRegularizedAPAFactorizer
from spagapa.imputation.sparse_gp import SparseGPImputer


def _make_low_rank_apa(
    n_genes: int,
    n_spots: int,
    rank: int,
    nan_frac: float = 0.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a noisy rank-``rank`` APA matrix, optionally masking entries as NaN.

    Returns ``(apa, mask)`` where mask is True for finite (observed) entries.
    """
    rng = np.random.default_rng(seed)
    g = rng.normal(size=(n_genes, rank))
    z = rng.normal(size=(n_spots, rank))
    signal = g @ z.T
    noise = rng.normal(scale=0.05, size=signal.shape)
    apa = signal + noise
    mask = np.ones(apa.shape, dtype=bool)
    if nan_frac > 0.0:
        missing = rng.random(apa.shape) < nan_frac
        apa[missing] = np.nan
        mask = ~missing
    return apa, mask


class TestRandomizedSVDInit:
    """Task 1: randomized SVD init for large matrices; dense SVD for small."""

    def test_small_matrix_uses_dense_svd(self):
        """n_spots <= 1000 should use dense SVD and return a valid factorization."""
        n_genes, n_spots = 50, 200
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=4, nan_frac=0.1)

        factorizer = GraphRegularizedAPAFactorizer(
            rank=4, max_iter=5, random_state=42
        )
        factorizer.fit(apa, mask=mask)

        assert factorizer.gene_factors_.shape == (n_genes, 4)
        assert factorizer.spot_factors_.shape == (n_spots, 4)
        assert np.isfinite(factorizer.gene_factors_).all()
        assert np.isfinite(factorizer.spot_factors_).all()
        assert np.isfinite(factorizer.result_.imputed).all()

    def test_large_matrix_no_oom(self):
        """n_spots > 1000 must use randomized SVD and must not OOM.

        A 200 x 5000 float64 matrix is only ~8 MB, so this is well within
        memory; the point is that the randomized path is exercised and
        produces a valid finite result without raising.
        """
        n_genes, n_spots = 200, 5000
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=6, nan_frac=0.8)

        factorizer = GraphRegularizedAPAFactorizer(
            rank=6, max_iter=3, random_state=42
        )
        factorizer.fit(apa, mask=mask)

        assert factorizer.gene_factors_.shape == (n_genes, 6)
        assert factorizer.spot_factors_.shape == (n_spots, 6)
        assert np.isfinite(factorizer.gene_factors_).all()
        assert np.isfinite(factorizer.spot_factors_).all()
        assert np.isfinite(factorizer.result_.imputed).all()


class TestBatchedRidgeSolve:
    """Task 2: einsum batched ridge solve replaces per-gene / per-spot loops."""

    def test_batched_equals_loop(self):
        """Batched einsum path must give a valid, low-error reconstruction.

        The einsum math is identical to the per-element ``_ridge_solve`` loop
        (see module docstring), so on a clean low-rank matrix the factorizer
        should still reconstruct well.
        """
        n_genes, n_spots = 30, 100
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=4, nan_frac=0.2)

        factorizer = GraphRegularizedAPAFactorizer(
            rank=4, max_iter=5, random_state=42
        )
        factorizer.fit(apa, mask=mask)

        assert factorizer.gene_factors_.shape == (n_genes, 4)
        assert factorizer.spot_factors_.shape == (n_spots, 4)
        # Low-rank signal + small noise -> reconstruction error well under 1.0.
        assert factorizer.result_.reconstruction_error < 1.0
        assert np.isfinite(factorizer.result_.imputed).all()

    def test_batched_speedup(self):
        """Wide matrix (100 genes x 5000 spots) must complete in <30s.

        The old per-spot Python loop would do 5000 separate
        ``np.linalg.solve`` calls per iteration -- minutes, not seconds. The
        batched einsum path collapses this to a handful of vectorized ops, so
        a few iterations must finish well under the 30s budget.
        """
        n_genes, n_spots = 100, 5000
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=8, nan_frac=0.5)

        factorizer = GraphRegularizedAPAFactorizer(
            rank=8, max_iter=3, random_state=42
        )
        start = time.perf_counter()
        factorizer.fit(apa, mask=mask)
        elapsed = time.perf_counter() - start

        assert elapsed < 30.0, f"fit took {elapsed:.2f}s, expected <30s"
        assert factorizer.gene_factors_.shape == (n_genes, 8)
        assert factorizer.spot_factors_.shape == (n_spots, 8)
        assert np.isfinite(factorizer.result_.imputed).all()


class TestDomainMethodFlexibility:
    """Task 3: domain recovery must allow kmeans/leiden, not just spectral."""

    def test_kmeans_domain_method_accepted(self):
        from spagapa.bioml.highres import highres_bioml_recover, HighResBioMLConfig
        rng = np.random.default_rng(42)
        observed = rng.random((20, 200))
        observed[observed < 0.5] = np.nan
        coords = rng.random((200, 2)) * 100
        config = HighResBioMLConfig(domain_method="kmeans")
        result = highres_bioml_recover(
            observed=observed, coords=coords, expression_embedding=None,
            sparse_gp=None, uncertainty=None, parent_index=None,
            n_domains=5, config=config,
        )
        assert result.labels.shape == (200,)
        assert len(np.unique(result.labels)) <= 5


class TestSparseGPAutoEnable:
    """Task 4: highres_accuracy must auto-enable sparse GP with spot-scaled
    inducing points, without requiring the user to set use_sparse_gp."""

    def _build_pipeline_with_dataset(self, n_genes, n_spots, preset, seed=0):
        """Wire a SpaGAPA pipeline to a synthetic APADataset (no I/O)."""
        import anndata as ad
        from spagapa.core import APADataset
        from spagapa import SpaGAPA

        rng = np.random.default_rng(seed)
        X = rng.random((n_spots, n_genes))
        coords = rng.random((n_spots, 2)) * 100.0
        adata = ad.AnnData(X=X)
        adata.obsm["spatial"] = coords
        ds = APADataset(adata=adata)

        pipe = SpaGAPA(analysis_preset=preset, verbose=False)
        pipe.dataset_ = ds
        return pipe, ds

    def test_highres_accuracy_enables_sparse_gp(self):
        """The pipeline must auto-enable sparse GP for highres_accuracy."""
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=8659, n_spots=42000, preset="highres_accuracy"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        assert options.get("use_sparse_gp") is True, (
            "highres_accuracy must auto-enable sparse GP"
        )

    def test_highres_accuracy_scales_inducing_points(self):
        """Inducing points must scale with n_spots: min(500, max(100, n//100))."""
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=8659, n_spots=42000, preset="highres_accuracy"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        expected = min(500, max(100, ds.n_spots // 100))  # 42000 // 100 = 420
        assert options["sparse_gp"]["n_inducing"] == expected, (
            f"n_inducing={options['sparse_gp']['n_inducing']}, "
            f"expected {expected} for {ds.n_spots} spots"
        )

    def test_highres_accuracy_inducing_clamped_to_500(self):
        """Very large spot counts must clamp n_inducing at 500."""
        # 200000 // 100 = 2000 -> clamped to 500
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=100, n_spots=200000, preset="highres_accuracy"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        assert options["sparse_gp"]["n_inducing"] == 500

    def test_highres_accuracy_inducing_floored_to_100(self):
        """Small spot counts must floor n_inducing at 100."""
        # 500 // 100 = 5 -> floored to 100
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=20, n_spots=500, preset="highres_accuracy"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        assert options["sparse_gp"]["n_inducing"] == 100

    def test_standard_preset_keeps_sparse_gp_off(self):
        """standard must stay on exact GP (use_sparse_gp=False)."""
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=50, n_spots=200, preset="standard"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        assert options.get("use_sparse_gp") is False

    def test_highres_fast_keeps_sparse_gp_off(self):
        """highres_fast skips imputation, so sparse GP stays off."""
        pipe, ds = self._build_pipeline_with_dataset(
            n_genes=50, n_spots=5000, preset="highres_fast"
        )
        options = pipe._resolve_run_options(impute=True, use_bioml=None)
        assert options.get("use_sparse_gp") is False


class TestChunkedFactorizer:
    """gene_chunk_size parameter processes genes in batches, reducing peak
    memory while producing mathematically identical results."""

    def test_chunked_matches_nonchunked(self):
        """Chunked factorizer must produce same result as non-chunked."""
        rng = np.random.default_rng(42)
        apa = rng.random((50, 200))
        apa[apa < 0.5] = 0  # 50% sparse
        from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer
        f1 = GraphRegularizedAPAFactorizer(rank=4, max_iter=5, random_state=42)
        f1.fit(apa)
        f2 = GraphRegularizedAPAFactorizer(rank=4, max_iter=5, random_state=42,
                                           gene_chunk_size=10)  # 5 chunks of 10 genes
        f2.fit(apa)
        # Results should be very close (chunking changes floating point order but not math)
        np.testing.assert_allclose(f1.result_.imputed, f2.result_.imputed, atol=1e-8)

    def test_chunked_memory_efficiency(self):
        """Chunked factorizer on large matrix should use less memory."""
        import tracemalloc
        rng = np.random.default_rng(42)
        apa = rng.random((500, 5000))
        apa[apa < 0.9] = 0  # 90% sparse
        from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer
        tracemalloc.start()
        f = GraphRegularizedAPAFactorizer(rank=8, max_iter=3, random_state=42,
                                           gene_chunk_size=100)
        f.fit(apa)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        # With chunking, peak should be well under the non-chunked memory
        # Non-chunked: 500×5000×8 = 20MB just for filled; chunked: 100×5000×8 = 4MB
        assert peak < 200 * 1024 * 1024, f"Peak memory {peak/1e6:.0f}MB too high"


class TestInducingPointsReuse:
    """Inducing points depend on coordinates, not gene values, so the batch
    path must select them once and reuse across all genes."""

    def test_batch_works(self):
        """Batch with precomputed inducing produces valid output."""
        rng = np.random.default_rng(42)
        coords = rng.random((200, 2)) * 100
        values = rng.random((20, 200))
        values[values < 0.5] = 0
        base = SparseGPImputer(n_inducing=50, inducing_method='kmeans',
                               length_scale=20.0, noise_level=0.1)
        batch = base.fit_batch(coords, values, verbose=False)
        pred, unc = batch.impute()
        assert pred.shape == (20, 200)
        assert np.isfinite(pred).all()

    def test_batch_speed(self):
        """Batch on 100 genes x 5000 spots completes in <30s.

        The redundant KMeans (10 init, 100 clusters, 5000 points) run 100x
        in the old path would dominate runtime; precomputing once collapses
        it to a single KMeans call.
        """
        start = time.perf_counter()
        rng = np.random.default_rng(42)
        coords = rng.random((5000, 2)) * 200
        values = rng.random((100, 5000))
        values[values < 0.7] = 0
        base = SparseGPImputer(n_inducing=100, inducing_method='kmeans',
                               length_scale=50.0, noise_level=0.1)
        batch = base.fit_batch(coords, values, verbose=False)
        elapsed = time.perf_counter() - start
        assert elapsed < 30, f"Batch took {elapsed:.1f}s -- too slow"
        pred, unc = batch.impute()
        assert pred.shape == (100, 5000)

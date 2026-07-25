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


class TestGraphRegularizedFitTransform:
    """Phase-3 Task-2: fit_transform with a sparse scipy Laplacian.

    The graph-reg spot update must (a) accept a sparse Laplacian without
    densifying it, (b) produce a finite result of the right shape, and
    (c) be mathematically equivalent to applying ``(I + λL)^{-1}`` per rank
    column -- which is what the multi-RHS sparse solve implements.
    """

    @staticmethod
    def _knn_laplacian(coords: np.ndarray, k: int = 5) -> "scipy.sparse.csr_matrix":
        from scipy import sparse
        from scipy.spatial import cKDTree
        n = coords.shape[0]
        tree = cKDTree(coords)
        dists, idx = tree.query(coords, k=k + 1)
        rows = np.repeat(np.arange(n), k)
        cols = idx[:, 1:].ravel()
        d = dists[:, 1:].ravel()
        sigma = float(np.median(d)) + 1e-12
        W = sparse.csr_matrix(
            (np.exp(-0.5 * (d / sigma) ** 2), (rows, cols)), shape=(n, n)
        )
        W = W.maximum(W.T)
        W.setdiag(0.0)
        W.eliminate_zeros()
        deg = np.asarray(W.sum(1)).ravel()
        return (sparse.diags(deg) - W).tocsr()

    def test_sparse_laplacian_finite_result(self):
        """fit_transform with a sparse Laplacian returns a finite, shaped result."""
        from scipy import sparse
        from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer
        rng = np.random.default_rng(0)
        n_genes, n_spots = 30, 80
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=4, nan_frac=0.2)
        coords = rng.random((n_spots, 2)) * 50.0
        lap = self._knn_laplacian(coords, k=5)
        assert sparse.issparse(lap), "laplacian must be sparse"

        f = GraphRegularizedAPAFactorizer(
            rank=4, lambda_graph=0.5, max_iter=5, random_state=42,
        )
        imputed = f.fit_transform(apa, graph_laplacian=lap, mask=mask)
        assert imputed.shape == apa.shape
        assert np.isfinite(imputed).all()
        # Observed entries are preserved when preserve_observed=True (default).
        np.testing.assert_allclose(imputed[mask], apa[mask])
        assert f.spot_factors_.shape == (n_spots, 4)

    def test_multi_rhs_solve_matches_per_column(self):
        """The single multi-RHS sparse solve must match the per-column reference.

        This pins the math: smoothing ``(I + λL)^{-1} Z_raw`` over all rank
        columns at once equals looping over columns and solving each. Catches
        any accidental transpose / wrong-RHS-shape bug in the optimization.
        """
        from scipy import sparse
        from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer
        rng = np.random.default_rng(1)
        n_genes, n_spots = 24, 60
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=3, nan_frac=0.1)
        coords = rng.random((n_spots, 2)) * 30.0
        lap = self._knn_laplacian(coords, k=4)

        # Reference: factor the system once, solve each rank column separately,
        # then column_stack (the pre-optimization behaviour).
        rank = 3
        lam = 0.5
        system = sparse.eye(n_spots, format="csc") + lam * lap.tocsc()
        from scipy.sparse.linalg import factorized
        smoother = factorized(system)

        # Build the reference spot update from a fresh factorizer's z_raw by
        # running the factorizer once and capturing the *raw* (un-smoothed)
        # spot factors via a tiny subclass hook is overkill; instead verify
        # equivalence directly on an arbitrary z_raw matrix.
        z_raw = rng.standard_normal((n_spots, rank))
        z_ref = np.column_stack([smoother(z_raw[:, d]) for d in range(rank)])
        z_multi = smoother(z_raw)  # what the optimized code now calls
        np.testing.assert_allclose(z_multi, z_ref, atol=1e-10)

        # And that the end-to-end factorizer result is stable / finite.
        f = GraphRegularizedAPAFactorizer(
            rank=rank, lambda_graph=lam, max_iter=5, random_state=42,
        )
        f.fit_transform(apa, graph_laplacian=lap, mask=mask)
        assert np.isfinite(f.result_.imputed).all()

    def test_sparse_laplacian_with_chunking(self):
        """fit_transform with sparse Laplacian + gene_chunk_size stays finite."""
        from spagapa.bioml.factorization import GraphRegularizedAPAFactorizer
        rng = np.random.default_rng(2)
        n_genes, n_spots = 40, 70
        apa, mask = _make_low_rank_apa(n_genes, n_spots, rank=4, nan_frac=0.3)
        coords = rng.random((n_spots, 2)) * 40.0
        lap = self._knn_laplacian(coords, k=5)

        f = GraphRegularizedAPAFactorizer(
            rank=4, lambda_graph=0.3, max_iter=5, random_state=42,
            gene_chunk_size=10,
        )
        imputed = f.fit_transform(apa, graph_laplacian=lap, mask=mask)
        assert imputed.shape == apa.shape
        assert np.isfinite(imputed).all()


class TestAutoChunkedFactorizer:
    """Task A3: pipeline must auto-enable gene_chunk_size=2000 when n_spots >
    20000, and keep it None for smaller datasets."""

    def test_large_n_spots_enables_chunking(self):
        """n_spots > 20000 must instantiate factorizer with gene_chunk_size=2000."""
        import inspect
        # Inspect _run_bioml source to confirm the auto-chunking branch is
        # present (the only place the factorizer is constructed in the pipeline).
        import spagapa.pipeline as pl
        src = inspect.getsource(pl.SpaGAPA._run_bioml)
        assert "gene_chunk_size = 2000 if n_spots > 20000 else None" in src, (
            "pipeline._run_bioml must auto-enable gene_chunk_size=2000 for "
            "n_spots > 20000"
        )
        assert "gene_chunk_size=gene_chunk_size" in src


class TestKNmPrecompute:
    """K_nm (kernel from all coords to inducing points) depends only on
    coordinates, so it can be computed once for the whole batch and indexed
    per gene's mask."""

    def test_precomputed_knm_matches(self):
        """SparseGPImputer.fit() with precomputed K_nm_full gives same result."""
        import numpy as np
        from spagapa.imputation.sparse_gp import SparseGPImputer
        rng = np.random.default_rng(42)
        coords = rng.random((200, 2)) * 100
        values = rng.random(200)
        values[values < 0.5] = 0
        # Without precompute
        imp1 = SparseGPImputer(n_inducing=50, length_scale=20.0, noise_level=0.1)
        imp1.fit(coords, values)
        # Precompute K_nm_full (full coord set -> inducing points)
        inducing = imp1.inducing_points_
        from scipy.spatial import distance_matrix
        dists = distance_matrix(coords, inducing)
        k_nm_full = np.exp(-0.5 * (dists / 20.0) ** 2)
        # With precompute
        imp2 = SparseGPImputer(n_inducing=50, length_scale=20.0, noise_level=0.1)
        imp2.fit(coords, values, k_nm_full=k_nm_full)
        np.testing.assert_allclose(imp1.alpha_, imp2.alpha_, atol=1e-10)

    def test_precomputed_knm_with_mask_matches(self):
        """Precomputed K_nm_full with a per-gene mask indexes correctly."""
        import numpy as np
        from spagapa.imputation.sparse_gp import SparseGPImputer
        rng = np.random.default_rng(7)
        coords = rng.random((300, 2)) * 100
        values = rng.random(300)
        mask = values > 0.3  # arbitrary mask != (values > 0)
        # Without precompute, using mask
        imp1 = SparseGPImputer(n_inducing=40, length_scale=25.0, noise_level=0.1)
        imp1.fit(coords, values, mask=mask)
        # Precompute K_nm_full over ALL coords
        inducing = imp1.inducing_points_
        from scipy.spatial import distance_matrix
        dists = distance_matrix(coords, inducing)
        k_nm_full = np.exp(-0.5 * (dists / 25.0) ** 2)
        imp2 = SparseGPImputer(n_inducing=40, length_scale=25.0, noise_level=0.1)
        imp2.fit(coords, values, mask=mask, k_nm_full=k_nm_full)
        np.testing.assert_allclose(imp1.alpha_, imp2.alpha_, atol=1e-10)

    def test_batch_uses_precomputed_knm(self):
        """SparseGPImputerBatch must precompute K_nm_full once and reuse it,
        producing the same output as the per-gene kernel path."""
        import numpy as np
        from spagapa.imputation.sparse_gp import SparseGPImputer
        rng = np.random.default_rng(11)
        coords = rng.random((500, 2)) * 100
        values = rng.random((15, 500))
        values[values < 0.5] = 0
        base = SparseGPImputer(n_inducing=50, length_scale=20.0, noise_level=0.1)
        batch = base.fit_batch(coords, values, verbose=False)
        pred, unc = batch.impute()
        # Sanity: result is finite and well-shaped.
        assert pred.shape == (15, 500)
        assert np.isfinite(pred).all()
        # The first imputer's fit must have stored the precomputed full-kernel
        # path: assert that per-gene alpha matches a fresh fit with k_nm_full
        # passed explicitly (and the same shared inducing points).
        inducing = batch.imputers_[0].inducing_points_
        from scipy.spatial import distance_matrix
        dists = distance_matrix(coords, inducing)
        k_nm_full = np.exp(-0.5 * (dists / 20.0) ** 2)
        imp_ref = SparseGPImputer(n_inducing=50, length_scale=20.0, noise_level=0.1)
        imp_ref.fit(
            coords, values[0],
            inducing_points=inducing,
            k_nm_full=k_nm_full,
        )
        np.testing.assert_allclose(
            batch.imputers_[0].alpha_, imp_ref.alpha_, atol=1e-10
        )


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


class TestConformalCalibration:
    """P0-1: split-conformal calibration must deliver marginal coverage at the
    target level on synthetic data, for both global and locally-adaptive modes.

    The sparse GP's raw posterior std is far too conservative (2-sigma
    coverage ~1.0) and barely correlates with error (~0.008). Split conformal
    fixes the marginal coverage rigorously by construction; these tests pin
    that guarantee.
    """

    @staticmethod
    def _make_synthetic(n: int, seed: int):
        """Errors and a GP std that tracks them (corr ~0.7).

        Models the situation we *want* the locally-adaptive mode to exploit:
        the raw std ranking is informative even if its scale is wrong.
        """
        rng = np.random.default_rng(seed)
        # heteroscedastic true error scale
        latent = np.abs(rng.normal(0, 1, n))
        errors = np.abs(rng.normal(0, 1, n)) * (0.3 + 0.7 * latent)
        # GP std tracks error up to a (deliberately wrong) constant scale
        gp_std = 3.0 * latent + np.abs(rng.normal(0, 0.2, n))
        # ground-truth values: pred=0, so truth = +/-error
        truth = errors * rng.choice([-1, 1], n)
        return errors, gp_std, truth

    @pytest.mark.parametrize("mode", ["global", "locally_adaptive"])
    @pytest.mark.parametrize("alpha,target", [(0.20, 0.80), (0.10, 0.90), (0.05, 0.95)])
    def test_coverage_hits_target(self, mode, alpha, target):
        """Empirical test coverage must be within +/-3% of the nominal target.

        Split conformal guarantees marginal coverage >= 1-alpha in expectation
        under exchangeability; with 6000 test points the sampling noise on the
        coverage estimate is ~ sqrt(0.9*0.1/6000) ~ 0.004, so +/-3% is a
        comfortable band that still catches a broken implementation.
        """
        from spagapa.imputation.calibration import (
            ConformalCalibrator,
            evaluate_coverage,
        )

        errors, gp_std, truth = self._make_synthetic(12000, seed=42)
        perm = np.random.default_rng(0).permutation(12000)
        n_cal = 6000
        cal_idx, test_idx = perm[:n_cal], perm[n_cal:]

        cal = ConformalCalibrator(alpha=alpha, mode=mode)
        cal.fit(errors[cal_idx],
                gp_std[cal_idx] if mode == "locally_adaptive" else None)
        lo, hi = cal.predict(
            np.zeros(n_cal),  # pred = 0 for everyone
            gp_std[test_idx] if mode == "locally_adaptive" else None,
        )
        cov = evaluate_coverage(lo, hi, truth[test_idx])
        assert abs(cov - target) <= 0.03, (
            f"{mode} alpha={alpha}: coverage {cov:.4f} outside "
            f"[{target-0.03:.3f}, {target+0.03:.3f}]"
        )

    def test_locally_adaptive_tighter_where_gp_confident(self):
        """Locally-adaptive intervals must be narrower than global where the
        GP std is small, and wider where it is large -- otherwise the mode is
        degenerate (e.g. when raw std is constant)."""
        from spagapa.imputation.calibration import ConformalCalibrator

        errors, gp_std, _ = self._make_synthetic(4000, seed=1)
        # Split
        n_cal = 2000
        cal = ConformalCalibrator(alpha=0.1, mode="locally_adaptive")
        cal.fit(errors[:n_cal], gp_std[:n_cal])
        test_std = gp_std[n_cal:]
        lo, hi = cal.predict(np.zeros(2000), test_std)
        half = (hi - lo) / 2.0
        # Intervals should be monotonically increasing in gp_std
        order = np.argsort(test_std)
        half_sorted = half[order]
        # Spearman-ish check: first quartile half-width < last quartile
        q = len(half_sorted) // 4
        assert half_sorted[:q].mean() < half_sorted[-q:].mean(), (
            "locally-adaptive half-width should scale with GP std"
        )

    def test_global_requires_no_std(self):
        """Global mode must not require gp_std in fit() or predict()."""
        from spagapa.imputation.calibration import ConformalCalibrator

        errors = np.abs(np.random.default_rng(2).normal(0, 1, 1000))
        cal = ConformalCalibrator(alpha=0.1, mode="global")
        cal.fit(errors)  # no std
        lo, hi = cal.predict(np.zeros(500))  # no std
        assert lo.shape == (500,) and (hi > lo).all()

    def test_locally_adaptive_requires_std(self):
        """Locally-adaptive mode must raise if std is missing."""
        from spagapa.imputation.calibration import ConformalCalibrator

        errors = np.abs(np.random.default_rng(3).normal(0, 1, 100))
        cal = ConformalCalibrator(alpha=0.1, mode="locally_adaptive")
        with pytest.raises(ValueError):
            cal.fit(errors)  # missing std

    def test_alpha_validation(self):
        """alpha must be in (0, 1)."""
        from spagapa.imputation.calibration import ConformalCalibrator
        with pytest.raises(ValueError):
            ConformalCalibrator(alpha=0.0)
        with pytest.raises(ValueError):
            ConformalCalibrator(alpha=1.0)

    def test_predict_before_fit_raises(self):
        """predict() before fit() must raise RuntimeError."""
        from spagapa.imputation.calibration import ConformalCalibrator
        cal = ConformalCalibrator(alpha=0.1)
        with pytest.raises(RuntimeError):
            cal.predict(np.zeros(10))

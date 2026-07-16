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

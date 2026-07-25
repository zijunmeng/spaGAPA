"""
Graph-regularized APA matrix factorization.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized

# Threshold on the number of spots (columns of the APA matrix) above which the
# SVD initialization switches from dense to randomized. Dense SVD on a matrix
# like 8659 x 42438 OOMs; randomized SVD only materializes ``rank`` components.
_SVD_RANDOMIZED_SPOT_THRESHOLD = 1000


def _randomized_svd_init(
    matrix: np.ndarray, n_components: int, n_iter: int, random_state: int
) -> np.ndarray:
    """Randomized SVD top components, transposed to (n_samples, n_components).

    Imported lazily so the factorizer stays importable on environments without
    scikit-learn.
    """
    from sklearn.utils.extmath import randomized_svd

    _, _, vt = randomized_svd(
        matrix,
        n_components=n_components,
        n_iter=n_iter,
        random_state=random_state,
    )
    return vt[:n_components].T


@dataclass
class FactorizationResult:
    """Result container for graph-regularized APA factorization."""

    imputed: np.ndarray
    gene_factors: np.ndarray
    spot_factors: np.ndarray
    n_iter: int
    reconstruction_error: float


class GraphRegularizedAPAFactorizer:
    """
    CPU-friendly graph-regularized low-rank APA factorization.

    The implementation uses alternating weighted ridge regression. Spot
    factors are smoothed on the supplied graph Laplacian after each spot update.
    """

    def __init__(
        self,
        rank: int = 8,
        lambda_graph: float = 1.0,
        lambda_l2: float = 1e-2,
        max_iter: int = 50,
        tol: float = 1e-4,
        random_state: int = 42,
        preserve_observed: bool = True,
        gene_chunk_size: Optional[int] = None,
    ):
        if rank <= 0:
            raise ValueError("rank must be positive")
        if lambda_graph < 0:
            raise ValueError("lambda_graph must be non-negative")
        if lambda_l2 < 0:
            raise ValueError("lambda_l2 must be non-negative")
        if max_iter <= 0:
            raise ValueError("max_iter must be positive")
        if gene_chunk_size is not None and gene_chunk_size <= 0:
            raise ValueError("gene_chunk_size must be positive or None")
        self.rank = int(rank)
        self.lambda_graph = float(lambda_graph)
        self.lambda_l2 = float(lambda_l2)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.random_state = int(random_state)
        self.preserve_observed = bool(preserve_observed)
        # When set, gene/spot update einsums are processed in batches of this
        # many genes. This caps peak memory on wide matrices (the einsum
        # otherwise materializes a full n_genes x n_spots intermediate).
        # None = no chunking (backward compatible).
        self.gene_chunk_size = (
            int(gene_chunk_size) if gene_chunk_size is not None else None
        )

        self.gene_factors_: Optional[np.ndarray] = None
        self.spot_factors_: Optional[np.ndarray] = None
        self.result_: Optional[FactorizationResult] = None

    @staticmethod
    def _fill_missing(values: np.ndarray) -> np.ndarray:
        filled = values.copy()
        gene_means = np.nanmean(filled, axis=1)
        gene_means = np.nan_to_num(gene_means, nan=0.0)
        missing_gene, missing_spot = np.where(~np.isfinite(filled))
        filled[missing_gene, missing_spot] = gene_means[missing_gene]
        return filled

    @staticmethod
    def _ridge_solve(design: np.ndarray, target: np.ndarray, weights: np.ndarray, lambda_l2: float) -> np.ndarray:
        valid = weights > 0
        if valid.sum() == 0:
            return np.zeros(design.shape[1], dtype=float)
        x = design[valid]
        y = target[valid]
        w = weights[valid]
        xw = x * np.sqrt(w)[:, None]
        yw = y * np.sqrt(w)
        lhs = xw.T @ xw + lambda_l2 * np.eye(design.shape[1])
        rhs = xw.T @ yw
        return np.linalg.solve(lhs, rhs)

    def _initial_spot_factors(self, filled: np.ndarray) -> np.ndarray:
        rng = np.random.default_rng(self.random_state)
        centered = filled - filled.mean(axis=1, keepdims=True)
        try:
            # Dense SVD OOMs on wide matrices (e.g. 8659 x 42438). For large
            # spot counts fall back to randomized SVD which only materializes
            # ``rank`` components. Small matrices keep the dense path: it is
            # exact and faster than the randomized estimator.
            if centered.shape[1] > _SVD_RANDOMIZED_SPOT_THRESHOLD:
                z = _randomized_svd_init(
                    centered,
                    n_components=self.rank,
                    n_iter=5,
                    random_state=self.random_state,
                )
            else:
                _, _, vt = np.linalg.svd(centered, full_matrices=False)
                z = vt[: self.rank, :].T
            if z.shape[1] < self.rank:
                pad = rng.normal(scale=0.01, size=(z.shape[0], self.rank - z.shape[1]))
                z = np.concatenate([z, pad], axis=1)
            return z
        except np.linalg.LinAlgError:
            return rng.normal(scale=0.01, size=(filled.shape[1], self.rank))

    def fit(
        self,
        apa_matrix: np.ndarray,
        graph_laplacian: Optional[sparse.spmatrix] = None,
        mask: Optional[np.ndarray] = None,
        confidence: Optional[np.ndarray] = None,
    ) -> "GraphRegularizedAPAFactorizer":
        y = np.asarray(apa_matrix, dtype=float)
        if y.ndim != 2:
            raise ValueError("apa_matrix must be 2-dimensional")
        n_genes, n_spots = y.shape

        if mask is None:
            mask = np.isfinite(y)
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != y.shape:
            raise ValueError("mask must have the same shape as apa_matrix")

        if confidence is None:
            weights = mask.astype(float)
        else:
            confidence = np.asarray(confidence, dtype=float)
            if confidence.shape != y.shape:
                raise ValueError("confidence must have the same shape as apa_matrix")
            weights = np.where(mask, np.nan_to_num(confidence, nan=0.0), 0.0)

        filled = self._fill_missing(y)
        z = self._initial_spot_factors(filled)
        g = np.zeros((n_genes, self.rank), dtype=float)

        smoother = None
        if graph_laplacian is not None and self.lambda_graph > 0:
            lap = graph_laplacian.tocsr()
            if lap.shape != (n_spots, n_spots):
                raise ValueError("graph_laplacian must have shape (n_spots, n_spots)")
            system = sparse.eye(n_spots, format="csc") + self.lambda_graph * lap.tocsc()
            smoother = factorized(system)

        last_error = np.inf
        n_iter = 0
        eye_k = self.lambda_l2 * np.eye(self.rank)
        for iteration in range(1, self.max_iter + 1):
            # Batched gene ridge solve: replace the per-gene Python loop with
            # a single vectorized einsum + batched solve. Mathematically
            # identical to calling ``_ridge_solve`` once per gene (rows of
            # ``weights`` that are all zero map to the zero RHS -> zero row,
            # matching ``_ridge_solve``'s ``valid.sum() == 0`` branch).
            #
            # When ``gene_chunk_size`` is set, the gene update is sliced over
            # chunks of genes. Each chunk materializes a (chunk, n_spots)
            # intermediate instead of (n_genes, n_spots), cutting peak memory
            # on wide APA matrices. The math is unchanged -- the sum over the
            # gene axis is independent per gene.
            if self.gene_chunk_size and self.gene_chunk_size < n_genes:
                chunk = self.gene_chunk_size
                g = np.zeros((n_genes, self.rank), dtype=float)
                for cs in range(0, n_genes, chunk):
                    ce = min(cs + chunk, n_genes)
                    w_chunk = weights[cs:ce, :]      # (chunk, n_spots)
                    f_chunk = filled[cs:ce, :]        # (chunk, n_spots)
                    lhs_g_c = np.einsum(
                        "sk,gs,sl->gkl", z, w_chunk, z, optimize=True
                    )
                    lhs_g_c += eye_k[None, :, :]
                    rhs_g_c = np.einsum(
                        "sk,gs,gs->gk", z, w_chunk, f_chunk, optimize=True
                    )
                    g[cs:ce] = np.linalg.solve(lhs_g_c, rhs_g_c[..., None])[..., 0]
            else:
                lhs_g = np.einsum("sk,gs,sl->gkl", z, weights, z, optimize=True)
                lhs_g += eye_k[None, :, :]
                rhs_g = np.einsum("sk,gs,gs->gk", z, weights, filled, optimize=True)
                g = np.linalg.solve(lhs_g, rhs_g[..., None])[..., 0]

            # Batched spot ridge solve: replace the per-spot Python loop with
            # one vectorized einsum + batched solve over all spots at once.
            # The spot lhs/rhs sum over the gene axis; when chunked, those
            # per-gene contributions are accumulated across chunks. Again the
            # math is unchanged -- summation is associative.
            if self.gene_chunk_size and self.gene_chunk_size < n_genes:
                chunk = self.gene_chunk_size
                lhs_s = np.zeros((n_spots, self.rank, self.rank))
                rhs_s = np.zeros((n_spots, self.rank))
                for cs in range(0, n_genes, chunk):
                    ce = min(cs + chunk, n_genes)
                    g_chunk = g[cs:ce, :]            # (chunk, rank)
                    w_chunk = weights[cs:ce, :]       # (chunk, n_spots)
                    f_chunk = filled[cs:ce, :]         # (chunk, n_spots)
                    lhs_s += np.einsum(
                        "gk,gs,gl->skl", g_chunk, w_chunk, g_chunk, optimize=True
                    )
                    rhs_s += np.einsum(
                        "gk,gs,gs->sk", g_chunk, w_chunk, f_chunk, optimize=True
                    )
                lhs_s += eye_k[None, :, :]
                z_raw = np.linalg.solve(lhs_s, rhs_s[..., None])[..., 0]
            else:
                lhs_s = np.einsum("gk,gs,gl->skl", g, weights, g, optimize=True)
                lhs_s += eye_k[None, :, :]
                rhs_s = np.einsum("gk,gs,gs->sk", g, weights, filled, optimize=True)
                z_raw = np.linalg.solve(lhs_s, rhs_s[..., None])[..., 0]

            if smoother is not None:
                # Apply the sparse graph smoothing operator ``(I + λL)^{-1}``
                # to all ``rank`` columns of ``z_raw`` in a single SuperLU
                # solve. ``scipy.sparse.linalg.factorized`` returns a callable
                # whose ``solve`` accepts a dense 2-D RHS and solves for every
                # column at once (one LU factorization, one triangular pass).
                # This is mathematically identical to looping over the rank
                # dimension and solving per column, but avoids ``rank`` Python
                # round-trips through the sparse solver each iteration. The
                # Laplacian stays sparse throughout (never densified to
                # (n_spots x n_spots)).
                z = smoother(z_raw)
            else:
                z = z_raw

            pred = g @ z.T
            residual = (filled - pred) * np.sqrt(weights)
            denom = max(float(weights.sum()), 1.0)
            error = float(np.sqrt(np.sum(residual * residual) / denom))
            n_iter = iteration
            if np.isfinite(last_error) and abs(last_error - error) / max(last_error, 1e-8) < self.tol:
                break
            last_error = error

        imputed = g @ z.T
        if self.preserve_observed:
            imputed[mask] = y[mask]

        self.gene_factors_ = g
        self.spot_factors_ = z
        self.result_ = FactorizationResult(
            imputed=imputed,
            gene_factors=g,
            spot_factors=z,
            n_iter=n_iter,
            reconstruction_error=last_error,
        )
        return self

    def fit_transform(
        self,
        apa_matrix: np.ndarray,
        graph_laplacian: Optional[sparse.spmatrix] = None,
        mask: Optional[np.ndarray] = None,
        confidence: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Fit the model and return the imputed APA matrix."""
        return self.fit(apa_matrix, graph_laplacian, mask, confidence).result_.imputed

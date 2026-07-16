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
    ):
        if rank <= 0:
            raise ValueError("rank must be positive")
        if lambda_graph < 0:
            raise ValueError("lambda_graph must be non-negative")
        if lambda_l2 < 0:
            raise ValueError("lambda_l2 must be non-negative")
        if max_iter <= 0:
            raise ValueError("max_iter must be positive")
        self.rank = int(rank)
        self.lambda_graph = float(lambda_graph)
        self.lambda_l2 = float(lambda_l2)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.random_state = int(random_state)
        self.preserve_observed = bool(preserve_observed)

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
        for iteration in range(1, self.max_iter + 1):
            for gene_idx in range(n_genes):
                g[gene_idx] = self._ridge_solve(
                    z,
                    filled[gene_idx],
                    weights[gene_idx],
                    self.lambda_l2,
                )

            z_raw = np.zeros_like(z)
            design = g
            for spot_idx in range(n_spots):
                z_raw[spot_idx] = self._ridge_solve(
                    design,
                    filled[:, spot_idx],
                    weights[:, spot_idx],
                    self.lambda_l2,
                )

            if smoother is not None:
                z = np.column_stack([smoother(z_raw[:, dim]) for dim in range(self.rank)])
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

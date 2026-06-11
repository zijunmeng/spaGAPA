"""
Multi-view graph construction for CPU-friendly BioML workflows.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


@dataclass
class MultiViewGraph:
    """Container for sparse spot-level multi-view graphs."""

    spatial: sparse.csr_matrix
    expression: Optional[sparse.csr_matrix]
    apa: Optional[sparse.csr_matrix]
    fused: sparse.csr_matrix
    weights: dict[str, float]

    def laplacian(self, normalized: bool = False) -> sparse.csr_matrix:
        """Return the fused graph Laplacian."""
        degrees = np.asarray(self.fused.sum(axis=1)).ravel()
        if not normalized:
            return sparse.diags(degrees) - self.fused

        inv_sqrt = np.zeros_like(degrees, dtype=float)
        valid = degrees > 0
        inv_sqrt[valid] = 1.0 / np.sqrt(degrees[valid])
        d_inv = sparse.diags(inv_sqrt)
        identity = sparse.eye(self.fused.shape[0], format="csr")
        return identity - d_inv @ self.fused @ d_inv


class MultiViewGraphBuilder:
    """
    Build sparse multi-view spot graphs from coordinates, expression, and APA.

    Parameters
    ----------
    n_neighbors : int, default=15
        KNN size for each view.
    spatial_weight, expression_weight, apa_weight : float
        Non-negative fusion weights. Available views are re-normalized.
    metric : str, default="euclidean"
        Distance metric used by nearest neighbors.
    """

    def __init__(
        self,
        n_neighbors: int = 15,
        spatial_weight: float = 0.5,
        expression_weight: float = 0.3,
        apa_weight: float = 0.2,
        metric: str = "euclidean",
    ):
        if n_neighbors <= 0:
            raise ValueError("n_neighbors must be positive")
        weights = [spatial_weight, expression_weight, apa_weight]
        if any(w < 0 for w in weights):
            raise ValueError("Graph fusion weights must be non-negative")
        if sum(weights) <= 0:
            raise ValueError("At least one graph fusion weight must be positive")

        self.n_neighbors = int(n_neighbors)
        self.spatial_weight = float(spatial_weight)
        self.expression_weight = float(expression_weight)
        self.apa_weight = float(apa_weight)
        self.metric = metric

    @staticmethod
    def _as_2d(name: str, values: np.ndarray) -> np.ndarray:
        arr = np.asarray(values, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"{name} must be 2-dimensional")
        return arr

    @staticmethod
    def _scale_features(values: np.ndarray) -> np.ndarray:
        scaler = StandardScaler()
        scaled = scaler.fit_transform(values)
        return np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)

    def _knn_rbf_graph(self, features: np.ndarray) -> sparse.csr_matrix:
        features = self._scale_features(features)
        n_samples = features.shape[0]
        if n_samples < 2:
            return sparse.csr_matrix((n_samples, n_samples), dtype=float)

        k = min(self.n_neighbors, n_samples - 1)
        nn = NearestNeighbors(n_neighbors=k + 1, metric=self.metric)
        nn.fit(features)
        distances, indices = nn.kneighbors(features)

        distances = distances[:, 1:]
        indices = indices[:, 1:]
        positive = distances[distances > 0]
        sigma = float(np.median(positive)) if positive.size else 1.0
        sigma = max(sigma, 1e-6)

        rows = np.repeat(np.arange(n_samples), k)
        cols = indices.ravel()
        weights = np.exp(-0.5 * (distances.ravel() / sigma) ** 2)
        graph = sparse.csr_matrix((weights, (rows, cols)), shape=(n_samples, n_samples))
        graph = graph.maximum(graph.T)
        graph.setdiag(0.0)
        graph.eliminate_zeros()
        return graph

    def _prepare_apa_features(
        self,
        apa_matrix: np.ndarray,
        uncertainty: Optional[np.ndarray],
    ) -> np.ndarray:
        apa = self._as_2d("apa_matrix", apa_matrix)
        gene_means = np.nanmean(apa, axis=1)
        gene_means = np.nan_to_num(gene_means, nan=0.0)
        missing_gene, missing_spot = np.where(~np.isfinite(apa))
        apa_filled = apa.copy()
        apa_filled[missing_gene, missing_spot] = gene_means[missing_gene]

        features = apa_filled.T
        if uncertainty is None:
            return features

        unc = self._as_2d("uncertainty", uncertainty)
        if unc.shape != apa.shape:
            raise ValueError("uncertainty must have the same shape as apa_matrix")
        confidence = 1.0 / (np.nan_to_num(unc, nan=np.nanmedian(unc)) + 1e-6)
        confidence = confidence / max(float(np.nanmedian(confidence)), 1e-6)
        confidence = np.clip(confidence, 0.0, 10.0)
        return (apa_filled * np.sqrt(confidence)).T

    @staticmethod
    def _normalize_weights(graphs: dict[str, sparse.csr_matrix], weights: dict[str, float]) -> dict[str, float]:
        active = {
            name: weight
            for name, weight in weights.items()
            if name in graphs and graphs[name] is not None and weight > 0
        }
        total = sum(active.values())
        if total <= 0:
            raise ValueError("No active graph view is available")
        return {name: weight / total for name, weight in active.items()}

    def build(
        self,
        coordinates: np.ndarray,
        expression_embedding: Optional[np.ndarray] = None,
        apa_matrix: Optional[np.ndarray] = None,
        uncertainty: Optional[np.ndarray] = None,
    ) -> MultiViewGraph:
        """Build a fused sparse graph from available views."""
        coords = self._as_2d("coordinates", coordinates)
        if coords.shape[1] != 2:
            raise ValueError("coordinates must have shape (n_spots, 2)")
        n_spots = coords.shape[0]

        spatial = self._knn_rbf_graph(coords)
        graphs: dict[str, sparse.csr_matrix] = {"spatial": spatial}
        expression = None
        apa = None

        if expression_embedding is not None:
            expr = self._as_2d("expression_embedding", expression_embedding)
            if expr.shape[0] != n_spots:
                raise ValueError("expression_embedding must have n_spots rows")
            expression = self._knn_rbf_graph(expr)
            graphs["expression"] = expression

        if apa_matrix is not None:
            apa_features = self._prepare_apa_features(apa_matrix, uncertainty)
            if apa_features.shape[0] != n_spots:
                raise ValueError("apa_matrix must have n_spots columns")
            apa = self._knn_rbf_graph(apa_features)
            graphs["apa"] = apa

        raw_weights = {
            "spatial": self.spatial_weight,
            "expression": self.expression_weight,
            "apa": self.apa_weight,
        }
        weights = self._normalize_weights(graphs, raw_weights)
        fused = sparse.csr_matrix((n_spots, n_spots), dtype=float)
        for name, weight in weights.items():
            fused = fused + weight * graphs[name]
        fused = fused.maximum(fused.T)
        fused.setdiag(0.0)
        fused.eliminate_zeros()

        return MultiViewGraph(
            spatial=spatial,
            expression=expression,
            apa=apa,
            fused=fused.tocsr(),
            weights=weights,
        )

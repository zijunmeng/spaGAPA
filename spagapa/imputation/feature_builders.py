"""
Feature builders for expression-informed spatial GP models.

This module provides reusable utilities for constructing spot-level
expression embeddings and geometry-aware spatial features.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


Orientation = Literal["genes_by_spots", "spots_by_genes"]


@dataclass
class ExpressionFeatureBuilder:
    """
    Build spot-level expression embeddings from expression matrices.

    Parameters
    ----------
    method : {"pca"}, default="pca"
        Embedding method. MVP currently supports PCA only.
    n_components : int, default=15
        Number of embedding dimensions.
    use_hvg : bool, default=False
        Whether to restrict to highly variable genes before PCA.
    n_top_genes : int, default=2000
        Number of genes to keep when `use_hvg=True`.
    scale : bool, default=True
        Whether to z-score genes before PCA.
    orientation : {"genes_by_spots", "spots_by_genes"}, default="genes_by_spots"
        Orientation of the input matrix.
    """

    method: str = "pca"
    n_components: int = 15
    use_hvg: bool = False
    n_top_genes: int = 2000
    scale: bool = True
    orientation: Orientation = "genes_by_spots"

    def _to_spot_by_gene(self, expression_matrix: pd.DataFrame | np.ndarray) -> np.ndarray:
        if isinstance(expression_matrix, pd.DataFrame):
            matrix = expression_matrix.to_numpy(dtype=float, copy=True)
        else:
            matrix = np.asarray(expression_matrix, dtype=float)

        if matrix.ndim != 2:
            raise ValueError("expression_matrix must be 2-dimensional")

        if self.orientation == "genes_by_spots":
            matrix = matrix.T
        elif self.orientation != "spots_by_genes":
            raise ValueError(f"Unknown orientation: {self.orientation}")

        return matrix

    def _select_hvg(self, matrix: np.ndarray) -> np.ndarray:
        if not self.use_hvg:
            return matrix

        variances = np.var(matrix, axis=0)
        n_keep = min(self.n_top_genes, matrix.shape[1])
        if n_keep <= 0:
            raise ValueError("n_top_genes must be positive when use_hvg=True")
        top_idx = np.argsort(variances)[::-1][:n_keep]
        return matrix[:, top_idx]

    def fit_transform(self, expression_matrix: pd.DataFrame | np.ndarray) -> np.ndarray:
        """
        Build a spot-level embedding matrix.

        Returns
        -------
        np.ndarray
            Shape `(n_spots, n_components_effective)`.
        """
        if self.method != "pca":
            raise ValueError(f"Unsupported method: {self.method}")

        matrix = self._to_spot_by_gene(expression_matrix)
        matrix = self._select_hvg(matrix)

        if self.scale:
            scaler = StandardScaler(with_mean=True, with_std=True)
            matrix = scaler.fit_transform(matrix)
            matrix = np.nan_to_num(matrix, nan=0.0)

        max_components = min(matrix.shape[0], matrix.shape[1])
        if max_components == 0:
            raise ValueError("expression_matrix has no usable rows or columns")
        n_components = min(self.n_components, max_components)
        if n_components < 1:
            raise ValueError("n_components must be at least 1")

        pca = PCA(n_components=n_components, random_state=42)
        return pca.fit_transform(matrix)


@dataclass
class GeometryFeatureBuilder:
    """
    Build geometry-aware spot features such as radius and polar angle.
    """

    use_theta: bool = True

    def build(
        self,
        coordinates: np.ndarray,
        layer_labels: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """
        Build geometry features from spatial coordinates.

        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates.
        layer_labels : np.ndarray, optional
            Optional categorical layer labels for pseudolayer generation.
        """
        coords = np.asarray(coordinates, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coordinates must have shape (n_spots, 2)")

        center = coords.mean(axis=0)
        shifted = coords - center
        radius = np.linalg.norm(shifted, axis=1)
        result = {
            "xy": coords,
            "radius": radius[:, None],
            "center": center,
        }

        if self.use_theta:
            theta = np.arctan2(shifted[:, 1], shifted[:, 0])
            result["theta"] = theta[:, None]

        if layer_labels is not None:
            labels = np.asarray(layer_labels).astype(str)
            layer_df = pd.DataFrame({"layer": labels, "radius": radius})
            ordered_layers = (
                layer_df.groupby("layer")["radius"]
                .mean()
                .sort_values()
                .index.astype(str)
                .tolist()
            )
            layer_to_order = {layer: i for i, layer in enumerate(ordered_layers)}
            pseudolayer = np.array([layer_to_order[layer] for layer in labels], dtype=float)
            result["pseudolayer"] = pseudolayer[:, None]
            result["layer_order"] = np.array(ordered_layers, dtype=object)

        return result

"""
Classical machine-learning domain recovery for BioML spot factors.
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.preprocessing import StandardScaler


class BioMLDomainDetector:
    """
    Recover biological domains from BioML spot factors or a fused graph.

    Parameters
    ----------
    method : {"kmeans", "spectral"}, default="kmeans"
        Clustering method.
    n_domains : int, default=5
        Number of domains.
    random_state : int, default=42
        Random seed for reproducibility.
    scale : bool, default=True
        Whether to z-score spot factors before KMeans.
    """

    def __init__(
        self,
        method: str = "kmeans",
        n_domains: int = 5,
        random_state: int = 42,
        scale: bool = True,
    ):
        if method not in {"kmeans", "spectral"}:
            raise ValueError(f"Unsupported domain detection method: {method}")
        if n_domains <= 0:
            raise ValueError("n_domains must be positive")
        self.method = method
        self.n_domains = int(n_domains)
        self.random_state = int(random_state)
        self.scale = bool(scale)
        self.labels_: Optional[np.ndarray] = None

    def fit_predict(
        self,
        spot_factors: Optional[np.ndarray] = None,
        graph: Optional[Union[sparse.spmatrix, np.ndarray]] = None,
    ) -> np.ndarray:
        """Cluster spots and return integer domain labels."""
        if self.method == "spectral":
            if graph is None:
                raise ValueError("graph is required for spectral domain detection")
            affinity = graph.toarray() if sparse.issparse(graph) else np.asarray(graph, dtype=float)
            if affinity.ndim != 2 or affinity.shape[0] != affinity.shape[1]:
                raise ValueError("graph must be a square affinity matrix")
            clustering = SpectralClustering(
                n_clusters=self.n_domains,
                affinity="precomputed",
                assign_labels="kmeans",
                random_state=self.random_state,
            )
            self.labels_ = clustering.fit_predict(affinity)
            return self.labels_

        if spot_factors is None:
            raise ValueError("spot_factors is required for kmeans domain detection")
        factors = np.asarray(spot_factors, dtype=float)
        if factors.ndim != 2:
            raise ValueError("spot_factors must be 2-dimensional")
        if self.scale:
            factors = StandardScaler().fit_transform(factors)
            factors = np.nan_to_num(factors, nan=0.0)
        kmeans = KMeans(n_clusters=self.n_domains, random_state=self.random_state, n_init=10)
        self.labels_ = kmeans.fit_predict(factors)
        return self.labels_

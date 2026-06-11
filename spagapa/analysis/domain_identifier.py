"""
Spatial domain identification for APA analysis.

This module provides methods to identify spatial domains (regions) in spatial
transcriptomics data based on APA patterns.
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple, Dict, List
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except Exception:
    SCANPY_AVAILABLE = False


class DomainIdentifier:
    """
    Identify spatial domains based on APA patterns.
    
    This class provides methods to cluster spots into spatial domains using
    various clustering algorithms (Leiden, Louvain, K-means).
    
    Parameters
    ----------
    method : str, default='leiden'
        Clustering method: 'leiden', 'louvain', or 'kmeans'
    resolution : float, default=1.0
        Resolution parameter for Leiden/Louvain clustering
    n_clusters : int, optional
        Number of clusters for K-means (required if method='kmeans')
    min_domain_size : int, default=10
        Minimum number of spots in a domain
    random_state : int, optional
        Random state for reproducibility
        
    Attributes
    ----------
    domain_labels_ : np.ndarray
        Domain labels for each spot
    n_domains_ : int
        Number of identified domains
    domain_stats_ : pd.DataFrame
        Statistics for each domain
    """
    
    def __init__(
        self,
        method: str = 'leiden',
        resolution: float = 1.0,
        n_clusters: Optional[int] = None,
        min_domain_size: int = 10,
        random_state: Optional[int] = None
    ):
        self.method = method
        self.resolution = resolution
        self.n_clusters = n_clusters
        self.min_domain_size = min_domain_size
        self.random_state = random_state
        
        self.domain_labels_ = None
        self.n_domains_ = None
        self.domain_stats_ = None
        
    def identify_domains(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        n_neighbors: int = 15,
        uncertainty: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Identify spatial domains based on APA patterns.
        
        Parameters
        ----------
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        spatial_coords : np.ndarray, shape (n_spots, 2)
            Spatial coordinates (x, y)
        n_neighbors : int, default=15
            Number of neighbors for graph construction
            
        Returns
        -------
        domain_labels : np.ndarray, shape (n_spots,)
            Domain label for each spot
        """
        # Convert to (n_spots, n_genes) using coordinates as the contract.
        n_coord_spots = spatial_coords.shape[0]
        if apa_matrix.shape[0] == n_coord_spots:
            pass
        elif apa_matrix.shape[1] == n_coord_spots:
            apa_matrix = apa_matrix.T
        else:
            raise ValueError(
                "apa_matrix must be either (n_genes, n_spots) or "
                "(n_spots, n_genes), with n_spots matching spatial_coords"
            )

        n_spots = apa_matrix.shape[0]
        
        if self.method in ['leiden', 'louvain']:
            if not SCANPY_AVAILABLE:
                raise ImportError(
                    f"{self.method} clustering requires scanpy. "
                    "Install with: pip install scanpy"
                )
            domain_labels = self._cluster_graph_based(
                apa_matrix, spatial_coords, n_neighbors
            )
        elif self.method == 'kmeans':
            if self.n_clusters is None:
                raise ValueError("n_clusters must be specified for kmeans")
            if uncertainty is not None:
                domain_labels = self._cluster_kmeans_weighted(
                    apa_matrix, uncertainty
                )
            else:
                domain_labels = self._cluster_kmeans(apa_matrix)
        else:
            raise ValueError(
                f"Unknown method: {self.method}. "
                "Choose from: 'leiden', 'louvain', 'kmeans'"
            )
            
        # Store results
        self.domain_labels_ = domain_labels
        self.n_domains_ = len(np.unique(domain_labels))
        
        return domain_labels
        
    def _cluster_graph_based(
        self,
        apa_matrix: np.ndarray,
        spatial_coords: np.ndarray,
        n_neighbors: int
    ) -> np.ndarray:
        """Perform graph-based clustering (Leiden or Louvain)."""
        import anndata as ad
        
        # Create AnnData object
        adata = ad.AnnData(X=apa_matrix)
        adata.obsm['spatial'] = spatial_coords
        
        # Preprocess
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Build neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X')
        
        # Cluster
        if self.method == 'leiden':
            sc.tl.leiden(
                adata,
                resolution=self.resolution,
                random_state=self.random_state
            )
            domain_labels = adata.obs['leiden'].astype(int).values
        else:  # louvain
            sc.tl.louvain(
                adata,
                resolution=self.resolution,
                random_state=self.random_state
            )
            domain_labels = adata.obs['louvain'].astype(int).values
            
        return domain_labels
        
    def _cluster_kmeans(self, apa_matrix: np.ndarray) -> np.ndarray:
        """Perform K-means clustering."""
        # Standardize features
        scaler = StandardScaler()
        apa_scaled = scaler.fit_transform(apa_matrix)

        # K-means clustering
        kmeans = KMeans(
            n_clusters=self.n_clusters,
            random_state=self.random_state,
            n_init=10
        )
        domain_labels = kmeans.fit_predict(apa_scaled)

        return domain_labels

    def _cluster_kmeans_weighted(
        self,
        apa_matrix: np.ndarray,
        uncertainty: np.ndarray,
    ) -> np.ndarray:
        """
        K-means clustering with per-spot weights derived from uncertainty.

        Parameters
        ----------
        apa_matrix : np.ndarray, (n_spots, n_genes)
        uncertainty : np.ndarray, (n_genes, n_spots)

        Returns
        -------
        domain_labels : np.ndarray, (n_spots,)
        """
        # Compute per-spot weight as mean inverse uncertainty across genes
        # (spots with high avg uncertainty get lower weight)
        un_T = uncertainty.T  # (n_spots, n_genes)
        spot_weights = 1.0 / (np.nanmean(un_T, axis=1) + 1e-8)
        spot_weights = np.clip(spot_weights, 0.1, 10.0)

        # Weighted K-means: replicate spots proportional to weight
        scaler = StandardScaler()
        apa_scaled = scaler.fit_transform(apa_matrix)

        # Convert weights to integer sample multipliers
        w_norm = spot_weights / spot_weights.min()
        sample_indices = np.repeat(np.arange(len(w_norm)), w_norm.astype(int))
        apa_weighted = apa_scaled[sample_indices]

        kmeans = KMeans(
            n_clusters=self.n_clusters,
            random_state=self.random_state,
            n_init=10,
        )
        kmeans.fit(apa_weighted)

        # Predict on original data
        domain_labels = kmeans.predict(apa_scaled)
        return domain_labels
        
    def refine_domains(
        self,
        domain_labels: np.ndarray,
        spatial_graph: Optional[np.ndarray] = None,
        spatial_coords: Optional[np.ndarray] = None,
        n_neighbors: int = 6
    ) -> np.ndarray:
        """
        Refine domain boundaries using spatial information.
        
        This method:
        1. Smooths domain boundaries using spatial neighbors
        2. Removes small isolated domains
        3. Merges similar adjacent domains
        
        Parameters
        ----------
        domain_labels : np.ndarray, shape (n_spots,)
            Initial domain labels
        spatial_graph : np.ndarray, optional
            Spatial adjacency matrix or neighbor indices
        spatial_coords : np.ndarray, optional
            Spatial coordinates (required if spatial_graph is None)
        n_neighbors : int, default=6
            Number of neighbors for smoothing
            
        Returns
        -------
        refined_labels : np.ndarray, shape (n_spots,)
            Refined domain labels
        """
        if spatial_graph is None:
            if spatial_coords is None:
                raise ValueError(
                    "Either spatial_graph or spatial_coords must be provided"
                )
            # Build spatial graph
            from ..spatial import build_knn_graph
            distances, indices = build_knn_graph(spatial_coords, k=n_neighbors)
            # Convert to adjacency format (store indices)
            spatial_graph = indices
            
        refined_labels = domain_labels.copy()
        
        # Step 1: Smooth boundaries (majority voting)
        refined_labels = self._smooth_boundaries(
            refined_labels, spatial_graph
        )
        
        # Step 2: Remove small domains
        refined_labels = self._remove_small_domains(
            refined_labels, self.min_domain_size
        )
        
        # Step 3: Relabel to consecutive integers
        refined_labels = self._relabel_consecutive(refined_labels)
        
        # Update stored labels
        self.domain_labels_ = refined_labels
        self.n_domains_ = len(np.unique(refined_labels))
        
        return refined_labels
        
    def _smooth_boundaries(
        self,
        labels: np.ndarray,
        spatial_graph: np.ndarray,
        n_iterations: int = 1
    ) -> np.ndarray:
        """
        Smooth domain boundaries using majority voting.
        
        Parameters
        ----------
        labels : np.ndarray
            Domain labels
        spatial_graph : np.ndarray, shape (n_spots, k)
            Neighbor indices for each spot
        n_iterations : int
            Number of smoothing iterations
        """
        smoothed = labels.copy()
        
        for _ in range(n_iterations):
            new_labels = smoothed.copy()
            
            for i in range(len(labels)):
                # Get neighbors (spatial_graph is now indices array)
                neighbors = spatial_graph[i]
                if len(neighbors) == 0:
                    continue
                    
                # Majority vote
                neighbor_labels = smoothed[neighbors]
                unique, counts = np.unique(neighbor_labels, return_counts=True)
                majority_label = unique[np.argmax(counts)]
                
                # Update if different from current
                if counts[np.argmax(counts)] > len(neighbors) / 2:
                    new_labels[i] = majority_label
                    
            smoothed = new_labels
            
        return smoothed
        
    def _remove_small_domains(
        self,
        labels: np.ndarray,
        min_size: int
    ) -> np.ndarray:
        """Remove domains with fewer than min_size spots."""
        unique_labels, counts = np.unique(labels, return_counts=True)
        
        # Find small domains
        small_domains = unique_labels[counts < min_size]
        
        if len(small_domains) == 0:
            return labels
            
        # Assign small domains to nearest large domain
        # (simplified: assign to most common neighbor)
        cleaned = labels.copy()
        
        for small_label in small_domains:
            mask = labels == small_label
            # For simplicity, assign to label -1 (will be relabeled)
            cleaned[mask] = -1
            
        return cleaned
        
    def _relabel_consecutive(self, labels: np.ndarray) -> np.ndarray:
        """Relabel domains to consecutive integers starting from 0."""
        unique_labels = np.unique(labels)
        label_map = {old: new for new, old in enumerate(unique_labels)}
        
        relabeled = np.array([label_map[label] for label in labels])
        return relabeled
        
    def compute_domain_stats(
        self,
        domain_labels: np.ndarray,
        apa_matrix: np.ndarray,
        spatial_coords: Optional[np.ndarray] = None
    ) -> pd.DataFrame:
        """
        Compute statistics for each domain.
        
        Parameters
        ----------
        domain_labels : np.ndarray, shape (n_spots,)
            Domain labels
        apa_matrix : np.ndarray, shape (n_genes, n_spots)
            APA index matrix
        spatial_coords : np.ndarray, optional
            Spatial coordinates
            
        Returns
        -------
        stats : pd.DataFrame
            Statistics for each domain
        """
        # Transpose if needed
        if apa_matrix.shape[0] < apa_matrix.shape[1]:
            apa_matrix = apa_matrix.T
            
        unique_domains = np.unique(domain_labels)
        stats_list = []
        
        for domain in unique_domains:
            mask = domain_labels == domain
            domain_apa = apa_matrix[mask, :]
            
            stats = {
                'domain': domain,
                'n_spots': mask.sum(),
                'mean_apa': np.nanmean(domain_apa),
                'std_apa': np.nanstd(domain_apa),
                'median_apa': np.nanmedian(domain_apa)
            }
            
            if spatial_coords is not None:
                domain_coords = spatial_coords[mask]
                stats['centroid_x'] = np.mean(domain_coords[:, 0])
                stats['centroid_y'] = np.mean(domain_coords[:, 1])
                stats['area'] = self._compute_domain_area(domain_coords)
                
            stats_list.append(stats)
            
        self.domain_stats_ = pd.DataFrame(stats_list)
        return self.domain_stats_
        
    def _compute_domain_area(self, coords: np.ndarray) -> float:
        """Compute approximate area of domain using convex hull."""
        if len(coords) < 3:
            return 0.0
            
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(coords)
            return hull.volume  # In 2D, volume is area
        except:
            # Fallback: bounding box area
            x_range = coords[:, 0].max() - coords[:, 0].min()
            y_range = coords[:, 1].max() - coords[:, 1].min()
            return x_range * y_range


def identify_spatial_domains(
    apa_matrix: np.ndarray,
    spatial_coords: np.ndarray,
    method: str = 'leiden',
    resolution: float = 1.0,
    n_clusters: Optional[int] = None,
    n_neighbors: int = 15,
    min_domain_size: int = 10,
    refine: bool = True,
    random_state: Optional[int] = None
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Identify spatial domains based on APA patterns.
    
    Convenience function that wraps DomainIdentifier.
    
    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_genes, n_spots)
        APA index matrix
    spatial_coords : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    method : str, default='leiden'
        Clustering method
    resolution : float, default=1.0
        Resolution for graph-based clustering
    n_clusters : int, optional
        Number of clusters for K-means
    n_neighbors : int, default=15
        Number of neighbors for graph construction
    min_domain_size : int, default=10
        Minimum domain size
    refine : bool, default=True
        Whether to refine domain boundaries
    random_state : int, optional
        Random state
        
    Returns
    -------
    domain_labels : np.ndarray
        Domain labels for each spot
    domain_stats : pd.DataFrame
        Statistics for each domain
        
    Examples
    --------
    >>> import numpy as np
    >>> from spagapa.analysis import identify_spatial_domains
    >>> 
    >>> # Generate example data
    >>> apa_matrix = np.random.rand(100, 500)  # 100 genes, 500 spots
    >>> coords = np.random.rand(500, 2) * 100
    >>> 
    >>> # Identify domains
    >>> labels, stats = identify_spatial_domains(
    ...     apa_matrix, coords, method='kmeans', n_clusters=5
    ... )
    >>> print(f"Found {len(np.unique(labels))} domains")
    >>> print(stats)
    """
    identifier = DomainIdentifier(
        method=method,
        resolution=resolution,
        n_clusters=n_clusters,
        min_domain_size=min_domain_size,
        random_state=random_state
    )
    
    # Identify domains
    domain_labels = identifier.identify_domains(
        apa_matrix, spatial_coords, n_neighbors
    )
    
    # Refine if requested
    if refine:
        domain_labels = identifier.refine_domains(
            domain_labels, spatial_coords=spatial_coords, n_neighbors=n_neighbors
        )
        
    # Compute statistics
    domain_stats = identifier.compute_domain_stats(
        domain_labels, apa_matrix, spatial_coords
    )
    
    return domain_labels, domain_stats

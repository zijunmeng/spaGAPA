"""
Spatial validation of APA sites.

This module implements spatial validation algorithms that use neighborhood
information to validate APA sites identified from sequencing data.
"""

import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict
from spagapa.spatial.neighbors import SpatialNeighbors
from spagapa.core import APASite, APASiteCollection
import logging

logger = logging.getLogger(__name__)


class SpatialValidator:
    """
    Validate APA sites using spatial neighborhood information.
    
    The key idea is that true APA sites should be consistently detected
    across spatially neighboring spots, while technical noise tends to
    be randomly distributed.
    
    Parameters
    ----------
    n_neighbors : int, default=6
        Number of spatial neighbors to consider
    support_threshold : float, default=0.3
        Minimum fraction of neighbors that must support a site
    method : str, default='knn'
        Method for finding neighbors: 'knn', 'radius', or 'delaunay'
    radius : float, optional
        Radius for radius-based neighbor finding
    
    Examples
    --------
    >>> validator = SpatialValidator(n_neighbors=6, support_threshold=0.3)
    >>> validator.fit(spatial_coords)
    >>> validated_sites = validator.validate_sites(candidate_sites, apa_counts)
    """
    
    def __init__(
        self,
        n_neighbors: int = 6,
        support_threshold: float = 0.3,
        method: str = 'knn',
        radius: Optional[float] = None
    ):
        self.n_neighbors = n_neighbors
        self.support_threshold = support_threshold
        self.method = method
        self.radius = radius
        self.spatial_neighbors = None
        self._neighbor_indices = None
        self._neighbor_distances = None
    
    def fit(self, coordinates: np.ndarray):
        """
        Fit the validator to spatial coordinates.
        
        Parameters
        ----------
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates (x, y)
        """
        self.spatial_neighbors = SpatialNeighbors(
            method=self.method,
            n_neighbors=self.n_neighbors,
            radius=self.radius
        )
        self.spatial_neighbors.fit(coordinates)
        
        # Pre-compute neighbors
        self._neighbor_distances, self._neighbor_indices = \
            self.spatial_neighbors.find_neighbors(coordinates)
        
        logger.info(
            f"Fitted spatial validator with {self.n_neighbors} neighbors, "
            f"threshold={self.support_threshold}"
        )
        
        return self
    
    def compute_spatial_support(
        self,
        apa_counts: np.ndarray,
        gene_idx: Optional[int] = None
    ) -> np.ndarray:
        """
        Compute spatial support score for each spot.
        
        The support score is the fraction of neighbors that also have
        non-zero counts for the same APA site.
        
        Parameters
        ----------
        apa_counts : np.ndarray, shape (n_genes, n_spots) or (n_spots,)
            APA count matrix or counts for single gene
        gene_idx : int, optional
            Gene index if apa_counts is 2D
        
        Returns
        -------
        np.ndarray, shape (n_spots,)
            Spatial support scores (0-1)
        """
        if apa_counts.ndim == 2:
            if gene_idx is None:
                raise ValueError("gene_idx required for 2D apa_counts")
            counts = apa_counts[gene_idx, :]
        else:
            counts = apa_counts
        
        n_spots = len(counts)
        support_scores = np.zeros(n_spots)
        
        # For each spot, check how many neighbors have non-zero counts
        for i in range(n_spots):
            if counts[i] == 0:
                support_scores[i] = 0.0
                continue
            
            neighbors = self._neighbor_indices[i]
            if len(neighbors) == 0:
                support_scores[i] = 0.0
                continue
            
            # Count neighbors with non-zero counts
            neighbor_counts = counts[neighbors]
            n_supporting = np.sum(neighbor_counts > 0)
            
            # Support score = fraction of neighbors with signal
            support_scores[i] = n_supporting / len(neighbors)
        
        return support_scores
    
    def compute_weighted_support(
        self,
        apa_counts: np.ndarray,
        gene_idx: Optional[int] = None
    ) -> np.ndarray:
        """
        Compute distance-weighted spatial support score.
        
        Closer neighbors have higher weight in the support calculation.
        
        Parameters
        ----------
        apa_counts : np.ndarray
            APA count matrix
        gene_idx : int, optional
            Gene index if apa_counts is 2D
        
        Returns
        -------
        np.ndarray
            Weighted support scores
        """
        if apa_counts.ndim == 2:
            if gene_idx is None:
                raise ValueError("gene_idx required for 2D apa_counts")
            counts = apa_counts[gene_idx, :]
        else:
            counts = apa_counts
        
        n_spots = len(counts)
        support_scores = np.zeros(n_spots)
        
        for i in range(n_spots):
            if counts[i] == 0:
                continue
            
            neighbors = self._neighbor_indices[i]
            distances = self._neighbor_distances[i]
            
            if len(neighbors) == 0:
                continue
            
            # Compute weights (inverse distance)
            weights = 1.0 / (distances + 1e-6)
            weights = weights / weights.sum()  # Normalize
            
            # Weighted support
            neighbor_counts = counts[neighbors]
            neighbor_support = (neighbor_counts > 0).astype(float)
            support_scores[i] = np.sum(weights * neighbor_support)
        
        return support_scores
    
    def validate_sites(
        self,
        candidate_sites: APASiteCollection,
        apa_counts: np.ndarray,
        use_weighted: bool = False
    ) -> Tuple[APASiteCollection, np.ndarray]:
        """
        Validate candidate APA sites using spatial information.
        
        Parameters
        ----------
        candidate_sites : APASiteCollection
            Collection of candidate APA sites
        apa_counts : np.ndarray, shape (n_genes, n_spots)
            APA count matrix
        use_weighted : bool, default=False
            Use distance-weighted support scores
        
        Returns
        -------
        validated_sites : APASiteCollection
            Sites that pass spatial validation
        support_scores : np.ndarray
            Support scores for all candidate sites
        """
        n_sites = len(candidate_sites)
        support_scores = np.zeros(n_sites)
        
        # Compute support for each site
        for i, site in enumerate(candidate_sites):
            # Find corresponding gene in count matrix
            # This assumes sites are ordered same as genes
            # In practice, you'd need to match by gene_id
            if use_weighted:
                support = self.compute_weighted_support(apa_counts, gene_idx=i)
            else:
                support = self.compute_spatial_support(apa_counts, gene_idx=i)
            
            # Average support across all spots
            support_scores[i] = np.mean(support[support > 0]) if np.any(support > 0) else 0.0
        
        # Filter sites by support threshold
        valid_mask = support_scores >= self.support_threshold
        
        validated_sites = APASiteCollection([
            site for site, valid in zip(candidate_sites, valid_mask) if valid
        ])
        
        logger.info(
            f"Validated {len(validated_sites)}/{n_sites} sites "
            f"(threshold={self.support_threshold})"
        )
        
        return validated_sites, support_scores
    
    def filter_by_spatial_consistency(
        self,
        apa_counts: np.ndarray,
        min_support: float = 0.3
    ) -> np.ndarray:
        """
        Filter APA count matrix by spatial consistency.
        
        Sets counts to zero for spots with low spatial support.
        
        Parameters
        ----------
        apa_counts : np.ndarray, shape (n_genes, n_spots)
            APA count matrix
        min_support : float, default=0.3
            Minimum support threshold
        
        Returns
        -------
        np.ndarray
            Filtered count matrix
        """
        n_genes, n_spots = apa_counts.shape
        filtered_counts = apa_counts.copy()
        
        for gene_idx in range(n_genes):
            support = self.compute_spatial_support(apa_counts, gene_idx)
            
            # Set low-support spots to zero
            low_support_mask = support < min_support
            filtered_counts[gene_idx, low_support_mask] = 0
        
        # Count how many values were filtered
        n_filtered = np.sum((apa_counts > 0) & (filtered_counts == 0))
        total_nonzero = np.sum(apa_counts > 0)
        
        logger.info(
            f"Filtered {n_filtered}/{total_nonzero} "
            f"({100*n_filtered/total_nonzero:.1f}%) low-support counts"
        )
        
        return filtered_counts
    
    def compute_spatial_autocorrelation(
        self,
        apa_counts: np.ndarray,
        gene_idx: int
    ) -> float:
        """
        Compute Moran's I spatial autocorrelation statistic.
        
        Parameters
        ----------
        apa_counts : np.ndarray
            APA count matrix
        gene_idx : int
            Gene index
        
        Returns
        -------
        float
            Moran's I statistic
        """
        counts = apa_counts[gene_idx, :]
        n = len(counts)
        
        # Compute spatial weights
        W = self.spatial_neighbors.compute_spatial_weights(
            self._neighbor_indices,
            weight_type='uniform'
        )
        
        # Moran's I formula
        mean_counts = np.mean(counts)
        numerator = 0
        denominator = np.sum((counts - mean_counts)**2)
        
        for i in range(n):
            for j in range(n):
                if W[i, j] > 0:
                    numerator += W[i, j] * (counts[i] - mean_counts) * (counts[j] - mean_counts)
        
        W_sum = W.sum()
        morans_i = (n / W_sum) * (numerator / denominator)
        
        return morans_i


def validate_apa_sites_spatial(
    candidate_sites: APASiteCollection,
    apa_counts: np.ndarray,
    spatial_coords: np.ndarray,
    n_neighbors: int = 6,
    support_threshold: float = 0.3
) -> Tuple[APASiteCollection, np.ndarray]:
    """
    Convenience function for spatial validation of APA sites.
    
    Parameters
    ----------
    candidate_sites : APASiteCollection
        Candidate APA sites
    apa_counts : np.ndarray
        APA count matrix
    spatial_coords : np.ndarray
        Spatial coordinates
    n_neighbors : int, default=6
        Number of neighbors
    support_threshold : float, default=0.3
        Support threshold
    
    Returns
    -------
    validated_sites : APASiteCollection
        Validated sites
    support_scores : np.ndarray
        Support scores
    """
    validator = SpatialValidator(
        n_neighbors=n_neighbors,
        support_threshold=support_threshold
    )
    validator.fit(spatial_coords)
    
    return validator.validate_sites(candidate_sites, apa_counts)

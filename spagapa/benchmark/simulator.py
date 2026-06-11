"""
Spatial APA Simulation Framework

This module provides tools for generating synthetic spatial transcriptomics
data with known APA patterns for benchmarking and validation.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans


class SpatialAPASimulator:
    """
    Generate synthetic spatial transcriptomics data with APA patterns.
    
    This simulator creates realistic spatial data with:
    - Spatial domains with distinct APA patterns
    - Dropout and technical noise
    - Spatial autocorrelation
    - Multiple genes with different patterns
    
    Parameters
    ----------
    n_spots : int
        Number of spatial spots to generate
    n_genes : int
        Number of genes to simulate
    n_domains : int
        Number of spatial domains
    spatial_pattern : str
        Type of spatial pattern ('domains', 'gradient', 'random')
    dropout_rate : float
        Proportion of values to set as missing (0-1)
    noise_level : float
        Standard deviation of Gaussian noise
    random_state : int, optional
        Random seed for reproducibility
    
    Attributes
    ----------
    coordinates : np.ndarray
        Spatial coordinates (n_spots, 2)
    domain_labels : np.ndarray
        Domain assignment for each spot (n_spots,)
    true_apa : np.ndarray
        Ground truth APA values (n_spots, n_genes)
    observed_apa : np.ndarray
        Observed APA with dropout and noise (n_spots, n_genes)
    """
    
    def __init__(
        self,
        n_spots: int = 500,
        n_genes: int = 100,
        n_domains: int = 3,
        spatial_pattern: str = 'domains',
        dropout_rate: float = 0.3,
        noise_level: float = 0.1,
        random_state: Optional[int] = None
    ):
        self.n_spots = n_spots
        self.n_genes = n_genes
        self.n_domains = n_domains
        self.spatial_pattern = spatial_pattern
        self.dropout_rate = dropout_rate
        self.noise_level = noise_level
        self.random_state = random_state
        
        if random_state is not None:
            np.random.seed(random_state)
        
        # Initialize attributes
        self.coordinates = None
        self.domain_labels = None
        self.true_apa = None
        self.observed_apa = None
        self.gene_patterns = None
    
    def generate(self) -> Dict[str, np.ndarray]:
        """
        Generate complete synthetic dataset.
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'coordinates': Spatial coordinates
            - 'domain_labels': Domain assignments
            - 'true_apa': Ground truth APA values
            - 'observed_apa': Observed APA with noise
            - 'gene_patterns': Pattern type for each gene
        """
        # Generate spatial coordinates
        self.coordinates = self._generate_coordinates()
        
        # Generate spatial domains
        self.domain_labels = self._generate_domains()
        
        # Generate APA patterns
        self.true_apa, self.gene_patterns = self._generate_apa_patterns()
        
        # Add dropout and noise
        self.observed_apa = self._add_noise_and_dropout()
        
        return {
            'coordinates': self.coordinates,
            'domain_labels': self.domain_labels,
            'true_apa': self.true_apa,
            'observed_apa': self.observed_apa,
            'gene_patterns': self.gene_patterns
        }
    
    def _generate_coordinates(self) -> np.ndarray:
        """Generate spatial coordinates."""
        if self.spatial_pattern == 'grid':
            # Regular grid
            side = int(np.ceil(np.sqrt(self.n_spots)))
            x = np.repeat(np.arange(side), side)[:self.n_spots]
            y = np.tile(np.arange(side), side)[:self.n_spots]
            coords = np.column_stack([x, y])
        else:
            # Random or clustered coordinates
            coords = np.random.rand(self.n_spots, 2) * 100
        
        return coords
    
    def _generate_domains(self) -> np.ndarray:
        """Generate spatial domain labels."""
        if self.spatial_pattern == 'domains':
            # Use K-means clustering on coordinates
            kmeans = KMeans(n_clusters=self.n_domains, random_state=self.random_state)
            labels = kmeans.fit_predict(self.coordinates)
        elif self.spatial_pattern == 'gradient':
            # Create gradient-based domains
            x_norm = (self.coordinates[:, 0] - self.coordinates[:, 0].min()) / \
                     (self.coordinates[:, 0].max() - self.coordinates[:, 0].min())
            labels = (x_norm * self.n_domains).astype(int)
            labels = np.clip(labels, 0, self.n_domains - 1)
        else:
            # Random assignment
            labels = np.random.randint(0, self.n_domains, self.n_spots)
        
        return labels

    
    def _generate_apa_patterns(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate APA patterns for genes.
        
        Returns
        -------
        apa_matrix : np.ndarray
            APA values (n_spots, n_genes)
        gene_patterns : np.ndarray
            Pattern type for each gene
        """
        apa_matrix = np.zeros((self.n_spots, self.n_genes))
        gene_patterns = np.zeros(self.n_genes, dtype=int)
        
        # Assign pattern types to genes
        # 0: domain-specific, 1: gradient, 2: random, 3: spatially autocorrelated
        pattern_types = np.random.choice([0, 1, 2, 3], size=self.n_genes, 
                                        p=[0.4, 0.3, 0.2, 0.1])
        
        for gene_idx in range(self.n_genes):
            pattern_type = pattern_types[gene_idx]
            gene_patterns[gene_idx] = pattern_type
            
            if pattern_type == 0:
                # Domain-specific pattern
                apa_matrix[:, gene_idx] = self._domain_specific_pattern()
            elif pattern_type == 1:
                # Gradient pattern
                apa_matrix[:, gene_idx] = self._gradient_pattern()
            elif pattern_type == 2:
                # Random pattern
                apa_matrix[:, gene_idx] = np.random.beta(2, 2, self.n_spots)
            else:
                # Spatially autocorrelated pattern
                apa_matrix[:, gene_idx] = self._spatially_autocorrelated_pattern()
        
        return apa_matrix, gene_patterns
    
    def _domain_specific_pattern(self) -> np.ndarray:
        """Generate domain-specific APA pattern."""
        # Each domain has a different mean APA value
        domain_means = np.random.beta(2, 2, self.n_domains)
        apa_values = domain_means[self.domain_labels]
        
        # Add within-domain variation
        apa_values += np.random.normal(0, 0.05, self.n_spots)
        apa_values = np.clip(apa_values, 0, 1)
        
        return apa_values
    
    def _gradient_pattern(self) -> np.ndarray:
        """Generate gradient APA pattern."""
        # Linear gradient along x-axis
        x_norm = (self.coordinates[:, 0] - self.coordinates[:, 0].min()) / \
                 (self.coordinates[:, 0].max() - self.coordinates[:, 0].min())
        
        # Add some curvature
        apa_values = 0.2 + 0.6 * x_norm + 0.2 * np.sin(x_norm * np.pi)
        apa_values = np.clip(apa_values, 0, 1)
        
        return apa_values
    
    def _spatially_autocorrelated_pattern(self) -> np.ndarray:
        """Generate spatially autocorrelated pattern using Gaussian process."""
        # Simple spatial autocorrelation using distance-based smoothing
        n_centers = 10
        center_indices = np.random.choice(self.n_spots, n_centers, replace=False)
        center_coords = self.coordinates[center_indices]
        center_values = np.random.beta(2, 2, n_centers)
        
        # Compute distances to centers
        distances = cdist(self.coordinates, center_coords)
        
        # Weight by inverse distance
        weights = 1.0 / (distances + 1.0)
        weights = weights / weights.sum(axis=1, keepdims=True)
        
        # Weighted average of center values
        apa_values = (weights @ center_values)
        apa_values = np.clip(apa_values, 0, 1)
        
        return apa_values
    
    def _add_noise_and_dropout(self) -> np.ndarray:
        """Add technical noise and dropout to true APA values."""
        observed = self.true_apa.copy()
        
        # Add Gaussian noise
        noise = np.random.normal(0, self.noise_level, observed.shape)
        observed += noise
        observed = np.clip(observed, 0, 1)
        
        # Add dropout
        dropout_mask = np.random.rand(*observed.shape) < self.dropout_rate
        observed[dropout_mask] = np.nan
        
        return observed
    
    def get_differential_genes(self, domain1: int, domain2: int, 
                              threshold: float = 0.2) -> np.ndarray:
        """
        Get genes with differential APA between two domains.
        
        Parameters
        ----------
        domain1 : int
            First domain ID
        domain2 : int
            Second domain ID
        threshold : float
            Minimum mean difference to consider differential
        
        Returns
        -------
        np.ndarray
            Boolean array indicating differential genes
        """
        mask1 = self.domain_labels == domain1
        mask2 = self.domain_labels == domain2
        
        mean1 = np.nanmean(self.true_apa[mask1], axis=0)
        mean2 = np.nanmean(self.true_apa[mask2], axis=0)
        
        diff = np.abs(mean1 - mean2)
        return diff >= threshold
    
    def get_svapa_genes(self, threshold: float = 0.1) -> np.ndarray:
        """
        Get spatially variable APA genes (ground truth).
        
        Parameters
        ----------
        threshold : float
            Minimum spatial variance to consider SVAPA
        
        Returns
        -------
        np.ndarray
            Boolean array indicating SVAPA genes
        """
        # Genes with domain-specific or gradient patterns are SVAPA
        return (self.gene_patterns == 0) | (self.gene_patterns == 1) | (self.gene_patterns == 3)
    
    def to_dataframe(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Convert simulation results to DataFrames.
        
        Returns
        -------
        apa_df : pd.DataFrame
            APA matrix (spots × genes)
        metadata_df : pd.DataFrame
            Spot metadata (coordinates, domains)
        """
        # APA matrix
        gene_names = [f"Gene_{i}" for i in range(self.n_genes)]
        spot_names = [f"Spot_{i}" for i in range(self.n_spots)]
        apa_df = pd.DataFrame(
            self.observed_apa,
            index=spot_names,
            columns=gene_names
        )
        
        # Metadata
        metadata_df = pd.DataFrame({
            'spot_id': spot_names,
            'x': self.coordinates[:, 0],
            'y': self.coordinates[:, 1],
            'domain': self.domain_labels
        })
        
        return apa_df, metadata_df


def simulate_spatial_apa(
    n_spots: int = 500,
    n_genes: int = 100,
    n_domains: int = 3,
    spatial_pattern: str = 'domains',
    dropout_rate: float = 0.3,
    noise_level: float = 0.1,
    random_state: Optional[int] = None
) -> Dict[str, np.ndarray]:
    """
    Convenience function to generate synthetic spatial APA data.
    
    Parameters
    ----------
    n_spots : int
        Number of spatial spots
    n_genes : int
        Number of genes
    n_domains : int
        Number of spatial domains
    spatial_pattern : str
        Type of spatial pattern ('domains', 'gradient', 'random')
    dropout_rate : float
        Proportion of missing values
    noise_level : float
        Standard deviation of noise
    random_state : int, optional
        Random seed
    
    Returns
    -------
    dict
        Simulation results
    
    Examples
    --------
    >>> data = simulate_spatial_apa(n_spots=200, n_genes=50, random_state=42)
    >>> print(data['observed_apa'].shape)
    (200, 50)
    """
    simulator = SpatialAPASimulator(
        n_spots=n_spots,
        n_genes=n_genes,
        n_domains=n_domains,
        spatial_pattern=spatial_pattern,
        dropout_rate=dropout_rate,
        noise_level=noise_level,
        random_state=random_state
    )
    return simulator.generate()

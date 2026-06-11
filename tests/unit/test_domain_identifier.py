"""
Unit tests for domain identification.
"""

import pytest
import numpy as np
import pandas as pd
from spagapa.analysis import DomainIdentifier, identify_spatial_domains


class TestDomainIdentifier:
    """Test DomainIdentifier class."""
    
    def test_init(self):
        """Test initialization."""
        identifier = DomainIdentifier(
            method='kmeans',
            n_clusters=5,
            min_domain_size=10
        )
        
        assert identifier.method == 'kmeans'
        assert identifier.n_clusters == 5
        assert identifier.min_domain_size == 10
        
    def test_kmeans_clustering(self):
        """Test K-means clustering."""
        # Generate data with clear clusters
        np.random.seed(42)
        n_spots = 200
        n_genes = 50
        
        # Create 4 clusters
        apa_matrix = np.zeros((n_genes, n_spots))
        for i in range(4):
            start = i * 50
            end = (i + 1) * 50
            apa_matrix[:, start:end] = np.random.rand(n_genes, 50) + i
            
        coords = np.random.rand(n_spots, 2) * 100
        
        identifier = DomainIdentifier(method='kmeans', n_clusters=4)
        labels = identifier.identify_domains(apa_matrix, coords)
        
        assert len(labels) == n_spots
        assert len(np.unique(labels)) == 4
        assert identifier.n_domains_ == 4
        
    def test_refine_domains(self):
        """Test domain refinement."""
        np.random.seed(42)
        n_spots = 100
        
        # Create initial labels with some noise
        labels = np.repeat([0, 1, 2, 3], 25)
        # Add some noise
        labels[10] = 2  # Outlier in domain 0
        labels[35] = 0  # Outlier in domain 1
        
        coords = np.random.rand(n_spots, 2) * 100
        
        identifier = DomainIdentifier(min_domain_size=5)
        refined = identifier.refine_domains(labels, spatial_coords=coords)
        
        assert len(refined) == n_spots
        assert len(np.unique(refined)) <= 4
        
    def test_compute_domain_stats(self):
        """Test domain statistics computation."""
        np.random.seed(42)
        n_spots = 100
        n_genes = 20
        
        labels = np.repeat([0, 1, 2], [30, 40, 30])
        apa_matrix = np.random.rand(n_genes, n_spots)
        coords = np.random.rand(n_spots, 2) * 100
        
        identifier = DomainIdentifier()
        stats = identifier.compute_domain_stats(labels, apa_matrix, coords)
        
        assert isinstance(stats, pd.DataFrame)
        assert len(stats) == 3
        assert 'domain' in stats.columns
        assert 'n_spots' in stats.columns
        assert 'mean_apa' in stats.columns
        assert 'centroid_x' in stats.columns
        
    def test_small_domain_removal(self):
        """Test removal of small domains."""
        np.random.seed(42)
        
        # Create labels with one small domain
        labels = np.array([0] * 50 + [1] * 45 + [2] * 5)  # Domain 2 is small
        
        identifier = DomainIdentifier(min_domain_size=10)
        cleaned = identifier._remove_small_domains(labels, min_size=10)
        
        # Small domain should be marked for removal
        assert np.sum(cleaned == 2) < 5 or np.sum(cleaned == -1) > 0
        
    def test_relabel_consecutive(self):
        """Test relabeling to consecutive integers."""
        labels = np.array([5, 5, 10, 10, 20, 20])
        
        identifier = DomainIdentifier()
        relabeled = identifier._relabel_consecutive(labels)
        
        assert set(relabeled) == {0, 1, 2}
        assert len(np.unique(relabeled)) == 3
        
    def test_invalid_method(self):
        """Test invalid clustering method."""
        identifier = DomainIdentifier(method='invalid')
        
        apa_matrix = np.random.rand(10, 50)
        coords = np.random.rand(50, 2)
        
        with pytest.raises(ValueError, match="Unknown method"):
            identifier.identify_domains(apa_matrix, coords)
            
    def test_kmeans_without_n_clusters(self):
        """Test K-means without specifying n_clusters."""
        identifier = DomainIdentifier(method='kmeans')
        
        apa_matrix = np.random.rand(10, 50)
        coords = np.random.rand(50, 2)
        
        with pytest.raises(ValueError, match="n_clusters must be specified"):
            identifier.identify_domains(apa_matrix, coords)


class TestIdentifySpatialDomains:
    """Test convenience function."""
    
    def test_basic_usage(self):
        """Test basic usage."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        coords = np.random.rand(100, 2) * 100
        
        labels, stats = identify_spatial_domains(
            apa_matrix, coords, method='kmeans', n_clusters=3
        )
        
        assert len(labels) == 100
        assert isinstance(stats, pd.DataFrame)
        assert len(stats) <= 3
        
    def test_with_refinement(self):
        """Test with domain refinement."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        coords = np.random.rand(100, 2) * 100
        
        labels, stats = identify_spatial_domains(
            apa_matrix, coords,
            method='kmeans',
            n_clusters=4,
            refine=True,
            min_domain_size=10
        )
        
        assert len(labels) == 100
        # After refinement, small domains should be removed
        for domain in np.unique(labels):
            assert np.sum(labels == domain) >= 10 or np.sum(labels == domain) == 0
            
    def test_without_refinement(self):
        """Test without domain refinement."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        coords = np.random.rand(100, 2) * 100
        
        labels, stats = identify_spatial_domains(
            apa_matrix, coords,
            method='kmeans',
            n_clusters=3,
            refine=False
        )
        
        assert len(labels) == 100
        assert len(np.unique(labels)) == 3

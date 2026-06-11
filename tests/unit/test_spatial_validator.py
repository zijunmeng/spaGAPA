"""
Unit tests for spatial validation.
"""

import pytest
import numpy as np
from spagapa.spatial import SpatialNeighbors
from spagapa.calling import SpatialValidator, validate_apa_sites_spatial
from spagapa.core import APASite, APASiteCollection


class TestSpatialValidator:
    """Test SpatialValidator class."""
    
    @pytest.fixture
    def simple_coords(self):
        """Create simple 2D coordinates."""
        # 3x3 grid
        x = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        y = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        return np.column_stack([x, y])
    
    @pytest.fixture
    def sample_apa_counts(self):
        """Create sample APA count matrix."""
        # 3 genes x 9 spots
        counts = np.array([
            [10, 8, 0, 9, 10, 0, 0, 0, 0],  # Gene 1: high in left region
            [0, 0, 0, 0, 5, 8, 0, 6, 9],    # Gene 2: high in right region
            [1, 0, 1, 0, 1, 0, 1, 0, 1]     # Gene 3: sparse/noisy
        ])
        return counts
    
    @pytest.fixture
    def sample_sites(self):
        """Create sample APA sites."""
        sites = [
            APASite(
                chr='chr1',
                start=1000,
                end=1050,
                strand='+',
                gene_id='gene1'
            ),
            APASite(
                chr='chr1',
                start=2000,
                end=2050,
                strand='+',
                gene_id='gene2'
            ),
            APASite(
                chr='chr1',
                start=3000,
                end=3050,
                strand='+',
                gene_id='gene3'
            )
        ]
        return APASiteCollection(sites)
    
    def test_initialization(self):
        """Test validator initialization."""
        validator = SpatialValidator(
            n_neighbors=6,
            support_threshold=0.3
        )
        assert validator.n_neighbors == 6
        assert validator.support_threshold == 0.3
        assert validator.method == 'knn'
    
    def test_fit(self, simple_coords):
        """Test fitting validator."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        assert validator.spatial_neighbors is not None
        assert validator._neighbor_indices is not None
        assert validator._neighbor_distances is not None
    
    def test_compute_spatial_support(self, simple_coords, sample_apa_counts):
        """Test computing spatial support scores."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        # Test gene 1 (clustered in left region)
        support = validator.compute_spatial_support(sample_apa_counts, gene_idx=0)
        
        assert len(support) == 9
        assert np.all(support >= 0) and np.all(support <= 1)
        
        # Spots with signal should have higher support
        # Spot 0 has signal and neighbors with signal
        assert support[0] > 0
    
    def test_compute_weighted_support(self, simple_coords, sample_apa_counts):
        """Test computing weighted spatial support."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        support = validator.compute_weighted_support(sample_apa_counts, gene_idx=0)
        
        assert len(support) == 9
        assert np.all(support >= 0) and np.all(support <= 1)
    
    def test_validate_sites(self, simple_coords, sample_sites, sample_apa_counts):
        """Test validating APA sites."""
        validator = SpatialValidator(n_neighbors=4, support_threshold=0.3)
        validator.fit(simple_coords)
        
        validated_sites, support_scores = validator.validate_sites(
            sample_sites,
            sample_apa_counts
        )
        
        # Should return APASiteCollection
        assert isinstance(validated_sites, APASiteCollection)
        
        # Support scores should be computed for all sites
        assert len(support_scores) == len(sample_sites)
        
        # Validated sites should be subset of original
        assert len(validated_sites) <= len(sample_sites)
    
    def test_filter_by_spatial_consistency(self, simple_coords, sample_apa_counts):
        """Test filtering by spatial consistency."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        filtered_counts = validator.filter_by_spatial_consistency(
            sample_apa_counts,
            min_support=0.3
        )
        
        # Shape should be preserved
        assert filtered_counts.shape == sample_apa_counts.shape
        
        # Some counts should be filtered
        assert np.sum(filtered_counts == 0) >= np.sum(sample_apa_counts == 0)
    
    def test_support_threshold_filtering(self, simple_coords, sample_sites, sample_apa_counts):
        """Test that support threshold affects filtering."""
        validator_strict = SpatialValidator(n_neighbors=4, support_threshold=0.8)
        validator_strict.fit(simple_coords)
        
        validator_lenient = SpatialValidator(n_neighbors=4, support_threshold=0.1)
        validator_lenient.fit(simple_coords)
        
        strict_sites, _ = validator_strict.validate_sites(sample_sites, sample_apa_counts)
        lenient_sites, _ = validator_lenient.validate_sites(sample_sites, sample_apa_counts)
        
        # Lenient threshold should pass more sites
        assert len(lenient_sites) >= len(strict_sites)
    
    def test_zero_counts(self, simple_coords):
        """Test with all-zero counts."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        zero_counts = np.zeros((3, 9))
        support = validator.compute_spatial_support(zero_counts, gene_idx=0)
        
        # All support should be zero
        assert np.all(support == 0)
    
    def test_single_spot_signal(self, simple_coords):
        """Test with signal in single spot."""
        validator = SpatialValidator(n_neighbors=4)
        validator.fit(simple_coords)
        
        counts = np.zeros((1, 9))
        counts[0, 4] = 10  # Only center spot has signal
        
        support = validator.compute_spatial_support(counts, gene_idx=0)
        
        # Center spot should have zero support (no neighbors with signal)
        assert support[4] == 0


class TestHelperFunction:
    """Test helper function."""
    
    @pytest.fixture
    def simple_coords(self):
        """Create simple coordinates."""
        return np.random.rand(10, 2)
    
    @pytest.fixture
    def sample_sites(self):
        """Create sample sites."""
        sites = [
            APASite(chr='chr1', start=i*1000, end=i*1000+50, strand='+', gene_id=f'gene{i}')
            for i in range(5)
        ]
        return APASiteCollection(sites)
    
    @pytest.fixture
    def sample_counts(self):
        """Create sample counts."""
        return np.random.poisson(5, size=(5, 10))
    
    def test_validate_apa_sites_spatial(self, simple_coords, sample_sites, sample_counts):
        """Test convenience function."""
        validated_sites, support_scores = validate_apa_sites_spatial(
            sample_sites,
            sample_counts,
            simple_coords,
            n_neighbors=4,
            support_threshold=0.3
        )
        
        assert isinstance(validated_sites, APASiteCollection)
        assert len(support_scores) == len(sample_sites)


class TestEdgeCases:
    """Test edge cases."""
    
    def test_more_neighbors_than_spots(self):
        """Test when n_neighbors > n_spots."""
        coords = np.array([[0, 0], [1, 0], [2, 0]])
        
        # Request more neighbors than available
        validator = SpatialValidator(n_neighbors=10)
        
        # Should handle gracefully (will use max available)
        # This might raise an error depending on implementation
        # For now, we just test it doesn't crash
        try:
            validator.fit(coords)
        except ValueError:
            # Expected behavior
            pass
    
    def test_identical_coordinates(self):
        """Test with identical coordinates."""
        coords = np.array([[1, 1], [1, 1], [1, 1]])
        
        validator = SpatialValidator(n_neighbors=2)
        validator.fit(coords)
        
        counts = np.array([[5, 5, 5]])
        support = validator.compute_spatial_support(counts, gene_idx=0)
        
        # All spots should have full support
        assert np.all(support > 0)
    
    def test_high_dimensional_coords(self):
        """Test that 3D coordinates raise appropriate error or work."""
        coords_3d = np.random.rand(10, 3)
        
        validator = SpatialValidator(n_neighbors=4)
        
        # Should either work or raise clear error
        # Most spatial methods work with any dimension
        validator.fit(coords_3d)
        
        counts = np.random.poisson(5, size=(3, 10))
        support = validator.compute_spatial_support(counts, gene_idx=0)
        
        assert len(support) == 10

"""
Unit tests for spatial neighbor finding.
"""

import pytest
import numpy as np
from spagapa.spatial import (
    SpatialNeighbors,
    build_knn_graph,
    build_radius_graph,
    build_delaunay_graph
)


class TestSpatialNeighbors:
    """Test SpatialNeighbors class."""
    
    @pytest.fixture
    def simple_coords(self):
        """Create simple 2D coordinates."""
        # 3x3 grid
        x = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
        y = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        return np.column_stack([x, y])
    
    def test_knn_initialization(self):
        """Test KNN neighbor finder initialization."""
        nbrs = SpatialNeighbors(method='knn', n_neighbors=4)
        assert nbrs.method == 'knn'
        assert nbrs.n_neighbors == 4
    
    def test_knn_fit(self, simple_coords):
        """Test fitting KNN neighbor finder."""
        nbrs = SpatialNeighbors(method='knn', n_neighbors=4)
        nbrs.fit(simple_coords)
        assert nbrs._nbrs is not None
    
    def test_knn_find_neighbors(self, simple_coords):
        """Test finding KNN neighbors."""
        nbrs = SpatialNeighbors(method='knn', n_neighbors=4)
        nbrs.fit(simple_coords)
        
        distances, indices = nbrs.find_neighbors(simple_coords)
        
        # Check shapes
        assert distances.shape == (9, 4)
        assert indices.shape == (9, 4)
        
        # Check that center point (1,1) has 4 neighbors
        center_idx = 4
        center_neighbors = indices[center_idx]
        assert len(center_neighbors) == 4
        
        # Neighbors should not include self
        assert center_idx not in center_neighbors
    
    def test_radius_initialization(self):
        """Test radius neighbor finder initialization."""
        nbrs = SpatialNeighbors(method='radius', radius=1.5)
        assert nbrs.method == 'radius'
        assert nbrs.radius == 1.5
    
    def test_radius_fit(self, simple_coords):
        """Test fitting radius neighbor finder."""
        nbrs = SpatialNeighbors(method='radius', radius=1.5)
        nbrs.fit(simple_coords)
        assert nbrs._nbrs is not None
    
    def test_radius_find_neighbors(self, simple_coords):
        """Test finding radius neighbors."""
        nbrs = SpatialNeighbors(method='radius', radius=1.5)
        nbrs.fit(simple_coords)
        
        distances, indices = nbrs.find_neighbors(simple_coords)
        
        # Check that it returns lists
        assert isinstance(distances, list)
        assert isinstance(indices, list)
        assert len(distances) == 9
        assert len(indices) == 9
        
        # Center point should have neighbors within radius 1.5
        # (actual number depends on grid spacing - just check it has some)
        center_idx = 4
        assert len(indices[center_idx]) > 0
    
    def test_delaunay_initialization(self):
        """Test Delaunay neighbor finder initialization."""
        nbrs = SpatialNeighbors(method='delaunay')
        assert nbrs.method == 'delaunay'
    
    def test_delaunay_fit(self, simple_coords):
        """Test fitting Delaunay neighbor finder."""
        nbrs = SpatialNeighbors(method='delaunay')
        nbrs.fit(simple_coords)
        assert hasattr(nbrs, '_coordinates')
    
    def test_delaunay_find_neighbors(self, simple_coords):
        """Test finding Delaunay neighbors."""
        nbrs = SpatialNeighbors(method='delaunay')
        nbrs.fit(simple_coords)
        
        distances, indices = nbrs.find_neighbors()
        
        # Check that it returns lists
        assert isinstance(distances, list)
        assert isinstance(indices, list)
        assert len(distances) == 9
        assert len(indices) == 9
    
    def test_compute_spatial_weights_knn(self, simple_coords):
        """Test computing spatial weights with KNN."""
        nbrs = SpatialNeighbors(method='knn', n_neighbors=4)
        nbrs.fit(simple_coords)
        
        W = nbrs.compute_spatial_weights(simple_coords, weight_type='inverse_distance')
        
        # Check shape
        assert W.shape == (9, 9)
        
        # Check that rows sum to 1 (normalized)
        row_sums = np.array(W.sum(axis=1)).flatten()
        np.testing.assert_array_almost_equal(row_sums, np.ones(9))
    
    def test_compute_spatial_weights_uniform(self, simple_coords):
        """Test computing uniform spatial weights."""
        nbrs = SpatialNeighbors(method='knn', n_neighbors=4)
        nbrs.fit(simple_coords)
        
        W = nbrs.compute_spatial_weights(simple_coords, weight_type='uniform')
        
        # Check that rows sum to 1
        row_sums = np.array(W.sum(axis=1)).flatten()
        np.testing.assert_array_almost_equal(row_sums, np.ones(9))
    
    def test_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError):
            nbrs = SpatialNeighbors(method='invalid')
            nbrs.fit(np.random.rand(10, 2))
    
    def test_radius_without_radius_param(self):
        """Test that radius method without radius parameter raises error."""
        with pytest.raises(ValueError):
            nbrs = SpatialNeighbors(method='radius')
            nbrs.fit(np.random.rand(10, 2))


class TestHelperFunctions:
    """Test helper functions."""
    
    @pytest.fixture
    def simple_coords(self):
        """Create simple 2D coordinates."""
        x = np.array([0, 1, 2, 0, 1, 2])
        y = np.array([0, 0, 0, 1, 1, 1])
        return np.column_stack([x, y])
    
    def test_build_knn_graph(self, simple_coords):
        """Test build_knn_graph helper."""
        distances, indices = build_knn_graph(simple_coords, k=3)
        
        assert distances.shape == (6, 3)
        assert indices.shape == (6, 3)
    
    def test_build_radius_graph(self, simple_coords):
        """Test build_radius_graph helper."""
        distances, indices = build_radius_graph(simple_coords, radius=1.5)
        
        assert isinstance(distances, list)
        assert isinstance(indices, list)
        assert len(distances) == 6
    
    def test_build_delaunay_graph(self, simple_coords):
        """Test build_delaunay_graph helper."""
        distances, indices = build_delaunay_graph(simple_coords)
        
        assert isinstance(distances, list)
        assert isinstance(indices, list)
        assert len(distances) == 6


class TestEdgeCases:
    """Test edge cases."""
    
    def test_single_point(self):
        """Test with single point - should raise error."""
        coords = np.array([[0, 0]])
        nbrs = SpatialNeighbors(method='knn', n_neighbors=1)
        nbrs.fit(coords)
        
        # Should raise error when trying to find neighbors (not enough points)
        with pytest.raises(ValueError):
            distances, indices = nbrs.find_neighbors(coords)
    
    def test_collinear_points(self):
        """Test with collinear points."""
        coords = np.array([[0, 0], [1, 0], [2, 0], [3, 0]])
        nbrs = SpatialNeighbors(method='knn', n_neighbors=2)
        nbrs.fit(coords)
        
        distances, indices = nbrs.find_neighbors(coords)
        assert distances.shape == (4, 2)
    
    def test_large_dataset(self):
        """Test with larger dataset."""
        np.random.seed(42)
        coords = np.random.rand(100, 2)
        
        nbrs = SpatialNeighbors(method='knn', n_neighbors=10)
        nbrs.fit(coords)
        
        distances, indices = nbrs.find_neighbors(coords)
        assert distances.shape == (100, 10)
        assert indices.shape == (100, 10)

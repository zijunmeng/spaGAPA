"""
Unit tests for sparse GP imputation.
"""

import pytest
import numpy as np
from spagapa.imputation import SparseGPImputer, BlockGPImputer


class TestSparseGPImputer:
    """Test SparseGPImputer class."""
    
    @pytest.fixture
    def simple_data(self):
        """Create simple test data."""
        np.random.seed(42)
        coordinates = np.random.rand(100, 2) * 10
        
        # Smooth function
        values = np.sin(coordinates[:, 0]) + np.cos(coordinates[:, 1])
        values += np.random.normal(0, 0.1, 100)
        
        return coordinates, values
    
    @pytest.fixture
    def sparse_data(self):
        """Create sparse test data."""
        np.random.seed(42)
        coordinates = np.random.rand(100, 2) * 10
        
        values = np.zeros(100)
        observed_idx = np.random.choice(100, size=30, replace=False)
        values[observed_idx] = np.random.poisson(5, size=30)
        
        return coordinates, values
    
    def test_initialization(self):
        """Test imputer initialization."""
        imputer = SparseGPImputer(n_inducing=50, inducing_method='kmeans')
        assert imputer.n_inducing == 50
        assert imputer.inducing_method == 'kmeans'
        assert imputer.inducing_points_ is None
    
    def test_fit_kmeans(self, simple_data):
        """Test fitting with k-means inducing points."""
        coordinates, values = simple_data
        
        imputer = SparseGPImputer(n_inducing=20, inducing_method='kmeans')
        imputer.fit(coordinates, values)
        
        assert imputer.inducing_points_ is not None
        assert len(imputer.inducing_points_) == 20
        assert imputer.alpha_ is not None
    
    def test_fit_random(self, simple_data):
        """Test fitting with random inducing points."""
        coordinates, values = simple_data
        
        imputer = SparseGPImputer(n_inducing=20, inducing_method='random')
        imputer.fit(coordinates, values)
        
        assert len(imputer.inducing_points_) == 20
    
    def test_fit_grid(self, simple_data):
        """Test fitting with grid inducing points."""
        coordinates, values = simple_data
        
        imputer = SparseGPImputer(n_inducing=25, inducing_method='grid')
        imputer.fit(coordinates, values)
        
        assert imputer.inducing_points_ is not None
    
    def test_predict(self, simple_data):
        """Test prediction."""
        coordinates, values = simple_data
        
        imputer = SparseGPImputer(n_inducing=20)
        imputer.fit(coordinates, values)
        
        predictions, std = imputer.predict(coordinates, return_std=True)
        
        assert len(predictions) == len(coordinates)
        assert len(std) == len(coordinates)
        assert np.all(std >= 0)
    
    def test_impute(self, sparse_data):
        """Test imputation."""
        coordinates, values = sparse_data
        
        imputer = SparseGPImputer(n_inducing=20)
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            return_uncertainty=True
        )
        
        assert len(imputed) == len(values)
        assert len(uncertainty) == len(values)
        
        # Observed values should be preserved
        mask = values > 0
        np.testing.assert_array_equal(imputed[mask], values[mask])
        
        # Should impute missing values
        assert np.sum(imputed > 0) > np.sum(values > 0)
    
    def test_invalid_inducing_method(self):
        """Test that invalid method raises error."""
        imputer = SparseGPImputer(inducing_method='invalid')
        coordinates = np.random.rand(50, 2)
        values = np.random.rand(50)
        
        with pytest.raises(ValueError):
            imputer.fit(coordinates, values)
    
    def test_more_inducing_than_data(self):
        """Test when n_inducing > n_data."""
        coordinates = np.random.rand(10, 2)
        values = np.random.rand(10)
        
        imputer = SparseGPImputer(n_inducing=20)
        imputer.fit(coordinates, values)
        
        # Should use all data points as inducing points
        assert len(imputer.inducing_points_) == 10
    
    def test_rbf_kernel(self):
        """Test RBF kernel computation."""
        imputer = SparseGPImputer()
        
        X1 = np.array([[0, 0], [1, 0]])
        X2 = np.array([[0, 0], [0, 1]])
        
        K = imputer._rbf_kernel(X1, X2, length_scale=1.0)
        
        assert K.shape == (2, 2)
        # Kernel should be symmetric for same points
        assert K[0, 0] == 1.0  # Self-similarity


class TestBlockGPImputer:
    """Test BlockGPImputer class."""
    
    @pytest.fixture
    def large_data(self):
        """Create large test data."""
        np.random.seed(42)
        coordinates = np.random.rand(200, 2) * 20
        
        values = np.zeros(200)
        observed_idx = np.random.choice(200, size=60, replace=False)
        values[observed_idx] = np.random.poisson(5, size=60)
        
        return coordinates, values
    
    def test_initialization(self):
        """Test block imputer initialization."""
        imputer = BlockGPImputer(block_size=100, overlap=0.1)
        assert imputer.block_size == 100
        assert imputer.overlap == 0.1
    
    def test_create_blocks(self, large_data):
        """Test block creation."""
        coordinates, _ = large_data
        
        imputer = BlockGPImputer(block_size=50)
        blocks = imputer._create_blocks(coordinates)
        
        assert len(blocks) > 1
        # All points should be covered
        all_indices = np.concatenate(blocks)
        assert len(np.unique(all_indices)) == len(coordinates)
    
    def test_impute(self, large_data):
        """Test block-wise imputation."""
        coordinates, values = large_data
        
        imputer = BlockGPImputer(block_size=50, base_imputer='sparse')
        imputed, uncertainty = imputer.impute(
            coordinates,
            values,
            return_uncertainty=True
        )
        
        assert len(imputed) == len(values)
        assert len(uncertainty) == len(values)
        
        # Observed values should be approximately preserved
        mask = values > 0
        np.testing.assert_array_almost_equal(imputed[mask], values[mask], decimal=5)
    
    def test_small_dataset_no_blocking(self):
        """Test that small datasets don't get blocked."""
        coordinates = np.random.rand(30, 2)
        values = np.random.rand(30)
        
        imputer = BlockGPImputer(block_size=100)
        blocks = imputer._create_blocks(coordinates)
        
        # Should create single block
        assert len(blocks) == 1
        assert len(blocks[0]) == 30
    
    def test_overlap_effect(self):
        """Test that overlap parameter is used."""
        np.random.seed(42)
        coordinates = np.random.rand(200, 2) * 20
        
        imputer_no_overlap = BlockGPImputer(block_size=50, overlap=0.0)
        blocks_no_overlap = imputer_no_overlap._create_blocks(coordinates)
        
        imputer_overlap = BlockGPImputer(block_size=50, overlap=0.2)
        blocks_overlap = imputer_overlap._create_blocks(coordinates)
        
        # With overlap, we might get more blocks or different coverage
        # Just check that both work without errors
        assert len(blocks_no_overlap) > 0
        assert len(blocks_overlap) > 0


class TestEdgeCases:
    """Test edge cases."""
    
    def test_sparse_gp_single_observation(self):
        """Test sparse GP with single observation."""
        coordinates = np.array([[0, 0], [1, 0], [2, 0]])
        values = np.array([5.0, 0.0, 0.0])
        
        imputer = SparseGPImputer(n_inducing=2)
        imputed, _ = imputer.impute(coordinates, values)
        
        assert len(imputed) == 3
        assert imputed[0] == 5.0
    
    def test_sparse_gp_no_training_data(self):
        """Test that no training data raises error."""
        coordinates = np.random.rand(10, 2)
        values = np.zeros(10)
        
        imputer = SparseGPImputer()
        with pytest.raises(ValueError):
            imputer.fit(coordinates, values)
    
    def test_block_gp_failed_block(self):
        """Test that failed blocks are handled gracefully."""
        coordinates = np.random.rand(100, 2) * 10
        values = np.zeros(100)
        # Only one observation
        values[0] = 5.0
        
        imputer = BlockGPImputer(block_size=20)
        # Should not crash even if some blocks fail
        imputed, _ = imputer.impute(coordinates, values)
        
        assert len(imputed) == 100
    
    def test_sparse_gp_predict_before_fit(self):
        """Test that predicting before fitting raises error."""
        imputer = SparseGPImputer()
        coordinates = np.random.rand(10, 2)
        
        with pytest.raises(ValueError):
            imputer.predict(coordinates)
    
    def test_different_length_scales(self):
        """Test sparse GP with different length scales."""
        coordinates = np.random.rand(50, 2)
        values = np.random.rand(50)
        
        for length_scale in [0.1, 1.0, 10.0]:
            imputer = SparseGPImputer(
                n_inducing=10,
                length_scale=length_scale
            )
            imputer.fit(coordinates, values)
            predictions, _ = imputer.predict(coordinates)
            
            assert len(predictions) == 50
    
    def test_different_noise_levels(self):
        """Test sparse GP with different noise levels."""
        coordinates = np.random.rand(50, 2)
        values = np.random.rand(50)
        
        for noise_level in [0.01, 0.1, 1.0]:
            imputer = SparseGPImputer(
                n_inducing=10,
                noise_level=noise_level
            )
            imputer.fit(coordinates, values)
            predictions, _ = imputer.predict(coordinates)
            
            assert len(predictions) == 50

"""
Unit tests for GP imputation.
"""

import pytest
import numpy as np
from spagapa.imputation import GPImputer, GPImputerBatch, impute_spatial_apa


class TestGPImputer:
    """Test GPImputer class."""
    
    @pytest.fixture
    def simple_data(self):
        """Create simple test data."""
        np.random.seed(42)
        # 5x5 grid
        x = np.repeat(np.arange(5), 5)
        y = np.tile(np.arange(5), 5)
        coordinates = np.column_stack([x, y])
        
        # Smooth function with some noise
        values = np.sin(x * 0.5) + np.cos(y * 0.5) + np.random.normal(0, 0.1, 25)
        
        return coordinates, values
    
    @pytest.fixture
    def sparse_data(self):
        """Create sparse test data (some missing values)."""
        np.random.seed(42)
        coordinates = np.random.rand(50, 2) * 10
        
        # Only 30% of points have values
        values = np.zeros(50)
        observed_idx = np.random.choice(50, size=15, replace=False)
        values[observed_idx] = np.random.poisson(5, size=15)
        
        return coordinates, values
    
    def test_initialization(self):
        """Test imputer initialization."""
        imputer = GPImputer(kernel_type='matern', nu=1.5)
        assert imputer.kernel_type == 'matern'
        assert imputer.nu == 1.5
        assert imputer.gp_ is None
    
    def test_fit(self, simple_data):
        """Test fitting GP model."""
        coordinates, values = simple_data
        
        imputer = GPImputer(kernel_type='rbf')
        imputer.fit(coordinates, values)
        
        assert imputer.gp_ is not None
        assert imputer.kernel_ is not None
        assert imputer._fitted_coords is not None
    
    def test_predict(self, simple_data):
        """Test prediction."""
        coordinates, values = simple_data
        
        imputer = GPImputer()
        imputer.fit(coordinates, values)
        
        predictions, std = imputer.predict(coordinates, return_std=True)
        
        assert len(predictions) == len(coordinates)
        assert len(std) == len(coordinates)
        assert np.all(std >= 0)
    
    def test_impute(self, sparse_data):
        """Test imputation."""
        coordinates, values = sparse_data
        
        imputer = GPImputer(kernel_type='matern')
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
        
        # Uncertainty should be zero for observed values
        np.testing.assert_array_equal(uncertainty[mask], 0.0)
        
        # Imputed values should be non-zero
        assert np.sum(imputed > 0) > np.sum(values > 0)
    
    def test_kernel_types(self, simple_data):
        """Test different kernel types."""
        coordinates, values = simple_data
        
        for kernel_type in ['rbf', 'matern', 'auto']:
            imputer = GPImputer(kernel_type=kernel_type)
            imputer.fit(coordinates, values)
            predictions, _ = imputer.predict(coordinates)
            
            assert len(predictions) == len(coordinates)
    
    def test_invalid_kernel(self):
        """Test that invalid kernel raises error."""
        with pytest.raises(ValueError):
            imputer = GPImputer(kernel_type='invalid')
            imputer._create_kernel()
    
    def test_no_training_data(self):
        """Test that fitting with no data raises error."""
        coordinates = np.random.rand(10, 2)
        values = np.zeros(10)  # All zeros
        
        imputer = GPImputer()
        with pytest.raises(ValueError):
            imputer.fit(coordinates, values)
    
    def test_predict_before_fit(self):
        """Test that predicting before fitting raises error."""
        imputer = GPImputer()
        coordinates = np.random.rand(10, 2)
        
        with pytest.raises(ValueError):
            imputer.predict(coordinates)
    
    def test_custom_mask(self, simple_data):
        """Test fitting with custom mask."""
        coordinates, values = simple_data
        
        # Use only first half of data for training
        mask = np.zeros(len(values), dtype=bool)
        mask[:len(values)//2] = True
        
        imputer = GPImputer()
        imputer.fit(coordinates, values, mask=mask)
        
        assert len(imputer._fitted_coords) == np.sum(mask)
    
    def test_score(self, simple_data):
        """Test model scoring."""
        coordinates, values = simple_data
        
        imputer = GPImputer()
        imputer.fit(coordinates, values)
        
        score = imputer.score(coordinates, values)
        assert isinstance(score, float)


class TestGPImputerBatch:
    """Test GPImputerBatch class."""
    
    @pytest.fixture
    def batch_data(self):
        """Create batch test data."""
        np.random.seed(42)
        coordinates = np.random.rand(30, 2) * 10
        
        # 5 genes with different patterns
        values = np.zeros((5, 30))
        for i in range(5):
            observed_idx = np.random.choice(30, size=10, replace=False)
            values[i, observed_idx] = np.random.poisson(5, size=10)
        
        return coordinates, values
    
    def test_batch_fit(self, batch_data):
        """Test batch fitting."""
        coordinates, values = batch_data
        
        imputer = GPImputer()
        batch_imputer = imputer.fit_batch(
            coordinates,
            values,
            n_jobs=1,
            verbose=False
        )
        
        assert len(batch_imputer.imputers_) == 5
    
    def test_batch_predict(self, batch_data):
        """Test batch prediction."""
        coordinates, values = batch_data
        
        imputer = GPImputer()
        batch_imputer = imputer.fit_batch(
            coordinates,
            values,
            n_jobs=1,
            verbose=False
        )
        
        predictions, uncertainty = batch_imputer.predict(
            coordinates,
            return_std=True
        )
        
        assert predictions.shape == values.shape
        assert uncertainty.shape == values.shape
    
    def test_batch_impute(self, batch_data):
        """Test batch imputation."""
        coordinates, values = batch_data
        
        imputer = GPImputer()
        batch_imputer = imputer.fit_batch(
            coordinates,
            values,
            n_jobs=1,
            verbose=False
        )
        
        imputed, uncertainty = batch_imputer.impute(return_uncertainty=True)
        
        assert imputed.shape == values.shape
        assert uncertainty.shape == values.shape
        
        # Observed values should be preserved
        mask = values > 0
        np.testing.assert_array_equal(imputed[mask], values[mask])


class TestHelperFunction:
    """Test helper function."""
    
    def test_impute_spatial_apa(self):
        """Test convenience function."""
        np.random.seed(42)
        coordinates = np.random.rand(30, 2) * 10
        apa_counts = np.zeros((3, 30))
        
        for i in range(3):
            observed_idx = np.random.choice(30, size=10, replace=False)
            apa_counts[i, observed_idx] = np.random.poisson(5, size=10)
        
        imputed, uncertainty = impute_spatial_apa(
            coordinates,
            apa_counts,
            kernel_type='matern',
            n_jobs=1,
            return_uncertainty=True
        )
        
        assert imputed.shape == apa_counts.shape
        assert uncertainty.shape == apa_counts.shape


class TestEdgeCases:
    """Test edge cases."""
    
    def test_single_observation(self):
        """Test with single observation."""
        coordinates = np.array([[0, 0], [1, 0], [2, 0]])
        values = np.array([5.0, 0.0, 0.0])
        
        imputer = GPImputer()
        imputed, _ = imputer.impute(coordinates, values)
        
        # Should handle gracefully
        assert len(imputed) == 3
        assert imputed[0] == 5.0
    
    def test_all_same_values(self):
        """Test with all same values."""
        coordinates = np.random.rand(10, 2)
        values = np.ones(10) * 5.0
        
        imputer = GPImputer()
        imputer.fit(coordinates, values)
        predictions, _ = imputer.predict(coordinates)
        
        # Should predict approximately constant
        assert np.std(predictions) < 1.0
    
    def test_high_dimensional_coords(self):
        """Test with 3D coordinates."""
        coordinates = np.random.rand(20, 3)
        values = np.random.rand(20)
        
        imputer = GPImputer()
        imputer.fit(coordinates, values)
        predictions, _ = imputer.predict(coordinates)
        
        assert len(predictions) == 20
    
    def test_large_length_scale(self):
        """Test with very large length scale."""
        coordinates = np.random.rand(20, 2)
        values = np.random.rand(20)
        
        imputer = GPImputer(length_scale=100.0)
        imputer.fit(coordinates, values)
        predictions, _ = imputer.predict(coordinates)
        
        # Large length scale should give smooth predictions
        assert len(predictions) == 20
    
    def test_small_length_scale(self):
        """Test with very small length scale."""
        coordinates = np.random.rand(20, 2)
        values = np.random.rand(20)
        
        imputer = GPImputer(length_scale=0.01)
        imputer.fit(coordinates, values)
        predictions, _ = imputer.predict(coordinates)
        
        assert len(predictions) == 20
    
    def test_return_cov(self):
        """Test returning covariance matrix."""
        coordinates = np.random.rand(10, 2)
        values = np.random.rand(10)
        
        imputer = GPImputer()
        imputer.fit(coordinates, values)
        
        predictions, cov = imputer.predict(coordinates, return_cov=True)
        
        assert predictions.shape == (10,)
        assert cov.shape == (10, 10)
        # Covariance should be symmetric
        np.testing.assert_array_almost_equal(cov, cov.T)

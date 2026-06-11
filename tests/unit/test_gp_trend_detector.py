"""
Tests for GP-based spatial trend detection
"""

import numpy as np
import pytest
from spagapa.analysis.gp_trend_detector import (
    GPTrendDetector,
    detect_svapa_genes_gp
)


@pytest.fixture
def spatial_data():
    """Create synthetic spatial data with known patterns."""
    np.random.seed(42)
    
    # Create grid coordinates
    x = np.linspace(0, 10, 20)
    y = np.linspace(0, 10, 20)
    xx, yy = np.meshgrid(x, y)
    coordinates = np.column_stack([xx.ravel(), yy.ravel()])
    
    # Gene 1: Strong spatial gradient
    values_spatial = xx.ravel() + 0.1 * np.random.randn(len(xx.ravel()))
    
    # Gene 2: Random noise (no spatial structure)
    values_random = np.random.randn(len(xx.ravel()))
    
    # Gene 3: Spatial hotspot
    dist_from_center = np.sqrt((xx.ravel() - 5)**2 + (yy.ravel() - 5)**2)
    values_hotspot = np.exp(-dist_from_center / 2) + 0.1 * np.random.randn(len(xx.ravel()))
    
    # Uncertainty (higher for some spots)
    uncertainty = 0.1 + 0.2 * np.random.rand(len(xx.ravel()))
    
    return {
        'coordinates': coordinates,
        'values_spatial': values_spatial,
        'values_random': values_random,
        'values_hotspot': values_hotspot,
        'uncertainty': uncertainty
    }


class TestGPTrendDetector:
    """Test GPTrendDetector class."""
    
    def test_initialization(self):
        """Test detector initialization."""
        detector = GPTrendDetector(kernel_type='matern')
        assert detector.kernel_type == 'matern'
        assert detector.alpha == 1e-10
        assert detector.n_restarts == 5
    
    def test_create_null_kernel(self):
        """Test null kernel creation."""
        detector = GPTrendDetector()
        kernel = detector._create_null_kernel()
        assert kernel is not None
    
    def test_create_spatial_kernel(self):
        """Test spatial kernel creation."""
        detector = GPTrendDetector(kernel_type='rbf')
        kernel = detector._create_spatial_kernel()
        assert kernel is not None
        
        detector = GPTrendDetector(kernel_type='matern')
        kernel = detector._create_spatial_kernel(length_scale=1.0)
        assert kernel is not None
    
    def test_likelihood_ratio_test_spatial(self, spatial_data):
        """Test LR test on data with spatial structure."""
        detector = GPTrendDetector(kernel_type='matern')
        
        result = detector.likelihood_ratio_test(
            spatial_data['values_spatial'],
            spatial_data['coordinates']
        )
        
        assert 'log_lik_null' in result
        assert 'log_lik_alt' in result
        assert 'lr_statistic' in result
        assert 'p_value' in result
        assert 'df' in result
        
        # Should detect spatial structure
        assert result['lr_statistic'] > 0
        assert result['p_value'] < 0.05
    
    def test_likelihood_ratio_test_random(self, spatial_data):
        """Test LR test on random data."""
        detector = GPTrendDetector(kernel_type='matern')
        
        result = detector.likelihood_ratio_test(
            spatial_data['values_random'],
            spatial_data['coordinates']
        )
        
        # Should not detect spatial structure
        assert result['p_value'] > 0.05 or result['lr_statistic'] < 10
    
    def test_likelihood_ratio_test_with_uncertainty(self, spatial_data):
        """Test LR test with uncertainty weighting."""
        detector = GPTrendDetector(kernel_type='matern')
        
        result = detector.likelihood_ratio_test(
            spatial_data['values_spatial'],
            spatial_data['coordinates'],
            uncertainty=spatial_data['uncertainty']
        )
        
        assert result['p_value'] < 0.1  # Should still detect
    
    def test_likelihood_ratio_test_small_sample(self):
        """Test LR test with small sample."""
        detector = GPTrendDetector()
        
        values = np.array([1, 2, 3])
        coords = np.array([[0, 0], [1, 1], [2, 2]])
        
        result = detector.likelihood_ratio_test(values, coords)
        
        assert np.isnan(result['lr_statistic'])
        assert result['p_value'] == 1.0
    
    def test_uncertainty_weighted_morans_i(self, spatial_data):
        """Test uncertainty-weighted Moran's I."""
        detector = GPTrendDetector()
        
        result = detector.uncertainty_weighted_morans_i(
            spatial_data['values_spatial'],
            spatial_data['coordinates'],
            spatial_data['uncertainty'],
            k=6
        )
        
        assert 'morans_i' in result
        assert 'expected_i' in result
        assert 'variance_i' in result
        assert 'z_score' in result
        assert 'p_value' in result
        
        # Should detect positive spatial autocorrelation
        assert result['morans_i'] > result['expected_i']
    
    def test_uncertainty_weighted_morans_i_random(self, spatial_data):
        """Test weighted Moran's I on random data."""
        detector = GPTrendDetector()
        
        result = detector.uncertainty_weighted_morans_i(
            spatial_data['values_random'],
            spatial_data['coordinates'],
            spatial_data['uncertainty'],
            k=6
        )
        
        # Should not detect strong autocorrelation
        assert abs(result['z_score']) < 3
    
    def test_spatial_variance_decomposition(self, spatial_data):
        """Test spatial variance decomposition."""
        detector = GPTrendDetector()
        
        result = detector.spatial_variance_decomposition(
            spatial_data['values_spatial'],
            spatial_data['coordinates']
        )
        
        assert 'total_variance' in result
        assert 'spatial_variance' in result
        assert 'random_variance' in result
        assert 'spatial_fraction' in result
        assert 'r_squared' in result
        
        # Should have high spatial fraction
        assert result['spatial_fraction'] > 0.5
        assert result['r_squared'] > 0.5
    
    def test_spatial_variance_decomposition_random(self, spatial_data):
        """Test variance decomposition on random data."""
        detector = GPTrendDetector()
        
        result = detector.spatial_variance_decomposition(
            spatial_data['values_random'],
            spatial_data['coordinates']
        )
        
        # Should have low spatial fraction
        assert result['spatial_fraction'] < 0.5
    
    def test_detect_gp_trends_likelihood_ratio(self, spatial_data):
        """Test SVAPA detection with likelihood ratio method."""
        detector = GPTrendDetector(kernel_type='matern')
        
        # Stack genes
        apa_values = np.vstack([
            spatial_data['values_spatial'],
            spatial_data['values_random'],
            spatial_data['values_hotspot']
        ])
        
        uncertainty = np.vstack([
            spatial_data['uncertainty'],
            spatial_data['uncertainty'],
            spatial_data['uncertainty']
        ])
        
        gene_names = ['Gene1_spatial', 'Gene2_random', 'Gene3_hotspot']
        
        results = detector.detect_gp_trends(
            apa_values,
            spatial_data['coordinates'],
            uncertainty,
            gene_names,
            method='likelihood_ratio',
            fdr_threshold=0.1
        )
        
        assert len(results) == 3
        assert 'gene' in results.columns
        assert 'lr_p_value' in results.columns
        assert 'q_value' in results.columns
        assert 'significant' in results.columns
        
        # Gene1 and Gene3 should be significant
        sig_genes = results[results['significant']]['gene'].tolist()
        assert 'Gene1_spatial' in sig_genes or 'Gene3_hotspot' in sig_genes
    
    def test_detect_gp_trends_weighted_morans(self, spatial_data):
        """Test SVAPA detection with weighted Moran's I."""
        detector = GPTrendDetector()
        
        apa_values = np.vstack([
            spatial_data['values_spatial'],
            spatial_data['values_random']
        ])
        
        uncertainty = np.vstack([
            spatial_data['uncertainty'],
            spatial_data['uncertainty']
        ])
        
        results = detector.detect_gp_trends(
            apa_values,
            spatial_data['coordinates'],
            uncertainty,
            gene_names=['Gene1', 'Gene2'],
            method='weighted_morans',
            fdr_threshold=0.1
        )
        
        assert 'morans_i' in results.columns
        assert 'morans_p_value' in results.columns
    
    def test_detect_gp_trends_all_methods(self, spatial_data):
        """Test SVAPA detection with all methods."""
        detector = GPTrendDetector()
        
        apa_values = spatial_data['values_spatial'].reshape(1, -1)
        uncertainty = spatial_data['uncertainty'].reshape(1, -1)
        
        results = detector.detect_gp_trends(
            apa_values,
            spatial_data['coordinates'],
            uncertainty,
            gene_names=['Gene1'],
            method='all',
            fdr_threshold=0.1
        )
        
        assert 'lr_p_value' in results.columns
        assert 'morans_p_value' in results.columns
        assert 'spatial_fraction' in results.columns
        assert 'r_squared' in results.columns
    
    def test_detect_gp_trends_no_uncertainty(self, spatial_data):
        """Test SVAPA detection without uncertainty."""
        detector = GPTrendDetector()
        
        apa_values = spatial_data['values_spatial'].reshape(1, -1)
        
        results = detector.detect_gp_trends(
            apa_values,
            spatial_data['coordinates'],
            uncertainty=None,
            gene_names=['Gene1'],
            method='likelihood_ratio'
        )
        
        assert len(results) == 1
        assert not np.isnan(results['lr_p_value'].iloc[0])


class TestConvenienceFunction:
    """Test convenience function."""
    
    def test_detect_svapa_genes_gp(self, spatial_data):
        """Test convenience function."""
        apa_values = np.vstack([
            spatial_data['values_spatial'],
            spatial_data['values_random']
        ])
        
        uncertainty = np.vstack([
            spatial_data['uncertainty'],
            spatial_data['uncertainty']
        ])
        
        results = detect_svapa_genes_gp(
            apa_values,
            spatial_data['coordinates'],
            uncertainty,
            gene_names=['Gene1', 'Gene2'],
            method='likelihood_ratio',
            kernel_type='matern',
            fdr_threshold=0.1
        )
        
        assert len(results) == 2
        assert 'significant' in results.columns
    
    def test_detect_svapa_genes_gp_minimal(self, spatial_data):
        """Test convenience function with minimal arguments."""
        apa_values = spatial_data['values_spatial'].reshape(1, -1)
        
        results = detect_svapa_genes_gp(
            apa_values,
            spatial_data['coordinates']
        )
        
        assert len(results) == 1
        assert 'gene' in results.columns


class TestEdgeCases:
    """Test edge cases."""
    
    def test_nan_values(self, spatial_data):
        """Test handling of NaN values."""
        detector = GPTrendDetector()
        
        values = spatial_data['values_spatial'].copy()
        values[:10] = np.nan
        
        result = detector.likelihood_ratio_test(
            values,
            spatial_data['coordinates']
        )
        
        assert not np.isnan(result['p_value'])
    
    def test_constant_values(self):
        """Test with constant values."""
        detector = GPTrendDetector()
        
        coords = np.random.rand(50, 2)
        values = np.ones(50)
        
        result = detector.spatial_variance_decomposition(values, coords)
        
        # Should handle gracefully
        assert result['total_variance'] < 1e-10
    
    def test_single_gene(self, spatial_data):
        """Test with single gene."""
        detector = GPTrendDetector()
        
        apa_values = spatial_data['values_spatial'].reshape(1, -1)
        
        results = detector.detect_gp_trends(
            apa_values,
            spatial_data['coordinates'],
            method='likelihood_ratio'
        )
        
        assert len(results) == 1

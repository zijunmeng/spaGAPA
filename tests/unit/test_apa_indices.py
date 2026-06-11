"""
Unit tests for APA index calculations.
"""

import pytest
import numpy as np
from spagapa.quantification.apa_indices import (
    calculate_rud,
    calculate_pdui,
    calculate_wul,
    calculate_pai,
    normalize_apa_index,
    APAIndexCalculator
)


class TestCalculateRUD:
    """Tests for RUD calculation."""
    
    def test_basic_rud(self):
        """Test basic RUD calculation."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        rud = calculate_rud(proximal, distal, pseudocount=0)
        
        # RUD = distal / (proximal + distal)
        expected = np.array([30/40, 10/30, 15/20])
        np.testing.assert_array_almost_equal(rud, expected)
    
    def test_rud_with_pseudocount(self):
        """Test RUD with pseudocount."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        rud = calculate_rud(proximal, distal, pseudocount=1.0)
        
        # RUD = distal / (proximal + distal + 1)
        expected = np.array([30/41, 10/31, 15/21])
        np.testing.assert_array_almost_equal(rud, expected)
    
    def test_rud_zeros(self):
        """Test RUD with zero counts."""
        proximal = np.array([0, 10, 0])
        distal = np.array([0, 20, 10])
        
        rud = calculate_rud(proximal, distal, pseudocount=1.0)
        
        # Should handle zeros gracefully
        assert np.all(rud >= 0) and np.all(rud <= 1)
    
    def test_rud_2d(self):
        """Test RUD with 2D arrays."""
        proximal = np.array([[10, 20], [5, 15]])
        distal = np.array([[30, 10], [15, 5]])
        
        rud = calculate_rud(proximal, distal, pseudocount=0)
        
        assert rud.shape == (2, 2)
        assert np.all(rud >= 0) and np.all(rud <= 1)
    
    def test_rud_shape_mismatch(self):
        """Test RUD with mismatched shapes."""
        proximal = np.array([10, 20])
        distal = np.array([30, 10, 15])
        
        with pytest.raises(ValueError):
            calculate_rud(proximal, distal)


class TestCalculatePDUI:
    """Tests for PDUI calculation."""
    
    def test_basic_pdui(self):
        """Test basic PDUI calculation."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        pdui = calculate_pdui(proximal, distal, pseudocount=0)
        
        # PDUI = 100 * distal / (proximal + distal)
        expected = 100 * np.array([30/40, 10/30, 15/20])
        np.testing.assert_array_almost_equal(pdui, expected)
    
    def test_pdui_with_long_form(self):
        """Test PDUI with long-form counts."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        long_form = np.array([5, 5, 5])
        
        pdui = calculate_pdui(proximal, distal, long_form, pseudocount=0)
        
        # PDUI = 100 * (distal + long_form) / (proximal + distal + long_form)
        expected = 100 * np.array([35/45, 15/35, 20/25])
        np.testing.assert_array_almost_equal(pdui, expected)
    
    def test_pdui_range(self):
        """Test PDUI is in [0, 100] range."""
        proximal = np.array([10, 0, 100])
        distal = np.array([0, 50, 10])
        
        pdui = calculate_pdui(proximal, distal, pseudocount=1.0)
        
        assert np.all(pdui >= 0) and np.all(pdui <= 100)


class TestCalculateWUL:
    """Tests for WUL calculation."""
    
    def test_basic_wul(self):
        """Test basic WUL calculation."""
        site_counts = {
            'proximal': np.array([10, 20, 5]),
            'distal': np.array([30, 10, 15])
        }
        site_positions = {'proximal': 1000, 'distal': 2000}
        
        wul = calculate_wul(site_counts, site_positions, normalize=False, pseudocount=0)
        
        # WUL = (10*1000 + 30*2000) / (10+30) for first spot
        expected = np.array([
            (10*1000 + 30*2000) / 40,
            (20*1000 + 10*2000) / 30,
            (5*1000 + 15*2000) / 20
        ])
        np.testing.assert_array_almost_equal(wul, expected)
    
    def test_wul_normalized(self):
        """Test normalized WUL."""
        site_counts = {
            'proximal': np.array([10, 20, 5]),
            'distal': np.array([30, 10, 15])
        }
        site_positions = {'proximal': 1000, 'distal': 2000}
        
        wul = calculate_wul(site_counts, site_positions, normalize=True, pseudocount=0)
        
        # Normalized WUL should be in [0, 1]
        assert np.all(wul >= 0) and np.all(wul <= 1)
    
    def test_wul_multiple_sites(self):
        """Test WUL with multiple sites."""
        site_counts = {
            'site1': np.array([10, 20]),
            'site2': np.array([20, 10]),
            'site3': np.array([30, 5])
        }
        site_positions = {'site1': 1000, 'site2': 1500, 'site3': 2000}
        
        wul = calculate_wul(site_counts, site_positions, normalize=True)
        
        assert wul.shape == (2,)
        assert np.all(wul >= 0) and np.all(wul <= 1)
    
    def test_wul_empty_input(self):
        """Test WUL with empty input."""
        with pytest.raises(ValueError):
            calculate_wul({}, {})
    
    def test_wul_mismatched_keys(self):
        """Test WUL with mismatched keys."""
        site_counts = {'site1': np.array([10, 20])}
        site_positions = {'site2': 1000}
        
        with pytest.raises(ValueError):
            calculate_wul(site_counts, site_positions)


class TestCalculatePAI:
    """Tests for PAI calculation."""
    
    def test_pai_ratio_method(self):
        """Test PAI with ratio method."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        pai = calculate_pai(proximal, distal, method='ratio', pseudocount=1.0)
        
        # PAI = log2((distal + 1) / (proximal + 1))
        expected = np.log2((distal + 1) / (proximal + 1))
        np.testing.assert_array_almost_equal(pai, expected)
    
    def test_pai_difference_method(self):
        """Test PAI with difference method."""
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        pai = calculate_pai(proximal, distal, method='difference', pseudocount=0)
        
        # PAI = (distal - proximal) / (distal + proximal)
        expected = (distal - proximal) / (distal + proximal)
        np.testing.assert_array_almost_equal(pai, expected)
    
    def test_pai_invalid_method(self):
        """Test PAI with invalid method."""
        proximal = np.array([10, 20])
        distal = np.array([30, 10])
        
        with pytest.raises(ValueError):
            calculate_pai(proximal, distal, method='invalid')


class TestNormalizeAPAIndex:
    """Tests for APA index normalization."""
    
    def test_zscore_normalization(self):
        """Test z-score normalization."""
        values = np.array([10, 20, 30, 40, 50])
        
        normalized = normalize_apa_index(values, method='zscore')
        
        # Should have mean ≈ 0 and std ≈ 1
        assert np.abs(np.mean(normalized)) < 1e-10
        assert np.abs(np.std(normalized) - 1.0) < 1e-10
    
    def test_minmax_normalization(self):
        """Test min-max normalization."""
        values = np.array([10, 20, 30, 40, 50])
        
        normalized = normalize_apa_index(values, method='minmax')
        
        # Should be in [0, 1] range
        assert np.abs(np.min(normalized) - 0.0) < 1e-10
        assert np.abs(np.max(normalized) - 1.0) < 1e-10
    
    def test_quantile_normalization(self):
        """Test quantile normalization."""
        values = np.array([10, 20, 30, 40, 50])
        
        normalized = normalize_apa_index(values, method='quantile')
        
        # Should be in [0, 1] range
        assert np.all(normalized >= 0) and np.all(normalized <= 1)
    
    def test_normalization_2d(self):
        """Test normalization with 2D array."""
        values = np.array([[10, 20, 30], [40, 50, 60]])
        
        # Normalize along axis 1
        normalized = normalize_apa_index(values, method='zscore', axis=1)
        
        assert normalized.shape == values.shape
    
    def test_invalid_method(self):
        """Test with invalid normalization method."""
        values = np.array([10, 20, 30])
        
        with pytest.raises(ValueError):
            normalize_apa_index(values, method='invalid')


class TestAPAIndexCalculator:
    """Tests for APAIndexCalculator class."""
    
    def test_calculator_initialization(self):
        """Test calculator initialization."""
        calc = APAIndexCalculator(pseudocount=1.0, normalize=False)
        
        assert calc.pseudocount == 1.0
        assert calc.normalize is False
    
    def test_calculate_all_basic(self):
        """Test calculating all indices."""
        calc = APAIndexCalculator(pseudocount=0, normalize=False)
        
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        indices = calc.calculate_all(proximal, distal)
        
        assert 'RUD' in indices
        assert 'PDUI' in indices
        assert 'PAI' in indices
        assert indices['RUD'].shape == (3,)
    
    def test_calculate_all_with_wul(self):
        """Test calculating all indices including WUL."""
        calc = APAIndexCalculator(pseudocount=0, normalize=False)
        
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        site_positions = {'proximal': 1000, 'distal': 2000}
        
        indices = calc.calculate_all(proximal, distal, site_positions=site_positions)
        
        assert 'WUL' in indices
        assert indices['WUL'].shape == (3,)
    
    def test_calculate_all_with_normalization(self):
        """Test calculating indices with normalization."""
        calc = APAIndexCalculator(pseudocount=0, normalize=True, normalization_method='zscore')
        
        proximal = np.array([10, 20, 5, 15, 25])
        distal = np.array([30, 10, 15, 20, 5])
        
        indices = calc.calculate_all(proximal, distal)
        
        # Check that values are normalized (mean ≈ 0)
        assert np.abs(np.mean(indices['RUD'])) < 0.5
    
    def test_to_dataframe_1d(self):
        """Test converting 1D indices to DataFrame."""
        calc = APAIndexCalculator()
        
        proximal = np.array([10, 20, 5])
        distal = np.array([30, 10, 15])
        
        indices = calc.calculate_all(proximal, distal)
        df = calc.to_dataframe(indices, spot_names=['spot1', 'spot2', 'spot3'])
        
        assert len(df) == 3
        assert 'RUD' in df.columns
        assert 'PDUI' in df.columns
    
    def test_to_dataframe_2d(self):
        """Test converting 2D indices to DataFrame."""
        calc = APAIndexCalculator()
        
        proximal = np.array([[10, 20], [5, 15]])
        distal = np.array([[30, 10], [15, 5]])
        
        indices = calc.calculate_all(proximal, distal)
        df = calc.to_dataframe(
            indices,
            gene_names=['gene1', 'gene2'],
            spot_names=['spot1', 'spot2']
        )
        
        assert 'index_type' in df.columns


class TestEdgeCases:
    """Tests for edge cases."""
    
    def test_all_zeros(self):
        """Test with all zero counts."""
        proximal = np.array([0, 0, 0])
        distal = np.array([0, 0, 0])
        
        rud = calculate_rud(proximal, distal, pseudocount=1.0)
        
        # Should not raise error
        assert rud.shape == (3,)
    
    def test_single_value(self):
        """Test with single value."""
        proximal = np.array([10])
        distal = np.array([30])
        
        rud = calculate_rud(proximal, distal)
        
        assert rud.shape == (1,)
        assert 0 <= rud[0] <= 1
    
    def test_large_values(self):
        """Test with large count values."""
        proximal = np.array([1e6, 2e6])
        distal = np.array([3e6, 1e6])
        
        rud = calculate_rud(proximal, distal)
        
        assert np.all(np.isfinite(rud))
        assert np.all(rud >= 0) and np.all(rud <= 1)

"""
Unit tests for QC metrics.
"""

import pytest
import numpy as np
import pandas as pd
from spagapa.quantification.qc_metrics import (
    calculate_rmse,
    calculate_mae,
    calculate_pearson,
    calculate_spearman,
    calculate_r2,
    cross_validate_imputation,
    evaluate_imputation_quality,
    QCReportGenerator,
    compare_imputation_methods
)


class TestBasicMetrics:
    """Tests for basic quality metrics."""
    
    def test_rmse_perfect(self):
        """Test RMSE with perfect predictions."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1, 2, 3, 4, 5])
        
        rmse = calculate_rmse(y_true, y_pred)
        
        assert rmse == 0.0
    
    def test_rmse_basic(self):
        """Test RMSE calculation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
        
        rmse = calculate_rmse(y_true, y_pred)
        
        # Manual calculation
        mse = np.mean((y_true - y_pred) ** 2)
        expected = np.sqrt(mse)
        
        assert np.abs(rmse - expected) < 1e-10
    
    def test_rmse_with_mask(self):
        """Test RMSE with mask."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
        mask = np.array([True, True, False, True, False])
        
        rmse = calculate_rmse(y_true, y_pred, mask)
        
        # Should only use masked values
        expected = np.sqrt(np.mean((y_true[mask] - y_pred[mask]) ** 2))
        assert np.abs(rmse - expected) < 1e-10
    
    def test_mae_perfect(self):
        """Test MAE with perfect predictions."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1, 2, 3, 4, 5])
        
        mae = calculate_mae(y_true, y_pred)
        
        assert mae == 0.0
    
    def test_mae_basic(self):
        """Test MAE calculation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1.5, 2.5, 2.5, 4.5, 4.5])
        
        mae = calculate_mae(y_true, y_pred)
        
        expected = np.mean(np.abs(y_true - y_pred))
        assert np.abs(mae - expected) < 1e-10
    
    def test_pearson_perfect(self):
        """Test Pearson correlation with perfect correlation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([2, 4, 6, 8, 10])  # Perfect linear relationship
        
        corr, pval = calculate_pearson(y_true, y_pred)
        
        assert np.abs(corr - 1.0) < 1e-10
        assert pval < 0.05
    
    def test_pearson_no_correlation(self):
        """Test Pearson with no correlation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([5, 4, 3, 2, 1])  # Negative correlation
        
        corr, pval = calculate_pearson(y_true, y_pred)
        
        assert corr < 0  # Negative correlation
    
    def test_pearson_insufficient_data(self):
        """Test Pearson with insufficient data."""
        y_true = np.array([1])
        y_pred = np.array([1])
        
        corr, pval = calculate_pearson(y_true, y_pred)
        
        assert np.isnan(corr)
        assert np.isnan(pval)
    
    def test_spearman_basic(self):
        """Test Spearman correlation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1.1, 2.2, 2.9, 4.1, 5.2])
        
        corr, pval = calculate_spearman(y_true, y_pred)
        
        assert corr > 0.9  # Should be highly correlated
        assert pval < 0.05
    
    def test_r2_perfect(self):
        """Test R² with perfect predictions."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1, 2, 3, 4, 5])
        
        r2 = calculate_r2(y_true, y_pred)
        
        assert np.abs(r2 - 1.0) < 1e-10
    
    def test_r2_basic(self):
        """Test R² calculation."""
        y_true = np.array([1, 2, 3, 4, 5])
        y_pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])
        
        r2 = calculate_r2(y_true, y_pred)
        
        # R² should be between 0 and 1 for reasonable predictions
        assert 0 <= r2 <= 1
    
    def test_r2_constant_predictions(self):
        """Test R² with constant true values."""
        y_true = np.array([5, 5, 5, 5, 5])
        y_pred = np.array([4, 5, 6, 5, 5])
        
        r2 = calculate_r2(y_true, y_pred)
        
        # Should return NaN when variance is zero
        assert np.isnan(r2)


class TestCrossValidation:
    """Tests for cross-validation."""
    
    def test_cross_validate_basic(self):
        """Test basic cross-validation."""
        # Create mock imputer
        class MockImputer:
            def fit(self, coords, values, mask):
                self.mean = np.mean(values[mask])
            
            def predict(self, coords):
                return np.full(len(coords), self.mean), None
        
        imputer = MockImputer()
        coordinates = np.random.rand(20, 2)
        values = np.random.rand(20) * 10
        
        metrics = cross_validate_imputation(
            imputer, coordinates, values, n_folds=3
        )
        
        assert 'rmse' in metrics
        assert 'mae' in metrics
        assert 'pearson' in metrics
        assert 'r2' in metrics
    
    def test_cross_validate_insufficient_data(self):
        """Test cross-validation with insufficient data."""
        class MockImputer:
            def fit(self, coords, values, mask):
                pass
            def predict(self, coords):
                return np.zeros(len(coords)), None
        
        imputer = MockImputer()
        coordinates = np.random.rand(3, 2)
        values = np.array([1, 2, 3])
        
        # Should adjust n_folds automatically
        metrics = cross_validate_imputation(
            imputer, coordinates, values, n_folds=5
        )
        
        assert 'rmse' in metrics


class TestEvaluateImputationQuality:
    """Tests for imputation quality evaluation."""
    
    def test_evaluate_basic(self):
        """Test basic quality evaluation."""
        true_values = np.array([1, 2, 0, 4, 0, 6])
        imputed_values = np.array([1, 2, 3, 4, 5, 6])
        mask = true_values == 0  # Imputed positions
        
        metrics = evaluate_imputation_quality(
            true_values, imputed_values, mask=mask
        )
        
        assert 'rmse' in metrics
        assert 'mae' in metrics
        assert 'pearson_r' in metrics
        assert 'r2' in metrics
        assert 'n_imputed' in metrics
        assert metrics['n_imputed'] == 2
    
    def test_evaluate_with_uncertainty(self):
        """Test evaluation with uncertainty estimates."""
        true_values = np.array([1, 2, 0, 4, 0, 6])
        imputed_values = np.array([1, 2, 3, 4, 5, 6])
        uncertainty = np.array([0, 0, 0.5, 0, 0.8, 0])
        mask = true_values == 0
        
        metrics = evaluate_imputation_quality(
            true_values, imputed_values, uncertainty, mask
        )
        
        assert 'mean_uncertainty' in metrics
        assert 'median_uncertainty' in metrics
        assert 'max_uncertainty' in metrics
        assert 'uncertainty_calibration' in metrics
    
    def test_evaluate_perfect_imputation(self):
        """Test evaluation with perfect imputation."""
        true_values = np.array([1, 2, 3, 4, 5, 6])
        imputed_values = np.array([1, 2, 3, 4, 5, 6])
        mask = np.array([False, False, True, True, False, False])
        
        metrics = evaluate_imputation_quality(
            true_values, imputed_values, mask=mask
        )
        
        assert metrics['rmse'] == 0.0
        assert metrics['mae'] == 0.0
        assert np.abs(metrics['r2'] - 1.0) < 1e-10


class TestQCReportGenerator:
    """Tests for QC report generator."""
    
    def test_generator_initialization(self):
        """Test QC report generator initialization."""
        qc = QCReportGenerator(dataset_name="Test Dataset")
        
        assert qc.dataset_name == "Test Dataset"
        assert len(qc.sections) == 0
    
    def test_add_imputation_metrics(self):
        """Test adding imputation metrics."""
        qc = QCReportGenerator()
        
        metrics = {
            'rmse': 0.5,
            'mae': 0.3,
            'pearson_r': 0.95
        }
        qc.add_imputation_metrics(metrics)
        
        assert 'imputation' in qc.sections
        assert qc.sections['imputation']['rmse'] == 0.5
    
    def test_add_quantification_metrics(self):
        """Test adding quantification metrics."""
        qc = QCReportGenerator()
        
        indices = {
            'RUD': np.array([0.5, 0.6, 0.7]),
            'PDUI': np.array([50, 60, 70])
        }
        qc.add_quantification_metrics(indices)
        
        assert 'quantification' in qc.sections
        assert 'RUD_mean' in qc.sections['quantification']
        assert 'PDUI_std' in qc.sections['quantification']
    
    def test_add_coverage_metrics(self):
        """Test adding coverage metrics."""
        qc = QCReportGenerator()
        
        qc.add_coverage_metrics(
            total_spots=100,
            observed_spots=30,
            imputed_spots=70
        )
        
        assert 'coverage' in qc.sections
        assert qc.sections['coverage']['total_spots'] == 100
        assert qc.sections['coverage']['observed_fraction'] == 0.3
    
    def test_generate_report_dict(self):
        """Test generating report as dict."""
        qc = QCReportGenerator(dataset_name="Test")
        qc.add_imputation_metrics({'rmse': 0.5})
        
        report = qc.generate_report(format='dict')
        
        assert isinstance(report, dict)
        assert report['dataset'] == "Test"
        assert 'imputation' in report
    
    def test_generate_report_dataframe(self):
        """Test generating report as DataFrame."""
        qc = QCReportGenerator(dataset_name="Test")
        qc.add_imputation_metrics({'rmse': 0.5, 'mae': 0.3})
        
        report = qc.generate_report(format='dataframe')
        
        assert isinstance(report, pd.DataFrame)
        assert len(report) == 1
        assert 'imputation_rmse' in report.columns
    
    def test_generate_report_text(self):
        """Test generating report as text."""
        qc = QCReportGenerator(dataset_name="Test")
        qc.add_imputation_metrics({'rmse': 0.5})
        
        report = qc.generate_report(format='text')
        
        assert isinstance(report, str)
        assert "Test" in report
        assert "IMPUTATION" in report
    
    def test_generate_report_invalid_format(self):
        """Test generating report with invalid format."""
        qc = QCReportGenerator()
        
        with pytest.raises(ValueError):
            qc.generate_report(format='invalid')
    
    def test_save_report_text(self, tmp_path):
        """Test saving report as text."""
        qc = QCReportGenerator(dataset_name="Test")
        qc.add_imputation_metrics({'rmse': 0.5})
        
        output_file = tmp_path / "report.txt"
        qc.save_report(str(output_file), format='text')
        
        assert output_file.exists()
        content = output_file.read_text()
        assert "Test" in content
    
    def test_save_report_dataframe(self, tmp_path):
        """Test saving report as CSV."""
        qc = QCReportGenerator(dataset_name="Test")
        qc.add_imputation_metrics({'rmse': 0.5})
        
        output_file = tmp_path / "report.csv"
        qc.save_report(str(output_file), format='dataframe')
        
        assert output_file.exists()
        df = pd.read_csv(output_file)
        assert len(df) == 1


class TestCompareImputationMethods:
    """Tests for comparing imputation methods."""
    
    def test_compare_methods_basic(self):
        """Test comparing multiple methods."""
        # Create mock imputers
        class MockImputer1:
            def fit(self, coords, values, mask):
                self.mean = np.mean(values[mask])
            def predict(self, coords):
                return np.full(len(coords), self.mean), None
        
        class MockImputer2:
            def fit(self, coords, values, mask):
                self.mean = np.mean(values[mask]) * 1.1
            def predict(self, coords):
                return np.full(len(coords), self.mean), None
        
        methods = {
            'method1': MockImputer1(),
            'method2': MockImputer2()
        }
        
        coordinates = np.random.rand(20, 2)
        values = np.random.rand(20) * 10
        
        results = compare_imputation_methods(
            methods, coordinates, values, n_folds=3
        )
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) == 2
        assert 'rmse' in results.columns
        assert 'mae' in results.columns


class TestEdgeCases:
    """Tests for edge cases."""
    
    def test_empty_arrays(self):
        """Test with empty arrays."""
        y_true = np.array([])
        y_pred = np.array([])
        
        # Should handle gracefully
        with pytest.raises((ValueError, IndexError)):
            calculate_rmse(y_true, y_pred)
    
    def test_single_value(self):
        """Test with single value."""
        y_true = np.array([5.0])
        y_pred = np.array([5.5])
        
        rmse = calculate_rmse(y_true, y_pred)
        mae = calculate_mae(y_true, y_pred)
        
        assert rmse == 0.5
        assert mae == 0.5
    
    def test_all_same_values(self):
        """Test with all same values."""
        y_true = np.array([5, 5, 5, 5, 5])
        y_pred = np.array([5, 5, 5, 5, 5])
        
        rmse = calculate_rmse(y_true, y_pred)
        r2 = calculate_r2(y_true, y_pred)
        
        assert rmse == 0.0
        # R² is undefined when variance is zero
        assert np.isnan(r2)

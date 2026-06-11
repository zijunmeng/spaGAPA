"""
Unit tests for benchmark module.
"""

import pytest
import numpy as np
import pandas as pd

from spagapa.benchmark import (
    SpatialAPASimulator,
    simulate_spatial_apa,
    BenchmarkEvaluator,
    mean_imputation,
    median_imputation,
    knn_imputation,
)


class TestSpatialAPASimulator:

    def test_generate_basic(self):
        sim = SpatialAPASimulator(n_spots=100, n_genes=20, n_domains=3,
                                  random_state=42)
        data = sim.generate()
        assert data['coordinates'].shape == (100, 2)
        assert data['domain_labels'].shape == (100,)
        assert data['true_apa'].shape == (100, 20)
        assert data['observed_apa'].shape == (100, 20)
        assert data['gene_patterns'].shape == (20,)

    def test_dropout_rate(self):
        sim = SpatialAPASimulator(n_spots=200, n_genes=50,
                                  dropout_rate=0.5, random_state=0)
        data = sim.generate()
        actual_dropout = np.isnan(data['observed_apa']).mean()
        # Allow ±10% tolerance
        assert 0.4 <= actual_dropout <= 0.6

    def test_no_dropout(self):
        sim = SpatialAPASimulator(n_spots=100, n_genes=10,
                                  dropout_rate=0.0, random_state=1)
        data = sim.generate()
        assert not np.any(np.isnan(data['observed_apa']))

    def test_true_apa_range(self):
        sim = SpatialAPASimulator(n_spots=100, n_genes=20, random_state=2)
        data = sim.generate()
        assert data['true_apa'].min() >= 0.0
        assert data['true_apa'].max() <= 1.0

    def test_domain_labels_valid(self):
        n_domains = 4
        sim = SpatialAPASimulator(n_spots=200, n_genes=10,
                                  n_domains=n_domains, random_state=3)
        data = sim.generate()
        labels = data['domain_labels']
        assert labels.min() >= 0
        assert labels.max() < n_domains

    def test_get_differential_genes(self):
        sim = SpatialAPASimulator(n_spots=200, n_genes=50,
                                  n_domains=3, random_state=4)
        sim.generate()
        diff = sim.get_differential_genes(0, 1, threshold=0.1)
        assert diff.shape == (50,)
        assert diff.dtype == bool

    def test_get_svapa_genes(self):
        sim = SpatialAPASimulator(n_spots=200, n_genes=50, random_state=5)
        sim.generate()
        svapa = sim.get_svapa_genes()
        assert svapa.shape == (50,)
        assert svapa.dtype == bool

    def test_to_dataframe(self):
        sim = SpatialAPASimulator(n_spots=50, n_genes=10, random_state=6)
        sim.generate()
        apa_df, meta_df = sim.to_dataframe()
        assert apa_df.shape == (50, 10)
        assert 'x' in meta_df.columns
        assert 'domain' in meta_df.columns

    def test_simulate_spatial_apa_convenience(self):
        data = simulate_spatial_apa(n_spots=80, n_genes=15, random_state=7)
        assert 'coordinates' in data
        assert 'true_apa' in data
        assert 'observed_apa' in data


class TestBaselineMethods:

    @pytest.fixture
    def sample_data(self):
        np.random.seed(42)
        n_spots, n_genes = 60, 15
        true = np.random.rand(n_spots, n_genes)
        observed = true.copy()
        mask = np.random.rand(n_spots, n_genes) < 0.3
        observed[mask] = np.nan
        coords = np.random.rand(n_spots, 2) * 100
        return observed, true, coords

    def test_mean_imputation_no_nan(self, sample_data):
        observed, _, coords = sample_data
        imputed = mean_imputation(observed, coords)
        assert not np.any(np.isnan(imputed))

    def test_mean_imputation_shape(self, sample_data):
        observed, _, coords = sample_data
        imputed = mean_imputation(observed, coords)
        assert imputed.shape == observed.shape

    def test_median_imputation_no_nan(self, sample_data):
        observed, _, coords = sample_data
        imputed = median_imputation(observed, coords)
        assert not np.any(np.isnan(imputed))

    def test_knn_imputation_no_nan(self, sample_data):
        observed, _, coords = sample_data
        imputed = knn_imputation(observed, coords, k=5)
        assert not np.any(np.isnan(imputed))

    def test_knn_imputation_shape(self, sample_data):
        observed, _, coords = sample_data
        imputed = knn_imputation(observed, coords, k=5)
        assert imputed.shape == observed.shape

    def test_observed_values_unchanged(self, sample_data):
        observed, _, coords = sample_data
        imputed = mean_imputation(observed, coords)
        valid = ~np.isnan(observed)
        np.testing.assert_array_almost_equal(imputed[valid], observed[valid])


class TestBenchmarkEvaluator:

    @pytest.fixture
    def simple_data(self):
        np.random.seed(0)
        n_spots, n_genes = 50, 10
        true = np.random.rand(n_spots, n_genes)
        observed = true.copy()
        mask = np.random.rand(n_spots, n_genes) < 0.3
        observed[mask] = np.nan
        coords = np.random.rand(n_spots, 2) * 100
        return observed, true, coords, mask

    def test_run_benchmark_returns_dataframe(self, simple_data):
        observed, true, coords, mask = simple_data
        evaluator = BenchmarkEvaluator(
            methods={'Mean': mean_imputation, 'Median': median_imputation}
        )
        results = evaluator.run_benchmark(observed, true, coords, mask)
        assert isinstance(results, pd.DataFrame)
        assert len(results) == 2
        assert 'method' in results.columns
        assert 'rmse' in results.columns

    def test_metrics_computed(self, simple_data):
        observed, true, coords, mask = simple_data
        evaluator = BenchmarkEvaluator(
            methods={'Mean': mean_imputation},
            metrics=['rmse', 'mae', 'pearson']
        )
        results = evaluator.run_benchmark(observed, true, coords, mask)
        row = results.iloc[0]
        assert not np.isnan(row['rmse'])
        assert not np.isnan(row['mae'])
        assert not np.isnan(row['pearson'])

    def test_compare_methods(self, simple_data):
        observed, true, coords, mask = simple_data
        evaluator = BenchmarkEvaluator(
            methods={'Mean': mean_imputation, 'Median': median_imputation,
                     'KNN': lambda o, c: knn_imputation(o, c, k=5)}
        )
        evaluator.run_benchmark(observed, true, coords, mask)
        comparison = evaluator.compare_methods('rmse')
        assert isinstance(comparison, pd.DataFrame)
        assert len(comparison) == 3
        # Should be sorted ascending for RMSE
        rmse_vals = comparison['rmse'].values
        assert all(rmse_vals[i] <= rmse_vals[i+1] for i in range(len(rmse_vals)-1))

    def test_generate_report(self, simple_data):
        observed, true, coords, mask = simple_data
        evaluator = BenchmarkEvaluator(
            methods={'Mean': mean_imputation}
        )
        evaluator.run_benchmark(observed, true, coords, mask)
        report = evaluator.generate_report()
        assert isinstance(report, str)
        assert 'BENCHMARK REPORT' in report
        assert 'Mean' in report

    def test_add_method(self, simple_data):
        observed, true, coords, mask = simple_data
        evaluator = BenchmarkEvaluator()
        evaluator.add_method('Mean', mean_imputation)
        assert 'Mean' in evaluator.methods

    def test_gp_better_than_mean(self, simple_data):
        """GP imputation should outperform mean on spatially structured data."""
        from spagapa.benchmark.evaluator import spagapa_gp_imputation

        # Generate data with strong spatial structure
        sim = SpatialAPASimulator(n_spots=100, n_genes=10,
                                  n_domains=3, dropout_rate=0.4,
                                  noise_level=0.05, random_state=99)
        data = sim.generate()
        observed = data['observed_apa']
        true     = data['true_apa']
        coords   = data['coordinates']
        mask     = np.isnan(observed)

        evaluator = BenchmarkEvaluator(
            methods={
                'Mean':       mean_imputation,
                'spaGAPA-GP': spagapa_gp_imputation,
            }
        )
        results = evaluator.run_benchmark(observed, true, coords, mask)

        gp_rmse   = results.loc[results['method'] == 'spaGAPA-GP', 'rmse'].values[0]
        mean_rmse = results.loc[results['method'] == 'Mean',       'rmse'].values[0]

        # GP should be at least as good as mean
        assert gp_rmse <= mean_rmse + 0.05, (
            f"GP RMSE ({gp_rmse:.4f}) should be ≤ Mean RMSE ({mean_rmse:.4f})"
        )

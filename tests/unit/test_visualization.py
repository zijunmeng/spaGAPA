"""
Unit tests for visualization module.
"""

import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing
import matplotlib.pyplot as plt

from spagapa.visualization import (
    SpatialPlotter,
    StatisticalPlotter,
    QCPlotter,
    plot_spatial_apa,
    plot_spatial_domains,
    plot_volcano,
    plot_heatmap,
    plot_imputation_quality
)


class TestSpatialPlotter:
    """Test SpatialPlotter class."""
    
    def test_init(self):
        """Test initialization."""
        plotter = SpatialPlotter(figsize=(10, 8), dpi=150)
        assert plotter.figsize == (10, 8)
        assert plotter.dpi == 150
        
    def test_plot_spatial_apa(self):
        """Test spatial APA plotting."""
        np.random.seed(42)
        coords = np.random.rand(100, 2) * 100
        apa_values = np.random.rand(100)
        
        plotter = SpatialPlotter()
        ax = plotter.plot_spatial_apa(coords, apa_values, gene_name="TEST1")
        
        assert ax is not None
        assert ax.get_title() == "TEST1 APA Usage"
        plt.close('all')
        
    def test_plot_spatial_apa_with_nan(self):
        """Test spatial APA plotting with NaN values."""
        np.random.seed(42)
        coords = np.random.rand(100, 2) * 100
        apa_values = np.random.rand(100)
        apa_values[:10] = np.nan
        
        plotter = SpatialPlotter()
        ax = plotter.plot_spatial_apa(coords, apa_values)
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_spatial_domains(self):
        """Test spatial domain plotting."""
        np.random.seed(42)
        coords = np.random.rand(100, 2) * 100
        domains = np.repeat([0, 1, 2, 3], 25)
        
        plotter = SpatialPlotter()
        ax = plotter.plot_spatial_domains(coords, domains)
        
        assert ax is not None
        assert ax.get_title() == "Spatial Domains"
        plt.close('all')
        
    def test_plot_apa_comparison(self):
        """Test multi-gene APA comparison."""
        np.random.seed(42)
        coords = np.random.rand(50, 2) * 100
        apa_matrix = np.random.rand(10, 50)
        genes = [0, 1, 2, 3]
        gene_names = [f"GENE{i}" for i in range(10)]
        
        plotter = SpatialPlotter()
        fig = plotter.plot_apa_comparison(
            coords, apa_matrix, genes, gene_names, ncols=2
        )
        
        assert fig is not None
        plt.close('all')


class TestStatisticalPlotter:
    """Test StatisticalPlotter class."""
    
    def test_init(self):
        """Test initialization."""
        plotter = StatisticalPlotter(figsize=(10, 8))
        assert plotter.figsize == (10, 8)
        
    def test_plot_volcano(self):
        """Test volcano plot."""
        np.random.seed(42)
        logfc = np.random.randn(100)
        pvalues = np.random.rand(100)
        gene_names = [f"GENE{i}" for i in range(100)]
        
        plotter = StatisticalPlotter()
        ax = plotter.plot_volcano(logfc, pvalues, gene_names)
        
        assert ax is not None
        assert ax.get_title() == "Volcano Plot"
        plt.close('all')
        
    def test_plot_heatmap(self):
        """Test heatmap plotting."""
        np.random.seed(42)
        matrix = np.random.rand(20, 10)
        row_labels = [f"Gene{i}" for i in range(20)]
        col_labels = [f"Spot{i}" for i in range(10)]
        
        plotter = StatisticalPlotter()
        fig, ax = plotter.plot_heatmap(
            matrix, row_labels, col_labels,
            cluster_rows=True, cluster_cols=True
        )
        
        assert fig is not None
        assert ax is not None
        plt.close('all')
        
    def test_plot_heatmap_no_clustering(self):
        """Test heatmap without clustering."""
        np.random.seed(42)
        matrix = np.random.rand(10, 5)
        
        plotter = StatisticalPlotter()
        fig, ax = plotter.plot_heatmap(
            matrix, cluster_rows=False, cluster_cols=False
        )
        
        assert fig is not None
        plt.close('all')
        
    def test_plot_boxplot(self):
        """Test box plot."""
        np.random.seed(42)
        apa_values = np.random.rand(100)
        groups = np.repeat([0, 1, 2, 3], 25)
        
        plotter = StatisticalPlotter()
        ax = plotter.plot_boxplot(apa_values, groups, gene_name="TEST1")
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_violin(self):
        """Test violin plot."""
        np.random.seed(42)
        apa_values = np.random.rand(100)
        groups = np.repeat(['A', 'B', 'C', 'D'], 25)
        
        plotter = StatisticalPlotter()
        ax = plotter.plot_violin(apa_values, groups, gene_name="TEST1")
        
        assert ax is not None
        plt.close('all')


class TestQCPlotter:
    """Test QCPlotter class."""
    
    def test_init(self):
        """Test initialization."""
        plotter = QCPlotter(figsize=(10, 8))
        assert plotter.figsize == (10, 8)
        
    def test_plot_imputation_quality(self):
        """Test imputation quality plot."""
        np.random.seed(42)
        observed = np.random.rand(100)
        imputed = observed + np.random.randn(100) * 0.1
        uncertainties = np.abs(np.random.randn(100) * 0.1)
        
        plotter = QCPlotter()
        fig = plotter.plot_imputation_quality(
            observed, imputed, uncertainties
        )
        
        assert fig is not None
        plt.close('all')
        
    def test_plot_imputation_quality_no_uncertainty(self):
        """Test imputation quality plot without uncertainty."""
        np.random.seed(42)
        observed = np.random.rand(100)
        imputed = observed + np.random.randn(100) * 0.1
        
        plotter = QCPlotter()
        fig = plotter.plot_imputation_quality(observed, imputed)
        
        assert fig is not None
        plt.close('all')
        
    def test_plot_spatial_support(self):
        """Test spatial support plot."""
        np.random.seed(42)
        coords = np.random.rand(100, 2) * 100
        support_scores = np.random.rand(100)
        
        plotter = QCPlotter()
        ax = plotter.plot_spatial_support(coords, support_scores)
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_dropout_stats(self):
        """Test dropout statistics plot."""
        np.random.seed(42)
        apa_matrix = np.random.rand(50, 100)
        # Add some NaN values
        apa_matrix[np.random.rand(50, 100) < 0.3] = np.nan
        coords = np.random.rand(100, 2) * 100
        
        plotter = QCPlotter()
        fig = plotter.plot_dropout_stats(apa_matrix, coords)
        
        assert fig is not None
        plt.close('all')


class TestConvenienceFunctions:
    """Test convenience functions."""
    
    def test_plot_spatial_apa(self):
        """Test plot_spatial_apa function."""
        np.random.seed(42)
        coords = np.random.rand(50, 2) * 100
        apa_values = np.random.rand(50)
        
        ax = plot_spatial_apa(coords, apa_values, gene_name="TEST1")
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_spatial_domains(self):
        """Test plot_spatial_domains function."""
        np.random.seed(42)
        coords = np.random.rand(100, 2) * 100
        domains = np.repeat([0, 1, 2], [30, 40, 30])
        
        ax = plot_spatial_domains(coords, domains)
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_volcano(self):
        """Test plot_volcano function."""
        np.random.seed(42)
        logfc = np.random.randn(50)
        pvalues = np.random.rand(50)
        
        ax = plot_volcano(logfc, pvalues)
        
        assert ax is not None
        plt.close('all')
        
    def test_plot_heatmap(self):
        """Test plot_heatmap function."""
        np.random.seed(42)
        matrix = np.random.rand(10, 8)
        
        fig, ax = plot_heatmap(matrix)
        
        assert fig is not None
        assert ax is not None
        plt.close('all')
        
    def test_plot_imputation_quality(self):
        """Test plot_imputation_quality function."""
        np.random.seed(42)
        observed = np.random.rand(50)
        imputed = observed + np.random.randn(50) * 0.1
        
        fig = plot_imputation_quality(observed, imputed)
        
        assert fig is not None
        plt.close('all')

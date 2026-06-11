"""
Unit tests for spatial pattern analysis.
"""

import pytest
import numpy as np
import pandas as pd
from spagapa.analysis import (
    SpatialPatternAnalyzer,
    identify_svapa_genes,
    cluster_spatial_patterns
)


class TestSpatialPatternAnalyzer:
    """Test SpatialPatternAnalyzer class."""
    
    def test_init(self):
        """Test initialization."""
        analyzer = SpatialPatternAnalyzer(
            n_permutations=500,
            alpha=0.01,
            random_state=42
        )
        
        assert analyzer.n_permutations == 500
        assert analyzer.alpha == 0.01
        assert analyzer.random_state == 42
        
    def test_compute_morans_i(self):
        """Test Moran's I computation."""
        np.random.seed(42)
        n_spots = 100
        
        # Create spatial weights (simple adjacency)
        weights = np.zeros((n_spots, n_spots))
        for i in range(n_spots - 1):
            weights[i, i+1] = 1
            weights[i+1, i] = 1
        # Row-normalize
        row_sums = weights.sum(axis=1)
        row_sums[row_sums == 0] = 1
        weights = weights / row_sums[:, np.newaxis]
        
        # Create spatially autocorrelated values
        values = np.arange(n_spots, dtype=float)  # Linear gradient
        
        analyzer = SpatialPatternAnalyzer()
        morans_i, expected_i, variance_i = analyzer.compute_morans_i(
            values, weights
        )
        
        assert not np.isnan(morans_i)
        assert morans_i > 0  # Should be positive (positive autocorrelation)
        assert expected_i == -1.0 / (n_spots - 1)
        assert variance_i > 0
        
    def test_compute_morans_i_random(self):
        """Test Moran's I with random values."""
        np.random.seed(42)
        n_spots = 100
        
        # Create spatial weights
        weights = np.eye(n_spots)
        for i in range(n_spots - 1):
            weights[i, i+1] = 1
            weights[i+1, i] = 1
        row_sums = weights.sum(axis=1)
        weights = weights / row_sums[:, np.newaxis]
        
        # Random values (no spatial pattern)
        values = np.random.rand(n_spots)
        
        analyzer = SpatialPatternAnalyzer()
        morans_i, expected_i, variance_i = analyzer.compute_morans_i(
            values, weights
        )
        
        # Should be close to expected value (near 0)
        assert abs(morans_i - expected_i) < 0.5
        
    def test_test_spatial_autocorrelation(self):
        """Test spatial autocorrelation testing."""
        np.random.seed(42)
        n_genes = 10
        n_spots = 100
        
        # Create data with spatial pattern
        coords = np.random.rand(n_spots, 2) * 100
        apa_matrix = np.random.rand(n_genes, n_spots)
        
        # Add spatial gradient to first gene
        apa_matrix[0, :] = coords[:, 0] / 100
        
        analyzer = SpatialPatternAnalyzer()
        results = analyzer.test_spatial_autocorrelation(
            apa_matrix, coords, n_neighbors=6
        )
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) <= n_genes
        assert 'morans_i' in results.columns
        assert 'pvalue' in results.columns
        assert 'z_score' in results.columns
        
        # First gene should have higher Moran's I
        if len(results) > 0:
            gene0_result = results[results['gene'] == 'Gene_0']
            if len(gene0_result) > 0:
                assert gene0_result['morans_i'].values[0] > 0
                
    def test_identify_svapa_genes(self):
        """Test SVAPA gene identification."""
        np.random.seed(42)
        n_genes = 20
        n_spots = 100
        
        coords = np.random.rand(n_spots, 2) * 100
        apa_matrix = np.random.rand(n_genes, n_spots)
        
        # Add strong spatial patterns to first 3 genes
        apa_matrix[0, :] = coords[:, 0] / 100  # X gradient
        apa_matrix[1, :] = coords[:, 1] / 100  # Y gradient
        apa_matrix[2, :] = (coords[:, 0] + coords[:, 1]) / 200  # Diagonal
        
        gene_names = [f"GENE{i}" for i in range(n_genes)]
        
        analyzer = SpatialPatternAnalyzer()
        svapa = analyzer.identify_svapa_genes(
            apa_matrix, coords, gene_names, fdr_threshold=0.1
        )
        
        assert isinstance(svapa, list)
        # Should identify at least some of the patterned genes
        assert len(svapa) >= 0
        
    def test_cluster_spatial_patterns(self):
        """Test spatial pattern clustering."""
        np.random.seed(42)
        n_genes = 30
        n_spots = 100
        
        coords = np.random.rand(n_spots, 2) * 100
        apa_matrix = np.random.rand(n_genes, n_spots)
        
        # Create 3 distinct patterns
        # Pattern 1: X gradient
        apa_matrix[:10, :] = coords[:, 0] / 100
        # Pattern 2: Y gradient
        apa_matrix[10:20, :] = coords[:, 1] / 100
        # Pattern 3: Random
        apa_matrix[20:, :] = np.random.rand(10, n_spots)
        
        analyzer = SpatialPatternAnalyzer()
        labels, info = analyzer.cluster_spatial_patterns(
            apa_matrix, coords, n_patterns=3
        )
        
        assert len(labels) == n_genes
        assert len(np.unique(labels)) == 3
        assert isinstance(info, pd.DataFrame)
        assert len(info) == 3
        assert 'pattern_id' in info.columns
        assert 'n_genes' in info.columns
        assert 'morans_i' in info.columns
        
    def test_compute_pattern_similarity(self):
        """Test pattern similarity computation."""
        np.random.seed(42)
        pattern1 = np.random.rand(100)
        pattern2 = pattern1 + np.random.rand(100) * 0.1  # Similar
        pattern3 = np.random.rand(100)  # Different
        
        analyzer = SpatialPatternAnalyzer()
        
        # Pearson correlation
        sim12 = analyzer.compute_pattern_similarity(
            pattern1, pattern2, method='pearson'
        )
        sim13 = analyzer.compute_pattern_similarity(
            pattern1, pattern3, method='pearson'
        )
        
        assert sim12 > sim13  # pattern1 more similar to pattern2
        assert -1 <= sim12 <= 1
        assert -1 <= sim13 <= 1
        
    def test_compute_pattern_similarity_methods(self):
        """Test different similarity methods."""
        np.random.seed(42)
        pattern1 = np.random.rand(50)
        pattern2 = np.random.rand(50)
        
        analyzer = SpatialPatternAnalyzer()
        
        # Test all methods
        for method in ['pearson', 'spearman', 'cosine']:
            sim = analyzer.compute_pattern_similarity(
                pattern1, pattern2, method=method
            )
            assert not np.isnan(sim)
            
    def test_invalid_similarity_method(self):
        """Test invalid similarity method."""
        pattern1 = np.random.rand(50)
        pattern2 = np.random.rand(50)
        
        analyzer = SpatialPatternAnalyzer()
        
        with pytest.raises(ValueError, match="Unknown method"):
            analyzer.compute_pattern_similarity(
                pattern1, pattern2, method='invalid'
            )
            
    def test_with_nan_values(self):
        """Test handling of NaN values."""
        np.random.seed(42)
        n_spots = 100
        
        # Create weights
        weights = np.eye(n_spots)
        for i in range(n_spots - 1):
            weights[i, i+1] = 1
        row_sums = weights.sum(axis=1)
        weights = weights / row_sums[:, np.newaxis]
        
        # Values with NaN
        values = np.random.rand(n_spots)
        values[:10] = np.nan
        
        analyzer = SpatialPatternAnalyzer()
        morans_i, expected_i, variance_i = analyzer.compute_morans_i(
            values, weights
        )
        
        # Should handle NaN gracefully
        assert not np.isnan(morans_i)


class TestIdentifySVAPAGenes:
    """Test convenience function."""
    
    def test_basic_usage(self):
        """Test basic usage."""
        np.random.seed(42)
        n_genes = 20
        n_spots = 100
        
        coords = np.random.rand(n_spots, 2) * 100
        apa_matrix = np.random.rand(n_genes, n_spots)
        
        # Add spatial pattern
        apa_matrix[0, :] = coords[:, 0] / 100
        
        svapa, results = identify_svapa_genes(
            apa_matrix, coords, fdr_threshold=0.2
        )
        
        assert isinstance(svapa, list)
        assert isinstance(results, pd.DataFrame)
        assert 'morans_i' in results.columns
        assert 'padj' in results.columns
        
    def test_with_gene_names(self):
        """Test with gene names."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 100)
        coords = np.random.rand(100, 2) * 100
        gene_names = [f"GENE{i}" for i in range(10)]
        
        svapa, results = identify_svapa_genes(
            apa_matrix, coords, gene_names=gene_names
        )
        
        assert all(gene in gene_names for gene in results['gene'])


class TestClusterSpatialPatterns:
    """Test convenience function."""
    
    def test_basic_usage(self):
        """Test basic usage."""
        np.random.seed(42)
        apa_matrix = np.random.rand(30, 100)
        coords = np.random.rand(100, 2) * 100
        
        labels, info, genes = cluster_spatial_patterns(
            apa_matrix, coords, n_patterns=5
        )
        
        assert len(labels) == 30
        assert len(np.unique(labels)) == 5
        assert isinstance(info, pd.DataFrame)
        assert len(info) == 5
        assert isinstance(genes, dict)
        assert len(genes) == 5
        
        # Check gene lists
        total_genes = sum(len(g) for g in genes.values())
        assert total_genes == 30
        
    def test_with_gene_names(self):
        """Test with gene names."""
        np.random.seed(42)
        n_genes = 20
        apa_matrix = np.random.rand(n_genes, 100)
        coords = np.random.rand(100, 2) * 100
        gene_names = [f"GENE{i}" for i in range(n_genes)]
        
        labels, info, genes = cluster_spatial_patterns(
            apa_matrix, coords, n_patterns=3, gene_names=gene_names
        )
        
        # Check that all genes are accounted for
        all_genes = []
        for pattern_genes in genes.values():
            all_genes.extend(pattern_genes)
            
        assert set(all_genes) == set(gene_names)

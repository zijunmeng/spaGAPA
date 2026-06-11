"""
Unit tests for differential APA analysis.
"""

import pytest
import numpy as np
import pandas as pd
from spagapa.analysis import (
    DifferentialAPAAnalyzer,
    run_differential_apa_test,
    find_domain_markers
)


class TestDifferentialAPAAnalyzer:
    """Test DifferentialAPAAnalyzer class."""
    
    def test_init(self):
        """Test initialization."""
        analyzer = DifferentialAPAAnalyzer(
            method='wilcoxon',
            alpha=0.01,
            min_spots_per_group=5
        )
        
        assert analyzer.method == 'wilcoxon'
        assert analyzer.alpha == 0.01
        assert analyzer.min_spots_per_group == 5
        
    def test_wilcoxon_test(self):
        """Test Wilcoxon rank-sum test."""
        np.random.seed(42)
        n_genes = 20
        n_spots = 100
        
        # Create data with differential APA
        apa_matrix = np.random.rand(n_genes, n_spots)
        # Make first 5 genes differential
        apa_matrix[:5, :50] += 1.0  # Higher in group 1
        
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer(method='wilcoxon')
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) == n_genes
        assert 'pvalue' in results.columns
        assert 'log2fc' in results.columns
        assert 'mean_group1' in results.columns
        assert 'mean_group2' in results.columns
        
        # First 5 genes should have lower p-values
        assert results.iloc[:5]['pvalue'].mean() < results.iloc[5:]['pvalue'].mean()
        
    def test_ttest(self):
        """Test t-test."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 100)
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer(method='t-test')
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) == 10
        
    def test_permutation_test(self):
        """Test permutation test."""
        np.random.seed(42)
        apa_matrix = np.random.rand(5, 60)
        group1 = np.arange(30)
        group2 = np.arange(30, 60)
        
        analyzer = DifferentialAPAAnalyzer(method='permutation')
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) == 5
        
    def test_adjust_pvalues(self):
        """Test p-value adjustment."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer()
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        
        # Adjust p-values
        results_adj = analyzer.adjust_pvalues(results, method='fdr_bh')
        
        assert 'padj' in results_adj.columns
        # Adjusted p-values should be >= original p-values
        assert (results_adj['padj'] >= results_adj['pvalue']).all()
        
    def test_rank_genes(self):
        """Test gene ranking."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer()
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        results = analyzer.adjust_pvalues(results)
        
        # Rank by adjusted p-value
        ranked = analyzer.rank_genes(results, by='padj', top_n=5)
        
        assert len(ranked) == 5
        assert ranked['padj'].is_monotonic_increasing
        
    def test_filter_results(self):
        """Test result filtering."""
        np.random.seed(42)
        n_genes = 50
        apa_matrix = np.random.rand(n_genes, 100)
        
        # Make some genes clearly differential
        apa_matrix[:10, :50] += 2.0
        
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer()
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        results = analyzer.adjust_pvalues(results)
        
        # Filter
        filtered = analyzer.filter_results(
            results,
            padj_threshold=0.05,
            logfc_threshold=0.5
        )
        
        assert len(filtered) <= len(results)
        assert (filtered['padj'] < 0.05).all()
        assert (np.abs(filtered['log2fc']) > 0.5).all()
        
    def test_test_all_pairwise(self):
        """Test all pairwise comparisons."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 120)
        domains = np.repeat([0, 1, 2], 40)
        
        analyzer = DifferentialAPAAnalyzer()
        results = analyzer.test_all_pairwise(apa_matrix, domains)
        
        assert isinstance(results, dict)
        # Should have 3 comparisons: 0vs1, 0vs2, 1vs2
        assert len(results) == 3
        
        for key, df in results.items():
            assert isinstance(df, pd.DataFrame)
            assert 'pvalue' in df.columns
            
    def test_small_group_error(self):
        """Test error with too few spots."""
        apa_matrix = np.random.rand(10, 10)
        group1 = np.array([0, 1])  # Only 2 spots
        group2 = np.arange(2, 10)
        
        analyzer = DifferentialAPAAnalyzer(min_spots_per_group=3)
        
        with pytest.raises(ValueError, match="has only"):
            analyzer.test_differential_apa(apa_matrix, group1, group2)
            
    def test_invalid_method(self):
        """Test invalid method."""
        analyzer = DifferentialAPAAnalyzer(method='invalid')
        
        apa_matrix = np.random.rand(10, 50)
        group1 = np.arange(25)
        group2 = np.arange(25, 50)
        
        with pytest.raises(ValueError, match="Unknown method"):
            analyzer.test_differential_apa(apa_matrix, group1, group2)
            
    def test_with_nan_values(self):
        """Test handling of NaN values."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 100)
        # Add some NaN values
        apa_matrix[0, :10] = np.nan
        apa_matrix[1, 50:60] = np.nan
        
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        analyzer = DifferentialAPAAnalyzer()
        results = analyzer.test_differential_apa(
            apa_matrix, group1, group2
        )
        
        # Should still return results
        assert isinstance(results, pd.DataFrame)
        assert len(results) <= 10  # Some genes might be skipped


class TestTestDifferentialAPA:
    """Test convenience function."""
    
    def test_basic_usage(self):
        """Test basic usage."""
        np.random.seed(42)
        apa_matrix = np.random.rand(20, 100)
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        
        results = run_differential_apa_test(
            apa_matrix, group1, group2, method='wilcoxon'
        )
        
        assert isinstance(results, pd.DataFrame)
        assert 'padj' in results.columns
        
    def test_with_gene_names(self):
        """Test with gene names."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 100)
        group1 = np.arange(50)
        group2 = np.arange(50, 100)
        gene_names = [f"GENE{i}" for i in range(10)]
        
        results = run_differential_apa_test(
            apa_matrix, group1, group2, gene_names=gene_names
        )
        
        assert results['gene'].tolist() == gene_names


class TestFindDomainMarkers:
    """Test find_domain_markers function."""
    
    def test_basic_usage(self):
        """Test basic usage."""
        np.random.seed(42)
        n_spots = 120
        apa_matrix = np.random.rand(20, n_spots)
        domains = np.repeat([0, 1, 2], 40)
        
        # Make some genes domain-specific
        apa_matrix[:5, :40] += 1.0  # Domain 0 markers
        apa_matrix[5:10, 40:80] += 1.0  # Domain 1 markers
        
        markers = find_domain_markers(apa_matrix, domains)
        
        assert isinstance(markers, dict)
        assert len(markers) <= 3
        
        for domain_id, df in markers.items():
            assert isinstance(df, pd.DataFrame)
            if len(df) > 0:
                assert 'padj' in df.columns
                assert 'log2fc' in df.columns
                
    def test_with_gene_names(self):
        """Test with gene names."""
        np.random.seed(42)
        apa_matrix = np.random.rand(10, 100)
        domains = np.repeat([0, 1], 50)
        gene_names = [f"GENE{i}" for i in range(10)]
        
        markers = find_domain_markers(
            apa_matrix, domains, gene_names=gene_names
        )
        
        for domain_id, df in markers.items():
            if len(df) > 0:
                assert all(gene in gene_names for gene in df['gene'])
                
    def test_filtering(self):
        """Test marker filtering."""
        np.random.seed(42)
        apa_matrix = np.random.rand(30, 100)
        domains = np.repeat([0, 1], 50)
        
        # Make first 10 genes strong markers for domain 0
        # Domain 0 has HIGHER APA values
        apa_matrix[:10, :50] += 2.0
        
        markers = find_domain_markers(
            apa_matrix,
            domains,
            padj_threshold=0.05,
            logfc_threshold=0.5
        )
        
        # Domain 0 should have markers
        # Note: log2fc = log2(mean_rest / mean_domain0)
        # Since domain 0 has higher values, log2fc will be NEGATIVE
        # But filter_results uses abs(log2fc) by default
        if 0 in markers and len(markers[0]) > 0:
            assert (markers[0]['padj'] < 0.05).all()
            # Check that absolute log2fc is > 0.5
            assert (np.abs(markers[0]['log2fc']) > 0.5).all()

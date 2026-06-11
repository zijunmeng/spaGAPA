"""
Unit tests for quality filtering.
"""

import pytest
import numpy as np
import pandas as pd
from spagapa.calling import QualityFilter, filter_apa_sites
from spagapa.core import APASite, APASiteCollection


class TestQualityFilter:
    """Test QualityFilter class."""
    
    @pytest.fixture
    def sample_apa_counts(self):
        """Create sample APA count matrix."""
        # 5 genes x 10 spots
        counts = np.array([
            [10, 8, 9, 10, 8, 9, 10, 8, 9, 10],  # Gene 1: high counts, many spots
            [5, 0, 0, 5, 0, 0, 5, 0, 0, 0],      # Gene 2: medium counts, few spots
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],      # Gene 3: low counts, many spots
            [100, 0, 0, 0, 0, 0, 0, 0, 0, 0],    # Gene 4: high count, single spot
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]       # Gene 5: no counts
        ])
        return counts
    
    @pytest.fixture
    def sample_sites(self):
        """Create sample APA sites."""
        sites = [
            APASite(chr='chr1', start=i*1000, end=i*1000+50, strand='+', gene_id=f'gene{i}')
            for i in range(5)
        ]
        return APASiteCollection(sites)
    
    def test_initialization(self):
        """Test filter initialization."""
        qf = QualityFilter(
            min_read_count=10,
            min_spots=3,
            min_spatial_support=0.3
        )
        assert qf.min_read_count == 10
        assert qf.min_spots == 3
        assert qf.min_spatial_support == 0.3
    
    def test_compute_site_metrics(self, sample_apa_counts):
        """Test computing metrics for a site."""
        qf = QualityFilter()
        
        metrics = qf.compute_site_metrics(sample_apa_counts, site_idx=0)
        
        assert 'total_count' in metrics
        assert 'n_spots' in metrics
        assert 'mean_count' in metrics
        assert 'median_count' in metrics
        assert 'max_count' in metrics
        assert 'cv' in metrics
        
        # Check values for gene 1
        assert metrics['total_count'] == 91
        assert metrics['n_spots'] == 10
        assert metrics['mean_count'] > 0
    
    def test_filter_by_read_count(self, sample_apa_counts):
        """Test filtering by read count."""
        qf = QualityFilter(min_read_count=20)
        
        mask = qf.filter_by_read_count(sample_apa_counts)
        
        assert len(mask) == 5
        assert mask[0] == True   # Gene 1: 91 reads
        assert mask[1] == False  # Gene 2: 15 reads
        assert mask[4] == False  # Gene 5: 0 reads
    
    def test_filter_by_spot_count(self, sample_apa_counts):
        """Test filtering by spot count."""
        qf = QualityFilter(min_spots=5)
        
        mask = qf.filter_by_spot_count(sample_apa_counts)
        
        assert len(mask) == 5
        assert mask[0] == True   # Gene 1: 10 spots
        assert mask[1] == False  # Gene 2: 3 spots
        assert mask[3] == False  # Gene 4: 1 spot
    
    def test_filter_by_mean_count(self, sample_apa_counts):
        """Test filtering by mean count."""
        qf = QualityFilter(min_mean_count=5.0)
        
        mask = qf.filter_by_mean_count(sample_apa_counts)
        
        assert len(mask) == 5
        assert mask[0] == True   # Gene 1: mean ~9
        assert mask[1] == True   # Gene 2: mean 5
        assert mask[2] == False  # Gene 3: mean 1
    
    def test_filter_by_cv(self, sample_apa_counts):
        """Test filtering by coefficient of variation."""
        qf = QualityFilter(max_cv=0.5)
        
        mask = qf.filter_by_cv(sample_apa_counts)
        
        assert len(mask) == 5
        # Gene 1 has low CV (consistent counts)
        # Gene 4 has zero CV (single value)
        assert mask[0] == True or mask[3] == True
    
    def test_filter_sites(self, sample_sites, sample_apa_counts):
        """Test filtering sites with all filters."""
        qf = QualityFilter(
            min_read_count=10,
            min_spots=3,
            min_mean_count=2.0
        )
        
        filtered_sites, pass_mask = qf.filter_sites(sample_sites, sample_apa_counts)
        
        assert isinstance(filtered_sites, APASiteCollection)
        assert len(pass_mask) == len(sample_sites)
        assert len(filtered_sites) <= len(sample_sites)
        
        # Gene 1 should pass all filters
        assert pass_mask[0] == True
        
        # Gene 5 (no counts) should fail
        assert pass_mask[4] == False
    
    def test_filter_sites_with_spatial_support(self, sample_sites, sample_apa_counts):
        """Test filtering with spatial support scores."""
        qf = QualityFilter(
            min_read_count=5,
            min_spots=1,
            min_spatial_support=0.5
        )
        
        # Create mock spatial support scores
        spatial_support = np.array([0.8, 0.6, 0.3, 0.9, 0.0])
        
        filtered_sites, pass_mask = qf.filter_sites(
            sample_sites,
            sample_apa_counts,
            spatial_support=spatial_support
        )
        
        # Sites with low spatial support should be filtered
        assert pass_mask[2] == False  # support = 0.3 < 0.5
        assert pass_mask[4] == False  # support = 0.0 < 0.5
    
    def test_generate_qc_report(self, sample_sites, sample_apa_counts):
        """Test generating QC report."""
        qf = QualityFilter(
            min_read_count=10,
            min_spots=3
        )
        
        report = qf.generate_qc_report(sample_sites, sample_apa_counts)
        
        assert isinstance(report, pd.DataFrame)
        assert len(report) == len(sample_sites)
        
        # Check required columns
        assert 'chr' in report.columns
        assert 'gene_id' in report.columns
        assert 'total_count' in report.columns
        assert 'n_spots' in report.columns
        assert 'pass_read_count' in report.columns
        assert 'pass_spot_count' in report.columns
        assert 'pass_all' in report.columns
    
    def test_generate_qc_report_with_spatial_support(self, sample_sites, sample_apa_counts):
        """Test QC report with spatial support."""
        qf = QualityFilter(min_spatial_support=0.5)
        
        spatial_support = np.array([0.8, 0.6, 0.3, 0.9, 0.0])
        
        report = qf.generate_qc_report(
            sample_sites,
            sample_apa_counts,
            spatial_support=spatial_support
        )
        
        assert 'spatial_support' in report.columns
        assert 'pass_spatial_support' in report.columns
        
        # Check spatial support values
        np.testing.assert_array_equal(
            report['spatial_support'].values,
            spatial_support
        )


class TestHelperFunction:
    """Test helper function."""
    
    @pytest.fixture
    def sample_sites(self):
        """Create sample sites."""
        sites = [
            APASite(chr='chr1', start=i*1000, end=i*1000+50, strand='+', gene_id=f'gene{i}')
            for i in range(3)
        ]
        return APASiteCollection(sites)
    
    @pytest.fixture
    def sample_counts(self):
        """Create sample counts."""
        return np.array([
            [10, 8, 9, 10, 8],
            [5, 0, 0, 5, 0],
            [1, 1, 1, 1, 1]
        ])
    
    def test_filter_apa_sites(self, sample_sites, sample_counts):
        """Test convenience function."""
        filtered_sites, pass_mask = filter_apa_sites(
            sample_sites,
            sample_counts,
            min_read_count=10,
            min_spots=3
        )
        
        assert isinstance(filtered_sites, APASiteCollection)
        assert len(pass_mask) == len(sample_sites)
    
    def test_filter_apa_sites_with_spatial_support(self, sample_sites, sample_counts):
        """Test convenience function with spatial support."""
        spatial_support = np.array([0.8, 0.2, 0.6])
        
        filtered_sites, pass_mask = filter_apa_sites(
            sample_sites,
            sample_counts,
            min_read_count=5,
            min_spots=2,
            min_spatial_support=0.5,
            spatial_support=spatial_support
        )
        
        # Site 2 should be filtered (support = 0.2 < 0.5)
        assert pass_mask[1] == False


class TestEdgeCases:
    """Test edge cases."""
    
    def test_all_zeros(self):
        """Test with all-zero count matrix."""
        qf = QualityFilter(min_read_count=1)
        
        counts = np.zeros((3, 5))
        mask = qf.filter_by_read_count(counts)
        
        # All should fail
        assert np.all(mask == False)
    
    def test_single_site(self):
        """Test with single site."""
        qf = QualityFilter()
        
        counts = np.array([[10, 8, 9]])
        mask = qf.filter_by_read_count(counts)
        
        assert len(mask) == 1
    
    def test_high_cv_filtering(self):
        """Test filtering highly variable sites."""
        qf = QualityFilter(max_cv=0.5)
        
        # Create counts with high CV
        counts = np.array([
            [1, 1, 1, 1, 1],      # Low CV
            [1, 10, 1, 10, 1]     # High CV
        ])
        
        mask = qf.filter_by_cv(counts)
        
        # Low CV should pass, high CV should fail
        assert mask[0] == True
        assert mask[1] == False
    
    def test_empty_site_collection(self):
        """Test with empty site collection."""
        qf = QualityFilter()
        
        empty_sites = APASiteCollection([])
        empty_counts = np.array([]).reshape(0, 5)
        
        filtered_sites, pass_mask = qf.filter_sites(empty_sites, empty_counts)
        
        assert len(filtered_sites) == 0
        assert len(pass_mask) == 0
    
    def test_metrics_with_single_nonzero(self):
        """Test metrics computation with single non-zero value."""
        qf = QualityFilter()
        
        counts = np.array([[10, 0, 0, 0, 0]])
        metrics = qf.compute_site_metrics(counts, 0)
        
        assert metrics['n_spots'] == 1
        assert metrics['mean_count'] == 10
        assert metrics['cv'] == 0  # Single value has zero variance

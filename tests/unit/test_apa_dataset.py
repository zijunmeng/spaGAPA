"""
Unit tests for APADataset class.
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch

# Mock anndata if not available
try:
    import anndata
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False


@pytest.mark.skipif(not ANNDATA_AVAILABLE, reason="anndata not installed")
class TestAPADataset:
    """Test suite for APADataset class."""
    
    def test_from_counts(self):
        """Test creating dataset from count matrix."""
        from spagapa.core import APADataset
        
        # Create test data
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        spot_names = [f"Spot_{i}" for i in range(n_spots)]
        
        # Create dataset
        dataset = APADataset.from_counts(
            apa_counts=apa_counts,
            spatial_coords=spatial_coords,
            gene_names=gene_names,
            spot_names=spot_names
        )
        
        # Check dimensions
        assert dataset.n_genes == n_genes
        assert dataset.n_spots == n_spots
        assert len(dataset.gene_names) == n_genes
        assert len(dataset.spot_names) == n_spots
    
    def test_from_dataframe(self):
        """Test creating dataset from DataFrame."""
        from spagapa.core import APADataset
        
        # Create test data
        n_genes, n_spots = 5, 10
        apa_counts = pd.DataFrame(
            np.random.poisson(5, (n_genes, n_spots)),
            index=[f"Gene_{i}" for i in range(n_genes)],
            columns=[f"Spot_{i}" for i in range(n_spots)]
        )
        spatial_coords = pd.DataFrame(
            np.random.rand(n_spots, 2) * 100,
            index=[f"Spot_{i}" for i in range(n_spots)],
            columns=['x', 'y']
        )
        
        # Create dataset
        dataset = APADataset.from_counts(
            apa_counts=apa_counts,
            spatial_coords=spatial_coords
        )
        
        assert dataset.n_genes == n_genes
        assert dataset.n_spots == n_spots
    
    def test_get_apa_counts(self):
        """Test retrieving APA counts."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        
        dataset = APADataset.from_counts(apa_counts, spatial_coords)
        
        # Get raw counts
        counts = dataset.get_apa_counts(imputed=False)
        assert counts.shape == (n_genes, n_spots)
        np.testing.assert_array_equal(counts, apa_counts)
    
    def test_add_imputation(self):
        """Test adding imputation results."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        
        dataset = APADataset.from_counts(apa_counts, spatial_coords)
        
        # Add imputation
        imputed = apa_counts + np.random.randn(n_genes, n_spots) * 0.1
        uncertainty = np.random.rand(n_genes, n_spots) * 0.5
        
        dataset.add_imputation(imputed, uncertainty)
        
        # Check imputed counts
        imputed_counts = dataset.get_apa_counts(imputed=True)
        np.testing.assert_array_almost_equal(imputed_counts, imputed)
    
    def test_subset_genes(self):
        """Test subsetting by genes."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        gene_names = [f"Gene_{i}" for i in range(n_genes)]
        
        dataset = APADataset.from_counts(
            apa_counts, spatial_coords, gene_names=gene_names
        )
        
        # Subset to first 5 genes
        subset = dataset.subset_genes(gene_names[:5])
        
        assert subset.n_genes == 5
        assert subset.n_spots == n_spots
    
    def test_subset_spots(self):
        """Test subsetting by spots."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        spot_names = [f"Spot_{i}" for i in range(n_spots)]
        
        dataset = APADataset.from_counts(
            apa_counts, spatial_coords, spot_names=spot_names
        )
        
        # Subset to first 10 spots
        subset = dataset.subset_spots(spot_names[:10])
        
        assert subset.n_genes == n_genes
        assert subset.n_spots == 10
    
    def test_save_load(self, tmp_path):
        """Test saving and loading dataset."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        
        dataset = APADataset.from_counts(apa_counts, spatial_coords)
        
        # Save
        filename = tmp_path / "test_dataset.h5ad"
        dataset.save(str(filename))
        
        # Load
        loaded = APADataset.load(str(filename))
        
        assert loaded.n_genes == n_genes
        assert loaded.n_spots == n_spots
    
    def test_repr(self):
        """Test string representation."""
        from spagapa.core import APADataset
        
        n_genes, n_spots = 10, 20
        apa_counts = np.random.poisson(5, (n_genes, n_spots))
        spatial_coords = np.random.rand(n_spots, 2) * 100
        
        dataset = APADataset.from_counts(apa_counts, spatial_coords)
        
        repr_str = repr(dataset)
        assert "APADataset" in repr_str
        assert f"n_spots={n_spots}" in repr_str
        assert f"n_genes={n_genes}" in repr_str

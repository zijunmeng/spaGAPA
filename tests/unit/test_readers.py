"""
Unit tests for readers module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path


class TestSpatialCoordinateReader:
    """Test suite for SpatialCoordinateReader."""
    
    def test_read_csv(self, tmp_path):
        """Test reading CSV coordinates."""
        from spagapa.io import SpatialCoordinateReader
        
        # Create test CSV
        csv_file = tmp_path / "coords.csv"
        df = pd.DataFrame({
            'barcode': ['SPOT_1', 'SPOT_2', 'SPOT_3'],
            'x': [100.0, 105.0, 110.0],
            'y': [200.0, 205.0, 210.0]
        })
        df.to_csv(csv_file, index=False)
        
        # Read coordinates
        coords = SpatialCoordinateReader.read_coordinates(str(csv_file))
        
        assert len(coords) == 3
        assert 'x' in coords.columns
        assert 'y' in coords.columns
        assert coords.index.name == 'barcode'
    
    def test_read_tsv(self, tmp_path):
        """Test reading TSV coordinates."""
        from spagapa.io import SpatialCoordinateReader
        
        # Create test TSV
        tsv_file = tmp_path / "coords.tsv"
        df = pd.DataFrame({
            'barcode': ['SPOT_1', 'SPOT_2'],
            'x': [100.0, 105.0],
            'y': [200.0, 205.0]
        })
        df.to_csv(tsv_file, sep='\t', index=False)
        
        # Read coordinates
        coords = SpatialCoordinateReader.read_coordinates(str(tsv_file))
        
        assert len(coords) == 2
    
    def test_read_10x_positions(self, tmp_path):
        """Test reading 10x Visium positions file."""
        from spagapa.io import SpatialCoordinateReader
        
        # Create test positions file
        positions_file = tmp_path / "tissue_positions.csv"
        df = pd.DataFrame({
            'barcode': ['SPOT_1', 'SPOT_2', 'SPOT_3'],
            'in_tissue': [1, 1, 0],
            'array_row': [0, 1, 2],
            'array_col': [0, 1, 2],
            'x': [100.0, 105.0, 110.0],
            'y': [200.0, 205.0, 210.0]
        })
        df.to_csv(positions_file, header=False, index=False)
        
        # Read positions
        coords = SpatialCoordinateReader.read_10x_positions(str(positions_file))
        
        assert len(coords) == 2  # Only in_tissue spots
        assert 'x' in coords.columns
        assert 'y' in coords.columns


class TestBEDReader:
    """Test suite for BEDReader."""
    
    def test_read_bed(self, tmp_path):
        """Test reading BED file."""
        from spagapa.io import BEDReader
        
        # Create test BED file
        bed_file = tmp_path / "test.bed"
        with open(bed_file, 'w') as f:
            f.write("chr1\t1000\t1050\tGENE1\t800\t+\n")
            f.write("chr1\t2000\t2050\tGENE2\t500\t+\n")
        
        # Read BED
        df = BEDReader.read_bed(str(bed_file))
        
        assert len(df) == 2
        assert df.iloc[0]['chr'] == 'chr1'
        assert df.iloc[0]['start'] == 1000
        assert df.iloc[0]['end'] == 1050


@pytest.mark.skipif(
    True,  # Skip by default as pysam may not be installed
    reason="pysam not available"
)
class TestBAMReader:
    """Test suite for BAMReader (requires pysam)."""
    
    def test_open_bam(self):
        """Test opening BAM file."""
        from spagapa.io import BAMReader
        
        # This would require a real BAM file
        # Skipped in unit tests
        pass

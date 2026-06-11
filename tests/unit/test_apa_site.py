"""
Unit tests for APASite and APASiteCollection classes.
"""

import pytest
import pandas as pd
from spagapa.core import APASite, APASiteCollection


class TestAPASite:
    """Test suite for APASite class."""
    
    def test_init(self):
        """Test site initialization."""
        site = APASite(
            chr="chr1",
            start=1000,
            end=1050,
            strand="+",
            gene_name="GENE1"
        )
        
        assert site.chr == "chr1"
        assert site.start == 1000
        assert site.end == 1050
        assert site.strand == "+"
        assert site.gene_name == "GENE1"
    
    def test_invalid_strand(self):
        """Test invalid strand raises error."""
        with pytest.raises(ValueError, match="Invalid strand"):
            APASite(chr="chr1", start=1000, end=1050, strand="X")
    
    def test_invalid_coordinates(self):
        """Test invalid coordinates raise error."""
        with pytest.raises(ValueError):
            APASite(chr="chr1", start=-1, end=1050, strand="+")
        
        with pytest.raises(ValueError):
            APASite(chr="chr1", start=1050, end=1000, strand="+")
    
    def test_length(self):
        """Test length calculation."""
        site = APASite(chr="chr1", start=1000, end=1050, strand="+")
        assert site.length == 50
    
    def test_coord(self):
        """Test coordinate calculation."""
        site = APASite(chr="chr1", start=1000, end=1050, strand="+")
        assert site.coord == 1025
    
    def test_overlaps(self):
        """Test overlap detection."""
        site1 = APASite(chr="chr1", start=1000, end=1050, strand="+")
        site2 = APASite(chr="chr1", start=1025, end=1075, strand="+")
        site3 = APASite(chr="chr1", start=1100, end=1150, strand="+")
        site4 = APASite(chr="chr2", start=1025, end=1075, strand="+")
        
        assert site1.overlaps(site2)
        assert not site1.overlaps(site3)
        assert not site1.overlaps(site4)  # Different chromosome
    
    def test_distance_to(self):
        """Test distance calculation."""
        site1 = APASite(chr="chr1", start=1000, end=1050, strand="+")
        site2 = APASite(chr="chr1", start=1025, end=1075, strand="+")
        site3 = APASite(chr="chr1", start=1100, end=1150, strand="+")
        
        assert site1.distance_to(site2) == 0  # Overlapping
        assert site1.distance_to(site3) == 50  # Gap of 50bp
    
    def test_to_bed_line(self):
        """Test BED format conversion."""
        site = APASite(
            chr="chr1",
            start=1000,
            end=1050,
            strand="+",
            gene_name="GENE1",
            support_score=0.8
        )
        
        bed_line = site.to_bed_line()
        fields = bed_line.split('\t')
        
        assert fields[0] == "chr1"
        assert fields[1] == "1000"
        assert fields[2] == "1050"
        assert fields[3] == "GENE1"
        assert fields[5] == "+"


class TestAPASiteCollection:
    """Test suite for APASiteCollection class."""
    
    def test_from_dataframe(self):
        """Test creating collection from DataFrame."""
        df = pd.DataFrame({
            'chr': ['chr1', 'chr1', 'chr2'],
            'start': [1000, 2000, 3000],
            'end': [1050, 2050, 3050],
            'strand': ['+', '+', '-'],
            'gene_name': ['GENE1', 'GENE2', 'GENE3']
        })
        
        collection = APASiteCollection.from_dataframe(df)
        
        assert len(collection) == 3
        assert collection[0].chr == "chr1"
        assert collection[0].gene_name == "GENE1"
    
    def test_to_dataframe(self):
        """Test converting collection to DataFrame."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", gene_name="GENE1"),
            APASite(chr="chr1", start=2000, end=2050, strand="+", gene_name="GENE2"),
        ]
        collection = APASiteCollection(sites)
        
        df = collection.to_dataframe()
        
        assert len(df) == 2
        assert df.iloc[0]['chr'] == "chr1"
        assert df.iloc[0]['gene_name'] == "GENE1"
    
    def test_filter_by_support(self):
        """Test filtering by support score."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", support_score=0.8),
            APASite(chr="chr1", start=2000, end=2050, strand="+", support_score=0.2),
            APASite(chr="chr1", start=3000, end=3050, strand="+", support_score=0.5),
        ]
        collection = APASiteCollection(sites)
        
        filtered = collection.filter_by_support(min_score=0.3)
        
        assert len(filtered) == 2  # Only sites with score >= 0.3
    
    def test_filter_by_reads(self):
        """Test filtering by read count."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", read_count=10),
            APASite(chr="chr1", start=2000, end=2050, strand="+", read_count=3),
            APASite(chr="chr1", start=3000, end=3050, strand="+", read_count=7),
        ]
        collection = APASiteCollection(sites)
        
        filtered = collection.filter_by_reads(min_reads=5)
        
        assert len(filtered) == 2  # Only sites with reads >= 5
    
    def test_filter_by_chromosome(self):
        """Test filtering by chromosome."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+"),
            APASite(chr="chr2", start=2000, end=2050, strand="+"),
            APASite(chr="chr3", start=3000, end=3050, strand="+"),
        ]
        collection = APASiteCollection(sites)
        
        filtered = collection.filter_by_chromosome(['chr1', 'chr2'])
        
        assert len(filtered) == 2
    
    def test_group_by_gene(self):
        """Test grouping sites by gene."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", gene_id="GENE1"),
            APASite(chr="chr1", start=2000, end=2050, strand="+", gene_id="GENE1"),
            APASite(chr="chr1", start=3000, end=3050, strand="+", gene_id="GENE2"),
        ]
        collection = APASiteCollection(sites)
        
        groups = collection.group_by_gene()
        
        assert len(groups) == 2
        assert len(groups['GENE1']) == 2
        assert len(groups['GENE2']) == 1
    
    def test_merge_overlapping(self):
        """Test merging overlapping sites."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", read_count=5),
            APASite(chr="chr1", start=1045, end=1095, strand="+", read_count=3),
            APASite(chr="chr1", start=2000, end=2050, strand="+", read_count=10),
        ]
        collection = APASiteCollection(sites)
        
        merged = collection.merge_overlapping(max_distance=10)
        
        assert len(merged) == 2  # First two should be merged
        assert merged[0].read_count == 8  # 5 + 3
    
    def test_to_bed(self, tmp_path):
        """Test saving to BED file."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+", gene_name="GENE1"),
            APASite(chr="chr1", start=2000, end=2050, strand="+", gene_name="GENE2"),
        ]
        collection = APASiteCollection(sites)
        
        filename = tmp_path / "test.bed"
        collection.to_bed(str(filename))
        
        assert filename.exists()
        
        # Read back
        with open(filename, 'r') as f:
            lines = f.readlines()
        assert len(lines) == 2
    
    def test_from_bed(self, tmp_path):
        """Test loading from BED file."""
        # Create BED file
        filename = tmp_path / "test.bed"
        with open(filename, 'w') as f:
            f.write("chr1\t1000\t1050\tGENE1\t800\t+\n")
            f.write("chr1\t2000\t2050\tGENE2\t500\t+\n")
        
        collection = APASiteCollection.from_bed(str(filename))
        
        assert len(collection) == 2
        assert collection[0].chr == "chr1"
        assert collection[0].gene_name == "GENE1"
    
    def test_iteration(self):
        """Test iterating over collection."""
        sites = [
            APASite(chr="chr1", start=1000, end=1050, strand="+"),
            APASite(chr="chr1", start=2000, end=2050, strand="+"),
        ]
        collection = APASiteCollection(sites)
        
        count = 0
        for site in collection:
            assert isinstance(site, APASite)
            count += 1
        
        assert count == 2

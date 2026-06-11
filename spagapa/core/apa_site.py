"""
APA site representation and annotation.

This module provides classes for representing individual APA sites
and collections of sites with genomic annotations.
"""

import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class APASite:
    """
    Represents a single alternative polyadenylation site.
    
    Attributes
    ----------
    chr : str
        Chromosome name
    start : int
        Start position (0-based)
    end : int
        End position (0-based, exclusive)
    strand : str
        Strand ('+' or '-')
    gene_id : str, optional
        Gene ID
    gene_name : str, optional
        Gene name
    site_type : str, optional
        Type of APA site (e.g., 'proximal', 'distal', '3UTR', 'intron')
    support_score : float, optional
        Spatial support score (0-1)
    read_count : int, optional
        Number of supporting reads
    """
    chr: str
    start: int
    end: int
    strand: str
    gene_id: Optional[str] = None
    gene_name: Optional[str] = None
    site_type: Optional[str] = None
    support_score: Optional[float] = None
    read_count: Optional[int] = None
    
    def __post_init__(self):
        """Validate site attributes."""
        if self.strand not in ['+', '-', '.']:
            raise ValueError(f"Invalid strand: {self.strand}")
        if self.start < 0 or self.end < 0:
            raise ValueError("Coordinates must be non-negative")
        if self.start >= self.end:
            raise ValueError("Start must be less than end")
    
    @property
    def length(self) -> int:
        """Length of the site."""
        return self.end - self.start
    
    @property
    def coord(self) -> int:
        """Representative coordinate (midpoint)."""
        return (self.start + self.end) // 2
    
    def to_bed_line(self) -> str:
        """
        Convert to BED format line.
        
        Returns
        -------
        str
            BED format string
        """
        name = self.gene_name or self.gene_id or "."
        score = int(self.support_score * 1000) if self.support_score else 0
        return f"{self.chr}\t{self.start}\t{self.end}\t{name}\t{score}\t{self.strand}"
    
    def overlaps(self, other: 'APASite') -> bool:
        """
        Check if this site overlaps with another site.
        
        Parameters
        ----------
        other : APASite
            Another APA site
        
        Returns
        -------
        bool
            True if sites overlap
        """
        if self.chr != other.chr or self.strand != other.strand:
            return False
        return not (self.end <= other.start or self.start >= other.end)
    
    def distance_to(self, other: 'APASite') -> int:
        """
        Calculate distance to another site.
        
        Parameters
        ----------
        other : APASite
            Another APA site
        
        Returns
        -------
        int
            Distance in base pairs (0 if overlapping)
        """
        if self.chr != other.chr or self.strand != other.strand:
            return float('inf')
        
        if self.overlaps(other):
            return 0
        
        if self.end <= other.start:
            return other.start - self.end
        else:
            return self.start - other.end
    
    def __repr__(self) -> str:
        """String representation."""
        gene = self.gene_name or self.gene_id or "unknown"
        return f"APASite({self.chr}:{self.start}-{self.end}:{self.strand}, {gene})"


class APASiteCollection:
    """
    Collection of APA sites with annotation and filtering capabilities.
    
    Parameters
    ----------
    sites : list of APASite or pd.DataFrame
        List of APASite objects or DataFrame with site information
    
    Examples
    --------
    >>> sites = APASiteCollection.from_dataframe(peaks_meta)
    >>> filtered = sites.filter_by_support(min_score=0.3)
    >>> sites.to_bed("apa_sites.bed")
    """
    
    def __init__(self, sites: List[APASite]):
        """Initialize collection."""
        self.sites = sites
    
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> 'APASiteCollection':
        """
        Create collection from DataFrame.
        
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with columns: chr, start, end, strand, etc.
        
        Returns
        -------
        APASiteCollection
            New collection
        """
        sites = []
        for idx, row in df.iterrows():
            site = APASite(
                chr=row.get('chr', row.get('seqnames', 'chr1')),
                start=int(row.get('start', 0)),
                end=int(row.get('end', 0)),
                strand=row.get('strand', '+'),
                gene_id=row.get('gene_id'),
                gene_name=row.get('gene_name'),
                site_type=row.get('site_type', row.get('ftr')),
                support_score=row.get('support_score'),
                read_count=row.get('read_count', row.get('count'))
            )
            sites.append(site)
        
        logger.info(f"Created collection with {len(sites)} sites")
        return cls(sites)
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert collection to DataFrame.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with site information
        """
        data = []
        for site in self.sites:
            data.append({
                'chr': site.chr,
                'start': site.start,
                'end': site.end,
                'strand': site.strand,
                'gene_id': site.gene_id,
                'gene_name': site.gene_name,
                'site_type': site.site_type,
                'support_score': site.support_score,
                'read_count': site.read_count,
                'length': site.length,
                'coord': site.coord
            })
        return pd.DataFrame(data)
    
    def filter_by_support(self, min_score: float = 0.3) -> 'APASiteCollection':
        """
        Filter sites by spatial support score.
        
        Parameters
        ----------
        min_score : float, default=0.3
            Minimum support score
        
        Returns
        -------
        APASiteCollection
            Filtered collection
        """
        filtered = [
            site for site in self.sites
            if site.support_score is not None and site.support_score >= min_score
        ]
        logger.info(f"Filtered by support: {len(filtered)}/{len(self.sites)} sites")
        return APASiteCollection(filtered)
    
    def filter_by_reads(self, min_reads: int = 5) -> 'APASiteCollection':
        """
        Filter sites by read count.
        
        Parameters
        ----------
        min_reads : int, default=5
            Minimum number of reads
        
        Returns
        -------
        APASiteCollection
            Filtered collection
        """
        filtered = [
            site for site in self.sites
            if site.read_count is not None and site.read_count >= min_reads
        ]
        logger.info(f"Filtered by reads: {len(filtered)}/{len(self.sites)} sites")
        return APASiteCollection(filtered)
    
    def filter_by_type(self, site_types: List[str]) -> 'APASiteCollection':
        """
        Filter sites by type.
        
        Parameters
        ----------
        site_types : list of str
            Allowed site types (e.g., ['3UTR', 'exon'])
        
        Returns
        -------
        APASiteCollection
            Filtered collection
        """
        filtered = [
            site for site in self.sites
            if site.site_type in site_types
        ]
        logger.info(f"Filtered by type: {len(filtered)}/{len(self.sites)} sites")
        return APASiteCollection(filtered)
    
    def filter_by_chromosome(self, chromosomes: List[str]) -> 'APASiteCollection':
        """
        Filter sites by chromosome.
        
        Parameters
        ----------
        chromosomes : list of str
            Allowed chromosomes
        
        Returns
        -------
        APASiteCollection
            Filtered collection
        """
        filtered = [site for site in self.sites if site.chr in chromosomes]
        logger.info(f"Filtered by chr: {len(filtered)}/{len(self.sites)} sites")
        return APASiteCollection(filtered)
    
    def group_by_gene(self) -> Dict[str, List[APASite]]:
        """
        Group sites by gene.
        
        Returns
        -------
        dict
            Dictionary mapping gene_id to list of sites
        """
        groups = {}
        for site in self.sites:
            gene_id = site.gene_id or site.gene_name or "unknown"
            if gene_id not in groups:
                groups[gene_id] = []
            groups[gene_id].append(site)
        return groups
    
    def merge_overlapping(self, max_distance: int = 10) -> 'APASiteCollection':
        """
        Merge overlapping or nearby sites.
        
        Parameters
        ----------
        max_distance : int, default=10
            Maximum distance to merge sites
        
        Returns
        -------
        APASiteCollection
            Collection with merged sites
        """
        # Sort sites by chromosome, strand, and position
        sorted_sites = sorted(
            self.sites,
            key=lambda s: (s.chr, s.strand, s.start)
        )
        
        merged = []
        current = None
        
        for site in sorted_sites:
            if current is None:
                current = site
            elif (current.chr == site.chr and 
                  current.strand == site.strand and
                  site.start - current.end <= max_distance):
                # Merge sites
                current = APASite(
                    chr=current.chr,
                    start=min(current.start, site.start),
                    end=max(current.end, site.end),
                    strand=current.strand,
                    gene_id=current.gene_id or site.gene_id,
                    gene_name=current.gene_name or site.gene_name,
                    site_type=current.site_type or site.site_type,
                    support_score=max(
                        current.support_score or 0,
                        site.support_score or 0
                    ),
                    read_count=(current.read_count or 0) + (site.read_count or 0)
                )
            else:
                merged.append(current)
                current = site
        
        if current is not None:
            merged.append(current)
        
        logger.info(f"Merged sites: {len(self.sites)} -> {len(merged)}")
        return APASiteCollection(merged)
    
    def to_bed(self, filename: str):
        """
        Save sites to BED file.
        
        Parameters
        ----------
        filename : str
            Output BED filename
        """
        with open(filename, 'w') as f:
            for site in self.sites:
                f.write(site.to_bed_line() + '\n')
        logger.info(f"Saved {len(self.sites)} sites to {filename}")
    
    @classmethod
    def from_bed(cls, filename: str) -> 'APASiteCollection':
        """
        Load sites from BED file.
        
        Parameters
        ----------
        filename : str
            Input BED filename
        
        Returns
        -------
        APASiteCollection
            Loaded collection
        """
        sites = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                site = APASite(
                    chr=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    strand=fields[5] if len(fields) > 5 else '+',
                    gene_name=fields[3] if len(fields) > 3 else None,
                    support_score=float(fields[4])/1000 if len(fields) > 4 else None
                )
                sites.append(site)
        
        logger.info(f"Loaded {len(sites)} sites from {filename}")
        return cls(sites)
    
    def __len__(self) -> int:
        """Number of sites."""
        return len(self.sites)
    
    def __getitem__(self, idx: int) -> APASite:
        """Get site by index."""
        return self.sites[idx]
    
    def __iter__(self):
        """Iterate over sites."""
        return iter(self.sites)
    
    def __repr__(self) -> str:
        """String representation."""
        n_genes = len(set(s.gene_id for s in self.sites if s.gene_id))
        return f"APASiteCollection({len(self.sites)} sites, {n_genes} genes)"

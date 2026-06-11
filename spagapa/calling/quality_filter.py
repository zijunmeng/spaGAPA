"""
Quality filtering for APA sites.

This module implements quality control filters for APA sites based on
read counts, spatial support, and other quality metrics.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple
from spagapa.core import APASite, APASiteCollection
import logging

logger = logging.getLogger(__name__)


class QualityFilter:
    """
    Filter APA sites based on quality metrics.
    
    Parameters
    ----------
    min_read_count : int, default=10
        Minimum total read count across all spots
    min_spots : int, default=3
        Minimum number of spots with non-zero counts
    min_spatial_support : float, default=0.3
        Minimum spatial support score (if available)
    min_mean_count : float, default=1.0
        Minimum mean count per spot (for non-zero spots)
    max_cv : float, optional
        Maximum coefficient of variation (filters noisy sites)
    
    Examples
    --------
    >>> qf = QualityFilter(min_read_count=10, min_spots=3)
    >>> filtered_sites = qf.filter_sites(sites, apa_counts)
    """
    
    def __init__(
        self,
        min_read_count: int = 10,
        min_spots: int = 3,
        min_spatial_support: Optional[float] = 0.3,
        min_mean_count: float = 1.0,
        max_cv: Optional[float] = None
    ):
        self.min_read_count = min_read_count
        self.min_spots = min_spots
        self.min_spatial_support = min_spatial_support
        self.min_mean_count = min_mean_count
        self.max_cv = max_cv
    
    def compute_site_metrics(
        self,
        apa_counts: np.ndarray,
        site_idx: int
    ) -> Dict[str, float]:
        """
        Compute quality metrics for a single APA site.
        
        Parameters
        ----------
        apa_counts : np.ndarray, shape (n_genes, n_spots)
            APA count matrix
        site_idx : int
            Index of the site
        
        Returns
        -------
        dict
            Dictionary of quality metrics
        """
        counts = apa_counts[site_idx, :]
        nonzero_counts = counts[counts > 0]
        
        metrics = {
            'total_count': np.sum(counts),
            'n_spots': len(nonzero_counts),
            'mean_count': np.mean(nonzero_counts) if len(nonzero_counts) > 0 else 0.0,
            'median_count': np.median(nonzero_counts) if len(nonzero_counts) > 0 else 0.0,
            'max_count': np.max(counts),
            'cv': np.std(nonzero_counts) / np.mean(nonzero_counts) if len(nonzero_counts) > 0 and np.mean(nonzero_counts) > 0 else np.inf
        }
        
        return metrics
    
    def filter_by_read_count(
        self,
        apa_counts: np.ndarray
    ) -> np.ndarray:
        """
        Filter sites by minimum read count.
        
        Parameters
        ----------
        apa_counts : np.ndarray, shape (n_genes, n_spots)
            APA count matrix
        
        Returns
        -------
        np.ndarray
            Boolean mask of sites passing filter
        """
        total_counts = np.sum(apa_counts, axis=1)
        mask = total_counts >= self.min_read_count
        
        logger.info(
            f"Read count filter: {np.sum(mask)}/{len(mask)} sites pass "
            f"(threshold={self.min_read_count})"
        )
        
        return mask
    
    def filter_by_spot_count(
        self,
        apa_counts: np.ndarray
    ) -> np.ndarray:
        """
        Filter sites by minimum number of spots.
        
        Parameters
        ----------
        apa_counts : np.ndarray
            APA count matrix
        
        Returns
        -------
        np.ndarray
            Boolean mask of sites passing filter
        """
        n_spots = np.sum(apa_counts > 0, axis=1)
        mask = n_spots >= self.min_spots
        
        logger.info(
            f"Spot count filter: {np.sum(mask)}/{len(mask)} sites pass "
            f"(threshold={self.min_spots})"
        )
        
        return mask
    
    def filter_by_mean_count(
        self,
        apa_counts: np.ndarray
    ) -> np.ndarray:
        """
        Filter sites by minimum mean count.
        
        Parameters
        ----------
        apa_counts : np.ndarray
            APA count matrix
        
        Returns
        -------
        np.ndarray
            Boolean mask of sites passing filter
        """
        n_genes = apa_counts.shape[0]
        mask = np.zeros(n_genes, dtype=bool)
        
        for i in range(n_genes):
            nonzero_counts = apa_counts[i, apa_counts[i, :] > 0]
            if len(nonzero_counts) > 0:
                mean_count = np.mean(nonzero_counts)
                mask[i] = mean_count >= self.min_mean_count
        
        logger.info(
            f"Mean count filter: {np.sum(mask)}/{len(mask)} sites pass "
            f"(threshold={self.min_mean_count})"
        )
        
        return mask
    
    def filter_by_cv(
        self,
        apa_counts: np.ndarray
    ) -> np.ndarray:
        """
        Filter sites by maximum coefficient of variation.
        
        High CV indicates noisy/inconsistent expression.
        
        Parameters
        ----------
        apa_counts : np.ndarray
            APA count matrix
        
        Returns
        -------
        np.ndarray
            Boolean mask of sites passing filter
        """
        if self.max_cv is None:
            return np.ones(apa_counts.shape[0], dtype=bool)
        
        n_genes = apa_counts.shape[0]
        mask = np.zeros(n_genes, dtype=bool)
        
        for i in range(n_genes):
            nonzero_counts = apa_counts[i, apa_counts[i, :] > 0]
            if len(nonzero_counts) > 0:
                mean_count = np.mean(nonzero_counts)
                if mean_count > 0:
                    cv = np.std(nonzero_counts) / mean_count
                    mask[i] = cv <= self.max_cv
                else:
                    mask[i] = False
            else:
                mask[i] = False
        
        logger.info(
            f"CV filter: {np.sum(mask)}/{len(mask)} sites pass "
            f"(threshold={self.max_cv})"
        )
        
        return mask
    
    def filter_sites(
        self,
        sites: APASiteCollection,
        apa_counts: np.ndarray,
        spatial_support: Optional[np.ndarray] = None
    ) -> Tuple[APASiteCollection, np.ndarray]:
        """
        Apply all quality filters to APA sites.
        
        Parameters
        ----------
        sites : APASiteCollection
            Collection of APA sites
        apa_counts : np.ndarray, shape (n_genes, n_spots)
            APA count matrix
        spatial_support : np.ndarray, optional
            Spatial support scores for each site
        
        Returns
        -------
        filtered_sites : APASiteCollection
            Sites passing all filters
        pass_mask : np.ndarray
            Boolean mask indicating which sites passed
        """
        n_sites = len(sites)
        
        # Initialize mask (all pass)
        pass_mask = np.ones(n_sites, dtype=bool)
        
        # Apply filters sequentially
        pass_mask &= self.filter_by_read_count(apa_counts)
        pass_mask &= self.filter_by_spot_count(apa_counts)
        pass_mask &= self.filter_by_mean_count(apa_counts)
        
        if self.max_cv is not None:
            pass_mask &= self.filter_by_cv(apa_counts)
        
        # Filter by spatial support if provided
        if spatial_support is not None and self.min_spatial_support is not None:
            support_mask = spatial_support >= self.min_spatial_support
            pass_mask &= support_mask
            logger.info(
                f"Spatial support filter: {np.sum(support_mask)}/{n_sites} sites pass "
                f"(threshold={self.min_spatial_support})"
            )
        
        # Create filtered collection
        filtered_sites = APASiteCollection([
            site for site, passed in zip(sites, pass_mask) if passed
        ])
        
        logger.info(
            f"Quality filtering: {len(filtered_sites)}/{n_sites} sites pass "
            f"({100*len(filtered_sites)/n_sites:.1f}%)" if n_sites > 0 else "Quality filtering: 0/0 sites"
        )
        
        return filtered_sites, pass_mask
    
    def generate_qc_report(
        self,
        sites: APASiteCollection,
        apa_counts: np.ndarray,
        spatial_support: Optional[np.ndarray] = None
    ) -> pd.DataFrame:
        """
        Generate quality control report for all sites.
        
        Parameters
        ----------
        sites : APASiteCollection
            Collection of APA sites
        apa_counts : np.ndarray
            APA count matrix
        spatial_support : np.ndarray, optional
            Spatial support scores
        
        Returns
        -------
        pd.DataFrame
            QC report with metrics for each site
        """
        n_sites = len(sites)
        
        # Compute metrics for all sites
        metrics_list = []
        for i in range(n_sites):
            metrics = self.compute_site_metrics(apa_counts, i)
            site = sites[i]
            metrics['chr'] = site.chr
            metrics['start'] = site.start
            metrics['end'] = site.end
            metrics['strand'] = site.strand
            metrics['gene_id'] = site.gene_id
            
            if spatial_support is not None:
                metrics['spatial_support'] = spatial_support[i]
            
            # Check which filters pass
            metrics['pass_read_count'] = metrics['total_count'] >= self.min_read_count
            metrics['pass_spot_count'] = metrics['n_spots'] >= self.min_spots
            metrics['pass_mean_count'] = metrics['mean_count'] >= self.min_mean_count
            
            if self.max_cv is not None:
                metrics['pass_cv'] = metrics['cv'] <= self.max_cv
            
            if spatial_support is not None and self.min_spatial_support is not None:
                metrics['pass_spatial_support'] = spatial_support[i] >= self.min_spatial_support
            
            metrics_list.append(metrics)
        
        df = pd.DataFrame(metrics_list)
        
        # Add overall pass/fail
        pass_cols = [col for col in df.columns if col.startswith('pass_')]
        df['pass_all'] = df[pass_cols].all(axis=1)
        
        logger.info(f"Generated QC report for {n_sites} sites")
        
        return df


def filter_apa_sites(
    sites: APASiteCollection,
    apa_counts: np.ndarray,
    min_read_count: int = 10,
    min_spots: int = 3,
    min_spatial_support: Optional[float] = 0.3,
    spatial_support: Optional[np.ndarray] = None
) -> Tuple[APASiteCollection, np.ndarray]:
    """
    Convenience function for quality filtering.
    
    Parameters
    ----------
    sites : APASiteCollection
        APA sites to filter
    apa_counts : np.ndarray
        APA count matrix
    min_read_count : int, default=10
        Minimum total read count
    min_spots : int, default=3
        Minimum number of spots
    min_spatial_support : float, default=0.3
        Minimum spatial support
    spatial_support : np.ndarray, optional
        Spatial support scores
    
    Returns
    -------
    filtered_sites : APASiteCollection
        Filtered sites
    pass_mask : np.ndarray
        Boolean mask
    """
    qf = QualityFilter(
        min_read_count=min_read_count,
        min_spots=min_spots,
        min_spatial_support=min_spatial_support
    )
    
    return qf.filter_sites(sites, apa_counts, spatial_support)

"""
File writers for various output formats.

This module provides writers for:
- CSV/TSV files
- BED files (APA sites)
- H5AD files (AnnData)
- Summary reports
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)

# Optional imports
try:
    import anndata
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False


class ResultWriter:
    """
    Writer for analysis results.
    
    Supports multiple output formats:
    - CSV/TSV for tables
    - BED for genomic regions
    - H5AD for AnnData objects
    """
    
    @staticmethod
    def write_apa_sites(
        sites_df: pd.DataFrame,
        filename: str,
        format: str = 'bed'
    ):
        """
        Write APA sites to file.
        
        Parameters
        ----------
        sites_df : pd.DataFrame
            DataFrame with APA site information
        filename : str
            Output filename
        format : str, default='bed'
            Output format ('bed', 'csv', 'tsv')
        """
        if format == 'bed':
            # Convert to BED format
            bed_df = sites_df[['chr', 'start', 'end', 'gene_name', 'score', 'strand']]
            bed_df.to_csv(filename, sep='\t', header=False, index=False)
        elif format == 'csv':
            sites_df.to_csv(filename, index=False)
        elif format == 'tsv':
            sites_df.to_csv(filename, sep='\t', index=False)
        else:
            raise ValueError(f"Unknown format: {format}")
        
        logger.info(f"Wrote {len(sites_df)} APA sites to {filename}")
    
    @staticmethod
    def write_matrix(
        matrix: np.ndarray,
        filename: str,
        row_names: Optional[list] = None,
        col_names: Optional[list] = None,
        format: str = 'csv'
    ):
        """
        Write matrix to file.
        
        Parameters
        ----------
        matrix : np.ndarray
            Matrix to write
        filename : str
            Output filename
        row_names : list, optional
            Row names
        col_names : list, optional
            Column names
        format : str, default='csv'
            Output format ('csv', 'tsv')
        """
        df = pd.DataFrame(matrix, index=row_names, columns=col_names)
        
        if format == 'csv':
            df.to_csv(filename)
        elif format == 'tsv':
            df.to_csv(filename, sep='\t')
        else:
            raise ValueError(f"Unknown format: {format}")
        
        logger.info(f"Wrote matrix {matrix.shape} to {filename}")
    
    @staticmethod
    def write_differential_results(
        results_df: pd.DataFrame,
        filename: str,
        format: str = 'csv'
    ):
        """
        Write differential APA results.
        
        Parameters
        ----------
        results_df : pd.DataFrame
            Differential APA results
        filename : str
            Output filename
        format : str, default='csv'
            Output format
        """
        if format == 'csv':
            results_df.to_csv(filename)
        elif format == 'tsv':
            results_df.to_csv(filename, sep='\t')
        else:
            raise ValueError(f"Unknown format: {format}")
        
        logger.info(f"Wrote differential results to {filename}")
    
    @staticmethod
    def write_h5ad(
        adata: 'anndata.AnnData',
        filename: str
    ):
        """
        Write AnnData to H5AD file.
        
        Parameters
        ----------
        adata : anndata.AnnData
            AnnData object
        filename : str
            Output filename
        """
        if not ANNDATA_AVAILABLE:
            raise ImportError("anndata required for H5AD writing")
        
        adata.write_h5ad(filename)
        logger.info(f"Wrote AnnData to {filename}")
    
    @staticmethod
    def write_summary_report(
        results: Dict[str, Any],
        filename: str
    ):
        """
        Write summary report.
        
        Parameters
        ----------
        results : dict
            Dictionary with analysis results
        filename : str
            Output filename (markdown format)
        """
        with open(filename, 'w') as f:
            f.write("# spaGAPA Analysis Report\n\n")
            
            # Dataset info
            if 'dataset_info' in results:
                f.write("## Dataset Information\n\n")
                for key, value in results['dataset_info'].items():
                    f.write(f"- {key}: {value}\n")
                f.write("\n")
            
            # APA sites
            if 'n_apa_sites' in results:
                f.write("## APA Sites\n\n")
                f.write(f"- Total sites identified: {results['n_apa_sites']}\n")
                if 'sites_by_type' in results:
                    f.write("- Sites by type:\n")
                    for site_type, count in results['sites_by_type'].items():
                        f.write(f"  - {site_type}: {count}\n")
                f.write("\n")
            
            # Imputation
            if 'imputation_metrics' in results:
                f.write("## Imputation Quality\n\n")
                for metric, value in results['imputation_metrics'].items():
                    f.write(f"- {metric}: {value:.4f}\n")
                f.write("\n")
            
            # Differential APA
            if 'n_differential_genes' in results:
                f.write("## Differential APA\n\n")
                f.write(f"- Significant genes: {results['n_differential_genes']}\n")
                f.write("\n")
        
        logger.info(f"Wrote summary report to {filename}")


class BEDWriter:
    """Writer for BED format files."""
    
    @staticmethod
    def write_bed(
        df: pd.DataFrame,
        filename: str,
        columns: Optional[list] = None
    ):
        """
        Write DataFrame to BED file.
        
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with BED columns
        filename : str
            Output filename
        columns : list, optional
            Column names to write (default: chr, start, end, name, score, strand)
        """
        if columns is None:
            columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
        
        # Select columns
        bed_df = df[columns]
        
        # Write without header
        bed_df.to_csv(filename, sep='\t', header=False, index=False)
        
        logger.info(f"Wrote {len(bed_df)} entries to BED file: {filename}")

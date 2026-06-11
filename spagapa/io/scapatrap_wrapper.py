"""
Python wrapper for scAPAtrap R package.

This module provides a Python interface to call scAPAtrap for APA site identification
from BAM files. scAPAtrap is an R package designed for single-cell RNA-seq data,
but we adapt it for spatial transcriptomics data.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class ScAPAtrapWrapper:
    """
    Wrapper class for calling scAPAtrap from Python.
    
    scAPAtrap identifies poly(A) sites through two modes:
    1. Peak calling based on change points
    2. Tail finding based on soft-clipping
    
    Parameters
    ----------
    scapatrap_path : str, optional
        Path to scAPAtrap R package. If None, assumes it's installed in R.
    r_executable : str, optional
        Path to R executable. Default is 'Rscript'.
    """
    
    def __init__(
        self,
        scapatrap_path: Optional[str] = None,
        r_executable: str = "Rscript"
    ):
        self.scapatrap_path = scapatrap_path
        self.r_executable = r_executable
        self._check_r_installation()
        self._check_scapatrap_installation()
    
    def _check_r_installation(self):
        """Check if R is installed and accessible."""
        try:
            result = subprocess.run(
                [self.r_executable, "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info(f"R version: {result.stdout.split()[2]}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                f"R executable not found at {self.r_executable}. "
                "Please install R or provide correct path."
            )
    
    def _check_scapatrap_installation(self):
        """Check if scAPAtrap is installed in R."""
        r_code = """
        if (!require("scAPAtrap", quietly = TRUE)) {
            cat("NOT_INSTALLED")
        } else {
            cat("INSTALLED")
        }
        """
        
        try:
            result = subprocess.run(
                [self.r_executable, "-e", r_code],
                capture_output=True,
                text=True,
                check=True
            )
            
            if "NOT_INSTALLED" in result.stdout:
                logger.warning(
                    "scAPAtrap not installed in R. "
                    "Install it using: devtools::install_github('BMILAB/scAPAtrap')"
                )
            else:
                logger.info("scAPAtrap is installed")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error checking scAPAtrap: {e}")
    
    def run_scapatrap(
        self,
        bam_file: str,
        output_dir: str,
        genome_fasta: str,
        gtf_file: str,
        barcode_file: Optional[str] = None,
        tails_search: str = "genome",
        n_cores: int = 4,
        **kwargs
    ) -> Dict[str, pd.DataFrame]:
        """
        Run scAPAtrap on a BAM file.
        
        Parameters
        ----------
        bam_file : str
            Path to input BAM file (sorted and indexed)
        output_dir : str
            Output directory for scAPAtrap results
        genome_fasta : str
            Path to genome FASTA file
        gtf_file : str
            Path to GTF annotation file
        barcode_file : str, optional
            Path to barcode file (for spatial data, this would be spot barcodes)
        tails_search : str, default='genome'
            Mode for tail search: 'genome', 'no', or 'both'
        n_cores : int, default=4
            Number of cores to use
        **kwargs
            Additional parameters for scAPAtrap
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'peaks': DataFrame of identified peaks
            - 'sites': DataFrame of poly(A) sites
            - 'counts': DataFrame of count matrix
        """
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate R script
        r_script = self._generate_r_script(
            bam_file=bam_file,
            output_dir=output_dir,
            genome_fasta=genome_fasta,
            gtf_file=gtf_file,
            barcode_file=barcode_file,
            tails_search=tails_search,
            n_cores=n_cores,
            **kwargs
        )
        
        # Write R script to temporary file
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.R',
            delete=False
        ) as f:
            f.write(r_script)
            script_path = f.name
        
        try:
            # Run R script
            logger.info(f"Running scAPAtrap on {bam_file}")
            result = subprocess.run(
                [self.r_executable, script_path],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("scAPAtrap completed successfully")
            logger.debug(f"R output: {result.stdout}")
            
            # Parse results
            results = self._parse_scapatrap_output(output_dir)
            return results
            
        except subprocess.CalledProcessError as e:
            logger.error(f"scAPAtrap failed: {e.stderr}")
            raise RuntimeError(f"scAPAtrap execution failed: {e.stderr}")
        finally:
            # Clean up temporary script
            if os.path.exists(script_path):
                os.remove(script_path)
    
    def _generate_r_script(
        self,
        bam_file: str,
        output_dir: str,
        genome_fasta: str,
        gtf_file: str,
        barcode_file: Optional[str],
        tails_search: str,
        n_cores: int,
        **kwargs
    ) -> str:
        """Generate R script for running scAPAtrap."""
        
        # Base R script template
        r_script = f"""
library(scAPAtrap)

# Set parameters
bam_file <- "{bam_file}"
output_dir <- "{output_dir}"
genome_fasta <- "{genome_fasta}"
gtf_file <- "{gtf_file}"
tails_search <- "{tails_search}"
n_cores <- {n_cores}

# Initialize scAPAtrap
trap_params <- setTrapParams(
    genome_fa = genome_fasta,
    gtf_file = gtf_file,
    tails.search = tails_search,
    n.cores = n_cores
)

# Run scAPAtrap
scAPAtrap(
    bam_file = bam_file,
    output_dir = output_dir,
    trap_params = trap_params
)

cat("scAPAtrap completed successfully\\n")
"""
        return r_script
    
    def _parse_scapatrap_output(
        self,
        output_dir: str
    ) -> Dict[str, pd.DataFrame]:
        """
        Parse scAPAtrap output files.
        
        Parameters
        ----------
        output_dir : str
            Directory containing scAPAtrap output
        
        Returns
        -------
        dict
            Dictionary with parsed results
        """
        results = {}
        
        # Look for scAPAtrapData.rda file
        rda_file = os.path.join(output_dir, "scAPAtrapData.rda")
        
        if not os.path.exists(rda_file):
            logger.warning(f"scAPAtrapData.rda not found in {output_dir}")
            return results
        
        # Use R to convert .rda to CSV for easier parsing in Python
        r_code = f"""
        load("{rda_file}")
        
        # Export peaks metadata
        write.csv(
            scAPAtrapData$peaks.meta,
            file.path("{output_dir}", "peaks_meta.csv"),
            row.names = TRUE
        )
        
        # Export peaks counts
        write.csv(
            as.matrix(scAPAtrapData$peaks.count),
            file.path("{output_dir}", "peaks_counts.csv"),
            row.names = TRUE
        )
        
        cat("Export completed\\n")
        """
        
        try:
            subprocess.run(
                [self.r_executable, "-e", r_code],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Read exported CSV files
            peaks_meta_file = os.path.join(output_dir, "peaks_meta.csv")
            peaks_counts_file = os.path.join(output_dir, "peaks_counts.csv")
            
            if os.path.exists(peaks_meta_file):
                results['peaks_meta'] = pd.read_csv(peaks_meta_file, index_col=0)
                logger.info(f"Loaded {len(results['peaks_meta'])} peaks")
            
            if os.path.exists(peaks_counts_file):
                results['peaks_counts'] = pd.read_csv(peaks_counts_file, index_col=0)
                logger.info(f"Loaded count matrix: {results['peaks_counts'].shape}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to parse scAPAtrap output: {e}")
        
        return results
    
    def filter_peaks_by_spatial_support(
        self,
        peaks_meta: pd.DataFrame,
        peaks_counts: pd.DataFrame,
        spatial_coords: pd.DataFrame,
        min_spots: int = 3
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Filter peaks that appear in multiple spatial spots.
        
        This is a preliminary spatial filter before our spatial validation.
        
        Parameters
        ----------
        peaks_meta : pd.DataFrame
            Peak metadata from scAPAtrap
        peaks_counts : pd.DataFrame
            Peak count matrix
        spatial_coords : pd.DataFrame
            Spatial coordinates of spots
        min_spots : int, default=3
            Minimum number of spots where peak must be detected
        
        Returns
        -------
        tuple
            Filtered (peaks_meta, peaks_counts)
        """
        # Count number of spots with non-zero counts for each peak
        n_spots_per_peak = (peaks_counts > 0).sum(axis=1)
        
        # Filter peaks
        valid_peaks = n_spots_per_peak >= min_spots
        
        filtered_meta = peaks_meta[valid_peaks]
        filtered_counts = peaks_counts[valid_peaks]
        
        logger.info(
            f"Filtered peaks: {valid_peaks.sum()}/{len(peaks_meta)} "
            f"({100*valid_peaks.sum()/len(peaks_meta):.1f}%) passed"
        )
        
        return filtered_meta, filtered_counts


def run_scapatrap_pipeline(
    bam_file: str,
    output_dir: str,
    genome_fasta: str,
    gtf_file: str,
    spatial_coords: Optional[pd.DataFrame] = None,
    **kwargs
) -> Dict[str, pd.DataFrame]:
    """
    Convenience function to run scAPAtrap pipeline.
    
    Parameters
    ----------
    bam_file : str
        Path to input BAM file
    output_dir : str
        Output directory
    genome_fasta : str
        Path to genome FASTA
    gtf_file : str
        Path to GTF annotation
    spatial_coords : pd.DataFrame, optional
        Spatial coordinates for filtering
    **kwargs
        Additional parameters for scAPAtrap
    
    Returns
    -------
    dict
        Dictionary with scAPAtrap results
    """
    wrapper = ScAPAtrapWrapper()
    results = wrapper.run_scapatrap(
        bam_file=bam_file,
        output_dir=output_dir,
        genome_fasta=genome_fasta,
        gtf_file=gtf_file,
        **kwargs
    )
    
    # Apply spatial filtering if coordinates provided
    if spatial_coords is not None and 'peaks_meta' in results:
        results['peaks_meta'], results['peaks_counts'] = \
            wrapper.filter_peaks_by_spatial_support(
                results['peaks_meta'],
                results['peaks_counts'],
                spatial_coords
            )
    
    return results

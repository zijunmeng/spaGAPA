"""
File readers for various input formats.

This module provides readers for:
- BAM files (aligned reads)
- Spatial coordinates (CSV, TSV)
- AnnData objects (H5AD)
- BED files (genomic regions)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Union, Literal
import logging

from spagapa.core import APADataset

logger = logging.getLogger(__name__)

# Optional imports
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    logger.warning("pysam not installed. BAM reading will not be available.")

try:
    import anndata
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    logger.warning("anndata not installed. H5AD reading will not be available.")


class BAMReader:
    """
    Reader for BAM files.
    
    Parameters
    ----------
    bam_file : str
        Path to BAM file (must be sorted and indexed)
    """
    
    def __init__(self, bam_file: str):
        """Initialize BAM reader."""
        if not PYSAM_AVAILABLE:
            raise ImportError(
                "pysam is required for BAM reading. "
                "Install with: pip install pysam"
            )
        
        self.bam_file = bam_file
        self.bam = None
        self._open()
    
    def _open(self):
        """Open BAM file."""
        try:
            self.bam = pysam.AlignmentFile(self.bam_file, "rb")
            logger.info(f"Opened BAM file: {self.bam_file}")
        except Exception as e:
            raise IOError(f"Failed to open BAM file: {e}")
    
    def get_barcodes(self) -> List[str]:
        """
        Extract unique cell/spot barcodes from BAM file.
        
        Returns
        -------
        list of str
            List of unique barcodes
        """
        barcodes = set()
        for read in self.bam:
            if read.has_tag('CB'):  # Cell barcode tag (10x format)
                barcodes.add(read.get_tag('CB'))
        
        logger.info(f"Found {len(barcodes)} unique barcodes")
        return sorted(list(barcodes))
    
    def get_coverage(
        self,
        region: str,
        barcode: Optional[str] = None
    ) -> np.ndarray:
        """
        Get read coverage for a genomic region.
        
        Parameters
        ----------
        region : str
            Genomic region (format: "chr:start-end")
        barcode : str, optional
            Filter reads by barcode
        
        Returns
        -------
        np.ndarray
            Coverage array
        """
        # Parse region
        chr_name, coords = region.split(':')
        start, end = map(int, coords.split('-'))
        
        # Initialize coverage array
        coverage = np.zeros(end - start)
        
        # Count reads
        for pileupcolumn in self.bam.pileup(chr_name, start, end):
            if start <= pileupcolumn.pos < end:
                if barcode is None:
                    coverage[pileupcolumn.pos - start] = pileupcolumn.n
                else:
                    # Filter by barcode
                    count = sum(
                        1 for read in pileupcolumn.pileups
                        if read.alignment.has_tag('CB') and
                        read.alignment.get_tag('CB') == barcode
                    )
                    coverage[pileupcolumn.pos - start] = count
        
        return coverage
    
    def close(self):
        """Close BAM file."""
        if self.bam is not None:
            self.bam.close()
            logger.info("Closed BAM file")
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


class SpatialCoordinateReader:
    """
    Reader for spatial coordinate files.
    
    Supports CSV and TSV formats with columns: barcode, x, y
    """
    
    @staticmethod
    def read_coordinates(
        coord_file: str,
        format: str = 'auto'
    ) -> pd.DataFrame:
        """
        Read spatial coordinates from file.
        
        Parameters
        ----------
        coord_file : str
            Path to coordinate file
        format : str, default='auto'
            File format ('csv', 'tsv', or 'auto' to detect)
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns: barcode (index), x, y
        """
        # Detect format
        if format == 'auto':
            if coord_file.endswith('.csv'):
                format = 'csv'
            elif coord_file.endswith('.tsv') or coord_file.endswith('.txt'):
                format = 'tsv'
            else:
                raise ValueError(f"Cannot detect format for {coord_file}")
        
        # Read file
        if format == 'csv':
            df = pd.read_csv(coord_file)
        elif format == 'tsv':
            df = pd.read_csv(coord_file, sep='\t')
        else:
            raise ValueError(f"Unknown format: {format}")
        
        # Validate columns
        required_cols = ['x', 'y']
        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"Coordinate file must have columns: {required_cols}")
        
        # Set spot/barcode ID as index if present
        if 'barcode' in df.columns:
            df = df.set_index('barcode')
        elif 'spot_id' in df.columns:
            df = df.set_index('spot_id')
        
        logger.info(f"Read {len(df)} spatial coordinates from {coord_file}")
        
        return df
    
    @staticmethod
    def read_10x_positions(positions_file: str) -> pd.DataFrame:
        """
        Read 10x Visium tissue_positions.csv file.
        
        Parameters
        ----------
        positions_file : str
            Path to tissue_positions.csv
        
        Returns
        -------
        pd.DataFrame
            DataFrame with spatial coordinates
        """
        # Read 10x format
        df = pd.read_csv(positions_file, header=None)
        
        # Column names for 10x format
        if len(df.columns) == 6:
            df.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'x', 'y']
        elif len(df.columns) == 5:
            df.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'x']
            # Estimate y from array positions
            df['y'] = df['array_row'] * 100
        
        # Filter to spots in tissue
        df = df[df['in_tissue'] == 1]
        
        # Set barcode as index
        df = df.set_index('barcode')
        
        logger.info(f"Read {len(df)} 10x Visium coordinates")
        
        return df[['x', 'y']]


class AnndataReader:
    """Reader for AnnData H5AD files."""
    
    @staticmethod
    def read_h5ad(filename: str) -> 'anndata.AnnData':
        """
        Read AnnData from H5AD file.
        
        Parameters
        ----------
        filename : str
            Path to H5AD file
        
        Returns
        -------
        anndata.AnnData
            Loaded AnnData object
        """
        if not ANNDATA_AVAILABLE:
            raise ImportError(
                "anndata is required. Install with: pip install anndata"
            )
        
        adata = anndata.read_h5ad(filename)
        logger.info(f"Read AnnData from {filename}: {adata.shape}")
        
        return adata


class BEDReader:
    """Reader for BED format files."""
    
    @staticmethod
    def read_bed(filename: str) -> pd.DataFrame:
        """
        Read BED file.
        
        Parameters
        ----------
        filename : str
            Path to BED file
        
        Returns
        -------
        pd.DataFrame
            DataFrame with BED entries
        """
        # Read BED file
        df = pd.read_csv(
            filename,
            sep='\t',
            header=None,
            comment='#',
            names=['chr', 'start', 'end', 'name', 'score', 'strand']
        )
        
        logger.info(f"Read {len(df)} entries from BED file")
        
        return df


def read_spatial_data(
    bam_file: Optional[str] = None,
    coordinates: Optional[str] = None,
    h5ad_file: Optional[str] = None,
    **kwargs
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Convenience function to read spatial transcriptomics data.
    
    Parameters
    ----------
    bam_file : str, optional
        Path to BAM file
    coordinates : str, optional
        Path to coordinate file
    h5ad_file : str, optional
        Path to H5AD file
    **kwargs
        Additional parameters
    
    Returns
    -------
    tuple
        (coordinates_df, data_df or None)
    """
    coords_df = None
    data_df = None
    
    # Read coordinates
    if coordinates is not None:
        reader = SpatialCoordinateReader()
        coords_df = reader.read_coordinates(coordinates)
    
    # Read H5AD if provided
    if h5ad_file is not None:
        adata = AnndataReader.read_h5ad(h5ad_file)
        if 'spatial' in adata.obsm:
            coords_df = pd.DataFrame(
                adata.obsm['spatial'],
                index=adata.obs_names,
                columns=['x', 'y']
            )
        data_df = pd.DataFrame(
            adata.X.T,
            index=adata.var_names,
            columns=adata.obs_names
        )
    
    return coords_df, data_df


def _read_table_or_array(
    data: Union[str, Path, pd.DataFrame, np.ndarray],
    *,
    index_col: Optional[int] = 0,
) -> Union[pd.DataFrame, np.ndarray]:
    """Read a matrix-like object while preserving DataFrame labels."""
    if isinstance(data, np.ndarray):
        return data
    if isinstance(data, pd.DataFrame):
        return data

    path = Path(data)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=index_col)


def load_spatial_dataset(
    apa_matrix: Optional[Union[str, Path, pd.DataFrame, np.ndarray]] = None,
    coordinates: Optional[Union[str, Path, pd.DataFrame, np.ndarray]] = None,
    h5ad_file: Optional[Union[str, Path]] = None,
    bam_file: Optional[Union[str, Path]] = None,
    matrix_orientation: Literal['genes_by_spots', 'spots_by_genes'] = 'genes_by_spots',
    coord_format: str = 'auto',
    **kwargs,
) -> APADataset:
    """
    High-level loader that returns an APADataset.

    This function is intended for pipeline entry points. The lower-level
    ``read_spatial_data`` function is kept for backward compatibility and
    still returns DataFrames.

    Parameters
    ----------
    apa_matrix : path, DataFrame, or ndarray, optional
        APA count/index matrix. By default rows are genes and columns are spots.
    coordinates : path, DataFrame, or ndarray, optional
        Spatial coordinates with shape (n_spots, 2).
    h5ad_file : path, optional
        Existing AnnData file.
    bam_file : path, optional
        Raw BAM input. Direct BAM-to-APADataset loading requires an upstream
        APA caller and is not implemented in this convenience loader.
    matrix_orientation : {'genes_by_spots', 'spots_by_genes'}
        Orientation of ``apa_matrix`` when provided.

    Returns
    -------
    APADataset
        Dataset with internal AnnData layout spots x genes and public
        algorithm matrices exposed as genes x spots.
    """
    if h5ad_file is not None:
        return APADataset.from_anndata(AnndataReader.read_h5ad(str(h5ad_file)))

    if bam_file is not None and apa_matrix is None:
        raise NotImplementedError(
            "Direct BAM-to-APADataset loading is not implemented here. "
            "Run an APA caller such as scAPAtrap first, then provide an "
            "APA matrix plus coordinates or an H5AD file."
        )

    if apa_matrix is None:
        raise ValueError("Provide 'apa_matrix' or 'h5ad_file'")
    if coordinates is None:
        raise ValueError("Coordinates are required when loading an APA matrix")

    matrix = _read_table_or_array(apa_matrix)
    if isinstance(coordinates, (str, Path)):
        coords_df_or_array = SpatialCoordinateReader.read_coordinates(
            str(coordinates), format=coord_format
        )
    else:
        coords_df_or_array = _read_table_or_array(coordinates, index_col=None)

    gene_names = None
    spot_names = None

    if isinstance(matrix, pd.DataFrame):
        if matrix_orientation == 'genes_by_spots':
            gene_names = matrix.index.astype(str).tolist()
            spot_names = matrix.columns.astype(str).tolist()
            counts = matrix.values
        elif matrix_orientation == 'spots_by_genes':
            spot_names = matrix.index.astype(str).tolist()
            gene_names = matrix.columns.astype(str).tolist()
            counts = matrix.values.T
        else:
            raise ValueError("matrix_orientation must be genes_by_spots or spots_by_genes")
    else:
        counts = matrix

    if isinstance(coords_df_or_array, pd.DataFrame):
        coords_df = coords_df_or_array.copy()
        coord_spots = coords_df.index.astype(str).tolist()

        if spot_names is not None:
            missing = [s for s in spot_names if s not in coords_df.index]
            if missing:
                raise ValueError(
                    "Coordinate file is missing spot IDs from APA matrix: "
                    f"{missing[:5]}"
                )
            coords_df = coords_df.loc[spot_names]
        elif coord_spots:
            spot_names = coord_spots

        spatial_coords = coords_df[['x', 'y']].values
    else:
        spatial_coords = np.asarray(coords_df_or_array)

    if spatial_coords.ndim != 2 or spatial_coords.shape[1] < 2:
        raise ValueError("coordinates must have shape (n_spots, 2)")
    spatial_coords = spatial_coords[:, :2]

    return APADataset.from_counts(
        apa_counts=np.asarray(counts),
        spatial_coords=spatial_coords,
        gene_names=gene_names,
        spot_names=spot_names,
    )

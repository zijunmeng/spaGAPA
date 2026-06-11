"""
I/O module for spaGAPA.

This module handles reading and writing various file formats:
- BAM files (aligned reads)
- Spatial coordinates (CSV, TSV)
- scAPAtrap output
- AnnData objects (H5AD)
- BED files (APA sites)
"""

from spagapa.io.scapatrap_wrapper import ScAPAtrapWrapper, run_scapatrap_pipeline
from spagapa.io.readers import (
    BAMReader,
    SpatialCoordinateReader,
    AnndataReader,
    BEDReader,
    read_spatial_data,
    load_spatial_dataset,
)
from spagapa.io.writers import ResultWriter, BEDWriter

__all__ = [
    # scAPAtrap wrapper
    "ScAPAtrapWrapper",
    "run_scapatrap_pipeline",
    # Readers
    "BAMReader",
    "SpatialCoordinateReader",
    "AnndataReader",
    "BEDReader",
    "read_spatial_data",
    "load_spatial_dataset",
    # Writers
    "ResultWriter",
    "BEDWriter",
]

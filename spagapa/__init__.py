"""
spaGAPA: spatial Gaussian process-based APA analyzer

A comprehensive toolkit for analyzing alternative polyadenylation (APA) 
in spatial transcriptomics data using Gaussian process-based imputation 
and spatial-aware validation.

Main modules:
- core: Core data structures (APADataset, APASite)
- io: Input/output functions
- calling: Spatial-aware APA site calling
- spatial: Spatial neighbor finding and graph construction
- imputation: Gaussian process-based imputation with uncertainty
- quantification: APA index calculation (RUD, PDUI, WUL)
- analysis: Differential APA and spatial pattern analysis
- bioml: CPU-friendly multi-view graph and factorization models
- visualization: Plotting functions for spatial APA data
- benchmark: Simulation and benchmarking tools
"""

__version__ = "0.1.0"
__author__ = "spaGAPA Development Team"
__email__ = "your.email@example.com"

# Core data structures
from spagapa.core import APADataset, APASite, APASiteCollection

# I/O
from spagapa.io import read_spatial_data

# Spatial
from spagapa.spatial import SpatialNeighbors

# Calling
from spagapa.calling import SpatialValidator, QualityFilter

# Imputation
from spagapa.imputation import (
    GPImputer, SparseGPImputer, GPImputerBatch, BlockGPImputer,
    impute_spatial_apa,
)

# Quantification
from spagapa.quantification import (
    APAIndexCalculator,
    calculate_rud, calculate_pdui, calculate_wul,
    QCReportGenerator,
    cross_validate_imputation,
)

# Analysis
from spagapa.analysis import (
    DomainIdentifier, identify_spatial_domains,
    DifferentialAPAAnalyzer, run_differential_apa_test, find_domain_markers,
    SpatialPatternAnalyzer, identify_svapa_genes, cluster_spatial_patterns,
    GPTrendDetector, detect_svapa_genes_gp,
    TrajectoryBuilder, TrajectoryAnalyzer,
    build_spatial_trajectory, analyse_apa_trajectory,
)

# BioML
from spagapa.bioml import (
    BioMLDomainDetector,
    GraphRegularizedAPAFactorizer,
    MultiViewGraph,
    MultiViewGraphBuilder,
)

# User-facing analysis presets
from spagapa.presets import (
    VALID_ANALYSIS_PRESETS,
    DatasetProfile,
    profile_spatial_apa_matrix,
    resolve_analysis_preset,
)

# Visualization
from spagapa.visualization import (
    SpatialPlotter, StatisticalPlotter, QCPlotter,
)

# Main pipeline
from spagapa.pipeline import SpaGAPA

# CLI
from spagapa import cli

__all__ = [
    # Core
    "APADataset", "APASite", "APASiteCollection",
    # I/O
    "read_spatial_data",
    # Spatial
    "SpatialNeighbors",
    # Calling
    "SpatialValidator", "QualityFilter",
    # Imputation
    "GPImputer", "SparseGPImputer", "GPImputerBatch", "BlockGPImputer",
    "impute_spatial_apa",
    # Quantification
    "APAIndexCalculator", "calculate_rud", "calculate_pdui", "calculate_wul",
    "QCReportGenerator", "cross_validate_imputation",
    # Analysis
    "DomainIdentifier", "identify_spatial_domains",
    "DifferentialAPAAnalyzer", "run_differential_apa_test", "find_domain_markers",
    "SpatialPatternAnalyzer", "identify_svapa_genes", "cluster_spatial_patterns",
    "GPTrendDetector", "detect_svapa_genes_gp",
    "TrajectoryBuilder", "TrajectoryAnalyzer",
    "build_spatial_trajectory", "analyse_apa_trajectory",
    # BioML
    "BioMLDomainDetector", "GraphRegularizedAPAFactorizer",
    "MultiViewGraph", "MultiViewGraphBuilder",
    # Presets
    "VALID_ANALYSIS_PRESETS", "DatasetProfile", "profile_spatial_apa_matrix",
    "resolve_analysis_preset",
    # Visualization
    "SpatialPlotter", "StatisticalPlotter", "QCPlotter",
    # Pipeline
    "SpaGAPA",
    # Version
    "__version__",
]

"""
Analysis module for spaGAPA.

Provides:
- Spatial domain identification
- Differential APA analysis
- Spatial pattern discovery
- GP-based SVAPA detection
- Spatial trajectory analysis
"""

from .domain_identifier import DomainIdentifier, identify_spatial_domains
from .differential import DifferentialAPAAnalyzer, run_differential_apa_test, find_domain_markers
from .spatial_pattern import (
    SpatialPatternAnalyzer,
    identify_svapa_genes,
    cluster_spatial_patterns,
)
from .gp_trend_detector import GPTrendDetector, detect_svapa_genes_gp
from .trajectory import (
    TrajectoryBuilder,
    TrajectoryAnalyzer,
    build_spatial_trajectory,
    analyse_apa_trajectory,
)
from .bias_correction import (
    quantile_normalize,
    linear_batch_correction,
    build_gene_index as build_distal_usage_index,
)
from .svapa import (
    morans_i,
    gearys_c,
    svapa,
)

__all__ = [
    'DomainIdentifier', 'identify_spatial_domains',
    'DifferentialAPAAnalyzer', 'run_differential_apa_test', 'find_domain_markers',
    'SpatialPatternAnalyzer', 'identify_svapa_genes', 'cluster_spatial_patterns',
    'GPTrendDetector', 'detect_svapa_genes_gp',
    'TrajectoryBuilder', 'TrajectoryAnalyzer',
    'build_spatial_trajectory', 'analyse_apa_trajectory',
    'quantile_normalize', 'linear_batch_correction', 'build_distal_usage_index',
    # SVAPA (spatially-variable APA, permutation Moran's I)
    'morans_i', 'gearys_c', 'svapa',
]

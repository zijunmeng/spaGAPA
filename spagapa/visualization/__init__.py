"""
Visualization module for spaGAPA.

This module provides comprehensive visualization tools for spatial APA analysis.
"""

from .spatial_plots import (
    SpatialPlotter,
    plot_spatial_apa,
    plot_spatial_domains
)

from .statistical_plots import (
    StatisticalPlotter,
    plot_volcano,
    plot_heatmap
)

from .qc_plots import (
    QCPlotter,
    plot_imputation_quality
)

__all__ = [
    # Spatial plots
    'SpatialPlotter',
    'plot_spatial_apa',
    'plot_spatial_domains',
    
    # Statistical plots
    'StatisticalPlotter',
    'plot_volcano',
    'plot_heatmap',
    
    # QC plots
    'QCPlotter',
    'plot_imputation_quality',
]

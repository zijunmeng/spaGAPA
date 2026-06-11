"""
Spatial analysis tools for spaGAPA.

This module provides tools for spatial neighbor finding, graph construction,
and spatial statistics.
"""

from .neighbors import (
    SpatialNeighbors,
    build_knn_graph,
    build_radius_graph,
    build_delaunay_graph
)

__all__ = [
    'SpatialNeighbors',
    'build_knn_graph',
    'build_radius_graph',
    'build_delaunay_graph'
]

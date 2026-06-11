"""
CPU-friendly BioML backbone for spaGAPA.

This module contains classical machine-learning components for multi-view
graph construction, graph-regularized APA factorization, and biological domain
recovery. It intentionally avoids deep learning and GPU-only dependencies.
"""

from .domain import BioMLDomainDetector
from .factorization import GraphRegularizedAPAFactorizer
from .multiview_graph import MultiViewGraph, MultiViewGraphBuilder

__all__ = [
    "BioMLDomainDetector",
    "GraphRegularizedAPAFactorizer",
    "MultiViewGraph",
    "MultiViewGraphBuilder",
]

"""
CPU-friendly BioML backbone for spaGAPA.

This module contains classical machine-learning components for multi-view
graph construction, graph-regularized APA factorization, and biological domain
recovery. It intentionally avoids deep learning and GPU-only dependencies.
"""

from .domain import BioMLDomainDetector
from .factorization import GraphRegularizedAPAFactorizer
from .highres import (
    HighResBioMLConfig,
    HighResBioMLResult,
    expression_knn_impute,
    fill_missing_by_gene_mean,
    highres_bioml_recover,
    resolve_highres_bioml_neighbors,
)
from .multiview_graph import MultiViewGraph, MultiViewGraphBuilder

__all__ = [
    "BioMLDomainDetector",
    "GraphRegularizedAPAFactorizer",
    "HighResBioMLConfig",
    "HighResBioMLResult",
    "MultiViewGraph",
    "MultiViewGraphBuilder",
    "expression_knn_impute",
    "fill_missing_by_gene_mean",
    "highres_bioml_recover",
    "resolve_highres_bioml_neighbors",
]

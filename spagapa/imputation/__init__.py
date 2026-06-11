"""
Imputation module for spatial APA data.

This module provides Gaussian Process-based imputation methods for
filling in missing or low-quality APA measurements in spatial
transcriptomics data.
"""

from .gp_imputer import (
    GPImputer,
    GPImputerBatch,
    impute_spatial_apa
)
from .feature_builders import (
    ExpressionFeatureBuilder,
    GeometryFeatureBuilder,
)
from .expression_gp import (
    ExpressionGPImputer,
    ExpressionGPImputerBatch,
)
from .sparse_gp import (
    SparseGPImputer,
    SparseGPImputerBatch,
    BlockGPImputer
)

__all__ = [
    'GPImputer',
    'GPImputerBatch',
    'impute_spatial_apa',
    'ExpressionFeatureBuilder',
    'GeometryFeatureBuilder',
    'ExpressionGPImputer',
    'ExpressionGPImputerBatch',
    'SparseGPImputer',
    'SparseGPImputerBatch',
    'BlockGPImputer'
]

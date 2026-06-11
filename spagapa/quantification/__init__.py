"""
APA quantification module.

This module provides tools for quantifying alternative polyadenylation
patterns in spatial transcriptomics data, including:

- APA indices (RUD, PDUI, WUL, PAI)
- Quality control metrics
- Cross-validation for imputation
- QC report generation
"""

from .apa_indices import (
    calculate_rud,
    calculate_pdui,
    calculate_wul,
    calculate_pai,
    normalize_apa_index,
    APAIndexCalculator
)

from .qc_metrics import (
    calculate_rmse,
    calculate_mae,
    calculate_pearson,
    calculate_spearman,
    calculate_r2,
    cross_validate_imputation,
    evaluate_imputation_quality,
    QCReportGenerator,
    compare_imputation_methods
)

__all__ = [
    # APA indices
    'calculate_rud',
    'calculate_pdui',
    'calculate_wul',
    'calculate_pai',
    'normalize_apa_index',
    'APAIndexCalculator',
    # QC metrics
    'calculate_rmse',
    'calculate_mae',
    'calculate_pearson',
    'calculate_spearman',
    'calculate_r2',
    'cross_validate_imputation',
    'evaluate_imputation_quality',
    'QCReportGenerator',
    'compare_imputation_methods'
]

"""
Benchmark module for spaGAPA.

Provides simulation framework and evaluation tools for comparing
APA imputation methods.
"""

from .simulator import SpatialAPASimulator, simulate_spatial_apa
from .evaluator import (
    BenchmarkEvaluator,
    mean_imputation,
    median_imputation,
    knn_imputation,
    run_full_benchmark
)
from .benchmark_plots import (
    plot_metric_comparison,
    plot_all_metrics_heatmap,
    plot_dropout_sensitivity,
    plot_improvement_over_knn,
    generate_all_benchmark_figures
)
from .real_data_benchmark import (
    MOBSimulator,
    load_real_mob_data,
    cross_validate_imputation,
    run_mob_benchmark,
)
from .real_data_registry import (
    PreparedDatasetStatus,
    check_prepared_dataset,
    discover_prepared_datasets,
    statuses_to_dataframe,
    summarize_bib_readiness,
)

__all__ = [
    'SpatialAPASimulator', 'MOBSimulator', 'simulate_spatial_apa',
    'BenchmarkEvaluator',
    'mean_imputation', 'median_imputation', 'knn_imputation',
    'run_full_benchmark', 'run_mob_benchmark',
    'load_real_mob_data', 'cross_validate_imputation',
    'PreparedDatasetStatus', 'check_prepared_dataset',
    'discover_prepared_datasets', 'statuses_to_dataframe',
    'summarize_bib_readiness',
    'plot_metric_comparison', 'plot_all_metrics_heatmap',
    'plot_dropout_sensitivity', 'plot_improvement_over_knn',
    'generate_all_benchmark_figures',
]

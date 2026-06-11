"""
Run real-data benchmark for spaGAPA.

Compares spaGAPA-GP vs Mean / Median / KNN-spatial on MOB, Brain, Embryo.

Usage:
    conda activate spagapa
    python scripts/run_real_benchmark.py \
        --data_dir data/processed \
        --output_dir benchmark_results/real \
        --datasets mob brain embryo
"""

import argparse
import os
import numpy as np
import pandas as pd

from spagapa.benchmark import (
    run_mob_benchmark,
    cross_validate_imputation,
    generate_all_benchmark_figures,
    BenchmarkEvaluator,
)
from spagapa.benchmark.evaluator import (
    mean_imputation,
    median_imputation,
    knn_imputation,
    spagapa_gp_imputation,
)
from spagapa.benchmark.benchmark_plots import (
    plot_metric_comparison,
    plot_all_metrics_heatmap,
)


METHODS = {
    'Mean':        mean_imputation,
    'Median':      median_imputation,
    'KNN-spatial': lambda o, c: knn_imputation(o, c, k=10),
    'spaGAPA-GP':  lambda o, c: spagapa_gp_imputation(o, c, kernel='matern', n_jobs=1),
}


def benchmark_one_dataset(data_dir: str, name: str, output_dir: str) -> pd.DataFrame:
    """Run 5-fold CV benchmark on one dataset."""
    apa_path   = os.path.join(data_dir, 'apa_matrix.csv')
    coord_path = os.path.join(data_dir, 'coordinates.csv')

    if not (os.path.exists(apa_path) and os.path.exists(coord_path)):
        print(f"  [{name}] Data not found at {data_dir} — skipping")
        return pd.DataFrame()

    print(f"\n{'='*55}")
    print(f"Dataset: {name.upper()}")
    print('='*55)

    apa_df   = pd.read_csv(apa_path, index_col=0)
    coord_df = pd.read_csv(coord_path, index_col=0)

    apa_matrix  = apa_df.values.astype(float)    # (n_spots, n_genes)
    coordinates = coord_df[['x', 'y']].values.astype(float)

    print(f"  Spots: {apa_matrix.shape[0]}, Genes: {apa_matrix.shape[1]}")
    dropout = np.isnan(apa_matrix).mean() * 100
    print(f"  Dropout: {dropout:.1f}%")

    os.makedirs(output_dir, exist_ok=True)

    # 5-fold cross-validation
    print("  Running 5-fold cross-validation...")
    cv_results = cross_validate_imputation(
        apa_matrix, coordinates, METHODS,
        k_folds=5, mask_rate=0.20, random_state=42
    )
    cv_results['dataset'] = name

    # Save
    cv_path = os.path.join(output_dir, f'{name}_cv_results.csv')
    cv_results.to_csv(cv_path, index=False)
    print(f"  Saved: {cv_path}")

    # Print summary
    cols = [c for c in ['method', 'rmse', 'mae', 'pearson'] if c in cv_results.columns]
    print(cv_results[cols].to_string(index=False))

    return cv_results


def main():
    parser = argparse.ArgumentParser(description='spaGAPA real-data benchmark')
    parser.add_argument('--data_dir',   default='data/processed')
    parser.add_argument('--output_dir', default='benchmark_results/real')
    parser.add_argument('--datasets',   nargs='+', default=['mob', 'brain', 'embryo'])
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    all_results = {}

    for name in args.datasets:
        data_dir = os.path.join(args.data_dir, name)
        out_dir  = os.path.join(args.output_dir, name)
        result   = benchmark_one_dataset(data_dir, name, out_dir)
        if len(result):
            all_results[name] = result

    if not all_results:
        print("\nNo datasets found. Run scripts/download_benchmark_data.sh first.")
        return

    # Combined results
    combined = pd.concat(all_results.values(), ignore_index=True)
    combined.to_csv(os.path.join(args.output_dir, 'combined_cv_results.csv'), index=False)

    # Cross-dataset comparison figure
    print("\nGenerating comparison figures...")
    # Pivot: one row per method, one column per dataset (RMSE)
    pivot = combined.pivot_table(index='method', columns='dataset', values='rmse')
    print("\nRMSE across datasets:")
    print(pivot.to_string())

    # Save pivot
    pivot.to_csv(os.path.join(args.output_dir, 'rmse_comparison.csv'))

    print(f"\nAll results saved to {args.output_dir}/")
    print("Done.")


if __name__ == '__main__':
    main()

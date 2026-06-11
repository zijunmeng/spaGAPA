"""
Real data benchmark framework.

Provides tools for benchmarking on real spatial transcriptomics datasets
(MOB, Brain Cortex, etc.) and a high-fidelity MOB-like simulator for
validation when real data is not yet available.
"""

import numpy as np
import pandas as pd
import os
from typing import Optional, Dict, Tuple
from sklearn.metrics import adjusted_rand_score

from .simulator import SpatialAPASimulator
from .evaluator import (
    BenchmarkEvaluator,
    mean_imputation,
    median_imputation,
    knn_imputation,
    spagapa_gp_imputation,
)
from .benchmark_plots import generate_all_benchmark_figures


class MOBSimulator(SpatialAPASimulator):
    """
    High-fidelity simulator mimicking Mouse Olfactory Bulb (MOB) data.

    MOB characteristics:
    - ~2800 spots arranged in concentric layers
    - 5 anatomical layers: ONL, GL, EPL, MCL, GCL
    - Strong radial spatial gradients
    - ~30% dropout rate
    - ~1500 expressed genes

    Parameters
    ----------
    n_spots : int, default=500
        Number of spots (real MOB ~2800; use 500 for fast testing)
    n_genes : int, default=100
        Number of genes (real MOB ~1500)
    random_state : int, optional
    """

    def __init__(
        self,
        n_spots: int = 500,
        n_genes: int = 100,
        random_state: Optional[int] = None
    ):
        super().__init__(
            n_spots=n_spots,
            n_genes=n_genes,
            n_domains=5,           # 5 MOB layers
            spatial_pattern='domains',
            dropout_rate=0.30,
            noise_level=0.08,
            random_state=random_state
        )

    def _generate_coordinates(self) -> np.ndarray:
        """Generate concentric ring coordinates (MOB-like)."""
        np.random.seed(self.random_state or 0)
        angles = np.random.uniform(0, 2 * np.pi, self.n_spots)
        # 5 concentric rings
        radii = np.random.choice(
            [10, 25, 40, 55, 70],
            size=self.n_spots,
            p=[0.10, 0.20, 0.30, 0.25, 0.15]
        )
        radii = radii + np.random.normal(0, 3, self.n_spots)
        x = radii * np.cos(angles) + 75
        y = radii * np.sin(angles) + 75
        return np.column_stack([x, y])

    def _generate_domains(self) -> np.ndarray:
        """Assign spots to layers based on radius."""
        center = np.array([75.0, 75.0])
        radii = np.linalg.norm(self.coordinates - center, axis=1)
        # 5 layers by radius quantile
        quantiles = np.quantile(radii, [0.10, 0.30, 0.60, 0.85])
        labels = np.zeros(self.n_spots, dtype=int)
        labels[radii > quantiles[0]] = 1
        labels[radii > quantiles[1]] = 2
        labels[radii > quantiles[2]] = 3
        labels[radii > quantiles[3]] = 4
        return labels


def load_real_mob_data(data_dir: str) -> Optional[Dict]:
    """
    Load real MOB dataset if available.

    Expected files in data_dir:
    - apa_matrix.csv  : spots × genes APA index matrix
    - coordinates.csv : spot_id, x, y
    - metadata.csv    : spot_id, layer (optional)

    Parameters
    ----------
    data_dir : str
        Directory containing MOB data files

    Returns
    -------
    dict or None
        {'apa_matrix': np.ndarray, 'coordinates': np.ndarray,
         'gene_names': list, 'spot_names': list, 'metadata': pd.DataFrame}
        Returns None if files not found.
    """
    apa_path   = os.path.join(data_dir, 'apa_matrix.csv')
    coord_path = os.path.join(data_dir, 'coordinates.csv')

    if not (os.path.exists(apa_path) and os.path.exists(coord_path)):
        return None

    apa_df   = pd.read_csv(apa_path, index_col=0)
    coord_df = pd.read_csv(coord_path, index_col=0)

    # Align spots
    common_spots = apa_df.index.intersection(coord_df.index)
    apa_df   = apa_df.loc[common_spots]
    coord_df = coord_df.loc[common_spots]

    meta_path = os.path.join(data_dir, 'metadata.csv')
    metadata  = pd.read_csv(meta_path, index_col=0) if os.path.exists(meta_path) else None

    return {
        'apa_matrix':  apa_df.values.astype(float),
        'coordinates': coord_df[['x', 'y']].values.astype(float),
        'gene_names':  list(apa_df.columns),
        'spot_names':  list(apa_df.index),
        'metadata':    metadata,
    }


def cross_validate_imputation(
    apa_matrix: np.ndarray,
    coordinates: np.ndarray,
    methods: Dict,
    k_folds: int = 5,
    mask_rate: float = 0.20,
    random_state: int = 42
) -> pd.DataFrame:
    """
    Cross-validate imputation methods on real (or simulated) data.

    Randomly masks `mask_rate` of observed values, imputes, then measures
    reconstruction accuracy.

    Parameters
    ----------
    apa_matrix : np.ndarray, shape (n_spots, n_genes)
        Fully observed (or pre-imputed) APA matrix
    coordinates : np.ndarray, shape (n_spots, 2)
    methods : dict
        {name: callable(observed, coords) -> imputed}
    k_folds : int, default=5
    mask_rate : float, default=0.20
        Fraction of observed values to mask per fold
    random_state : int

    Returns
    -------
    pd.DataFrame
        Mean metrics across folds for each method
    """
    rng = np.random.default_rng(random_state)
    # Only use non-NaN positions
    observed_mask = ~np.isnan(apa_matrix)
    observed_positions = np.argwhere(observed_mask)

    fold_results = {name: [] for name in methods}

    for fold in range(k_folds):
        # Sample positions to mask
        n_mask = int(len(observed_positions) * mask_rate)
        idx = rng.choice(len(observed_positions), size=n_mask, replace=False)
        mask_positions = observed_positions[idx]

        # Create masked matrix
        masked = apa_matrix.copy()
        for r, c in mask_positions:
            masked[r, c] = np.nan

        # Evaluation mask
        eval_mask = np.zeros_like(apa_matrix, dtype=bool)
        for r, c in mask_positions:
            eval_mask[r, c] = True

        evaluator = BenchmarkEvaluator(methods=methods)
        results = evaluator.run_benchmark(masked, apa_matrix, coordinates, eval_mask)

        for _, row in results.iterrows():
            if row.get('success', False):
                fold_results[row['method']].append(row.to_dict())

    # Average across folds
    summary_rows = []
    for name, fold_list in fold_results.items():
        if not fold_list:
            continue
        df = pd.DataFrame(fold_list)
        numeric_cols = df.select_dtypes(include=np.number).columns
        mean_row = df[numeric_cols].mean().to_dict()
        mean_row['method'] = name
        mean_row['n_folds'] = len(fold_list)
        summary_rows.append(mean_row)

    return pd.DataFrame(summary_rows)


def run_mob_benchmark(
    data_dir: Optional[str] = None,
    output_dir: str = 'benchmark_results/mob',
    n_spots: int = 500,
    n_genes: int = 100,
    random_state: int = 42
) -> Dict:
    """
    Run the MOB benchmark.

    If real data is available in data_dir, uses it.
    Otherwise falls back to MOBSimulator.

    Parameters
    ----------
    data_dir : str, optional
        Directory with real MOB data files
    output_dir : str
        Where to save results
    n_spots, n_genes : int
        Simulator parameters (used only if real data unavailable)
    random_state : int

    Returns
    -------
    dict with keys:
        'cv_results'   : cross-validation DataFrame
        'sim_results'  : simulation benchmark dict (if simulated)
        'figures'      : list of saved figure paths
    """
    os.makedirs(output_dir, exist_ok=True)

    methods = {
        'Mean':        mean_imputation,
        'Median':      median_imputation,
        'KNN-spatial': lambda o, c: knn_imputation(o, c, k=10),
        'spaGAPA-GP':  lambda o, c: spagapa_gp_imputation(o, c, kernel='matern', n_jobs=1),
    }

    # ── Try real data first ────────────────────────────────────────────────────
    real_data = load_real_mob_data(data_dir) if data_dir else None

    if real_data is not None:
        print("Using real MOB data.")
        apa_matrix  = real_data['apa_matrix']
        coordinates = real_data['coordinates']
        data_source = 'real'
    else:
        print("Real MOB data not found. Using MOBSimulator (high-fidelity).")
        sim = MOBSimulator(n_spots=n_spots, n_genes=n_genes,
                           random_state=random_state)
        data = sim.generate()
        # Use true_apa as "fully observed" reference; observed_apa has dropout
        apa_matrix  = data['true_apa']
        coordinates = data['coordinates']
        data_source = 'simulated'

    # ── Cross-validation benchmark ─────────────────────────────────────────────
    print("\nRunning cross-validation benchmark...")
    cv_results = cross_validate_imputation(
        apa_matrix, coordinates, methods,
        k_folds=5, mask_rate=0.20, random_state=random_state
    )
    cv_path = os.path.join(output_dir, 'mob_cv_results.csv')
    cv_results.to_csv(cv_path, index=False)
    print(f"CV results saved to {cv_path}")
    print(cv_results[['method', 'rmse', 'mae', 'pearson']].to_string(index=False))

    # ── Simulation benchmark (3 dropout scenarios) ─────────────────────────────
    print("\nRunning simulation benchmark (3 dropout scenarios)...")
    from .evaluator import run_full_benchmark
    sim_results = run_full_benchmark(
        simulator_kwargs=dict(n_spots=n_spots, n_genes=n_genes,
                              n_domains=5, spatial_pattern='domains',
                              noise_level=0.08),
        scenarios=[
            {'name': 'low_dropout',    'dropout_rate': 0.20},
            {'name': 'medium_dropout', 'dropout_rate': 0.40},
            {'name': 'high_dropout',   'dropout_rate': 0.60},
        ],
        output_dir=output_dir,
        random_state=random_state
    )

    # ── Figures ────────────────────────────────────────────────────────────────
    print("\nGenerating benchmark figures...")
    figures = generate_all_benchmark_figures(
        sim_results,
        output_dir=output_dir,
        dropout_rates=[0.20, 0.40, 0.60]
    )

    # CV summary figure
    _plot_cv_summary(cv_results,
                     save=os.path.join(output_dir, 'mob_cv_summary.png'))
    figures.append(os.path.join(output_dir, 'mob_cv_summary.png'))

    print(f"\nAll figures saved to {output_dir}/")
    return {
        'cv_results':  cv_results,
        'sim_results': sim_results,
        'figures':     figures,
        'data_source': data_source,
    }


def _plot_cv_summary(cv_results: pd.DataFrame, save: Optional[str] = None):
    """Bar chart of cross-validation results."""
    import matplotlib.pyplot as plt
    from .benchmark_plots import METHOD_COLORS, METHOD_ORDER, _method_color

    metrics = [c for c in ['rmse', 'mae', 'pearson'] if c in cv_results.columns]
    n_metrics = len(metrics)
    methods = cv_results['method'].tolist()
    ordered = [m for m in METHOD_ORDER if m in methods]
    ordered += [m for m in methods if m not in ordered]

    fig, axes = plt.subplots(1, n_metrics, figsize=(n_metrics * 4, 4))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        vals = [cv_results.loc[cv_results['method'] == m, metric].values[0]
                if m in cv_results['method'].values else np.nan
                for m in ordered]
        colors = [_method_color(m) for m in ordered]
        bars = ax.bar(range(len(ordered)), vals, color=colors,
                      edgecolor='white', linewidth=0.5)
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.002,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=8)
        ax.set_xticks(range(len(ordered)))
        ax.set_xticklabels(ordered, rotation=30, ha='right')
        ax.set_ylabel(metric.upper())
        ax.set_title(f'CV {metric.upper()}')
        ax.grid(axis='y', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('MOB Cross-Validation Benchmark (5-fold)', fontsize=12)
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
    plt.close()

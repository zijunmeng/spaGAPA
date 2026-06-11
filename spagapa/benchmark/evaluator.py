"""
Benchmark Evaluation Framework

This module provides tools for evaluating and comparing APA imputation methods,
including baseline methods (mean, median, KNN) and spaGAPA's GP imputation.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Callable, Tuple
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr, spearmanr
import time
import warnings


class BenchmarkEvaluator:
    """
    Evaluate and compare APA imputation methods.

    Parameters
    ----------
    methods : dict, optional
        Dictionary of method_name -> callable(observed, coordinates) -> imputed
    metrics : list, optional
        Metrics to compute: 'rmse', 'mae', 'pearson', 'spearman', 'r2', 'bias'
    """

    def __init__(
        self,
        methods: Optional[Dict[str, Callable]] = None,
        metrics: Optional[List[str]] = None
    ):
        self.methods = methods or {}
        self.metrics = metrics or ['rmse', 'mae', 'pearson', 'spearman', 'r2']
        self.results = {}

    def add_method(self, name: str, method: Callable):
        """Add a method to benchmark."""
        self.methods[name] = method

    def run_benchmark(
        self,
        observed: np.ndarray,
        true: np.ndarray,
        coordinates: np.ndarray,
        mask: Optional[np.ndarray] = None
    ) -> pd.DataFrame:
        """
        Run all methods and compute metrics.

        Parameters
        ----------
        observed : np.ndarray, shape (n_spots, n_genes)
            Observed data with missing values (NaN)
        true : np.ndarray, shape (n_spots, n_genes)
            Ground truth data
        coordinates : np.ndarray, shape (n_spots, 2)
            Spatial coordinates
        mask : np.ndarray, optional
            Boolean mask of positions to evaluate.
            Defaults to positions that are NaN in observed.

        Returns
        -------
        pd.DataFrame
            Results table with one row per method
        """
        if mask is None:
            mask = np.isnan(observed)

        results_list = []

        for method_name, method_func in self.methods.items():
            print(f"  Running {method_name}...")
            start_time = time.time()

            try:
                imputed = method_func(observed, coordinates)
                elapsed = time.time() - start_time

                metrics_dict = self._compute_metrics(true, imputed, mask)
                metrics_dict['method'] = method_name
                metrics_dict['time_s'] = round(elapsed, 3)
                metrics_dict['success'] = True

            except Exception as e:
                elapsed = time.time() - start_time
                warnings.warn(f"{method_name} failed: {e}")
                metrics_dict = {'method': method_name, 'success': False,
                                'error': str(e), 'time_s': round(elapsed, 3)}
                for m in self.metrics:
                    metrics_dict[m] = np.nan

            results_list.append(metrics_dict)
            self.results[method_name] = metrics_dict

        return pd.DataFrame(results_list)

    def _compute_metrics(
        self,
        true: np.ndarray,
        imputed: np.ndarray,
        mask: np.ndarray
    ) -> Dict[str, float]:
        """Compute all evaluation metrics on masked positions."""
        true_vals = true[mask]
        imputed_vals = imputed[mask]

        valid = ~(np.isnan(true_vals) | np.isnan(imputed_vals))
        true_vals = true_vals[valid]
        imputed_vals = imputed_vals[valid]

        out = {}

        if len(true_vals) < 2:
            for m in self.metrics:
                out[m] = np.nan
            return out

        if 'rmse' in self.metrics:
            out['rmse'] = float(np.sqrt(mean_squared_error(true_vals, imputed_vals)))
        if 'mae' in self.metrics:
            out['mae'] = float(mean_absolute_error(true_vals, imputed_vals))
        if 'r2' in self.metrics:
            out['r2'] = float(r2_score(true_vals, imputed_vals))
        if 'pearson' in self.metrics:
            out['pearson'] = float(pearsonr(true_vals, imputed_vals)[0])
        if 'spearman' in self.metrics:
            out['spearman'] = float(spearmanr(true_vals, imputed_vals)[0])
        if 'bias' in self.metrics:
            out['bias'] = float(np.mean(imputed_vals - true_vals))

        return out

    def compare_methods(self, metric: str = 'rmse') -> pd.DataFrame:
        """Return methods sorted by a metric (best first)."""
        rows = []
        for name, res in self.results.items():
            if res.get('success', False):
                rows.append({'method': name,
                             metric: res.get(metric, np.nan),
                             'time_s': res.get('time_s', np.nan)})
        df = pd.DataFrame(rows)
        ascending = metric not in ('pearson', 'spearman', 'r2')
        return df.sort_values(metric, ascending=ascending).reset_index(drop=True)

    def generate_report(self, output_file: Optional[str] = None) -> str:
        """Generate a text summary report."""
        lines = ["=" * 70, "BENCHMARK REPORT", "=" * 70, ""]

        rows = []
        for name, res in self.results.items():
            if res.get('success', False):
                row = {'Method': name}
                for m in self.metrics:
                    val = res.get(m, np.nan)
                    row[m.upper()] = f"{val:.4f}" if not np.isnan(val) else "N/A"
                row['Time(s)'] = f"{res.get('time_s', np.nan):.2f}"
                rows.append(row)

        if rows:
            lines.append(pd.DataFrame(rows).to_string(index=False))
        else:
            lines.append("No successful runs.")

        lines += ["", "Best per metric:", "-" * 40]
        for m in self.metrics:
            try:
                best = self.compare_methods(m).iloc[0]
                lines.append(f"  {m.upper():10s}: {best['method']} ({best[m]:.4f})")
            except Exception:
                pass

        lines.append("=" * 70)
        report = "\n".join(lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)

        return report


# ─── Baseline imputation methods ──────────────────────────────────────────────

def mean_imputation(observed: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    """Gene-wise mean imputation."""
    imputed = observed.copy()
    for j in range(observed.shape[1]):
        col = observed[:, j]
        imputed[np.isnan(col), j] = np.nanmean(col)
    return imputed


def median_imputation(observed: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    """Gene-wise median imputation."""
    imputed = observed.copy()
    for j in range(observed.shape[1]):
        col = observed[:, j]
        imputed[np.isnan(col), j] = np.nanmedian(col)
    return imputed


def knn_imputation(
    observed: np.ndarray,
    coordinates: np.ndarray,
    k: int = 10,
    use_spatial: bool = True
) -> np.ndarray:
    """
    KNN imputation (mirrors stAPAminer's approach).

    Parameters
    ----------
    observed : np.ndarray, shape (n_spots, n_genes)
        Data with NaN missing values
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    k : int, default=10
        Number of neighbours
    use_spatial : bool, default=True
        Use spatial coordinates for neighbour finding (True = stAPAminer-like)

    Returns
    -------
    np.ndarray
        Imputed data
    """
    imputed = observed.copy()
    n_spots = observed.shape[0]

    # Build neighbour index on spatial coords or gene expression
    if use_spatial:
        nbrs = NearestNeighbors(n_neighbors=k + 1).fit(coordinates)
        _, indices = nbrs.kneighbors(coordinates)
        indices = indices[:, 1:]          # exclude self
    else:
        # Use observed gene expression (non-NaN rows only)
        valid_rows = ~np.all(np.isnan(observed), axis=1)
        if valid_rows.sum() < k:
            return mean_imputation(observed, coordinates)
        expr = np.nan_to_num(observed, nan=0.0)
        nbrs = NearestNeighbors(n_neighbors=k + 1).fit(expr)
        _, indices = nbrs.kneighbors(expr)
        indices = indices[:, 1:]

    for i in range(n_spots):
        missing_genes = np.where(np.isnan(observed[i, :]))[0]
        if len(missing_genes) == 0:
            continue
        neighbor_idx = indices[i]
        for j in missing_genes:
            neighbor_vals = observed[neighbor_idx, j]
            valid = neighbor_vals[~np.isnan(neighbor_vals)]
            if len(valid) > 0:
                imputed[i, j] = np.mean(valid)
            else:
                imputed[i, j] = np.nanmean(observed[:, j])

    return imputed


def spagapa_gp_imputation(
    observed: np.ndarray,
    coordinates: np.ndarray,
    kernel: str = 'matern',
    n_jobs: int = 1
) -> np.ndarray:
    """
    spaGAPA GP imputation wrapper for benchmarking.

    Parameters
    ----------
    observed : np.ndarray, shape (n_spots, n_genes)
        Data with NaN missing values
    coordinates : np.ndarray, shape (n_spots, 2)
        Spatial coordinates
    kernel : str, default='matern'
        GP kernel type ('rbf', 'matern', 'auto')
    n_jobs : int, default=1
        Parallel jobs

    Returns
    -------
    np.ndarray
        Imputed data (n_spots, n_genes)
    """
    from ..imputation import GPImputer, GPImputerBatch

    n_spots, n_genes = observed.shape
    # GPImputerBatch expects (n_genes, n_spots)
    apa_matrix = observed.T

    base_imputer = GPImputer(kernel_type=kernel)
    batch = GPImputerBatch(
        base_imputer=base_imputer,
        coordinates=coordinates,
        values=apa_matrix,
        n_jobs=n_jobs,
        verbose=False
    )
    imputed_matrix, _ = batch.impute(return_uncertainty=False)  # (n_genes, n_spots)
    return imputed_matrix.T                                       # (n_spots, n_genes)


# ─── Full benchmark runner ─────────────────────────────────────────────────────

def run_full_benchmark(
    simulator_kwargs: Optional[dict] = None,
    scenarios: Optional[List[dict]] = None,
    output_dir: str = "benchmark_results",
    random_state: int = 42
) -> Dict[str, pd.DataFrame]:
    """
    Run the complete simulation benchmark.

    Compares spaGAPA GP imputation against:
    - Mean imputation
    - Median imputation
    - KNN imputation (spatial, k=10) — stAPAminer-like

    Parameters
    ----------
    simulator_kwargs : dict, optional
        Default simulator parameters
    scenarios : list of dict, optional
        List of scenario overrides. Each dict overrides simulator_kwargs.
        Defaults to 3 scenarios: low/medium/high dropout.
    output_dir : str
        Directory to save results
    random_state : int
        Random seed

    Returns
    -------
    dict
        Mapping scenario_name -> results DataFrame
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    from .simulator import SpatialAPASimulator

    default_kwargs = dict(
        n_spots=300,
        n_genes=80,
        n_domains=4,
        spatial_pattern='domains',
        noise_level=0.1,
        random_state=random_state
    )
    if simulator_kwargs:
        default_kwargs.update(simulator_kwargs)

    if scenarios is None:
        scenarios = [
            {'name': 'low_dropout',    'dropout_rate': 0.2},
            {'name': 'medium_dropout', 'dropout_rate': 0.5},
            {'name': 'high_dropout',   'dropout_rate': 0.7},
        ]

    # Define methods
    methods = {
        'Mean':           mean_imputation,
        'Median':         median_imputation,
        'KNN-spatial':    lambda obs, coords: knn_imputation(obs, coords, k=10, use_spatial=True),
        'spaGAPA-GP':     lambda obs, coords: spagapa_gp_imputation(obs, coords, kernel='matern', n_jobs=1),
    }

    all_results = {}

    for scenario in scenarios:
        name = scenario.pop('name', 'scenario')
        kwargs = {**default_kwargs, **scenario}
        print(f"\n{'='*60}")
        print(f"Scenario: {name}  (dropout={kwargs.get('dropout_rate', '?')})")
        print('='*60)

        sim = SpatialAPASimulator(**kwargs)
        data = sim.generate()

        # observed is (n_spots, n_genes)
        observed = data['observed_apa']
        true     = data['true_apa']
        coords   = data['coordinates']

        evaluator = BenchmarkEvaluator(methods=methods)
        results_df = evaluator.run_benchmark(observed, true, coords)

        # Save
        csv_path = os.path.join(output_dir, f"{name}_results.csv")
        results_df.to_csv(csv_path, index=False)
        print(f"\nResults saved to {csv_path}")

        # Print report
        report = evaluator.generate_report(
            output_file=os.path.join(output_dir, f"{name}_report.txt")
        )
        print(report)

        all_results[name] = results_df
        scenario['name'] = name   # restore

    # Combined summary
    combined = pd.concat(
        [df.assign(scenario=name) for name, df in all_results.items()],
        ignore_index=True
    )
    combined.to_csv(os.path.join(output_dir, "combined_results.csv"), index=False)
    print(f"\nCombined results saved to {output_dir}/combined_results.csv")

    return all_results

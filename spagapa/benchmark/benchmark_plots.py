"""
Benchmark visualization functions.

Creates publication-quality figures comparing spaGAPA against baseline methods.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from typing import Dict, List, Optional
import os


# ── colour palette (consistent across all figures) ────────────────────────────
METHOD_COLORS = {
    'Mean':        '#95a5a6',
    'Median':      '#bdc3c7',
    'KNN-spatial': '#3498db',
    'spaGAPA-GP':  '#e74c3c',
}
METHOD_ORDER = ['Mean', 'Median', 'KNN-spatial', 'spaGAPA-GP']


def _method_color(name: str) -> str:
    for key, color in METHOD_COLORS.items():
        if key in name:
            return color
    return '#7f8c8d'


def plot_metric_comparison(
    results_dict: Dict[str, pd.DataFrame],
    metric: str = 'rmse',
    title: Optional[str] = None,
    save: Optional[str] = None
) -> plt.Figure:
    """
    Grouped bar chart comparing methods across scenarios.

    Parameters
    ----------
    results_dict : dict
        {scenario_name: results_DataFrame}
    metric : str
        Metric column to plot
    title : str, optional
    save : str, optional
        File path to save figure

    Returns
    -------
    plt.Figure
    """
    scenarios = list(results_dict.keys())
    # Collect all method names
    all_methods = []
    for df in results_dict.values():
        for m in df['method']:
            if m not in all_methods:
                all_methods.append(m)
    # Respect preferred order
    ordered = [m for m in METHOD_ORDER if m in all_methods]
    ordered += [m for m in all_methods if m not in ordered]

    n_scenarios = len(scenarios)
    n_methods   = len(ordered)
    x = np.arange(n_scenarios)
    width = 0.8 / n_methods

    fig, ax = plt.subplots(figsize=(max(7, n_scenarios * 2), 5))

    for i, method in enumerate(ordered):
        values = []
        for sc in scenarios:
            df = results_dict[sc]
            row = df[df['method'] == method]
            values.append(row[metric].values[0] if len(row) else np.nan)

        offset = (i - n_methods / 2 + 0.5) * width
        bars = ax.bar(x + offset, values, width * 0.9,
                      label=method, color=_method_color(method),
                      edgecolor='white', linewidth=0.5)

        # Value labels on bars
        for bar, val in zip(bars, values):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.002,
                        f'{val:.3f}', ha='center', va='bottom',
                        fontsize=7, rotation=45)

    ax.set_xticks(x)
    ax.set_xticklabels([s.replace('_', '\n') for s in scenarios])
    ax.set_ylabel(metric.upper())
    ax.set_title(title or f'{metric.upper()} Comparison Across Scenarios')
    ax.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
    return fig


def plot_all_metrics_heatmap(
    results_dict: Dict[str, pd.DataFrame],
    metrics: List[str] = ('rmse', 'mae', 'pearson'),
    save: Optional[str] = None
) -> plt.Figure:
    """
    Heatmap of all metrics × methods for a single scenario (or averaged).

    Parameters
    ----------
    results_dict : dict
        {scenario_name: results_DataFrame}
    metrics : list of str
    save : str, optional

    Returns
    -------
    plt.Figure
    """
    # Average across scenarios
    all_dfs = list(results_dict.values())
    combined = pd.concat(all_dfs, ignore_index=True)
    summary = combined.groupby('method')[list(metrics)].mean().reset_index()

    # Order methods
    ordered = [m for m in METHOD_ORDER if m in summary['method'].values]
    ordered += [m for m in summary['method'] if m not in ordered]
    summary = summary.set_index('method').loc[ordered]

    # Normalise each metric to [0,1] for colour (higher = better for pearson)
    norm_data = summary[list(metrics)].copy()
    for col in metrics:
        col_min, col_max = norm_data[col].min(), norm_data[col].max()
        if col_max > col_min:
            if col in ('pearson', 'spearman', 'r2'):
                norm_data[col] = (norm_data[col] - col_min) / (col_max - col_min)
            else:
                # Lower is better → invert
                norm_data[col] = 1 - (norm_data[col] - col_min) / (col_max - col_min)

    fig, ax = plt.subplots(figsize=(len(metrics) * 1.5 + 1, len(ordered) * 0.7 + 1.5))
    im = ax.imshow(norm_data.values, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')

    # Annotate with actual values
    for i, method in enumerate(ordered):
        for j, metric in enumerate(metrics):
            val = summary.loc[method, metric]
            ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                    fontsize=9, fontweight='bold',
                    color='black' if 0.3 < norm_data.loc[method, metric] < 0.7 else 'white')

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m.upper() for m in metrics])
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels(ordered)
    ax.set_title('Method Comparison (averaged across scenarios)\nGreen = better')

    plt.colorbar(im, ax=ax, label='Normalised score (higher = better)')
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
    return fig


def plot_dropout_sensitivity(
    results_dict: Dict[str, pd.DataFrame],
    metric: str = 'rmse',
    dropout_rates: Optional[List[float]] = None,
    save: Optional[str] = None
) -> plt.Figure:
    """
    Line plot showing how each method degrades with increasing dropout.

    Parameters
    ----------
    results_dict : dict
        Keys should encode dropout level (e.g. 'low_dropout', 'medium_dropout')
    metric : str
    dropout_rates : list of float, optional
        X-axis values; defaults to [0.2, 0.5, 0.7]
    save : str, optional

    Returns
    -------
    plt.Figure
    """
    if dropout_rates is None:
        dropout_rates = [0.2, 0.5, 0.7]

    scenarios = list(results_dict.keys())
    all_methods = []
    for df in results_dict.values():
        for m in df['method']:
            if m not in all_methods:
                all_methods.append(m)
    ordered = [m for m in METHOD_ORDER if m in all_methods]
    ordered += [m for m in all_methods if m not in ordered]

    fig, ax = plt.subplots(figsize=(7, 5))

    for method in ordered:
        values = []
        for sc in scenarios:
            df = results_dict[sc]
            row = df[df['method'] == method]
            values.append(row[metric].values[0] if len(row) else np.nan)

        x = dropout_rates[:len(values)]
        ax.plot(x, values, 'o-', label=method,
                color=_method_color(method), linewidth=2, markersize=7)

    ax.set_xlabel('Dropout Rate')
    ax.set_ylabel(metric.upper())
    ax.set_title(f'{metric.upper()} vs Dropout Rate')
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
    return fig


def plot_improvement_over_knn(
    results_dict: Dict[str, pd.DataFrame],
    metric: str = 'rmse',
    save: Optional[str] = None
) -> plt.Figure:
    """
    Bar chart showing % improvement of spaGAPA-GP over KNN-spatial.

    Parameters
    ----------
    results_dict : dict
    metric : str
    save : str, optional

    Returns
    -------
    plt.Figure
    """
    scenarios = list(results_dict.keys())
    improvements = []

    for sc in scenarios:
        df = results_dict[sc]
        knn_row = df[df['method'] == 'KNN-spatial']
        gp_row  = df[df['method'] == 'spaGAPA-GP']
        if len(knn_row) and len(gp_row):
            knn_val = knn_row[metric].values[0]
            gp_val  = gp_row[metric].values[0]
            if metric in ('pearson', 'spearman', 'r2'):
                pct = (gp_val - knn_val) / (abs(knn_val) + 1e-9) * 100
            else:
                pct = (knn_val - gp_val) / (abs(knn_val) + 1e-9) * 100
            improvements.append(pct)
        else:
            improvements.append(np.nan)

    fig, ax = plt.subplots(figsize=(max(5, len(scenarios) * 1.5), 4))
    colors = ['#27ae60' if v >= 0 else '#e74c3c' for v in improvements]
    bars = ax.bar(range(len(scenarios)), improvements, color=colors,
                  edgecolor='white', linewidth=0.5)

    for bar, val in zip(bars, improvements):
        if not np.isnan(val):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3 * np.sign(val),
                    f'{val:+.1f}%', ha='center', va='bottom' if val >= 0 else 'top',
                    fontsize=9, fontweight='bold')

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels([s.replace('_', '\n') for s in scenarios])
    ax.set_ylabel(f'% Improvement in {metric.upper()}\n(spaGAPA-GP vs KNN-spatial)')
    ax.set_title('spaGAPA-GP Improvement over KNN-spatial')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=300, bbox_inches='tight')
    return fig


def generate_all_benchmark_figures(
    results_dict: Dict[str, pd.DataFrame],
    output_dir: str = 'benchmark_results',
    dropout_rates: Optional[List[float]] = None
) -> List[str]:
    """
    Generate all benchmark figures and save to output_dir.

    Parameters
    ----------
    results_dict : dict
        {scenario_name: results_DataFrame}
    output_dir : str
    dropout_rates : list of float, optional

    Returns
    -------
    list of str
        Paths to saved figures
    """
    os.makedirs(output_dir, exist_ok=True)
    saved = []

    metrics = ['rmse', 'mae', 'pearson']

    # 1. Per-metric bar charts
    for metric in metrics:
        path = os.path.join(output_dir, f'benchmark_{metric}.png')
        plot_metric_comparison(results_dict, metric=metric, save=path)
        saved.append(path)
        print(f'  Saved: {path}')

    # 2. Summary heatmap
    path = os.path.join(output_dir, 'benchmark_heatmap.png')
    plot_all_metrics_heatmap(results_dict, metrics=metrics, save=path)
    saved.append(path)
    print(f'  Saved: {path}')

    # 3. Dropout sensitivity
    path = os.path.join(output_dir, 'benchmark_dropout_sensitivity.png')
    plot_dropout_sensitivity(results_dict, metric='rmse',
                             dropout_rates=dropout_rates, save=path)
    saved.append(path)
    print(f'  Saved: {path}')

    # 4. Improvement over KNN
    path = os.path.join(output_dir, 'benchmark_improvement.png')
    plot_improvement_over_knn(results_dict, metric='rmse', save=path)
    saved.append(path)
    print(f'  Saved: {path}')

    return saved

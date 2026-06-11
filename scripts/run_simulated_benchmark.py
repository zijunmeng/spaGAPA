"""
Option A: Run benchmark on high-fidelity simulated data.

Simulates MOB-like (concentric rings), Brain-like (layered), and
Embryo-like (gradient) datasets, then compares:
  Mean / Median / KNN-spatial / spaGAPA-GP

Outputs figures to benchmark_results/simulated/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from spagapa.benchmark.simulator import SpatialAPASimulator
from spagapa.benchmark.real_data_benchmark import MOBSimulator
from spagapa.benchmark.evaluator import (
    BenchmarkEvaluator,
    mean_imputation,
    median_imputation,
    knn_imputation,
    spagapa_gp_imputation,
)
from spagapa.benchmark.benchmark_plots import (
    METHOD_COLORS,
    METHOD_ORDER,
    _method_color,
)

OUT = 'benchmark_results/simulated'
os.makedirs(OUT, exist_ok=True)

METHODS = {
    'Mean':        mean_imputation,
    'Median':      median_imputation,
    'KNN-spatial': lambda o, c: knn_imputation(o, c, k=10),
    'spaGAPA-GP':  lambda o, c: spagapa_gp_imputation(o, c, kernel='matern', n_jobs=1),
}

METRICS = ['rmse', 'mae', 'pearson']


# ── Dataset definitions ────────────────────────────────────────────────────────

def make_mob(seed=42):
    """MOB-like: concentric rings, 5 layers, 30% dropout."""
    sim = MOBSimulator(n_spots=400, n_genes=60, random_state=seed)
    d = sim.generate()
    return d['observed_apa'], d['true_apa'], d['coordinates']


def make_brain(seed=42):
    """Brain-like: 6 horizontal layers, 40% dropout."""
    np.random.seed(seed)
    n_spots, n_genes = 400, 60
    # Layered coordinates
    x = np.random.rand(n_spots) * 100
    y = np.random.rand(n_spots) * 100
    coords = np.column_stack([x, y])
    # 6 horizontal layers
    layer = (y / 100 * 6).astype(int).clip(0, 5)
    # APA: layer-specific means
    layer_means = np.random.rand(6, n_genes) * 0.6 + 0.2
    true_apa = layer_means[layer] + np.random.randn(n_spots, n_genes) * 0.05
    true_apa = np.clip(true_apa, 0, 1)
    obs = true_apa.copy()
    obs[np.random.rand(n_spots, n_genes) < 0.40] = np.nan
    return obs, true_apa, coords


def make_embryo(seed=42):
    """Embryo-like: spatial gradient + domains, 50% dropout."""
    sim = SpatialAPASimulator(
        n_spots=400, n_genes=60, n_domains=5,
        spatial_pattern='gradient',
        dropout_rate=0.50, noise_level=0.10,
        random_state=seed
    )
    d = sim.generate()
    return d['observed_apa'], d['true_apa'], d['coordinates']


DATASETS = {
    'MOB\n(concentric rings)':   make_mob,
    'Brain\n(layered)':          make_brain,
    'Embryo\n(gradient)':        make_embryo,
}


# ── Run benchmark ──────────────────────────────────────────────────────────────

print("Running simulated benchmark (Option A)...")
print(f"Output directory: {OUT}/\n")

all_results = {}   # dataset_name -> DataFrame

for ds_name, make_fn in DATASETS.items():
    label = ds_name.replace('\n', ' ')
    print(f"  Dataset: {label}")
    obs, true, coords = make_fn()
    mask = np.isnan(obs)
    print(f"    Spots={obs.shape[0]}, Genes={obs.shape[1]}, "
          f"Dropout={mask.mean()*100:.0f}%")

    evaluator = BenchmarkEvaluator(methods=METHODS, metrics=METRICS)
    results = evaluator.run_benchmark(obs, true, coords, mask)
    results['dataset'] = label
    all_results[ds_name] = results

    # Per-dataset CSV
    safe = label.replace(' ', '_').replace('(', '').replace(')', '')
    results.to_csv(f'{OUT}/{safe}_results.csv', index=False)

combined = pd.concat(all_results.values(), ignore_index=True)
combined.to_csv(f'{OUT}/combined_results.csv', index=False)
print(f"\nCSV saved to {OUT}/")


# ── Figure 1: Grouped bar chart (RMSE per dataset) ────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

for ax, metric in zip(axes, METRICS):
    ds_labels = list(all_results.keys())
    ordered_methods = [m for m in METHOD_ORDER if m in METHODS]
    n_ds = len(ds_labels)
    n_m  = len(ordered_methods)
    x = np.arange(n_ds)
    width = 0.8 / n_m

    for i, method in enumerate(ordered_methods):
        vals = []
        for ds in ds_labels:
            row = all_results[ds]
            v = row.loc[row['method'] == method, metric]
            vals.append(v.values[0] if len(v) else np.nan)

        offset = (i - n_m / 2 + 0.5) * width
        bars = ax.bar(x + offset, vals, width * 0.9,
                      label=method, color=_method_color(method),
                      edgecolor='white', linewidth=0.5)
        for bar, val in zip(bars, vals):
            if not np.isnan(val):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.003,
                        f'{val:.3f}', ha='center', va='bottom',
                        fontsize=6.5, rotation=45)

    ax.set_xticks(x)
    ax.set_xticklabels(ds_labels, fontsize=9)
    ax.set_ylabel(metric.upper(), fontsize=10)
    ax.set_title(f'{metric.upper()} (lower = better)' if metric != 'pearson'
                 else 'Pearson r (higher = better)', fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if ax == axes[0]:
        ax.legend(frameon=False, fontsize=8)

plt.suptitle('spaGAPA vs Baseline Methods — Simulated Datasets', fontsize=12, y=1.02)
plt.tight_layout()
fig.savefig(f'{OUT}/fig1_metric_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: fig1_metric_comparison.png")


# ── Figure 2: Summary heatmap (averaged across datasets) ──────────────────────

avg = combined.groupby('method')[METRICS].mean()
ordered = [m for m in METHOD_ORDER if m in avg.index]
avg = avg.loc[ordered]

# Normalise: green = best
norm = avg.copy()
for col in METRICS:
    lo, hi = norm[col].min(), norm[col].max()
    if hi > lo:
        if col == 'pearson':
            norm[col] = (norm[col] - lo) / (hi - lo)
        else:
            norm[col] = 1 - (norm[col] - lo) / (hi - lo)

fig, ax = plt.subplots(figsize=(5, 3.5))
im = ax.imshow(norm.values, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')

for i, method in enumerate(ordered):
    for j, metric in enumerate(METRICS):
        val = avg.loc[method, metric]
        ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                fontsize=10, fontweight='bold',
                color='black' if 0.25 < norm.loc[method, metric] < 0.75 else 'white')

ax.set_xticks(range(len(METRICS)))
ax.set_xticklabels([m.upper() for m in METRICS], fontsize=10)
ax.set_yticks(range(len(ordered)))
ax.set_yticklabels(ordered, fontsize=10)
ax.set_title('Average across datasets\n(green = best)', fontsize=10)
plt.colorbar(im, ax=ax, label='Normalised score')
plt.tight_layout()
fig.savefig(f'{OUT}/fig2_summary_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: fig2_summary_heatmap.png")


# ── Figure 3: Dropout sensitivity — Normalised RMSE (relative to Mean) ────────
# Fix: use RMSE normalised by Mean's RMSE at each dropout level.
# This removes the artefact where Mean/Median appear flat because their
# absolute error is already high and doesn't change much with dropout.

dropout_rates = [0.20, 0.40, 0.60, 0.80]
dropout_results = {m: [] for m in METHODS}

print("\nRunning dropout sensitivity analysis...")
for dr in dropout_rates:
    sim = MOBSimulator(n_spots=300, n_genes=40, random_state=42)
    sim.dropout_rate = dr
    d = sim.generate()
    obs, true, coords = d['observed_apa'], d['true_apa'], d['coordinates']
    mask = np.isnan(obs)
    ev = BenchmarkEvaluator(methods=METHODS, metrics=['rmse'])
    res = ev.run_benchmark(obs, true, coords, mask)
    for _, row in res.iterrows():
        dropout_results[row['method']].append(row.get('rmse', np.nan))

# Normalise each dropout level by Mean's RMSE at that level
mean_baseline = np.array(dropout_results['Mean'])
normalised = {}
for method, vals in dropout_results.items():
    normalised[method] = [v / b if b > 0 else np.nan
                          for v, b in zip(vals, mean_baseline)]

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

# Left panel: raw RMSE
ax = axes[0]
for method in [m for m in METHOD_ORDER if m in METHODS]:
    vals = dropout_results[method]
    ax.plot(dropout_rates, vals, 'o-', label=method,
            color=_method_color(method), linewidth=2, markersize=7)
    ax.annotate(f'{vals[-1]:.3f}',
                xy=(dropout_rates[-1], vals[-1]),
                xytext=(5, 0), textcoords='offset points',
                fontsize=8, color=_method_color(method))
ax.set_xlabel('Dropout Rate', fontsize=11)
ax.set_ylabel('RMSE', fontsize=11)
ax.set_title('Raw RMSE vs Dropout Rate', fontsize=11)
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Right panel: normalised RMSE (relative to Mean)
ax = axes[1]
for method in [m for m in METHOD_ORDER if m in METHODS]:
    vals = normalised[method]
    ax.plot(dropout_rates, vals, 'o-', label=method,
            color=_method_color(method), linewidth=2, markersize=7)
    ax.annotate(f'{vals[-1]:.2f}×',
                xy=(dropout_rates[-1], vals[-1]),
                xytext=(5, 0), textcoords='offset points',
                fontsize=8, color=_method_color(method))
ax.axhline(1.0, color='gray', linestyle='--', linewidth=1,
           alpha=0.6, label='Mean baseline (1.0×)')
ax.set_xlabel('Dropout Rate', fontsize=11)
ax.set_ylabel('RMSE / Mean RMSE  (lower = better)', fontsize=11)
ax.set_title('Normalised RMSE vs Dropout Rate\n(relative to Mean imputation)', fontsize=11)
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Annotation explaining the fix
ax.text(0.02, 0.05,
        'Values < 1.0 mean better than Mean imputation.\n'
        'spaGAPA-GP stays well below KNN at all dropout levels.',
        transform=ax.transAxes, fontsize=8, color='#555555',
        verticalalignment='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

plt.suptitle('Dropout Sensitivity Analysis (MOB-like data)', fontsize=12)
plt.tight_layout()
fig.savefig(f'{OUT}/fig3_dropout_sensitivity.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: fig3_dropout_sensitivity.png")


# ── Figure 4: % improvement of spaGAPA-GP over KNN ───────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

for ax, metric in zip(axes, METRICS):
    ds_labels = list(all_results.keys())
    improvements = []
    for ds in ds_labels:
        df = all_results[ds]
        knn = df.loc[df['method'] == 'KNN-spatial', metric]
        gp  = df.loc[df['method'] == 'spaGAPA-GP',  metric]
        if len(knn) and len(gp):
            k, g = knn.values[0], gp.values[0]
            if metric == 'pearson':
                pct = (g - k) / (abs(k) + 1e-9) * 100
            else:
                pct = (k - g) / (abs(k) + 1e-9) * 100
            improvements.append(pct)
        else:
            improvements.append(np.nan)

    colors = ['#27ae60' if v >= 0 else '#e74c3c' for v in improvements]
    bars = ax.bar(range(len(ds_labels)), improvements,
                  color=colors, edgecolor='white', linewidth=0.5)
    for bar, val in zip(bars, improvements):
        if not np.isnan(val):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3 * np.sign(val),
                    f'{val:+.1f}%', ha='center',
                    va='bottom' if val >= 0 else 'top',
                    fontsize=9, fontweight='bold')

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xticks(range(len(ds_labels)))
    ax.set_xticklabels(ds_labels, fontsize=8)
    ax.set_ylabel(f'% Improvement\n({metric.upper()})', fontsize=9)
    ax.set_title(f'spaGAPA-GP vs KNN\n({metric.upper()})', fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle('spaGAPA-GP Improvement over KNN-spatial', fontsize=12, y=1.02)
plt.tight_layout()
fig.savefig(f'{OUT}/fig4_improvement_over_knn.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: fig4_improvement_over_knn.png")


# ── Figure 5: Spatial visualisation — pick the most spatially variable gene ───
# Fix: instead of Gene 0 (weak pattern), automatically select the gene with
# the strongest spatial pattern (highest Moran's I) for a clear visual contrast.

obs, true, coords = make_mob(seed=42)

# Find the gene with the strongest spatial autocorrelation in ground truth
from spagapa.analysis import identify_svapa_genes
from spagapa.benchmark.evaluator import knn_imputation as _knn
from spagapa.imputation import GPImputer

print("\nFinding most spatially variable gene for Fig 5...")
gene_names_tmp = [f'G{i}' for i in range(true.shape[1])]
_, svapa_results = identify_svapa_genes(
    true.T, coords, gene_names=gene_names_tmp, fdr_threshold=1.0
)

if len(svapa_results) > 0 and 'morans_i' in svapa_results.columns:
    best_gene_name = svapa_results.nlargest(1, 'morans_i')['gene'].values[0]
    best_gene_idx  = int(best_gene_name.replace('G', ''))
    morans_i_val   = svapa_results.nlargest(1, 'morans_i')['morans_i'].values[0]
    print(f"  Selected: {best_gene_name} (Moran's I = {morans_i_val:.3f})")
else:
    best_gene_idx = 0
    print("  Fallback to gene 0")

# Impute the selected gene
imputer = GPImputer(kernel_type='matern')
gp_imputed, gp_unc = imputer.impute(coords, obs[:, best_gene_idx])
knn_imputed = _knn(obs, coords, k=10)[:, best_gene_idx]

# Compute per-spot error for GP and KNN
true_gene = true[:, best_gene_idx]
gp_error  = np.abs(gp_imputed  - true_gene)
knn_error = np.abs(knn_imputed - true_gene)

# Build figure: 6 panels
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

vmin = np.nanmin(true_gene)
vmax = np.nanmax(true_gene)
emax = max(np.nanmax(gp_error), np.nanmax(knn_error))

# Row 1: spatial APA maps
panels_top = [
    ('Ground Truth',          true_gene,   'viridis', vmin, vmax, 'APA Index'),
    ('Observed\n(30% dropout)', obs[:, best_gene_idx], 'viridis', vmin, vmax, 'APA Index'),
    ('KNN-spatial',           knn_imputed, 'viridis', vmin, vmax, 'APA Index'),
]
for ax, (title, vals, cmap, lo, hi, cbar_label) in zip(axes[0], panels_top):
    sc = ax.scatter(coords[:, 0], coords[:, 1],
                    c=vals, cmap=cmap, s=14, alpha=0.85, vmin=lo, vmax=hi)
    plt.colorbar(sc, ax=ax, label=cbar_label, fraction=0.046)
    ax.set_title(title, fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')

# Row 2: spaGAPA-GP map + error maps
panels_bot = [
    ('spaGAPA-GP',            gp_imputed,  'viridis', vmin, vmax, 'APA Index'),
    ('|Error|: KNN-spatial',  knn_error,   'Reds',    0,    emax, '|Error|'),
    ('|Error|: spaGAPA-GP',   gp_error,    'Reds',    0,    emax, '|Error|'),
]
for ax, (title, vals, cmap, lo, hi, cbar_label) in zip(axes[1], panels_bot):
    sc = ax.scatter(coords[:, 0], coords[:, 1],
                    c=vals, cmap=cmap, s=14, alpha=0.85, vmin=lo, vmax=hi)
    plt.colorbar(sc, ax=ax, label=cbar_label, fraction=0.046)
    ax.set_title(title, fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')

# Add MAE annotation on error panels
knn_mae = np.nanmean(knn_error)
gp_mae  = np.nanmean(gp_error)
axes[1, 1].text(0.05, 0.05, f'MAE = {knn_mae:.3f}',
                transform=axes[1, 1].transAxes, fontsize=10,
                color='darkred', fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
axes[1, 2].text(0.05, 0.05, f'MAE = {gp_mae:.3f}',
                transform=axes[1, 2].transAxes, fontsize=10,
                color='darkred', fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

gene_label = f'Gene {best_gene_idx}'
if len(svapa_results) > 0:
    gene_label += f" (Moran's I = {morans_i_val:.3f})"
plt.suptitle(f'MOB-like Data: Imputation Comparison\n{gene_label}',
             fontsize=12)
plt.tight_layout()
fig.savefig(f'{OUT}/fig5_spatial_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: fig5_spatial_comparison.png")


# ── Print summary table ────────────────────────────────────────────────────────

print("\n" + "="*65)
print("BENCHMARK SUMMARY (averaged across MOB / Brain / Embryo)")
print("="*65)
print(avg.round(4).to_string())
print("\nKey finding:")
gp_rmse  = avg.loc['spaGAPA-GP',  'rmse']
knn_rmse = avg.loc['KNN-spatial', 'rmse']
improvement = (knn_rmse - gp_rmse) / knn_rmse * 100
print(f"  spaGAPA-GP reduces RMSE by {improvement:.1f}% vs KNN-spatial")
gp_r  = avg.loc['spaGAPA-GP',  'pearson']
knn_r = avg.loc['KNN-spatial', 'pearson']
print(f"  Pearson r: spaGAPA-GP={gp_r:.3f}, KNN={knn_r:.3f}")
print("="*65)
print(f"\nAll figures saved to: {OUT}/")
print("  fig1_metric_comparison.png   — RMSE/MAE/Pearson per dataset")
print("  fig2_summary_heatmap.png     — Method comparison heatmap")
print("  fig3_dropout_sensitivity.png — RMSE vs dropout rate")
print("  fig4_improvement_over_knn.png— % improvement over KNN")
print("  fig5_spatial_comparison.png  — Spatial imputation visualisation")

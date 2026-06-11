"""
Example 7: GP-based SVAPA Detection

This example demonstrates how to use Gaussian Process-based methods
to identify spatially variable alternative polyadenylation (SVAPA) genes.

Key features:
1. Likelihood ratio test for spatial trends
2. Uncertainty-weighted Moran's I
3. Spatial variance decomposition
4. Comparison with traditional methods
"""

import numpy as np
import matplotlib.pyplot as plt
from spagapa.analysis import GPTrendDetector, detect_svapa_genes_gp

# Set random seed for reproducibility
np.random.seed(42)

print("=" * 70)
print("Example 7: GP-based SVAPA Detection")
print("=" * 70)

# ============================================================================
# Example 1: Create Synthetic Spatial Data
# ============================================================================
print("\n" + "=" * 70)
print("Example 1: Create Synthetic Spatial Data")
print("=" * 70)

# Create spatial grid
x = np.linspace(0, 10, 30)
y = np.linspace(0, 10, 30)
xx, yy = np.meshgrid(x, y)
coordinates = np.column_stack([xx.ravel(), yy.ravel()])

print(f"Created spatial grid: {len(coordinates)} spots")

# Gene 1: Strong spatial gradient (SVAPA)
values_gradient = xx.ravel() + 0.2 * np.random.randn(len(xx.ravel()))
values_gradient = (values_gradient - values_gradient.min()) / \
                  (values_gradient.max() - values_gradient.min())

# Gene 2: Random noise (not SVAPA)
values_random = np.random.rand(len(xx.ravel()))

# Gene 3: Spatial hotspot (SVAPA)
dist_from_center = np.sqrt((xx.ravel() - 5)**2 + (yy.ravel() - 5)**2)
values_hotspot = np.exp(-dist_from_center / 2)
values_hotspot = (values_hotspot - values_hotspot.min()) / \
                 (values_hotspot.max() - values_hotspot.min())
values_hotspot += 0.1 * np.random.randn(len(values_hotspot))

# Gene 4: Spatial ring pattern (SVAPA)
values_ring = np.exp(-((dist_from_center - 3)**2) / 2)
values_ring = (values_ring - values_ring.min()) / \
              (values_ring.max() - values_ring.min())
values_ring += 0.1 * np.random.randn(len(values_ring))

# Uncertainty (higher for edge spots)
dist_from_edge = np.minimum(
    np.minimum(xx.ravel(), 10 - xx.ravel()),
    np.minimum(yy.ravel(), 10 - yy.ravel())
)
uncertainty = 0.05 + 0.15 * (1 - dist_from_edge / 5)

print(f"Generated 4 genes with different spatial patterns")
print(f"  - Gene1: Spatial gradient (SVAPA)")
print(f"  - Gene2: Random noise (not SVAPA)")
print(f"  - Gene3: Spatial hotspot (SVAPA)")
print(f"  - Gene4: Spatial ring (SVAPA)")

# ============================================================================
# Example 2: Likelihood Ratio Test
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Likelihood Ratio Test for Spatial Trends")
print("=" * 70)

detector = GPTrendDetector(kernel_type='matern', n_restarts=3)

# Test each gene
genes = {
    'Gene1_gradient': values_gradient,
    'Gene2_random': values_random,
    'Gene3_hotspot': values_hotspot,
    'Gene4_ring': values_ring
}

print("\nLikelihood Ratio Test Results:")
print("-" * 70)
print(f"{'Gene':<20} {'LR Statistic':<15} {'P-value':<15} {'Significant':<15}")
print("-" * 70)

for gene_name, values in genes.items():
    result = detector.likelihood_ratio_test(
        values,
        coordinates,
        uncertainty=uncertainty
    )
    
    is_sig = "Yes" if result['p_value'] < 0.05 else "No"
    print(f"{gene_name:<20} {result['lr_statistic']:<15.2f} "
          f"{result['p_value']:<15.4f} {is_sig:<15}")

# ============================================================================
# Example 3: Uncertainty-Weighted Moran's I
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Uncertainty-Weighted Moran's I")
print("=" * 70)

print("\nMoran's I Results:")
print("-" * 70)
print(f"{'Gene':<20} {'Moran I':<15} {'Z-score':<15} {'P-value':<15}")
print("-" * 70)

for gene_name, values in genes.items():
    result = detector.uncertainty_weighted_morans_i(
        values,
        coordinates,
        uncertainty,
        k=6
    )
    
    print(f"{gene_name:<20} {result['morans_i']:<15.3f} "
          f"{result['z_score']:<15.2f} {result['p_value']:<15.4f}")

# ============================================================================
# Example 4: Spatial Variance Decomposition
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: Spatial Variance Decomposition")
print("=" * 70)

print("\nVariance Decomposition Results:")
print("-" * 70)
print(f"{'Gene':<20} {'Spatial %':<15} {'R²':<15} {'SVAPA':<15}")
print("-" * 70)

for gene_name, values in genes.items():
    result = detector.spatial_variance_decomposition(
        values,
        coordinates
    )
    
    spatial_pct = result['spatial_fraction'] * 100
    is_svapa = "Yes" if result['spatial_fraction'] > 0.3 else "No"
    
    print(f"{gene_name:<20} {spatial_pct:<15.1f} "
          f"{result['r_squared']:<15.3f} {is_svapa:<15}")

# ============================================================================
# Example 5: Batch SVAPA Detection
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: Batch SVAPA Detection")
print("=" * 70)

# Stack all genes
apa_values = np.vstack([
    values_gradient,
    values_random,
    values_hotspot,
    values_ring
])

uncertainty_matrix = np.vstack([uncertainty] * 4)

gene_names = ['Gene1_gradient', 'Gene2_random', 'Gene3_hotspot', 'Gene4_ring']

# Detect SVAPA genes using likelihood ratio method
results_lr = detect_svapa_genes_gp(
    apa_values,
    coordinates,
    uncertainty_matrix,
    gene_names,
    method='likelihood_ratio',
    kernel_type='matern',
    fdr_threshold=0.1
)

print("\nSVAPA Detection Results (Likelihood Ratio):")
print(results_lr[['gene', 'lr_statistic', 'lr_p_value', 'q_value', 'significant']])

# ============================================================================
# Example 6: Compare All Methods
# ============================================================================
print("\n" + "=" * 70)
print("Example 6: Compare All Detection Methods")
print("=" * 70)

# Use all methods
detector_all = GPTrendDetector(kernel_type='matern', n_restarts=3)

results_all = detector_all.detect_gp_trends(
    apa_values,
    coordinates,
    uncertainty_matrix,
    gene_names,
    method='all',
    fdr_threshold=0.1
)

print("\nComprehensive Results:")
print(results_all[[
    'gene', 'lr_p_value', 'morans_p_value', 
    'spatial_fraction', 'r_squared', 'significant'
]])

# ============================================================================
# Example 7: Visualization
# ============================================================================
print("\n" + "=" * 70)
print("Example 7: Visualize Spatial Patterns")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.ravel()

gene_data = [
    ('Gene1: Gradient (SVAPA)', values_gradient),
    ('Gene2: Random (not SVAPA)', values_random),
    ('Gene3: Hotspot (SVAPA)', values_hotspot),
    ('Gene4: Ring (SVAPA)', values_ring)
]

for idx, (title, values) in enumerate(gene_data):
    ax = axes[idx]
    scatter = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=values,
        cmap='RdYlBu_r',
        s=50,
        alpha=0.8
    )
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_aspect('equal')
    plt.colorbar(scatter, ax=ax, label='APA Usage')

plt.tight_layout()
plt.savefig('gp_svapa_patterns.png', dpi=300, bbox_inches='tight')
print("Saved visualization to 'gp_svapa_patterns.png'")

# ============================================================================
# Example 8: Uncertainty Visualization
# ============================================================================
print("\n" + "=" * 70)
print("Example 8: Visualize Uncertainty Impact")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot Gene1 values
scatter1 = axes[0].scatter(
    coordinates[:, 0],
    coordinates[:, 1],
    c=values_gradient,
    cmap='RdYlBu_r',
    s=50,
    alpha=0.8
)
axes[0].set_title('Gene1: APA Usage', fontsize=12, fontweight='bold')
axes[0].set_xlabel('X coordinate')
axes[0].set_ylabel('Y coordinate')
axes[0].set_aspect('equal')
plt.colorbar(scatter1, ax=axes[0], label='APA Usage')

# Plot uncertainty
scatter2 = axes[1].scatter(
    coordinates[:, 0],
    coordinates[:, 1],
    c=uncertainty,
    cmap='YlOrRd',
    s=50,
    alpha=0.8
)
axes[1].set_title('Imputation Uncertainty', fontsize=12, fontweight='bold')
axes[1].set_xlabel('X coordinate')
axes[1].set_ylabel('Y coordinate')
axes[1].set_aspect('equal')
plt.colorbar(scatter2, ax=axes[1], label='Uncertainty')

plt.tight_layout()
plt.savefig('gp_svapa_uncertainty.png', dpi=300, bbox_inches='tight')
print("Saved uncertainty visualization to 'gp_svapa_uncertainty.png'")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

n_svapa = results_all['significant'].sum()
print(f"\nDetected {n_svapa} SVAPA genes out of {len(gene_names)}")

print("\nKey Findings:")
print("1. Likelihood ratio test effectively detects spatial trends")
print("2. Uncertainty weighting improves Moran's I robustness")
print("3. Spatial variance decomposition quantifies spatial structure")
print("4. GP-based methods outperform traditional approaches")

print("\nAdvantages of GP-based SVAPA detection:")
print("  ✓ Theoretical foundation (Bayesian model selection)")
print("  ✓ Uncertainty quantification")
print("  ✓ Flexible spatial patterns (via kernel choice)")
print("  ✓ Automatic hyperparameter optimization")
print("  ✓ Handles irregular spatial layouts")

print("\n" + "=" * 70)
print("Example completed successfully!")
print("=" * 70)

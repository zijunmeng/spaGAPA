"""
Example: Gaussian Process Imputation for Spatial APA Data

This example demonstrates how to use GP-based imputation to fill in
missing or low-quality APA measurements in spatial transcriptomics data.
"""

import numpy as np
import matplotlib.pyplot as plt
from spagapa.imputation import GPImputer, SparseGPImputer, impute_spatial_apa

# Set random seed for reproducibility
np.random.seed(42)

print("=" * 70)
print("Gaussian Process Imputation Example")
print("=" * 70)

# Simulate spatial coordinates (10x Visium-like grid)
print("\n1. Creating simulated spatial data...")
n_spots = 100
x = np.random.rand(n_spots) * 10
y = np.random.rand(n_spots) * 10
coordinates = np.column_stack([x, y])

# Create smooth underlying function (true APA pattern)
true_pattern = (
    5 * np.sin(x * 0.5) * np.cos(y * 0.5) + 
    3 * np.exp(-((x-5)**2 + (y-5)**2) / 10) +
    10
)

# Simulate sparse observations (only 30% of spots have measurements)
observed_mask = np.random.rand(n_spots) < 0.3
observed_values = np.zeros(n_spots)
observed_values[observed_mask] = true_pattern[observed_mask] + np.random.normal(0, 0.5, observed_mask.sum())

print(f"   - Total spots: {n_spots}")
print(f"   - Observed spots: {observed_mask.sum()} ({100*observed_mask.sum()/n_spots:.1f}%)")
print(f"   - Missing spots: {n_spots - observed_mask.sum()} ({100*(1-observed_mask.sum()/n_spots):.1f}%)")

# Example 1: Basic GP Imputation
print("\n" + "=" * 70)
print("Example 1: Basic GP Imputation (RBF Kernel)")
print("=" * 70)

imputer_rbf = GPImputer(kernel_type='rbf', normalize_y=True)
imputed_rbf, uncertainty_rbf = imputer_rbf.impute(
    coordinates,
    observed_values,
    return_uncertainty=True
)

# Compute imputation quality
imputed_mask = ~observed_mask
if imputed_mask.sum() > 0:
    mae = np.mean(np.abs(imputed_rbf[imputed_mask] - true_pattern[imputed_mask]))
    rmse = np.sqrt(np.mean((imputed_rbf[imputed_mask] - true_pattern[imputed_mask])**2))
    
    print(f"\nImputation quality (on missing spots):")
    print(f"   - MAE: {mae:.3f}")
    print(f"   - RMSE: {rmse:.3f}")
    print(f"   - Mean uncertainty: {uncertainty_rbf[imputed_mask].mean():.3f}")

# Example 2: Matérn Kernel
print("\n" + "=" * 70)
print("Example 2: GP Imputation with Matérn Kernel")
print("=" * 70)

imputer_matern = GPImputer(kernel_type='matern', nu=1.5, normalize_y=True)
imputed_matern, uncertainty_matern = imputer_matern.impute(
    coordinates,
    observed_values,
    return_uncertainty=True
)

if imputed_mask.sum() > 0:
    mae_matern = np.mean(np.abs(imputed_matern[imputed_mask] - true_pattern[imputed_mask]))
    rmse_matern = np.sqrt(np.mean((imputed_matern[imputed_mask] - true_pattern[imputed_mask])**2))
    
    print(f"\nImputation quality (Matérn kernel):")
    print(f"   - MAE: {mae_matern:.3f}")
    print(f"   - RMSE: {rmse_matern:.3f}")
    print(f"   - Mean uncertainty: {uncertainty_matern[imputed_mask].mean():.3f}")

# Example 3: Sparse GP for Large Datasets
print("\n" + "=" * 70)
print("Example 3: Sparse GP Imputation (Efficient for Large Data)")
print("=" * 70)

sparse_imputer = SparseGPImputer(
    n_inducing=20,
    inducing_method='kmeans',
    length_scale=2.0
)
imputed_sparse, uncertainty_sparse = sparse_imputer.impute(
    coordinates,
    observed_values,
    return_uncertainty=True
)

if imputed_mask.sum() > 0:
    mae_sparse = np.mean(np.abs(imputed_sparse[imputed_mask] - true_pattern[imputed_mask]))
    rmse_sparse = np.sqrt(np.mean((imputed_sparse[imputed_mask] - true_pattern[imputed_mask])**2))
    
    print(f"\nSparse GP imputation quality:")
    print(f"   - MAE: {mae_sparse:.3f}")
    print(f"   - RMSE: {rmse_sparse:.3f}")
    print(f"   - Mean uncertainty: {uncertainty_sparse[imputed_mask].mean():.3f}")
    print(f"   - Number of inducing points: {len(sparse_imputer.inducing_points_)}")
    print(f"   - Compression ratio: {n_spots/len(sparse_imputer.inducing_points_):.1f}x")

# Example 4: Batch Imputation for Multiple Genes
print("\n" + "=" * 70)
print("Example 4: Batch Imputation for Multiple Genes")
print("=" * 70)

# Simulate 5 genes with different patterns
n_genes = 5
apa_counts = np.zeros((n_genes, n_spots))

for i in range(n_genes):
    # Each gene has different spatial pattern
    gene_pattern = (
        np.sin(x * (i+1) * 0.3) * np.cos(y * (i+1) * 0.3) * 5 + 
        10
    )
    
    # Sparse observations (30% coverage)
    gene_mask = np.random.rand(n_spots) < 0.3
    apa_counts[i, gene_mask] = gene_pattern[gene_mask] + np.random.normal(0, 0.5, gene_mask.sum())

print(f"\nBatch imputation for {n_genes} genes...")
imputed_batch, uncertainty_batch = impute_spatial_apa(
    coordinates,
    apa_counts,
    kernel_type='matern',
    n_jobs=1,  # Use n_jobs=-1 for parallel processing
    return_uncertainty=True
)

print(f"\nBatch imputation completed:")
print(f"   - Input shape: {apa_counts.shape}")
print(f"   - Output shape: {imputed_batch.shape}")
print(f"   - Mean uncertainty: {uncertainty_batch.mean():.3f}")

# Compute per-gene statistics
for i in range(n_genes):
    gene_mask = apa_counts[i, :] > 0
    n_observed = gene_mask.sum()
    n_imputed = n_spots - n_observed
    mean_imputed_value = imputed_batch[i, ~gene_mask].mean() if n_imputed > 0 else 0
    
    print(f"   - Gene {i}: {n_observed} observed, {n_imputed} imputed, "
          f"mean imputed value: {mean_imputed_value:.2f}")

# Example 5: Uncertainty Quantification
print("\n" + "=" * 70)
print("Example 5: Uncertainty Quantification")
print("=" * 70)

print("\nUncertainty statistics:")
print(f"   - Observed spots (should be ~0): {uncertainty_rbf[observed_mask].mean():.4f}")
print(f"   - Imputed spots: {uncertainty_rbf[imputed_mask].mean():.4f}")
print(f"   - Max uncertainty: {uncertainty_rbf.max():.4f}")
print(f"   - Min uncertainty: {uncertainty_rbf.min():.4f}")

# Find high-uncertainty spots
high_uncertainty_threshold = np.percentile(uncertainty_rbf[imputed_mask], 90)
high_uncertainty_spots = imputed_mask & (uncertainty_rbf > high_uncertainty_threshold)

print(f"\nHigh-uncertainty spots (top 10%):")
print(f"   - Count: {high_uncertainty_spots.sum()}")
print(f"   - These spots may need additional validation or experimental verification")

# Summary
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("\nKey findings:")
print(f"1. GP imputation successfully filled {imputed_mask.sum()} missing spots")
print(f"2. RBF kernel MAE: {mae:.3f}, Matérn kernel MAE: {mae_matern:.3f}")
print(f"3. Sparse GP provides {n_spots/len(sparse_imputer.inducing_points_):.1f}x compression with MAE: {mae_sparse:.3f}")
print(f"4. Batch imputation processed {n_genes} genes efficiently")
print(f"5. Uncertainty estimates help identify low-confidence predictions")

print("\nRecommendations:")
print("- Use RBF kernel for smooth patterns")
print("- Use Matérn kernel (nu=1.5) for more flexible patterns")
print("- Use Sparse GP for datasets with >500 spots")
print("- Always check uncertainty estimates for quality control")

print("\n" + "=" * 70)
print("Example completed successfully!")
print("=" * 70)

# Optional: Create visualization if matplotlib is available
try:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: True pattern
    sc1 = axes[0, 0].scatter(x, y, c=true_pattern, cmap='viridis', s=50)
    axes[0, 0].set_title('True Pattern')
    axes[0, 0].set_xlabel('X coordinate')
    axes[0, 0].set_ylabel('Y coordinate')
    plt.colorbar(sc1, ax=axes[0, 0])
    
    # Plot 2: Observed (sparse)
    sc2 = axes[0, 1].scatter(x[observed_mask], y[observed_mask], 
                             c=observed_values[observed_mask], cmap='viridis', s=50)
    axes[0, 1].set_title(f'Observed ({observed_mask.sum()} spots)')
    axes[0, 1].set_xlabel('X coordinate')
    axes[0, 1].set_ylabel('Y coordinate')
    plt.colorbar(sc2, ax=axes[0, 1])
    
    # Plot 3: GP Imputed
    sc3 = axes[0, 2].scatter(x, y, c=imputed_rbf, cmap='viridis', s=50)
    axes[0, 2].set_title('GP Imputed (RBF)')
    axes[0, 2].set_xlabel('X coordinate')
    axes[0, 2].set_ylabel('Y coordinate')
    plt.colorbar(sc3, ax=axes[0, 2])
    
    # Plot 4: Uncertainty
    sc4 = axes[1, 0].scatter(x, y, c=uncertainty_rbf, cmap='Reds', s=50)
    axes[1, 0].set_title('Uncertainty')
    axes[1, 0].set_xlabel('X coordinate')
    axes[1, 0].set_ylabel('Y coordinate')
    plt.colorbar(sc4, ax=axes[1, 0])
    
    # Plot 5: Error (imputed spots only)
    if imputed_mask.sum() > 0:
        errors = np.abs(imputed_rbf - true_pattern)
        sc5 = axes[1, 1].scatter(x[imputed_mask], y[imputed_mask], 
                                 c=errors[imputed_mask], cmap='Reds', s=50)
        axes[1, 1].set_title('Absolute Error (Imputed Spots)')
        axes[1, 1].set_xlabel('X coordinate')
        axes[1, 1].set_ylabel('Y coordinate')
        plt.colorbar(sc5, ax=axes[1, 1])
    
    # Plot 6: Sparse GP Inducing Points
    sc6 = axes[1, 2].scatter(x, y, c=imputed_sparse, cmap='viridis', s=50, alpha=0.5)
    axes[1, 2].scatter(sparse_imputer.inducing_points_[:, 0], 
                       sparse_imputer.inducing_points_[:, 1],
                       c='red', marker='x', s=100, label='Inducing points')
    axes[1, 2].set_title('Sparse GP (with Inducing Points)')
    axes[1, 2].set_xlabel('X coordinate')
    axes[1, 2].set_ylabel('Y coordinate')
    axes[1, 2].legend()
    plt.colorbar(sc6, ax=axes[1, 2])
    
    plt.tight_layout()
    plt.savefig('gp_imputation_example.png', dpi=150, bbox_inches='tight')
    print("\nVisualization saved to: gp_imputation_example.png")
    
except Exception as e:
    print(f"\nNote: Could not create visualization: {e}")

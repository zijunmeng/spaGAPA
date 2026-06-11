"""
Example 8: Complete spaGAPA Pipeline

This example demonstrates the end-to-end spaGAPA pipeline.
"""

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA')

from spagapa import SpaGAPA, APADataset
import anndata as ad

np.random.seed(42)

print("=" * 70)
print("Example 8: Complete spaGAPA Pipeline")
print("=" * 70)

# Create synthetic dataset
print("\n[Setup] Creating synthetic spatial APA dataset...")

# Spatial coordinates (10x10 grid)
x = np.linspace(0, 10, 10)
y = np.linspace(0, 10, 10)
xx, yy = np.meshgrid(x, y)
coordinates = np.column_stack([xx.ravel(), yy.ravel()])

# APA counts for 5 genes
n_spots = len(coordinates)
n_genes = 5

# Gene 1: Spatial gradient (SVAPA)
gene1 = xx.ravel() + 0.2 * np.random.randn(n_spots)
gene1 = (gene1 - gene1.min()) / (gene1.max() - gene1.min()) * 10

# Gene 2: Random (not SVAPA)
gene2 = np.random.rand(n_spots) * 10

# Gene 3: Hotspot (SVAPA)
dist = np.sqrt((xx.ravel() - 5)**2 + (yy.ravel() - 5)**2)
gene3 = np.exp(-dist / 2) * 10

# Gene 4: Ring pattern (SVAPA)
gene4 = np.exp(-((dist - 3)**2) / 2) * 10

# Gene 5: Low expression
gene5 = np.random.rand(n_spots) * 2

# Stack genes
apa_counts = np.vstack([gene1, gene2, gene3, gene4, gene5]).T

# Add some dropout (30%)
dropout_mask = np.random.rand(*apa_counts.shape) < 0.3
apa_counts[dropout_mask] = 0

# Create AnnData object
adata = ad.AnnData(
    X=apa_counts,
    obs=pd.DataFrame({
        'x': coordinates[:, 0],
        'y': coordinates[:, 1]
    }),
    var=pd.DataFrame({
        'gene_name': [f'Gene{i+1}' for i in range(n_genes)]
    })
)

# Create APADataset
dataset = APADataset(adata)
dataset.coordinates = coordinates

print(f"  ✓ Created dataset: {n_genes} genes, {n_spots} spots")
print(f"  ✓ Dropout rate: {(apa_counts == 0).mean():.1%}")

# Run complete pipeline
print("\n" + "=" * 70)
print("Running spaGAPA Pipeline")
print("=" * 70)

spa = SpaGAPA(
    n_neighbors=6,
    kernel_type='matern',
    use_sparse_gp=False,
    min_spatial_support=0.3,
    verbose=True
)

results = spa.run(
    dataset=dataset,
    impute=True,
    quantify=True,
    identify_domains=True,
    differential_analysis=False,
    detect_svapa=True,
    n_domains=3,
    fdr_threshold=0.1
)

# Display results
print("\n" + "=" * 70)
print("Results Summary")
print("=" * 70)

print("\n1. Imputation Results:")
imputed, uncertainty = spa.get_imputed_values()
print(f"   - Imputed shape: {imputed.shape}")
print(f"   - Mean uncertainty: {uncertainty.mean():.3f}")
print(f"   - Imputation coverage: {(imputed > 0).mean():.1%}")

print("\n2. SVAPA Detection:")
svapa_genes = spa.get_svapa_genes()
print(f"   - Total genes tested: {len(results['svapa_genes'])}")
print(f"   - SVAPA genes found: {len(svapa_genes)}")
if len(svapa_genes) > 0:
    print("\n   Significant SVAPA genes:")
    for _, row in svapa_genes.iterrows():
        print(f"     - {row['gene']}: p={row['lr_p_value']:.2e}, "
              f"q={row['q_value']:.2e}")

print("\n3. Spatial Domains:")
domains = spa.get_domains()
print(f"   - Number of domains: {len(np.unique(domains['labels']))}")
print(f"   - Domain sizes: {np.bincount(domains['labels'])}")

print("\n4. Quality Control:")
qc = results['qc_report']
print(f"   - QC report generated: {qc is not None}")

print("\n" + "=" * 70)
print("Pipeline completed successfully!")
print("=" * 70)

# Save results
print("\nSaving results...")
spa.save_results("pipeline_output")
print("  ✓ Results saved to 'pipeline_output/'")

print("\n✓ Example completed!")

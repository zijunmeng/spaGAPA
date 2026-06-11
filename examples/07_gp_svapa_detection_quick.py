"""
Quick test of GP-based SVAPA Detection
"""

import numpy as np
import sys
sys.path.insert(0, '/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA')

from spagapa.analysis import GPTrendDetector, detect_svapa_genes_gp

np.random.seed(42)

print("Testing GP-based SVAPA Detection...")

# Create small spatial grid
x = np.linspace(0, 10, 10)
y = np.linspace(0, 10, 10)
xx, yy = np.meshgrid(x, y)
coordinates = np.column_stack([xx.ravel(), yy.ravel()])

print(f"Created spatial grid: {len(coordinates)} spots")

# Gene 1: Spatial gradient (SVAPA)
values_gradient = xx.ravel() + 0.2 * np.random.randn(len(xx.ravel()))
values_gradient = (values_gradient - values_gradient.min()) / \
                  (values_gradient.max() - values_gradient.min())

# Gene 2: Random noise (not SVAPA)
values_random = np.random.rand(len(xx.ravel()))

# Uncertainty
uncertainty = 0.1 + 0.1 * np.random.rand(len(xx.ravel()))

print("\n1. Testing Likelihood Ratio Test...")
detector = GPTrendDetector(kernel_type='matern', n_restarts=2)

result1 = detector.likelihood_ratio_test(values_gradient, coordinates, uncertainty)
print(f"   Gene1 (gradient): LR={result1['lr_statistic']:.2f}, p={result1['p_value']:.4f}")

result2 = detector.likelihood_ratio_test(values_random, coordinates, uncertainty)
print(f"   Gene2 (random): LR={result2['lr_statistic']:.2f}, p={result2['p_value']:.4f}")

print("\n2. Testing Uncertainty-Weighted Moran's I...")
morans1 = detector.uncertainty_weighted_morans_i(values_gradient, coordinates, uncertainty)
print(f"   Gene1: I={morans1['morans_i']:.3f}, p={morans1['p_value']:.4f}")

morans2 = detector.uncertainty_weighted_morans_i(values_random, coordinates, uncertainty)
print(f"   Gene2: I={morans2['morans_i']:.3f}, p={morans2['p_value']:.4f}")

print("\n3. Testing Spatial Variance Decomposition...")
var1 = detector.spatial_variance_decomposition(values_gradient, coordinates)
print(f"   Gene1: Spatial fraction={var1['spatial_fraction']:.2%}, R²={var1['r_squared']:.3f}")

var2 = detector.spatial_variance_decomposition(values_random, coordinates)
print(f"   Gene2: Spatial fraction={var2['spatial_fraction']:.2%}, R²={var2['r_squared']:.3f}")

print("\n4. Testing Batch SVAPA Detection...")
apa_values = np.vstack([values_gradient, values_random])
uncertainty_matrix = np.vstack([uncertainty, uncertainty])
gene_names = ['Gene1_gradient', 'Gene2_random']

results = detect_svapa_genes_gp(
    apa_values,
    coordinates,
    uncertainty_matrix,
    gene_names,
    method='likelihood_ratio',
    kernel_type='matern',
    fdr_threshold=0.1
)

print("\nResults:")
print(results[['gene', 'lr_statistic', 'lr_p_value', 'q_value', 'significant']])

print("\n✓ All tests completed successfully!")
print(f"✓ Detected {results['significant'].sum()} SVAPA gene(s)")

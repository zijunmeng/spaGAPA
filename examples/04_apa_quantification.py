"""
Example 4: APA Quantification and Quality Control

This example demonstrates:
1. Calculating APA indices (RUD, PDUI, WUL, PAI)
2. Cross-validation for imputation quality
3. Quality control metrics
4. Generating QC reports

Author: spaGAPA Development Team
"""

import numpy as np
import matplotlib.pyplot as plt
from spagapa.quantification import (
    APAIndexCalculator,
    calculate_rud,
    calculate_pdui,
    calculate_wul,
    cross_validate_imputation,
    evaluate_imputation_quality,
    QCReportGenerator
)
from spagapa.imputation import GPImputer, SparseGPImputer


def generate_example_data(n_spots=100, n_genes=5, seed=42):
    """Generate synthetic APA data for demonstration."""
    np.random.seed(seed)
    
    # Generate spatial coordinates
    coordinates = np.random.rand(n_spots, 2) * 10
    
    # Generate proximal and distal site counts
    # Simulate spatial pattern: distal usage increases with x-coordinate
    x_coords = coordinates[:, 0]
    
    proximal_counts = np.zeros((n_genes, n_spots))
    distal_counts = np.zeros((n_genes, n_spots))
    
    for i in range(n_genes):
        # Base expression level
        base_expr = np.random.randint(20, 100)
        
        # Proximal usage decreases with x
        proximal_fraction = 0.8 - 0.6 * (x_coords / 10.0)
        proximal_fraction = np.clip(proximal_fraction, 0.1, 0.9)
        
        # Add noise
        proximal_fraction += np.random.normal(0, 0.1, n_spots)
        proximal_fraction = np.clip(proximal_fraction, 0.1, 0.9)
        
        # Generate counts
        total_counts = np.random.poisson(base_expr, n_spots)
        proximal_counts[i, :] = np.random.binomial(total_counts, proximal_fraction)
        distal_counts[i, :] = total_counts - proximal_counts[i, :]
    
    # Create sparse data (30% observed)
    observed_mask = np.random.rand(n_genes, n_spots) < 0.3
    proximal_sparse = proximal_counts.copy()
    distal_sparse = distal_counts.copy()
    proximal_sparse[~observed_mask] = 0
    distal_sparse[~observed_mask] = 0
    
    return {
        'coordinates': coordinates,
        'proximal_counts': proximal_counts,
        'distal_counts': distal_counts,
        'proximal_sparse': proximal_sparse,
        'distal_sparse': distal_sparse,
        'observed_mask': observed_mask
    }


def example1_basic_indices():
    """Example 1: Calculate basic APA indices."""
    print("=" * 60)
    print("Example 1: Basic APA Index Calculation")
    print("=" * 60)
    
    # Generate data
    data = generate_example_data(n_spots=50, n_genes=3)
    proximal = data['proximal_counts']
    distal = data['distal_counts']
    
    # Calculate individual indices
    print("\n1. Individual Index Calculations:")
    
    # RUD (Relative Usage of Distal site)
    rud = calculate_rud(proximal, distal)
    print(f"   RUD: mean={np.mean(rud):.3f}, std={np.std(rud):.3f}")
    print(f"        range=[{np.min(rud):.3f}, {np.max(rud):.3f}]")
    
    # PDUI (Percentage of Distal Usage Index)
    pdui = calculate_pdui(proximal, distal)
    print(f"   PDUI: mean={np.mean(pdui):.2f}%, std={np.std(pdui):.2f}%")
    print(f"         range=[{np.min(pdui):.2f}%, {np.max(pdui):.2f}%]")
    
    # WUL (Weighted 3' UTR Length)
    site_counts = {'proximal': proximal, 'distal': distal}
    site_positions = {'proximal': 1000, 'distal': 2500}
    wul = calculate_wul(site_counts, site_positions, normalize=True)
    print(f"   WUL: mean={np.mean(wul):.3f}, std={np.std(wul):.3f}")
    print(f"        range=[{np.min(wul):.3f}, {np.max(wul):.3f}]")
    
    # Use calculator for all indices
    print("\n2. Using APAIndexCalculator:")
    calculator = APAIndexCalculator(pseudocount=1.0, normalize=False)
    indices = calculator.calculate_all(
        proximal, distal,
        site_positions=site_positions
    )
    
    for index_name, values in indices.items():
        print(f"   {index_name}: shape={values.shape}, mean={np.mean(values):.3f}")
    
    return indices


def example2_spatial_pattern():
    """Example 2: Analyze spatial patterns in APA indices."""
    print("\n" + "=" * 60)
    print("Example 2: Spatial Patterns in APA Indices")
    print("=" * 60)
    
    # Generate data with spatial pattern
    data = generate_example_data(n_spots=100, n_genes=1)
    coordinates = data['coordinates']
    proximal = data['proximal_counts'][0, :]
    distal = data['distal_counts'][0, :]
    
    # Calculate RUD
    rud = calculate_rud(proximal, distal)
    
    # Analyze correlation with spatial coordinates
    from scipy.stats import pearsonr
    
    x_corr, x_pval = pearsonr(coordinates[:, 0], rud)
    y_corr, y_pval = pearsonr(coordinates[:, 1], rud)
    
    print(f"\nSpatial correlation with RUD:")
    print(f"  X-coordinate: r={x_corr:.3f}, p={x_pval:.4f}")
    print(f"  Y-coordinate: r={y_corr:.3f}, p={y_pval:.4f}")
    
    # Visualize
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Plot 1: Proximal counts
    sc1 = axes[0].scatter(
        coordinates[:, 0], coordinates[:, 1],
        c=proximal, cmap='Reds', s=50
    )
    axes[0].set_title('Proximal Site Counts')
    axes[0].set_xlabel('X coordinate')
    axes[0].set_ylabel('Y coordinate')
    plt.colorbar(sc1, ax=axes[0])
    
    # Plot 2: Distal counts
    sc2 = axes[1].scatter(
        coordinates[:, 0], coordinates[:, 1],
        c=distal, cmap='Blues', s=50
    )
    axes[1].set_title('Distal Site Counts')
    axes[1].set_xlabel('X coordinate')
    axes[1].set_ylabel('Y coordinate')
    plt.colorbar(sc2, ax=axes[1])
    
    # Plot 3: RUD
    sc3 = axes[2].scatter(
        coordinates[:, 0], coordinates[:, 1],
        c=rud, cmap='viridis', s=50, vmin=0, vmax=1
    )
    axes[2].set_title(f'RUD (r={x_corr:.3f} with X)')
    axes[2].set_xlabel('X coordinate')
    axes[2].set_ylabel('Y coordinate')
    plt.colorbar(sc3, ax=axes[2])
    
    plt.tight_layout()
    plt.savefig('apa_quantification_spatial.png', dpi=150, bbox_inches='tight')
    print("\nSaved spatial pattern visualization to 'apa_quantification_spatial.png'")
    
    return rud


def example3_imputation_qc():
    """Example 3: Quality control for imputation."""
    print("\n" + "=" * 60)
    print("Example 3: Imputation Quality Control")
    print("=" * 60)
    
    # Generate sparse data
    data = generate_example_data(n_spots=80, n_genes=1)
    coordinates = data['coordinates']
    true_values = data['distal_counts'][0, :]
    sparse_values = data['distal_sparse'][0, :]
    observed_mask = data['observed_mask'][0, :]
    
    print(f"\nData statistics:")
    print(f"  Total spots: {len(true_values)}")
    print(f"  Observed spots: {np.sum(observed_mask)} ({100*np.mean(observed_mask):.1f}%)")
    print(f"  Missing spots: {np.sum(~observed_mask)} ({100*np.mean(~observed_mask):.1f}%)")
    
    # Perform imputation
    print("\n1. Performing GP imputation...")
    imputer = GPImputer(kernel_type='matern', nu=1.5)
    imputed_values, uncertainty = imputer.impute(
        coordinates, sparse_values, observed_mask
    )
    
    # Evaluate imputation quality
    print("\n2. Evaluating imputation quality...")
    from spagapa.quantification import evaluate_imputation_quality
    
    metrics = evaluate_imputation_quality(
        true_values, imputed_values, uncertainty, ~observed_mask
    )
    
    print("\nImputation Quality Metrics:")
    print(f"  RMSE: {metrics['rmse']:.3f}")
    print(f"  MAE: {metrics['mae']:.3f}")
    print(f"  Pearson r: {metrics['pearson_r']:.3f} (p={metrics['pearson_pval']:.4f})")
    print(f"  Spearman r: {metrics['spearman_r']:.3f} (p={metrics['spearman_pval']:.4f})")
    print(f"  R²: {metrics['r2']:.3f}")
    print(f"  Mean uncertainty: {metrics['mean_uncertainty']:.3f}")
    print(f"  Uncertainty calibration: {metrics['uncertainty_calibration']:.3f}")
    
    # Cross-validation
    print("\n3. Performing cross-validation...")
    cv_metrics = cross_validate_imputation(
        imputer, coordinates, true_values,
        n_folds=5, mask=observed_mask
    )
    
    print("\nCross-Validation Results:")
    print(f"  RMSE: {cv_metrics['rmse']:.3f} ± {cv_metrics['rmse_std']:.3f}")
    print(f"  MAE: {cv_metrics['mae']:.3f} ± {cv_metrics['mae_std']:.3f}")
    print(f"  Pearson r: {cv_metrics['pearson']:.3f} ± {cv_metrics['pearson_std']:.3f}")
    print(f"  R²: {cv_metrics['r2']:.3f} ± {cv_metrics['r2_std']:.3f}")
    
    return metrics, cv_metrics


def example4_qc_report():
    """Example 4: Generate comprehensive QC report."""
    print("\n" + "=" * 60)
    print("Example 4: Comprehensive QC Report")
    print("=" * 60)
    
    # Generate data
    data = generate_example_data(n_spots=100, n_genes=5)
    
    # Calculate APA indices
    calculator = APAIndexCalculator(pseudocount=1.0)
    site_positions = {'proximal': 1000, 'distal': 2500}
    indices = calculator.calculate_all(
        data['proximal_counts'],
        data['distal_counts'],
        site_positions=site_positions
    )
    
    # Create QC report
    qc = QCReportGenerator(dataset_name="Example APA Dataset")
    
    # Add quantification metrics
    qc.add_quantification_metrics(indices)
    
    # Add coverage metrics
    qc.add_coverage_metrics(
        total_spots=100,
        observed_spots=30,
        imputed_spots=70
    )
    
    # Add imputation metrics (from previous example)
    qc.add_imputation_metrics({
        'rmse': 2.5,
        'mae': 1.8,
        'pearson_r': 0.92,
        'r2': 0.85,
        'mean_uncertainty': 0.65
    })
    
    # Generate report
    print("\n1. Report as dictionary:")
    report_dict = qc.generate_report(format='dict')
    print(f"   Sections: {list(report_dict.keys())}")
    
    print("\n2. Report as text:")
    report_text = qc.generate_report(format='text')
    print(report_text)
    
    # Save report
    qc.save_report('apa_qc_report.txt', format='text')
    print("\nSaved QC report to 'apa_qc_report.txt'")
    
    return qc


def example5_compare_methods():
    """Example 5: Compare imputation methods."""
    print("\n" + "=" * 60)
    print("Example 5: Compare Imputation Methods")
    print("=" * 60)
    
    # Generate data
    data = generate_example_data(n_spots=60, n_genes=1)
    coordinates = data['coordinates']
    values = data['distal_counts'][0, :]
    observed_mask = data['observed_mask'][0, :]
    
    # Define methods to compare
    methods = {
        'GP-RBF': GPImputer(kernel_type='rbf'),
        'GP-Matern': GPImputer(kernel_type='matern', nu=1.5),
        'Sparse-GP': SparseGPImputer(n_inducing=15, inducing_method='kmeans')
    }
    
    print(f"\nComparing {len(methods)} imputation methods...")
    print(f"Data: {len(values)} spots, {np.sum(observed_mask)} observed")
    
    # Compare methods
    from spagapa.quantification import compare_imputation_methods
    
    results = compare_imputation_methods(
        methods, coordinates, values,
        mask=observed_mask, n_folds=3
    )
    
    print("\nComparison Results:")
    print(results.to_string())
    
    # Find best method
    best_method = results['rmse'].idxmin()
    print(f"\nBest method (lowest RMSE): {best_method}")
    print(f"  RMSE: {results.loc[best_method, 'rmse']:.3f}")
    print(f"  MAE: {results.loc[best_method, 'mae']:.3f}")
    print(f"  Pearson r: {results.loc[best_method, 'pearson']:.3f}")
    
    return results


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("APA Quantification and Quality Control Examples")
    print("=" * 60)
    
    # Run examples
    indices = example1_basic_indices()
    rud = example2_spatial_pattern()
    metrics, cv_metrics = example3_imputation_qc()
    qc = example4_qc_report()
    results = example5_compare_methods()
    
    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  - apa_quantification_spatial.png")
    print("  - apa_qc_report.txt")
    print("\nKey findings:")
    print(f"  - RUD values range from {np.min(rud):.3f} to {np.max(rud):.3f}")
    print(f"  - Imputation RMSE: {metrics['rmse']:.3f}")
    print(f"  - Cross-validation R²: {cv_metrics['r2']:.3f}")
    print(f"  - Best imputation method: {results['rmse'].idxmin()}")


if __name__ == "__main__":
    main()

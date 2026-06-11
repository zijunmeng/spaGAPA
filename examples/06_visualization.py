"""
Example 6: Visualization

This example demonstrates all visualization capabilities of spaGAPA.
"""

import numpy as np
import matplotlib.pyplot as plt

from spagapa.visualization import (
    plot_spatial_apa,
    plot_spatial_domains,
    plot_volcano,
    plot_heatmap,
    plot_imputation_quality,
    SpatialPlotter,
    StatisticalPlotter,
    QCPlotter
)


def generate_example_data():
    """Generate synthetic data for visualization."""
    np.random.seed(42)
    
    # Spatial coordinates
    n_spots = 200
    coords = np.random.rand(n_spots, 2) * 100
    
    # APA matrix
    n_genes = 50
    apa_matrix = np.random.rand(n_genes, n_spots) * 0.5 + 0.25
    
    # Add spatial pattern to some genes
    for i in range(5):
        apa_matrix[i, :] += coords[:, 0] / 200
    
    # Domain labels
    domains = np.repeat([0, 1, 2, 3], 50)
    
    # Differential analysis results
    logfc = np.random.randn(n_genes)
    pvalues = np.random.rand(n_genes)
    
    # Imputation data
    observed = np.random.rand(100)
    imputed = observed + np.random.randn(100) * 0.1
    uncertainties = np.abs(np.random.randn(100) * 0.1)
    
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    
    return {
        'coords': coords,
        'apa_matrix': apa_matrix,
        'domains': domains,
        'logfc': logfc,
        'pvalues': pvalues,
        'observed': observed,
        'imputed': imputed,
        'uncertainties': uncertainties,
        'gene_names': gene_names
    }


def example1_spatial_plots(data):
    """Example 1: Spatial visualization."""
    print("=" * 70)
    print("Example 1: Spatial Plots")
    print("=" * 70)
    
    # Single gene spatial plot
    print("\n1. Plotting single gene APA distribution...")
    ax = plot_spatial_apa(
        data['coords'],
        data['apa_matrix'][0, :],
        gene_name="GENE000",
        cmap='viridis',
        save='spatial_apa_single.png'
    )
    plt.close()
    print("   Saved: spatial_apa_single.png")
    
    # Domain visualization
    print("\n2. Plotting spatial domains...")
    domain_names = {0: 'Domain A', 1: 'Domain B', 2: 'Domain C', 3: 'Domain D'}
    ax = plot_spatial_domains(
        data['coords'],
        data['domains'],
        domain_names=domain_names,
        show_boundaries=False,
        save='spatial_domains.png'
    )
    plt.close()
    print("   Saved: spatial_domains.png")
    
    # Multi-gene comparison
    print("\n3. Plotting multi-gene comparison...")
    plotter = SpatialPlotter()
    fig = plotter.plot_apa_comparison(
        data['coords'],
        data['apa_matrix'],
        genes=[0, 1, 2, 3, 4, 5],
        gene_names=data['gene_names'],
        ncols=3,
        save='spatial_comparison.png'
    )
    plt.close()
    print("   Saved: spatial_comparison.png")


def example2_statistical_plots(data):
    """Example 2: Statistical visualization."""
    print("\n" + "=" * 70)
    print("Example 2: Statistical Plots")
    print("=" * 70)
    
    # Volcano plot
    print("\n1. Creating volcano plot...")
    ax = plot_volcano(
        data['logfc'],
        data['pvalues'],
        gene_names=data['gene_names'],
        threshold_fc=0.5,
        threshold_p=0.05,
        label_top=5,
        save='volcano_plot.png'
    )
    plt.close()
    print("   Saved: volcano_plot.png")
    
    # Heatmap
    print("\n2. Creating clustered heatmap...")
    # Use subset for better visualization
    subset_matrix = data['apa_matrix'][:20, :30]
    subset_genes = data['gene_names'][:20]
    subset_spots = [f"Spot{i}" for i in range(30)]
    
    fig, ax = plot_heatmap(
        subset_matrix,
        row_labels=subset_genes,
        col_labels=subset_spots,
        cluster_rows=True,
        cluster_cols=True,
        cmap='RdBu_r',
        save='heatmap.png'
    )
    plt.close()
    print("   Saved: heatmap.png")
    
    # Box plot
    print("\n3. Creating box plot...")
    plotter = StatisticalPlotter()
    
    # APA values by domain
    apa_values = data['apa_matrix'][0, :]
    ax = plotter.plot_boxplot(
        apa_values,
        data['domains'],
        gene_name="GENE000",
        show_points=True,
        save='boxplot.png'
    )
    plt.close()
    print("   Saved: boxplot.png")
    
    # Violin plot
    print("\n4. Creating violin plot...")
    ax = plotter.plot_violin(
        apa_values,
        data['domains'],
        gene_name="GENE000",
        save='violin_plot.png'
    )
    plt.close()
    print("   Saved: violin_plot.png")


def example3_qc_plots(data):
    """Example 3: Quality control visualization."""
    print("\n" + "=" * 70)
    print("Example 3: Quality Control Plots")
    print("=" * 70)
    
    # Imputation quality
    print("\n1. Plotting imputation quality...")
    fig = plot_imputation_quality(
        data['observed'],
        data['imputed'],
        data['uncertainties'],
        save='qc_imputation.png'
    )
    plt.close()
    print("   Saved: qc_imputation.png")
    
    # Spatial support
    print("\n2. Plotting spatial support...")
    plotter = QCPlotter()
    support_scores = np.random.rand(len(data['coords']))
    ax = plotter.plot_spatial_support(
        data['coords'],
        support_scores,
        save='qc_spatial_support.png'
    )
    plt.close()
    print("   Saved: qc_spatial_support.png")
    
    # Dropout statistics
    print("\n3. Plotting dropout statistics...")
    # Add some NaN values
    apa_with_dropout = data['apa_matrix'].copy()
    apa_with_dropout[np.random.rand(*apa_with_dropout.shape) < 0.3] = np.nan
    
    fig = plotter.plot_dropout_stats(
        apa_with_dropout,
        data['coords'],
        save='qc_dropout.png'
    )
    plt.close()
    print("   Saved: qc_dropout.png")


def example4_custom_styling():
    """Example 4: Custom styling and advanced features."""
    print("\n" + "=" * 70)
    print("Example 4: Custom Styling")
    print("=" * 70)
    
    np.random.seed(42)
    coords = np.random.rand(100, 2) * 100
    apa_values = np.random.rand(100)
    
    # Custom colormap and styling
    print("\n1. Creating custom styled plot...")
    plotter = SpatialPlotter(figsize=(10, 8), dpi=150)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Different colormaps
    cmaps = ['viridis', 'plasma', 'RdYlBu_r', 'coolwarm']
    titles = ['Viridis', 'Plasma', 'RdYlBu_r', 'Coolwarm']
    
    for ax, cmap, title in zip(axes.flat, cmaps, titles):
        plotter.plot_spatial_apa(
            coords,
            apa_values,
            gene_name=title,
            cmap=cmap,
            size=30,
            alpha=0.8,
            ax=ax
        )
    
    plt.tight_layout()
    plt.savefig('custom_styling.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("   Saved: custom_styling.png")


def main():
    """Run all visualization examples."""
    print("\n" + "=" * 70)
    print("spaGAPA Visualization Examples")
    print("=" * 70)
    
    # Generate data
    print("\nGenerating example data...")
    data = generate_example_data()
    print(f"  Spots: {len(data['coords'])}")
    print(f"  Genes: {len(data['gene_names'])}")
    print(f"  Domains: {len(np.unique(data['domains']))}")
    
    # Run examples
    example1_spatial_plots(data)
    example2_statistical_plots(data)
    example3_qc_plots(data)
    example4_custom_styling()
    
    print("\n" + "=" * 70)
    print("All visualization examples completed!")
    print("Generated files:")
    print("  - spatial_apa_single.png")
    print("  - spatial_domains.png")
    print("  - spatial_comparison.png")
    print("  - volcano_plot.png")
    print("  - heatmap.png")
    print("  - boxplot.png")
    print("  - violin_plot.png")
    print("  - qc_imputation.png")
    print("  - qc_spatial_support.png")
    print("  - qc_dropout.png")
    print("  - custom_styling.png")
    print("=" * 70)


if __name__ == '__main__':
    main()

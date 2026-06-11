"""
Example 5: Differential APA Analysis

This example demonstrates:
1. Spatial domain identification
2. Differential APA testing between domains
3. Finding domain-specific marker genes
4. Identifying spatially variable APA genes (SVAPA)
5. Clustering genes by spatial patterns
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import spaGAPA analysis functions
from spagapa.analysis import (
    identify_spatial_domains,
    test_differential_apa,
    find_domain_markers,
    identify_svapa_genes,
    cluster_spatial_patterns
)


def generate_example_data(n_genes=50, n_spots=200, n_domains=4):
    """Generate synthetic spatial APA data with domains."""
    np.random.seed(42)
    
    # Generate spatial coordinates (grid-like)
    grid_size = int(np.sqrt(n_spots))
    x = np.repeat(np.arange(grid_size), grid_size)[:n_spots]
    y = np.tile(np.arange(grid_size), grid_size)[:n_spots]
    coords = np.column_stack([x, y])
    
    # Create domains (spatial regions)
    true_domains = np.zeros(n_spots, dtype=int)
    for i in range(n_domains):
        start = i * (n_spots // n_domains)
        end = (i + 1) * (n_spots // n_domains)
        true_domains[start:end] = i
    
    # Generate APA matrix with domain-specific patterns
    apa_matrix = np.random.rand(n_genes, n_spots) * 0.5 + 0.25
    
    # Add domain-specific APA patterns
    for domain in range(n_domains):
        domain_mask = true_domains == domain
        # Make some genes domain-specific
        domain_genes = range(domain * 10, (domain + 1) * 10)
        for gene in domain_genes:
            if gene < n_genes:
                apa_matrix[gene, domain_mask] += 0.3
    
    # Add spatial gradients to some genes (SVAPA genes)
    for i in range(5):
        apa_matrix[i, :] += coords[:, 0] / (2 * grid_size)
    
    # Generate gene names
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    
    return apa_matrix, coords, gene_names, true_domains


def example1_identify_domains():
    """Example 1: Identify spatial domains."""
    print("=" * 70)
    print("Example 1: Spatial Domain Identification")
    print("=" * 70)
    
    # Generate data
    apa_matrix, coords, gene_names, true_domains = generate_example_data()
    print(f"Data: {apa_matrix.shape[0]} genes × {apa_matrix.shape[1]} spots")
    
    # Identify domains using K-means
    print("\n1. Identifying domains using K-means...")
    labels, stats = identify_spatial_domains(
        apa_matrix,
        coords,
        method='kmeans',
        n_clusters=4,
        refine=True,
        min_domain_size=10
    )
    
    print(f"\nFound {len(stats)} domains:")
    print(stats)
    
    # Compare with true domains
    from sklearn.metrics import adjusted_rand_score
    ari = adjusted_rand_score(true_domains, labels)
    print(f"\nAdjusted Rand Index (vs true domains): {ari:.3f}")
    
    return apa_matrix, coords, gene_names, labels


def example2_differential_apa(apa_matrix, coords, gene_names, labels):
    """Example 2: Test for differential APA between domains."""
    print("\n" + "=" * 70)
    print("Example 2: Differential APA Testing")
    print("=" * 70)
    
    # Compare domain 0 vs domain 1
    domain0_spots = np.where(labels == 0)[0]
    domain1_spots = np.where(labels == 1)[0]
    
    print(f"\nComparing Domain 0 ({len(domain0_spots)} spots) vs "
          f"Domain 1 ({len(domain1_spots)} spots)")
    
    # Test differential APA
    print("\n1. Running Wilcoxon test...")
    results = test_differential_apa(
        apa_matrix,
        domain0_spots,
        domain1_spots,
        gene_names=gene_names,
        method='wilcoxon'
    )
    
    # Filter significant genes
    sig_genes = results[results['padj'] < 0.05]
    print(f"\nFound {len(sig_genes)} significant genes (FDR < 0.05)")
    
    # Show top differential genes
    print("\nTop 10 differential genes:")
    top_genes = sig_genes.nsmallest(10, 'padj')
    print(top_genes[['gene', 'mean_group1', 'mean_group2', 'log2fc', 'padj']])
    
    return results


def example3_find_markers(apa_matrix, coords, gene_names, labels):
    """Example 3: Find domain-specific marker genes."""
    print("\n" + "=" * 70)
    print("Example 3: Domain Marker Identification")
    print("=" * 70)
    
    print("\n1. Finding markers for each domain (one-vs-rest)...")
    markers = find_domain_markers(
        apa_matrix,
        labels,
        gene_names=gene_names,
        padj_threshold=0.05,
        logfc_threshold=0.5
    )
    
    print(f"\nFound markers for {len(markers)} domains:")
    for domain_id, marker_df in markers.items():
        print(f"\nDomain {domain_id}: {len(marker_df)} marker genes")
        if len(marker_df) > 0:
            print(f"  Top 5 markers:")
            top5 = marker_df.head(5)
            for _, row in top5.iterrows():
                print(f"    {row['gene']}: log2FC={row['log2fc']:.2f}, "
                      f"padj={row['padj']:.2e}")
    
    return markers


def example4_svapa_genes(apa_matrix, coords, gene_names):
    """Example 4: Identify spatially variable APA genes."""
    print("\n" + "=" * 70)
    print("Example 4: Spatially Variable APA Genes (SVAPA)")
    print("=" * 70)
    
    print("\n1. Computing Moran's I for spatial autocorrelation...")
    svapa_genes, results = identify_svapa_genes(
        apa_matrix,
        coords,
        gene_names=gene_names,
        fdr_threshold=0.1,
        n_neighbors=6
    )
    
    print(f"\nFound {len(svapa_genes)} SVAPA genes (FDR < 0.1)")
    
    # Show top SVAPA genes
    print("\nTop 10 SVAPA genes by Moran's I:")
    top_svapa = results.nlargest(10, 'morans_i')
    print(top_svapa[['gene', 'morans_i', 'z_score', 'pvalue']])
    
    return svapa_genes, results


def example5_pattern_clustering(apa_matrix, coords, gene_names):
    """Example 5: Cluster genes by spatial patterns."""
    print("\n" + "=" * 70)
    print("Example 5: Spatial Pattern Clustering")
    print("=" * 70)
    
    print("\n1. Clustering genes by spatial expression patterns...")
    labels, info, genes = cluster_spatial_patterns(
        apa_matrix,
        coords,
        n_patterns=5,
        gene_names=gene_names
    )
    
    print(f"\nIdentified {len(info)} spatial patterns:")
    print(info)
    
    # Show genes in each pattern
    print("\nGenes per pattern:")
    for pattern_id, pattern_genes in genes.items():
        print(f"  Pattern {pattern_id}: {len(pattern_genes)} genes")
        if len(pattern_genes) <= 5:
            print(f"    {', '.join(pattern_genes)}")
        else:
            print(f"    {', '.join(pattern_genes[:5])}, ...")
    
    return labels, info, genes


def visualize_results(apa_matrix, coords, gene_names, domain_labels, 
                     svapa_genes, pattern_labels):
    """Visualize analysis results."""
    print("\n" + "=" * 70)
    print("Visualization")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Spatial domains
    ax = axes[0, 0]
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=domain_labels, 
                        cmap='tab10', s=50, alpha=0.7)
    ax.set_title('Spatial Domains')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    plt.colorbar(scatter, ax=ax, label='Domain')
    
    # 2. Example gene APA distribution
    ax = axes[0, 1]
    gene_idx = 0
    scatter = ax.scatter(coords[:, 0], coords[:, 1], 
                        c=apa_matrix[gene_idx, :],
                        cmap='viridis', s=50, alpha=0.7)
    ax.set_title(f'{gene_names[gene_idx]} APA Usage')
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    plt.colorbar(scatter, ax=ax, label='APA Index')
    
    # 3. SVAPA gene example
    ax = axes[0, 2]
    if len(svapa_genes) > 0:
        svapa_idx = gene_names.index(svapa_genes[0])
        scatter = ax.scatter(coords[:, 0], coords[:, 1],
                           c=apa_matrix[svapa_idx, :],
                           cmap='RdYlBu_r', s=50, alpha=0.7)
        ax.set_title(f'SVAPA Gene: {svapa_genes[0]}')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        plt.colorbar(scatter, ax=ax, label='APA Index')
    else:
        ax.text(0.5, 0.5, 'No SVAPA genes found', 
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('SVAPA Gene')
    
    # 4. Domain sizes
    ax = axes[1, 0]
    unique_domains, counts = np.unique(domain_labels, return_counts=True)
    ax.bar(unique_domains, counts, color='steelblue', alpha=0.7)
    ax.set_xlabel('Domain')
    ax.set_ylabel('Number of spots')
    ax.set_title('Domain Sizes')
    ax.set_xticks(unique_domains)
    
    # 5. APA distribution by domain
    ax = axes[1, 1]
    domain_means = []
    for domain in unique_domains:
        domain_mask = domain_labels == domain
        domain_means.append(np.mean(apa_matrix[:, domain_mask]))
    ax.bar(unique_domains, domain_means, color='coral', alpha=0.7)
    ax.set_xlabel('Domain')
    ax.set_ylabel('Mean APA Index')
    ax.set_title('Mean APA by Domain')
    ax.set_xticks(unique_domains)
    
    # 6. Pattern distribution
    ax = axes[1, 2]
    unique_patterns, counts = np.unique(pattern_labels, return_counts=True)
    ax.bar(unique_patterns, counts, color='mediumseagreen', alpha=0.7)
    ax.set_xlabel('Pattern')
    ax.set_ylabel('Number of genes')
    ax.set_title('Genes per Spatial Pattern')
    ax.set_xticks(unique_patterns)
    
    plt.tight_layout()
    plt.savefig('differential_analysis_results.png', dpi=300, bbox_inches='tight')
    print("\nSaved visualization to: differential_analysis_results.png")
    plt.close()


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("spaGAPA Differential APA Analysis Examples")
    print("=" * 70)
    
    # Example 1: Identify domains
    apa_matrix, coords, gene_names, labels = example1_identify_domains()
    
    # Example 2: Differential APA testing
    diff_results = example2_differential_apa(apa_matrix, coords, gene_names, labels)
    
    # Example 3: Find domain markers
    markers = example3_find_markers(apa_matrix, coords, gene_names, labels)
    
    # Example 4: Identify SVAPA genes
    svapa_genes, svapa_results = example4_svapa_genes(apa_matrix, coords, gene_names)
    
    # Example 5: Pattern clustering
    pattern_labels, pattern_info, pattern_genes = example5_pattern_clustering(
        apa_matrix, coords, gene_names
    )
    
    # Visualize results
    visualize_results(apa_matrix, coords, gene_names, labels, 
                     svapa_genes, pattern_labels)
    
    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)


if __name__ == '__main__':
    main()

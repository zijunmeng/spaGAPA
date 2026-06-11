"""
Example: Using scAPAtrap wrapper to identify APA sites from spatial transcriptomics data.

This example demonstrates how to:
1. Run scAPAtrap on a BAM file
2. Parse the results
3. Apply basic spatial filtering
"""

import pandas as pd
from spagapa.io import run_scapatrap_pipeline

# Example usage
if __name__ == "__main__":
    # Input files
    bam_file = "data/sample_spatial.bam"  # Your spatial transcriptomics BAM file
    genome_fasta = "data/genome.fa"        # Reference genome
    gtf_file = "data/genes.gtf"            # Gene annotation
    output_dir = "results/scapatrap"       # Output directory
    
    # Optional: Load spatial coordinates for filtering
    spatial_coords = pd.read_csv("data/spatial_coordinates.csv", index_col=0)
    # Expected format:
    # barcode,x,y
    # SPOT_1,100.5,200.3
    # SPOT_2,105.2,198.7
    # ...
    
    # Run scAPAtrap
    print("Running scAPAtrap...")
    results = run_scapatrap_pipeline(
        bam_file=bam_file,
        output_dir=output_dir,
        genome_fasta=genome_fasta,
        gtf_file=gtf_file,
        spatial_coords=spatial_coords,
        tails_search="genome",  # Search for poly(A) tails
        n_cores=4               # Use 4 CPU cores
    )
    
    # Examine results
    if 'peaks_meta' in results:
        peaks_meta = results['peaks_meta']
        print(f"\nIdentified {len(peaks_meta)} APA peaks")
        print("\nPeak metadata columns:")
        print(peaks_meta.columns.tolist())
        print("\nFirst few peaks:")
        print(peaks_meta.head())
    
    if 'peaks_counts' in results:
        peaks_counts = results['peaks_counts']
        print(f"\nCount matrix shape: {peaks_counts.shape}")
        print(f"  - {peaks_counts.shape[0]} peaks")
        print(f"  - {peaks_counts.shape[1]} spots")
        
        # Basic statistics
        print(f"\nBasic statistics:")
        print(f"  - Total counts: {peaks_counts.sum().sum()}")
        print(f"  - Mean counts per peak: {peaks_counts.sum(axis=1).mean():.2f}")
        print(f"  - Sparsity: {(peaks_counts == 0).sum().sum() / peaks_counts.size * 100:.1f}%")
    
    print("\nscAPAtrap analysis completed!")
    print(f"Results saved to: {output_dir}")

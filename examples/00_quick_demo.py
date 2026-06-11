"""
Quick demo of spaGAPA core functionality.

This script demonstrates the basic usage of spaGAPA's core data structures
and I/O functions.
"""

import numpy as np
import pandas as pd
import spagapa
from spagapa.core import APADataset, APASite, APASiteCollection

print("=" * 60)
print("spaGAPA Quick Demo")
print("=" * 60)
print(f"Version: {spagapa.__version__}\n")

# 1. Create APADataset
print("1. Creating APADataset...")
n_genes, n_spots = 10, 20
apa_counts = np.random.poisson(5, (n_genes, n_spots))
spatial_coords = np.random.rand(n_spots, 2) * 100
gene_names = [f"Gene_{i}" for i in range(n_genes)]
spot_names = [f"Spot_{i}" for i in range(n_spots)]

dataset = APADataset.from_counts(
    apa_counts=apa_counts,
    spatial_coords=spatial_coords,
    gene_names=gene_names,
    spot_names=spot_names
)

print(f"   Created dataset: {dataset.n_genes} genes, {dataset.n_spots} spots")
print(f"   Spatial coords shape: {dataset.get_spatial_coords().shape}")

# 2. Add imputation results
print("\n2. Adding imputation results...")
imputed = apa_counts + np.random.randn(n_genes, n_spots) * 0.1
uncertainty = np.random.rand(n_genes, n_spots) * 0.5
dataset.add_imputation(imputed, uncertainty)
print("   Imputation added successfully")

# 3. Create APA sites
print("\n3. Creating APA sites...")
sites = [
    APASite(
        chr="chr1",
        start=1000 + i*1000,
        end=1050 + i*1000,
        strand="+",
        gene_name=f"Gene_{i}",
        support_score=np.random.rand(),
        read_count=np.random.randint(5, 50)
    )
    for i in range(5)
]

collection = APASiteCollection(sites)
print(f"   Created collection: {len(collection)} sites")

# 4. Filter sites
print("\n4. Filtering sites...")
filtered = collection.filter_by_support(min_score=0.3)
print(f"   After filtering: {len(filtered)} sites (support >= 0.3)")

# 5. Group by gene
print("\n5. Grouping sites by gene...")
groups = collection.group_by_gene()
print(f"   Found {len(groups)} genes")

# 6. Subset dataset
print("\n6. Subsetting dataset...")
subset = dataset.subset_genes(gene_names[:5])
print(f"   Subset: {subset.n_genes} genes, {subset.n_spots} spots")

# 7. Display dataset info
print("\n7. Dataset summary:")
print(dataset)

print("\n" + "=" * 60)
print("Demo completed successfully! ✅")
print("=" * 60)

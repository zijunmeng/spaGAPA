# spaGAPA User Guide

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Data Input](#data-input)
4. [APA Calling](#apa-calling)
5. [GP Imputation](#gp-imputation)
6. [APA Quantification](#apa-quantification)
7. [Differential Analysis](#differential-analysis)
8. [Visualization](#visualization)
9. [Benchmarking](#benchmarking)
10. [Complete Pipeline](#complete-pipeline)
11. [FAQ](#faq)

---

## Installation

```bash
pip install -e .          # from source
pip install -e ".[dev]"   # with dev dependencies
```

Conda environment (recommended):

```bash
conda create -n spagapa python=3.10
conda activate spagapa
pip install -e ".[dev]"
```

---

## Quick Start

```python
from spagapa.benchmark import simulate_spatial_apa
from spagapa.imputation import GPImputer
from spagapa.analysis import identify_spatial_domains, test_differential_apa
from spagapa.visualization import plot_spatial_apa

# Simulate data
data = simulate_spatial_apa(n_spots=300, n_genes=50, dropout_rate=0.4)
coords   = data['coordinates']
observed = data['observed_apa'].T   # (n_genes, n_spots)

# Impute one gene
imputer = GPImputer(kernel_type='matern')
imputed, uncertainty = imputer.impute(observed[0], coords)

# Identify domains
labels, stats = identify_spatial_domains(observed, coords,
                                         method='kmeans', n_clusters=4)

# Differential APA
import numpy as np
g1 = np.where(labels == 0)[0]
g2 = np.where(labels == 1)[0]
results = test_differential_apa(observed, g1, g2)

# Plot
ax = plot_spatial_apa(coords, imputed, gene_name='Gene_0')
```

---

## Data Input

### From CSV files

```python
from spagapa.io import SpatialCoordinateReader
import pandas as pd
import numpy as np

# APA matrix: rows = spots, columns = genes
apa_df = pd.read_csv('apa_matrix.csv', index_col=0)
apa_matrix = apa_df.values.T   # (n_genes, n_spots)

# Coordinates
reader = SpatialCoordinateReader()
coords = reader.read_coordinates('coordinates.csv')   # (n_spots, 2)
```

### From 10x Visium

```python
from spagapa.io import SpatialCoordinateReader
coords = reader.read_visium('spatial/tissue_positions_list.csv')
```

### From AnnData / H5AD

```python
import anndata as ad
adata = ad.read_h5ad('data.h5ad')
apa_matrix = adata.X.T          # (n_genes, n_spots)
coords     = adata.obsm['spatial']
```

---

## APA Calling

spaGAPA uses scAPAtrap for initial APA site detection, then applies spatial validation.

```python
from spagapa.io import ScAPAtrapWrapper
from spagapa.calling import SpatialValidator, QualityFilter

# Step 1: Run scAPAtrap (requires R)
wrapper = ScAPAtrapWrapper()
sites = wrapper.run(bam_file='sample.bam', gtf_file='genes.gtf')

# Step 2: Spatial validation
validator = SpatialValidator(n_neighbors=6, support_threshold=0.3)
validated = validator.validate(sites, coords)

# Step 3: Quality filtering
qfilter = QualityFilter(min_reads=5, min_spots=3)
filtered = qfilter.filter(validated)
```

---

## GP Imputation

### Single gene

```python
from spagapa.imputation import GPImputer
import numpy as np

imputer = GPImputer(kernel_type='matern')   # 'rbf', 'matern', 'auto'

gene_apa = apa_matrix[0, :]                 # (n_spots,) with NaN
imputed, uncertainty = imputer.impute(gene_apa, coords)

print(f"Imputed {np.isnan(gene_apa).sum()} missing values")
print(f"Mean uncertainty: {uncertainty.mean():.3f}")
```

### Batch (all genes)

```python
from spagapa.imputation import GPImputer, GPImputerBatch

base = GPImputer(kernel_type='matern')
batch = GPImputerBatch(base_imputer=base, coordinates=coords,
                       values=apa_matrix, n_jobs=4)
imputed_matrix, uncertainty_matrix = batch.impute()
# imputed_matrix: (n_genes, n_spots)
```

### Large datasets (sparse GP)

```python
from spagapa.imputation import SparseGPImputer

imputer = SparseGPImputer(n_inducing=50, method='kmeans')
imputed, uncertainty = imputer.impute(gene_apa, coords)
```

---

## APA Quantification

```python
from spagapa.quantification import APAIndexCalculator

calc = APAIndexCalculator()

# RUD: Relative Usage of Distal site
rud = calc.calculate_rud(proximal_counts, distal_counts)

# PDUI: Percentage of Distal Usage Index
pdui = calc.calculate_pdui(proximal_counts, distal_counts)

# WUL: Weighted 3' UTR Length
wul = calc.calculate_wul(site_counts, site_positions)

# All indices at once
indices = calc.compute_all(proximal_counts, distal_counts, site_positions)
```

---

## Differential Analysis

### Identify spatial domains

```python
from spagapa.analysis import identify_spatial_domains

labels, stats = identify_spatial_domains(
    apa_matrix, coords,
    method='kmeans',    # or 'leiden', 'louvain'
    n_clusters=5,
    refine=True,
    min_domain_size=10
)
print(stats)
```

### Test differential APA

```python
from spagapa.analysis import test_differential_apa, find_domain_markers
import numpy as np

# Two-group comparison
g1 = np.where(labels == 0)[0]
g2 = np.where(labels == 1)[0]
results = test_differential_apa(apa_matrix, g1, g2,
                                 gene_names=gene_names,
                                 method='wilcoxon')
sig = results[results['padj'] < 0.05]
print(f"{len(sig)} significant genes")

# Domain markers (one-vs-rest)
markers = find_domain_markers(apa_matrix, labels, gene_names=gene_names)
print(markers[0].head())   # markers for domain 0
```

### Identify SVAPA genes

```python
from spagapa.analysis import identify_svapa_genes

svapa, results = identify_svapa_genes(apa_matrix, coords,
                                       gene_names=gene_names,
                                       fdr_threshold=0.05)
print(f"Found {len(svapa)} SVAPA genes")
```

---

## Visualization

```python
from spagapa.visualization import (
    plot_spatial_apa, plot_spatial_domains,
    plot_volcano, plot_heatmap,
    plot_imputation_quality
)

# Spatial APA map
ax = plot_spatial_apa(coords, apa_matrix[0], gene_name='GENE1',
                      cmap='viridis', save='gene1.png')

# Domain map
ax = plot_spatial_domains(coords, labels, save='domains.png')

# Volcano plot
ax = plot_volcano(results['log2fc'].values, results['padj'].values,
                  gene_names=gene_names, save='volcano.png')

# Heatmap
fig, ax = plot_heatmap(apa_matrix[:20, :50],
                       row_labels=gene_names[:20],
                       save='heatmap.png')

# Imputation QC
fig = plot_imputation_quality(observed_vals, imputed_vals, uncertainty_vals,
                               save='qc.png')
```

---

## Benchmarking

```python
from spagapa.benchmark import run_full_benchmark, generate_all_benchmark_figures

# Run simulation benchmark
results = run_full_benchmark(
    simulator_kwargs=dict(n_spots=300, n_genes=80),
    scenarios=[
        {'name': 'low_dropout',    'dropout_rate': 0.2},
        {'name': 'medium_dropout', 'dropout_rate': 0.5},
        {'name': 'high_dropout',   'dropout_rate': 0.7},
    ],
    output_dir='benchmark_results'
)

# Generate figures
figs = generate_all_benchmark_figures(results, output_dir='benchmark_results')
```

### Real data benchmark (MOB)

```python
from spagapa.benchmark import run_mob_benchmark

out = run_mob_benchmark(
    data_dir='path/to/mob_data',   # None = use simulator
    output_dir='benchmark_results/mob'
)
print(out['cv_results'])
```

---

## Complete Pipeline

```python
from spagapa import SpaGAPA

spa = SpaGAPA(
    n_neighbors=6,
    gp_kernel='matern',
    n_domains=5,
    diff_method='wilcoxon',
    svapa_fdr=0.05
)

results = spa.run(
    apa_matrix=apa_matrix,
    coordinates=coords,
    gene_names=gene_names
)

print(results['domain_labels'])
print(results['differential_apa'].head())
print(results['svapa_genes'][:10])
```

---

## FAQ

**Q: What input format does spaGAPA expect?**  
A: APA matrix as `np.ndarray` of shape `(n_genes, n_spots)` with NaN for missing values, and spatial coordinates as `(n_spots, 2)`.

**Q: How do I choose the GP kernel?**  
A: Use `'matern'` (default) for most biological data. `'rbf'` assumes smoother patterns. `'auto'` selects by cross-validation.

**Q: How long does GP imputation take?**  
A: ~1–5 s per gene on 500 spots. Use `SparseGPImputer` for >2000 spots or `BlockGPImputer` for very large datasets.

**Q: What is SVAPA?**  
A: Spatially Variable APA genes — genes whose APA usage pattern is non-random in space (detected by Moran's I or GP likelihood ratio test).

**Q: Can I use real MOB data?**  
A: Yes. Place `apa_matrix.csv`, `coordinates.csv` (and optionally `metadata.csv`) in a directory and pass it to `run_mob_benchmark(data_dir=...)`.

**Q: How do I cite spaGAPA?**  
A: Manuscript in preparation. Please check the README for the latest citation.

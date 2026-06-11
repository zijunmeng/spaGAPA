# spaGAPA

**spatial Gaussian process-based APA analyzer**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-268%20passing-brightgreen)]()
[![Status](https://img.shields.io/badge/status-beta-orange)]()

spaGAPA is a Python toolkit for analysing **alternative polyadenylation (APA)** in spatial transcriptomics data. It combines spatial-aware APA site validation, Gaussian process imputation, uncertainty quantification, and CPU-friendly multi-view graph machine learning for biologically coherent APA domain recovery.

---

## Key innovations

| Feature | spaGAPA | stAPAminer |
|---------|---------|-----------|
| APA calling | Spatial-aware validation | scAPAtrap only |
| Imputation | GP with uncertainty | Simple KNN |
| Spatial patterns | Moran's I + GP trend | Basic |
| Uncertainty estimates | ✅ | ❌ |
| Domain identification | GP + BioML multi-view graph | Manual / KNN |

---

## Installation

```bash
# From source (recommended during development)
git clone https://github.com/yourlab/spaGAPA.git
cd spaGAPA
pip install -e .

# With development dependencies
pip install -e ".[dev]"
```

**Requirements**: Python ≥ 3.8, numpy, scipy, scikit-learn, anndata, matplotlib, seaborn, statsmodels

---

## Quick start

```python
import numpy as np
from spagapa.benchmark import simulate_spatial_apa
from spagapa.imputation import GPImputer
from spagapa.analysis import identify_spatial_domains, test_differential_apa
from spagapa.visualization import plot_spatial_apa, plot_volcano

# 1. Generate / load data
data = simulate_spatial_apa(n_spots=300, n_genes=50, dropout_rate=0.4, random_state=42)
coords   = data['coordinates']          # (n_spots, 2)
observed = data['observed_apa'].T       # (n_genes, n_spots)

# 2. GP imputation with uncertainty
imputer = GPImputer(kernel_type='matern')
imputed, uncertainty = imputer.impute(observed[0], coords)

# 3. Identify spatial domains
labels, stats = identify_spatial_domains(observed, coords, method='kmeans', n_clusters=4)

# 4. Differential APA
group1 = np.where(labels == 0)[0]
group2 = np.where(labels == 1)[0]
results = test_differential_apa(observed, group1, group2)
print(results.head())

# 5. Visualise
ax = plot_spatial_apa(coords, observed[0], gene_name='Gene_0')
ax = plot_volcano(results['log2fc'].values, results['padj'].values)
```

### CLI with analysis presets

```bash
spagapa run \
  --apa-matrix apa_matrix.csv \
  --coordinates coordinates.csv \
  --expression-matrix expression_matrix.csv \
  --analysis-preset auto \
  --n-domains 5 \
  --output spagapa_results
```

`--analysis-preset auto` inspects the APA matrix shape and sparsity, then resolves
to either `standard` or `highres_accuracy`. Advanced users can force
`standard`, `highres_accuracy`, or `highres_fast`; they can also override the
automatic BioML decision with `--enable-bioml` or `--disable-bioml`.

BioML/preset outputs include `analysis_preset.json`, `domains.csv`,
`bioml_metadata.json`, `bioml_spot_factors.npy`, and
`bioml_imputed_values.npy`.

---

## Analysis workflow

```
BAM + Coordinates
      │
      ▼
Spatial-Aware APA Calling   ← spatial validation, quality filtering
      │
      ▼
GP Imputation               ← Matérn kernel, uncertainty estimates
      │
      ▼
APA Quantification          ← RUD, PDUI, WUL indices
      │
      ▼
Domain Identification       ← BioML multi-view graph / K-means
      │
      ├─► Differential APA  ← Wilcoxon / t-test, FDR correction
      │
      └─► SVAPA Detection   ← Moran's I + GP likelihood ratio test
              │
              ▼
         Visualisation      ← spatial plots, volcano, heatmap, QC
```

---

## Module overview

| Module | Description |
|--------|-------------|
| `spagapa.core` | `APADataset`, `APASite` data structures |
| `spagapa.io` | BAM reader, coordinate reader, result writers |
| `spagapa.spatial` | KNN / radius / Delaunay graph construction |
| `spagapa.calling` | Spatial validation, quality filtering |
| `spagapa.imputation` | `GPImputer`, `SparseGPImputer`, `BlockGPImputer` |
| `spagapa.quantification` | RUD, PDUI, WUL, PAI indices; QC metrics |
| `spagapa.analysis` | Domain ID, differential APA, SVAPA, spatial patterns |
| `spagapa.visualization` | Spatial, statistical, and QC plots |
| `spagapa.benchmark` | Simulator, evaluator, baseline methods |
| `spagapa.pipeline` | End-to-end `SpaGAPA` pipeline class |

---

## Examples

See the [`examples/`](examples/) directory:

| Script | Description |
|--------|-------------|
| `00_quick_demo.py` | 5-minute end-to-end demo |
| `03_gp_imputation.py` | GP imputation with uncertainty |
| `04_apa_quantification.py` | RUD / PDUI / WUL calculation |
| `05_differential_analysis.py` | Domain identification + differential APA |
| `06_visualization.py` | All plot types |
| `08_complete_pipeline.py` | Full pipeline |
| `09_benchmark.py` | Simulation benchmark |

---

## Benchmark results (simulation)

Comparison on simulated data (500 spots, 100 genes, 50% dropout):

| Method | RMSE ↓ | MAE ↓ | Pearson ↑ |
|--------|--------|-------|-----------|
| Mean | ~0.18 | ~0.14 | ~0.65 |
| Median | ~0.18 | ~0.14 | ~0.65 |
| KNN-spatial | ~0.15 | ~0.12 | ~0.75 |
| **spaGAPA-GP** | **~0.12** | **~0.09** | **~0.85** |

---

## Citation

If you use spaGAPA in your research, please cite:

> [Manuscript in preparation]  
> spaGAPA: spatial Gaussian process-based APA analysis for spatial transcriptomics

---

## License

MIT License — see [LICENSE](LICENSE).

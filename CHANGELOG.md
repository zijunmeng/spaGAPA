# Changelog

All notable changes to spaGAPA are documented here.

---

## [0.1.0] — 2026-03-11 (Beta)

### Added

**Core**
- `APADataset` class wrapping AnnData for APA-specific data
- `APASite` and `APASiteCollection` for APA site management

**I/O**
- `BAMReader`, `SpatialCoordinateReader`, `AnndataReader`, `BEDReader`
- `ResultWriter`, `BEDWriter` for output
- `ScAPAtrapWrapper` for R/scAPAtrap integration

**Spatial**
- `SpatialNeighbors`: KNN, radius, and Delaunay graph construction
- `build_knn_graph`, `build_radius_graph`, `build_delaunay_graph`

**Calling**
- `SpatialValidator`: spatial support scoring and consistency filtering
- `QualityFilter`: read count, spot count, CV-based filtering

**Imputation**
- `GPImputer`: Gaussian process imputation with Matérn/RBF kernels
- `GPImputerBatch`: parallel multi-gene imputation
- `SparseGPImputer`: inducing-point approximation (O(nm²))
- `BlockGPImputer`: block processing for large datasets
- Uncertainty quantification (standard deviation) for all methods

**Quantification**
- APA indices: RUD, PDUI, WUL, PAI
- Normalisation: z-score, min-max, quantile
- `QCReportGenerator` with cross-validation metrics

**Analysis**
- `DomainIdentifier`: K-means, Leiden, Louvain clustering with spatial refinement
- `DifferentialAPAAnalyzer`: Wilcoxon, t-test, permutation; FDR/Bonferroni correction
- `SpatialPatternAnalyzer`: Moran's I, SVAPA identification, pattern clustering
- `GPTrendDetector`: GP likelihood ratio test for SVAPA (core innovation)
- `find_domain_markers`: one-vs-rest marker identification

**Visualization**
- `SpatialPlotter`: spatial APA maps, domain plots, multi-gene grids
- `StatisticalPlotter`: volcano plots, heatmaps, box/violin plots
- `QCPlotter`: imputation quality, dropout statistics, spatial support

**Benchmark**
- `SpatialAPASimulator`: 4 spatial pattern types, configurable dropout/noise
- `MOBSimulator`: high-fidelity Mouse Olfactory Bulb simulator
- `BenchmarkEvaluator`: RMSE, MAE, Pearson, Spearman, R² metrics
- Baseline methods: mean, median, KNN-spatial (stAPAminer-like)
- `run_full_benchmark`, `run_mob_benchmark` convenience functions
- Publication-quality benchmark figures

**Pipeline**
- `SpaGAPA` end-to-end pipeline class integrating all modules

### Tests
- 268 unit tests, 100% passing
- Coverage: 80–99% for core modules

---

## Deferred to v1.1

- Spatial trajectory analysis (TrajectoryBuilder, TrajectoryAnalyzer)
- Interactive Plotly visualizations
- Command-line interface (CLI)
- Local Moran's I (LISA)
- Real MOB / Brain / Embryo dataset benchmarks (pending data download)
- PyPI release

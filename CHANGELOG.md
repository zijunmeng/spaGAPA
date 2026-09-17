# Changelog

All notable changes to spaGAPA are documented here.

---

## [Unreleased]

### Added

- `SparseGPImputer(transform='logit')`: latent logit-space sparse GP for
  bounded compositional data (APA usage fractions). Training values are
  mapped through `logit(clip(y, epsilon, 1 - epsilon))` and predictions
  back through the logistic sigmoid, guaranteeing output in (0, 1)
  without post-hoc clipping. `epsilon` (default 0.05) acts as
  pseudo-count shrinkage of boundary proportions. Posterior std is
  reported in latent space; conformal calibration stays an
  original-space procedure (`|y_true - y_pred|` scores). Default
  `transform=None` preserves the historical raw-value GP.
- `highres_accuracy` preset now enables `transform='logit'` for sparse
  GP imputation; other presets keep the raw-value GP.
- `scripts/test_logit_gp.py`: MOB st11 benchmark (top-100 multi-PAS
  genes, 20% masking, 5 seeds) comparing raw vs logit GP on RMSE,
  out-of-bound counts, conformal coverage/width, spatial gradient
  recovery, and runtime. Results in `benchmark_results/logit_gp/`.

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

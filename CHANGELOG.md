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
- `scripts/hd_gp_benchmark.py`: HD-scale sparse-GP + conformal benchmark
  on the Visium HD adrenal scAPAtrap counts (887 peaks x 2.0M barcodes),
  superseding the per-peak-mean-only `hd_conformal_quick.json`. Samples
  5k/20k/50k/100k barcodes (seed 42), top-300 peaks, 20% masking with
  50/50 cal/eval split, and compares the logit sparse GP
  (n_inducing=100) at the pipeline-default length scale vs an
  inducing-spacing-adapted scale vs the per-peak mean on RMSE, pooled
  split-conformal coverage/width at 80/90/95%, runtime, and peak RSS.
  Headline: coverage is exact at every scale (max deviation ~0.015);
  mean RMSE 0.273-0.283 vs GP 0.301-0.335 under pseudo-coordinates (no
  spatial-skill claim); 95% intervals ~full [0,1] range for every
  method (heavy-tailed single-read usage); fit 15-37s at 100k while
  full-surface predict is dominated by a 300x-redundant per-gene test
  kernel recomputation (cacheable); peak RSS 2.7 GB at 100k. Runs on
  random 2D pseudo-coordinates (HD barcode bridging unresolved).
  Results in `pipeline_output/hd_gp_benchmark/`.

### Fixed

- `morans_i_recovery` in `scripts/benchmark_stapaminer_headtohead.py`
  compared Moran's I computed on different spot supports: mean /
  spatial-KNN keep the sparse truth support (~363 spots/gene) while the
  GP batch `impute()` (and the R tools) return dense matrices, so the
  GP's per-gene Moran's I was computed over ALL spots on a different
  kNN graph — driving its recovery to -0.28 with no modelling failure.
  Predictions are now masked to the ground-truth observed support before
  the imputed statistic; corrected recovery is +0.65/+0.67 on
  gse183456/gse220442. Published benchmark tables regenerated from the
  corrected values; audit (reproduction, support table, mask-fraction
  sweep showing retention dominance, oracle check) in
  `pipeline_output/morans_i_audit/`.

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

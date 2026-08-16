# Supplementary Figure Legends (S1–S16)

Panel-level legends for spaGAPA's supplementary figures, matched to the
manuscript's *Supplementary Figures* section. All figures share the main-figure
formatting: vector PDF + 300-DPI PNG preview at print width (~178 mm),
Okabe–Ito palette, DejaVu Sans. Generation scripts live in this directory
(`supp_figNN_*.py`; see `README.md` for the script/data index and reproduce
commands); cached source data under `pipeline_output/supplementary_figures/`.

---

- **Supplementary Figure S1 — Dataset overview.** (A) PAS (scAPAtrap-called
  site) count per sample across the 32-sample Visium benchmark, grouped by GSE
  series and coloured by tissue. (B) Per-series summary table: species,
  platform, tissue, sample count, and mean/total spot and PAS counts.

- **Supplementary Figure S2 — Sparse-GP inducing-point sensitivity
  (GSE183456).** (A) Held-out RMSE (overall and median-per-gene) versus number
  of inducing points m ∈ {50, 100, 200, 500}; both are essentially flat, so the
  posterior mean is robust to m. (B) Wall time versus m, scaling near-linearly.

- **Supplementary Figure S3 — Masking-level sensitivity (GSE183456).**
  Imputation benchmark re-run at 10/20/30/50% masking for the three
  Python-native methods (spaGAPA-GP, mean, spatial-KNN): median per-gene
  held-out RMSE by method × masking level. stAPAminer/spvAPA (R-based) are
  reported at the 20% level used in the main benchmark.

- **Supplementary Figure S4 — Mean-baseline stratification.** Box plot of
  per-gene (GP RMSE − mean RMSE) stratified by gene-level spatial signal
  (Moran's I quintile). The per-gene mean is a strong RMSE baseline: GP beats
  it on only ~15% of genes, and the GP−mean gap does not shrink with spatial
  signal — the highest-spatial-signal quintile shows the largest positive
  median Δ (GP worse).

- **Supplementary Figure S5 — Benchmark parameter table.** The five benchmarked
  methods: implementation language, neighbour basis, key hyperparameters,
  free-parameter count, seed, and single-run wall time (GSE183456, same
  hardware). (Presentation note: a parameter table is inherently tabular and is
  best converted to a Supplementary Table.)

- **Supplementary Figure S6 — Full conformal coverage across 11 samples.**
  (A) Empirical vs nominal coverage at 80/90/95% (per sample). (B) Interval
  width 2·q̂ at each level. (C) Exact per-observation Winkler interval score
  (Gneiting & Raftery 2007), lower = better — computed from the persisted
  per-observation bounds (an earlier q̂+RMSE approximation had overstated the
  score by ~2.3× and is replaced). (D) Empirical subgroup coverage by
  uncertainty quintile (deviation from the 90% target) — a marginal, not
  conditional, diagnostic; the highest-uncertainty quintile over-covers.

- **Supplementary Figure S7 — Leakage audit and LOOCV method-selection
  stability.** (A) Clean vs leaked per-sample uncertainty/error correlation
  for the local-gene method across five datasets (≈ identical ⇒ no train/test
  leakage). (B) 5-fold leave-one-dataset-out: each held-out dataset selects
  the same method (5/5 stable), with held-out r per method per fold.

- **Supplementary Figure S8 — Per-gene uncertainty–error correlation.**
  Faceted per-dataset histograms of per-gene Pearson r (uncertainty ↔ absolute
  error) with the pooled across-gene r marked. Per-gene r is modest on average
  but centred above zero — uncertainty is informative within genes, not only
  globally.

- **Supplementary Figure S9 — Spatial-block vs random-split conformal
  coverage.** Empirical coverage under random vs spatially-disjoint block
  splitting, per dataset, at 80/90/95%. Intervals remain valid (within ~2 pp
  of nominal) under the harder block split — coverage is not an artefact of
  adjacent spots leaking across the train/test boundary.

- **Supplementary Figure S10 — Domain-recovery weight-configuration robustness
  on MOB.** (A) Leiden-resolution sweep (ARI/NMI vs resolution) for each weight
  configuration (mean_apa, apa_dominant, balanced, spatial_apa, expression_apa)
  with the best resolution per config marked. (B) ARI/NMI per weight regime:
  a weight-configuration robustness test (random-seed stability was not
  assessed). The global-best expression_apa configuration (ARI 0.597, k=5) and
  the mean_apa baseline (0.576) are reported alongside the apa_dominant
  configuration shown in main Figure 5.

- **Supplementary Figure S11 — Pseudoreplication analysis (GSE220442).**
  (A) Spot-level vs donor-level significant gene counts: the pseudoreplicated
  spot-level test calls 91 genes significant whereas the correct donor-level
  test (n = 3 ctrl vs n = 3 AD) calls 0; direction-agreement genes = 29.
  (B) Effect-size agreement for the 29 direction-agreement genes.

- **Supplementary Figure S12 — Stereo-seq QC and binning sensitivity
  (GSE263789).** (A) PAS count distribution by chromosome from scAPAtrap peak
  calls on the full 21,455-PAS call set (consistent with main Figure 7).
  (B) Read 3′-end-to-TES signed-distance histogram (500k gene-annotated reads),
  showing the 3′-enrichment (polyA-capture) signature of Stereo-seq.

- **Supplementary Figure S13 — Batch correction: QN vs Harmony (APA
  matrices).** (A) Mean pairwise PCC (biological signal recovery) vs residual
  batch signal per correction method. (B) Per-gene variance ratio
  (corrected / before). Caveat: Harmony is designed for expression counts, so
  this APA-matrix comparison is APA-specific and preliminary; batch correction
  remains an optional, still-under-evaluation module.

- **Supplementary Figure S14 — AD vs WT Stereo-seq effect sizes (n=1 per
  condition, descriptive).** Top effect sizes by |Δ| among the 1,090 genes
  observed in both conditions, plus the direction split of the 607 candidates
  exceeding |Δ| > 0.10. Effect sizes only — no p-values, no FDR (biological
  n=1 per condition); hypothesis-generating candidates.

- **Supplementary Figure S15 — Conformal coverage deviation forest.**
  Per-sample empirical−nominal deviation at 80/90/95% with binomial CI
  envelopes, sorted; moved out of main Figure 3 to keep it at print width.

- **Supplementary Figure S16 — MAE risk-coverage (paired unit for Figure 4F).**
  (A) Pooled MAE vs retention fraction. (B) Per-dataset MAE improvement at 80%
  retention (paired dots + mean, n = 5 datasets — the same pairing unit as
  main Figure 4F).

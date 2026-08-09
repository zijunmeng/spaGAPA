# spaGAPA BIB Paper — Main Figure Index

All figures: **300 DPI**, colorblind-friendly Okabe-Ito palette
(`#0072B2` blue / `#D55E00` orange / `#009E73` green),
DejaVu Sans, bold panel letters (A, B, C…), labeled axes with units, paired
per-dataset dots (not only mean bars). Generated 2026-07-28 on S91.

| # | File | Title | Panels |
|---|------|-------|--------|
| 1 | `fig1_framework_overview.png` | Framework overview | A: three gaps · B: input pipeline · C: sparse GP · D: split-conformal · E: outputs · F: dataset atlas |
| 2 | `fig2_spatial_vs_mean_benchmark.png` | Spatial vs mean benchmark | A: masking · B: per-dataset RMSE/corr · C: spatial fidelity · D: accuracy-spatial 2D · E: gene maps · F: stratified by spatial signal |
| 3 | `fig3_conformal_marginal_coverage.png` | Conformal marginal coverage | A: flow · B: per-sample 80/90/95 · C: calibration curve · D: deviation forest · E: width A vs D · F: spatial instance |
| 4 | `fig4_noise_and_risk_coverage.png` | Noise methods + risk-coverage | A: four noise models · B: multi-objective · C: pooled r · D: per-gene vs pooled · E: quintile coverage · F: subgroup heatmap · **G: risk-coverage (key)** |
| 5 | `fig5_domain_recovery.png` | Domain recovery | A: MOB layers · B: mean vs spaGAPA maps · C: ARI/NMI · D: Moran's-I recovery · E: gradient genes · F: cross-sample note |
| 6 | `fig6_scalability.png` | Scalability | A: complexity · B: runtime log-log · C: memory · D: completion matrix · E: accuracy-runtime Pareto · F: ablation → Supplementary |
| 7 | `fig7_stereo_seq.png` | Stereo-seq pilot | A: workflow · B: 3'-bias QC · C: scale & sparsity · D: atlas · E: PAS-by-gene · F: binning · (AD-WT descriptive -> Supplementary S14) |

---

## Figure 1 — Framework overview (`fig1_framework_overview.png`)
Conceptual/schematic figure (no quantitative data).
- **A. Three gaps in spatial APA.** Sparse PAS matrix (<12% observed), no calibrated
  uncertainty intervals, and O(N²) methods that cannot scale past 15k spots.
- **B. Input pipeline.** Spatial reads → PAS peak counts → APA usage matrix
  [gene × spot] → spaGAPA. Sparsity grid visualizes the <12% observed fraction.
- **C. Sparse Gaussian Process.** Posterior mean + 68/95% CI with M=150 inducing
  points; O(N·M²) complexity enables 100k-spot scaling.
- **D. Split-conformal calibration.** Train / calibrate / test partition;
  nonconformity score sᵢ=|yᵢ−ŷᵢ|; quantile qₜ gives ŷ±qₜ interval with finite-sample
  coverage guarantee.
- **E. spaGAPA outputs.** Posterior mean, interval width, uncertainty map,
  tissue domains.
- **F. Dataset atlas.** 8 tissues · 2 species (human/mouse) · 2 platforms
  (Visium/Stereo-seq); 11 calibration samples; up to 100k spots.

## Figure 2 — Spatial vs mean benchmark (`fig2_spatial_vs_mean_benchmark.png`)
Head-to-head, 5 methods × 2 datasets (GSE183456 kidney, GSE220442 brain), same 20%
per-gene mask (seed 42). **Honest**: mean is shown to beat GP on entry-wise RMSE.
- **A. Masking design.** 20%-of-observed per-gene holdout; same entries masked for
  every method.
- **B. Per-dataset metrics (paired dots).** RMSE (left) and Pearson/Spearman
  (right) with per-dataset dots + mean bars. **Mean wins RMSE/Pearson/Spearman;
  this is shown openly.**
- **C. Spatial fidelity.** GP 0.42 vs Mean 0.00 — mean imputes a per-gene constant
  and recovers no spatial gradient by construction.
- **D. Accuracy-spatial 2D scatter.** One point per method per dataset; mean is
  high-accuracy/zero-gradient, GP moderate-accuracy/strong-gradient.
- **E. Representative APA spatial maps (GSE183456).** Two high-spatial-variance
  genes: truth / mean (flat) / spaGAPA-GP (spatially smoothed) / |error|.
- **F. GP advantage grows with gene spatial signal.** Low/med/high variance strata;
  the GP-vs-mean RMSE gap narrows as spatial signal rises.

## Figure 3 — Conformal marginal coverage (`fig3_conformal_marginal_coverage.png`)
11 samples, 523,174 test points, 80/90/95% nominal.
- **A. Split-conformal flow.** Calibration set → nonconformity |yᵢ−ŷᵢ| → quantile qₜ
  → interval ŷ±qₜ → test coverage ≥ 1−α.
- **B. Per-sample coverage at 80/90/95%.** Dots per sample + nominal lines + 90%
  binomial CI band. Empirical lands within ~0.2% of nominal.
- **C. Calibration curve.** Nominal vs empirical; thin per-sample lines + bold
  pooled mean; y=x reference. Mean abs deviation 0.21% (80%), 0.16% (90%),
  0.10% (95%).
- **D. Coverage deviation forest.** Per-sample empirical−nominal at all three α.
- **E. Interval width: Constant (A) vs Residual-spot (D).** D ≈ same width as A →
  adaptive intervals without inflation.
- **F. Spatial instance.** One GSE183456 gene: posterior mean / interval width /
  |error| / covered-vs-uncovered (green/red).

## Figure 4 — Noise methods + risk-coverage curve (`fig4_noise_and_risk_coverage.png`)
4 noise methods (A constant, B local-gene, C spatial-spot, D residual-spot) × 5 datasets.
- **A. Four uncertainty models.** Schematic of each noise definition.
- **B. Multi-objective method selection.** x=uncertainty-error r, y=interval width,
  bubble size=|coverage deviation|. D (spaGAPA) is high-r, near-minimal-width.
- **C. Pooled r per dataset: A vs D.** A near-zero baseline; D clears the r≥0.3 target.
- **D. Per-gene vs pooled r.** Boxplot of within-gene r distribution + pooled-r star.
  Pooled r partly driven by cross-gene heterogeneity; within-gene discrimination is
  modest (median r ≈ 0.10–0.19).
- **E. Uncertainty quintile coverage.** Coverage climbs Q1→Q5 (high-uncertainty
  points correctly over-covered).
- **F. Conditional coverage heatmap.** Tissue × uncertainty quintile; green≈nominal.
- **G. Risk-coverage curve (KEY NEW PANEL).** Retention fraction vs RMSE for
  uncertainty-ranked / oracle / random retention. Uncertainty triage gives −66% RMSE
  at 20% retention and −23% at 80% (paired t-test p=0.004 vs random).

## Figure 5 — Domain recovery (`fig5_domain_recovery.png`)
MOB (mouse olfactory bulb) spvAPA ST array, 5 annotated layers, 260 spots.
- **A. MOB tissue: 5 olfactory bulb layers** (GCL/GL/MCL/ONL/OPL).
- **B. Domain maps.** Mean impute (no APA signal) vs spaGAPA (APA-dominant).
- **C. ARI/NMI bars.** All run configs; best ARI=0.60, NMI=0.68
  (expression_apa / leiden res0.8).
- **D. Spatial signal recovery (Moran's-I).** Per-method across the 2 benchmark
  datasets; GP reconstructs autocorrelation that mean cannot.
- **E. Representative APA gradient genes.** Four layer-stratified gradients.
- **F. Cross-sample note.** MOB is single-sample (n=1); multi-section domain
  stability reported in Supplementary (GSE237183).

## Figure 6 — Scalability (`fig6_scalability.png`)
Synthetic grids 1k–100k spots, 4 methods, 1200s wall cap.
- **A. Algorithmic complexity.** O(N²) KNN vs O(N·M²) sparse GP vs O(N) mean.
- **B. Runtime scaling (log-log).** spaGAPA-fast completes 100k in 511s;
  stAPAminer/spvAPA timeout/fail at 42k+.
- **C. Peak memory scaling.** spaGAPA-fast 8.8 GB at 100k; R-tools exceed limits.
- **D. Completion matrix.** Method × scale: green=completed, yellow=skipped
  (out of scope), red=timeout/OOM.
- **E. Accuracy-runtime Pareto.** At matched accuracy spaGAPA-GP is 5–7× faster
  than stAPAminer/spvAPA; mean/spatial-KNN are fast but recover no gradient.
- **F. Ablation** (inducing-point sweep) → moved to Supplementary per reviewer spec.

## Figure 7 — Stereo-seq pilot (`fig7_stereo_seq.png`)
GSE263789 mouse AD brain (Stereo-seq + scAPAtrap). Six panels (A–F); the
n=1-vs-n=1 AD-WT descriptive effect-size panel was moved to **Supplementary
Figure S14** to keep the main figure free of pseudoreplicated-looking claims.
- **A. Stereo-seq workflow.** DNB array → SAW count → scAPAtrap PAS peaks →
  bin 50/100/200 → spaGAPA impute.
- **B. 3'-end bias QC.** Reuses existing `fig1_3prime_enrichment_near_TES.png`
  (read coverage enrichment near TES confirms APA-relevant 3' capture).
- **C. Data scale & sparsity.** 20.7M DNB, 21,455 PAS, 15,235 bin200 spots.
  Sparsity is reported with a single consistent definition: the
  peak × spot **usage matrix** is **~10.3% observed (~90% empty)** = nnz 33.8M /
  326.9M entries (from `binned_200/qc_summary.json`). The raw peak × spot
  record count (97.8M) is a separate count and is **not** conflated with the
  observed fraction.
- **D. Tissue APA atlas.** GP-imputed spatial domains (n=41) + per-spot
  posterior uncertainty across bin200 coordinates. **Descriptive bin-level
  domains** computed on a 1,500 high-variance-gene subset; at this gene count
  the partition likely over-segments relative to a canonical tissue taxonomy,
  so it is shown as an exploratory spatial summary, not a ground-truth cell-type
  map.
- **E. Representative PAS, gene-annotated.** Three real high-spatial-variance
  PAS maps (bin200 usage), each annotated with its overlapping gene on mm10:
  **Cdk8** (peak_69378, chr5:146.26M), **Apoe** (peak_312125, chr7:19.70M —
  Alzheimer's APOE), **Gnb1l** (peak_415285, chr16:18.53M).
- **F. Binning consistency.** Cross-bin Pearson r for 50→100 and 50→200 over
  8 high-variance PAS (16 pairs); mean r=0.84, all pairs >0.7. Moran's I
  decreases with coarser binning as expected. The 8-PAS / 16-pair panel is a
  **descriptive cross-bin agreement check**, not a population-level claim.

---

## Data sources (all under `pipeline_output/`)
- Benchmark: `benchmark_mean_transparent/transparent_comparison.csv`,
  `benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json`
- Conformal: `conformal_validation/all_samples_coverage.csv` + `summary.json`
- Conditional coverage: `conformal_conditional_coverage/*.csv`
- Noise methods: `uncertainty_corr_improvement/per_dataset_corr.csv`
- Per-gene: `uncertainty_within_gene_audit/per_gene_corr_distribution.csv`
- Risk-coverage: `risk_coverage_curve/risk_coverage_data.csv`
- MOB: `mob_domain_recovery/spagapa_metrics.json` + `spagapa_domains.csv`
- Scalability: `benchmark_runtime/runtime_table.csv`
- Stereo-seq: `stereo_binning_consistency/binning_correlation.csv`,
  `gse263789_stereo_pilot/figures/fig1_3prime_enrichment_near_TES.png`

## Generation scripts
`fig1_overview.py`, `fig2_benchmark.py`, `fig3_conformal.py`, `fig4_noise.py`,
`fig5_domain.py`, `fig6_scalability.py`, `fig7_stereo_seq.py` (shared `"_style.py`).
Reproduce with:
```bash
OPENBLAS_NUM_THREADS=8 TMPDIR=/s3/mengzijun/tmp \
  ~/anaconda3/envs/spagapa/bin/python figN_*.py
```

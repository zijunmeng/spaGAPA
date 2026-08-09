# spaGAPA BIB Paper — Main Figure Index

All figures: **vector PDF + 300-DPI PNG preview**, colorblind-friendly Okabe-Ito
palette (`#0072B2` blue / `#D55E00` orange / `#009E73` green), DejaVu Sans, bold
panel letters. Every figure is sized at `PAGE_WIDTH_IN` (~7 in / ~178 mm) so text
stays ≥7–8 pt at print with **no downscaling**. Panels show data only;
explanatory prose lives in the figure caption. Regenerated 2026-08-10 on S91.

| # | File | Title | Panels |
|---|------|-------|--------|
| 1 | `fig1_framework_overview` | Framework overview | A: three gaps · B: input pipeline · C: sparse GP · D: split-conformal |
| 2 | `fig2_spatial_vs_mean_benchmark` | Spatial vs mean benchmark | A: masking · B: per-dataset RMSE/corr · C: spatial fidelity · D: accuracy-spatial 2D · E: representative spatial reconstruction · F: ΔRMSE by spatial-signal quintile |
| 3 | `fig3_conformal_marginal_coverage` | Conformal marginal coverage | A: flow · B: per-sample 80/90/95 · C: calibration curve · D: interval width A vs D · E: spatial instance |
| 4 | `fig4_noise_and_risk_coverage` | Uncertainty + risk-coverage | A: four noise models · B: multi-objective · C: pooled r · D: per-gene vs pooled · E: empirical subgroup coverage · F: risk-coverage (RMSE) |
| 5 | `fig5_domain_recovery` | Domain recovery | A: MOB layers · B: mean vs spaGAPA domains · C: ARI/NMI all configs · D: domain-layer confusion matrix · E: real gradient genes |
| 6 | `fig6_scalability` | Scalability | A: runtime scaling + empirical slopes · B: memory · C: completion matrix · D: spatial-fidelity–runtime Pareto |
| 7 | `fig7_stereo_seq` | Stereo-seq pilot | A: workflow · B: APA QC · C: scale & sparsity · D: real spatial APA examples · E: binning robustness |

> Deviation forest → **Supp. Fig. S15**; MAE risk-coverage → **Supp. Fig. S16**;
> n=1 AD-WT effect sizes → **Supp. Fig. S14**; global-best MOB config
> (expression_apa) → **Supp. Fig. S10**.

---

## Figure 1 — Framework overview (`fig1_framework_overview`)
Conceptual/schematic figure (no quantitative data). 4 panels.
- **A. Three gaps in spatial APA.** Sparse APA-usage matrix, no calibrated
  uncertainty intervals, and neighbor-based comparators that did not complete at
  ≥42k locations under benchmark limits.
- **B. Input pipeline.** Spatial reads → PAS peak counts → APA-usage matrix
  [gene × spot] → spaGAPA.
- **C. Sparse Gaussian Process.** Posterior mean + 68/95% CI with M inducing
  points; O(N·M²) complexity.
- **D. Split-conformal calibration.** Train / calibrate / test partition;
  nonconformity sᵢ=|yᵢ−ŷᵢ|; quantile q̂ gives ŷ±q̂ with finite-sample marginal
  coverage guarantee.

## Figure 2 — Spatial vs mean benchmark (`fig2_spatial_vs_mean_benchmark`)
Head-to-head, 5 methods × 2 datasets (GSE183456 kidney, GSE220442 brain), same
20% per-gene mask (seed 42). **Honest framing: per-gene mean is a strong
pointwise RMSE baseline; spaGAPA's value is spatial reconstruction + calibrated
uncertainty + scalable probabilistic inference.**
- **A. Masking design.** 20%-of-observed per-gene holdout; same entries masked
  for every method.
- **B. Per-dataset metrics (paired dots).** RMSE (left), Pearson/Spearman
  (right). Mean wins pointwise RMSE/Pearson/Spearman — shown openly.
- **C. Spatial fidelity.** GP 0.42 vs Mean 0.00 (mean imputes a per-gene
  constant, recovering no spatial gradient by construction).
- **D. Accuracy–spatial 2D scatter.** One point per method per dataset.
- **E. Representative spatial reconstruction.** Held-out truth / mean /
  **real spaGAPA-GP** (`SparseGPImputer`, peak_94938) / |error|. (Gene/Moran's-I
  in caption; data: `_data/fig2E_gse183456_peak_94938_gp_posterior.csv`.)
- **F. ΔRMSE by spatial-signal quintile.** Real per-gene Δ=GP−mean stratified by
  Moran's-I quintile. Δ>0 = mean performs better; the gap does **not** shrink
  with spatial signal. (Data: `_cache/s4_stratification.csv`.)

## Figure 3 — Conformal marginal coverage (`fig3_conformal_marginal_coverage`)
11 samples, 523,174 test points, 80/90/95% nominal. 5 panels.
- **A. Split-conformal flow.** Calibrate → |yᵢ−ŷᵢ| → q̂ → ŷ±q̂ → coverage ≥ 1−α.
- **B. Per-sample coverage at 80/90/95%.** Dots per sample + nominal lines + 90%
  binomial CI band.
- **C. Calibration curve.** Nominal vs empirical; y=x reference. **Mean abs
  deviation 0.21% (80%) / 0.16% (90%) / 0.10% (95%)** — reported precisely (80%
  is 0.21%, so "within 0.2%" is not strictly stated; max deviation also reported).
- **D. Interval width: Constant (A) vs Residual-spot (D).**
- **E. Spatial instance.** GSE183456 peak_20919: real posterior mean μ / σ /
  |error| / covered-vs-uncovered at conformal q̂₉₀=0.284. (Per-sample deviation
  forest moved to **Supp. Fig. S15**.)

## Figure 4 — Uncertainty + risk-coverage (`fig4_noise_and_risk_coverage`)
4 noise methods (A constant, B local-gene, C spatial-spot, D residual-spot) × 5
datasets. 6 panels.
- **A. Four uncertainty models.** Schematic.
- **B. Multi-objective method selection.** x=uncertainty–error r, y=interval
  width, bubble=|coverage deviation|.
- **C. Pooled r per dataset: A vs D.** A near-zero; D clears r≥0.3.
- **D. Per-gene vs pooled r.** Within-gene median r ≈ 0.10–0.19; pooled r ≈ 0.5
  partly driven by cross-gene heterogeneity (reported transparently).
- **E. Empirical subgroup coverage.** Left: tissue × uncertainty-quintile
  heatmap; **right: expression bins — high-expression bin undercovers (~0.82)**,
  shown as an honest limitation (guarantee is marginal, not per-subgroup).
- **F. Risk-coverage curve (RMSE).** Retention vs RMSE; uncertainty triage
  −23% RMSE at 80% retention (**paired across 5 datasets**, p=0.004 read from
  `risk_coverage_curve/summary.json`). MAE variant → **Supp. Fig. S16**.

## Figure 5 — Domain recovery (`fig5_domain_recovery`)
MOB (mouse olfactory bulb), 5 annotated layers, 260 spots. 5 panels.
- **A. MOB anatomy.** Ground-truth 5 layers (GCL/GL/MCL/ONL/OPL).
- **B. Domain maps (fair comparison).** Same Leiden pipeline (spatial graph,
  expression component, weights, resolution); **only the APA source differs** —
  per-gene **mean**-imputed vs spaGAPA **GP**-imputed. Both are real runs.
- **C. ARI/NMI across configs.** mean_apa, apa_dominant (displayed), spatial_apa
  (pure APA), balanced, expression_apa (global best, → S10). mean-impute
  recovers substantial structure on MOB — reported honestly.
- **D. Domain–layer confusion matrix.** Real correspondence of the displayed
  config's domains to the 5 anatomical layers (replaces the prior Moran's-I
  panel, whose caption had contradicted its data).
- **E. Representative real gradient genes.** Trim16 / Gdap2 / Tor3a / Pcdhb21 —
  sparse input (top) vs real spaGAPA GP reconstruction (bottom), color scaled
  per gene.

## Figure 6 — Scalability (`fig6_scalability`)
Synthetic grids 1k–100k spots, 4 methods, 1200 s wall cap. 4 panels.
- **A. Runtime scaling + empirical slopes.** Measured runtime points + fitted
  power-law slopes (spaGAPA-fast 0.84, spaGAPA-accuracy 1.01, stAPAminer 1.01,
  spvAPA 0.49). Slopes reflect deployed implementations and may differ from
  textbook complexity (approximate-neighbor structures can flatten observed
  scaling). No bare O(N²) assertion.
- **B. Peak memory scaling.** spaGAPA-fast 8.8 GB at 100k; R-tools exceed limits.
- **C. Completion matrix.** Method × scale: green=completed, red=timeout/OOM.
- **D. Spatial-fidelity–runtime Pareto.** Spatial-KNN has the highest fidelity
  (0.91) but is non-probabilistic and does not scale; spaGAPA-GP gives
  moderate fidelity + uncertainty + scalability at 6–7× the speed of
  stAPAminer/spvAPA; mean is fastest with ~0 fidelity.

## Figure 7 — Stereo-seq pilot (`fig7_stereo_seq`)
GSE263789 mouse AD brain (Stereo-seq + scAPAtrap). 5 panels.
- **A. Workflow.** FASTQ+mask → SAW 8.2.2 → retag → scAPAtrap → PAS matrix →
  spaGAPA.
- **B. Stereo-seq APA QC (native redraw).** TES-distance distribution +
  cumulative 3′-enrichment (53.4% within 500 bp, 73.5% within 2 kb).
- **C. Data scale & sparsity (single definition).** 20.7M DNB · 21,455 PAS ·
  15,235 bin200 spots · **~10.3% observed** (peak×spot usage matrix nnz fraction;
  not conflated with raw record counts).
- **D. Real spatial APA examples (main visual, enlarged).** Three real
  high-spatial-variance PAS maps (bin200), gene-annotated on mm10: **Cdk8**
  (peak_69378), **Apoe** (peak_312125, Alzheimer's APOE), **Gnb1l** (peak_415285).
- **E. Binning robustness.** Cross-bin Pearson r (50→100, 50→200) + Moran's I
  vs bin size. (The earlier 41-domain + raw-σ uncertainty panel was removed:
  its uncertainty carried 200×200 bin-grid × inducing-point geometry artifacts,
  and the domain count over-segmented.)

---

## Data sources (all under `pipeline_output/`)
- Benchmark: `benchmark_mean_transparent/transparent_comparison.csv`,
  `benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json`
- Conformal: `conformal_validation/all_samples_coverage.csv` (+ `winkler_*`
  columns) + `summary.json`; per-observation bounds `per_observation_bounds.npz`
- Conditional coverage: `conformal_conditional_coverage/*.csv`
- Noise methods: `uncertainty_corr_improvement/per_dataset_corr.csv`
- Per-gene: `uncertainty_within_gene_audit/per_gene_corr_distribution.csv`
- Risk-coverage: `risk_coverage_curve/{risk_coverage_data.csv,summary.json}`
- MOB: `mob_domain_recovery/{spagapa_metrics.json,spagapa_domains.csv}`
- Scalability: `benchmark_runtime/runtime_table.csv`
- Stereo-seq: `gse263789_stereo_pilot/spagapa_downstream_full/`,
  `stereo_binning_consistency/binning_correlation.csv`,
  `gse263789_ad_vs_wt_differential/sample_level_effect_sizes.csv`

## Generation scripts
`fig1_overview.py`, `fig2_benchmark.py`, `fig3_conformal.py`, `fig4_noise.py`,
`fig5_domain.py`, `fig6_scalability.py`, `fig7_stereo_seq.py` (shared `_style.py`).
Reproduce with:
```bash
OPENBLAS_NUM_THREADS=8 TMPDIR=/s3/mengzijun/tmp \
  ~/anaconda3/envs/spagapa/bin/python figN_*.py
```

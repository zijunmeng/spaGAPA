# spaGAPA BIB Paper — Supplementary Figure Generation Scripts

Canonical, version-controlled generation scripts for **Supplementary Figures S1–S16**.
Mirrors `scripts/main_figures/` (main figures 1–7); the rendered PNGs and cached
intermediate data are written under `pipeline_output/supplementary_figures/`
(gitignored — treated as generated outputs).

All figures: **300 DPI**, colorblind-friendly Okabe–Ito palette, DejaVu Sans,
bold panel letters. Shared formatting lives in `_style.py` (which in turn loads
`../main_figures/_style.py`).

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=8 TMPDIR=/s3/mengzijun/tmp   # S91; see CLAUDE.md per-host
~/anaconda3/envs/spagapa/bin/python supp_figNN_*.py
```

Two scripts pre-compute cached data used by downstream figure scripts — run them
first if the cache is stale or missing (they write to
`pipeline_output/supplementary_figures/_cache/`):

- `_gp_sweeps.py`  → `s2_inducing_sweep.csv`, `s3_masking_sweep.csv` (S2, S3)
- `_stereo_3end.py` → `s12_tes_distances.npy` (S12)

## Index

| # | Script | Title | Source data |
|---|--------|-------|-------------|
| S1 | `supp_fig01_dataset_overview.py` | Dataset overview (8 GSE, 32 samples) | `data/processed/*_scapatrap/qc_summary.json` |
| S2 | `supp_fig02_inducing_points.py` | Sparse-GP inducing-point sensitivity | `_cache/s2_inducing_sweep.csv` (`_gp_sweeps.py`) |
| S3 | `supp_fig03_masking_sensitivity.py` | Masking-level sensitivity (10–50%) | `_cache/s3_masking_sweep.csv` (`_gp_sweeps.py`) |
| S4 | `supp_fig04_mean_stratification.py` | Mean-baseline stratification by spatial signal | `_cache/s4_stratification.csv` (inline; cached) |
| S5 | `supp_fig05_benchmark_parameters.py` | Benchmark parameter table (5 methods) | `benchmark_fairness/parameter_table.csv` |
| S6 | `supp_fig06_conformal_coverage.py` | Full conformal coverage (11 samples) | `conformal_validation/all_samples_coverage.csv` |
| S7 | `supp_fig07_leakage_loocv.py` | Leakage audit + LOOCV method-selection stability | uncertainty audit outputs |
| S8 | `supp_fig08_per_gene_uncertainty.py` | Per-gene uncertainty–error correlation | `uncertainty_within_gene_audit/per_gene_corr_distribution.csv` |
| S9 | `supp_fig09_spatial_block_conformal.py` | Spatial-block vs random-split coverage | `conformal_conditional_coverage/*.csv` |
| S10 | `supp_fig10_domain_recovery.py` | Domain-recovery details on MOB | `mob_domain_recovery/` |
| S11 | `supp_fig11_pseudoreplication.py` | Pseudoreplication analysis (GSE220442) | donor-level diff-APA outputs |
| S12 | `supp_fig12_stereo_qc_binning.py` | Stereo-seq QC + binning sensitivity | `_cache/s12_tes_distances.npy` (`_stereo_3end.py`) |
| S13 | `supp_fig13_batch_correction.py` | Batch correction: QN vs Harmony | `apa_bias_correction/` |
| S14 | `supp_fig14_stereo_adwt_descriptive.py` | AD vs WT Stereo-seq APA effect sizes (n=1, descriptive) | `gse263789_ad_vs_wt_differential/sample_level_effect_sizes.csv` |
| S15 | `supp_fig15_conformal_deviation_forest.py` | Per-sample coverage deviation forest (80/90/95%) | `conformal_validation/all_samples_coverage.csv` |
| S16 | `supp_fig16_mae_risk_coverage.py` | MAE risk-coverage (paired unit for Fig 4F) | `risk_coverage_curve/risk_coverage_data.csv` |

Panel-level figure descriptions are documented in each `supp_figNN_*.py`
docstring (top-of-file); main figures 1–7 are indexed in
`scripts/main_figures/figure_index.md`. S14 houses the n=1-vs-n=1 AD-WT
effect-size analysis moved out of main Figure 7 — effect sizes only, **no
p-value / no FDR** (biological n=1 per condition). S15 (conformal deviation
forest) and S16 (MAE risk-coverage) were moved out of main Figures 3 and 4
respectively to keep the main figures at print width.

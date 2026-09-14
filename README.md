# spaGAPA: Spatial Gaussian-Process and Graph-Aware APA Analyzer

**A statistical framework for spatial alternative polyadenylation (APA) analysis with calibrated uncertainty quantification.**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests: 439](https://img.shields.io/badge/tests-439-green.svg)]()

---

## Overview

spaGAPA is a Python toolkit for spatial transcriptomics APA analysis. It addresses a critical gap in the field: **existing spatial APA tools (stAPAminer, spvAPA) lack uncertainty quantification, cannot scale to high-resolution platforms (Stereo-seq), and provide no mechanism for cross-study APA integration.**

spaGAPA solves these problems through four core innovations:

2. **Sparse Gaussian Process Framework** — O(n·m²) probabilistic imputation with heteroscedastic noise estimation; scales to 100k spots (competitors fail at 42k). Conformal coverage validated across **16 samples / 3 species / 2 platforms / 2 callers** (5.3M+ test points; mean |deviation| 0.21/0.16/0.10 pp at 80/90/95% on the 11-sample frozen set, ≤0.07 pp on 5 new samples). Cross-sample transfer decay <0.4 pp on average.
3. **Subcellular Stereo-seq Support** — the only tool validated on subcellular-resolution Stereo-seq APA data (21,455 PAS × 20.7M DNBs)
4. **Caller-Agnostic Design** — statistical guarantees hold unchanged when the PAS caller is swapped (scAPAtrap → Sierra on 3 datasets, 2 species × 3 tissues; max coverage deviation 0.5 pp)

APA batch correction (quantile normalization + linear removal) is provided as an optional module.

### Pipeline

```
Spatial Validation → Sparse GP Imputation (+ Conformal UQ) → APA Quantification
    → Leiden Domain Detection → Differential APA → SVAPA → Bias Correction
    → sAPA-RegNet (miRNA/RBP regulatory annotation)
```

---

## Key Results

### Head-to-Head Benchmark (vs stAPAminer + spvAPA)

| Method | RMSE | Pearson r | Spearman ρ | Wall Time |
|--------|------|-----------|------------|-----------|
| **spaGAPA-GP** | **0.122** | **0.942** | **0.842** | **19s** |
| stAPAminer | 0.208 | 0.852 | 0.778 | 147s |
| spvAPA | 0.272 | 0.770 | 0.697 | 189s |
| spatial-KNN | 0.147 | 0.913 | 0.849 | 1s |
| mean | 0.080 | 0.975 | 0.894 | 0.1s |

GP outperforms stAPAminer and spvAPA on spatial fidelity (0.42 vs 0.00 for the mean, which recovers no spatial gradient by construction) and runtime (**5.7–7.2× faster**), while the per-gene mean remains a strong RMSE baseline on the bimodal APA index — reported openly (see Limitations).

### Scalability (Competitors Fail at High Resolution)

| Spots | spaGAPA-fast | stAPAminer | spvAPA |
|-------|-------------|-----------|--------|
| 1k–15k | ✅ COMPLETED | ✅ COMPLETED | ✅ COMPLETED |
| 42k | ✅ **162s** | ❌ TIMEOUT (1200s) | ❌ FAILED (rc=1) |
| 100k | ✅ **511s**, 8.8 GB peak | ❌ TIMEOUT | ❌ FAILED |

Measured power-law slopes (deployed implementations): spaGAPA-fast 0.84, spaGAPA-accuracy 1.01, stAPAminer 1.01, spvAPA 0.49.

### Conformal Uncertainty (Universal Coverage)

**16 samples, 5.3M+ test points, 7 GSE datasets (Visium) + 3 GSE (Stereo-seq), 4 tissue types, 3 species (human/mouse/rat), 2 callers (scAPAtrap/Sierra):**

| Target | Frozen 11-sample set (mean |dev| / max) | New 5-sample set (mean |dev| / max) | Cross-sample transfer (mean decay / max) |
|--------|---------------------------------------|-------------------------------------|----------------------------------------|
| 80% | 0.21 pp / 0.4 pp | ≤0.05 pp / 0.07 pp | +0.19 pp / 11.77 pp |
| 90% | 0.16 pp / 0.5 pp | ≤0.07 pp / 0.07 pp | +0.33 pp / 9.34 pp |
| 95% | 0.10 pp / 0.2 pp | ≤0.04 pp / 0.04 pp | +0.38 pp / 6.75 pp |

Coverage also holds under spatial-block splits and locally adaptive intervals. Uncertainty-guided triage removes 23% of RMSE at 80% retention (paired across 5 datasets, p = 0.004). Cross-sample transfer (leave-one-out, 110 pairs): average decay <0.4 pp; 95% level is most robust.
### Multi-Caller Robustness (Sierra × spaGAPA)

Same Space Ranger BAMs re-called with Sierra 0.99.27, converted to spaGAPA input, **no re-tuning**:

| Dataset | Species / tissue | Sierra sites (genes) | Coverage 80/90/95 (global) | Gene-level distal-usage r | PAS Jaccard (±50 bp) |
|---------|-----------------|---------------------|---------------------------|---------------------------|----------------------|
| GSE183456 | human / kidney | 10,456 (4,258) | 0.801 / 0.900 / 0.950 | 0.51 (n=814) | 0.07 |
| GSE220442 | human / brain | 8,400 (3,384) | 0.796 / 0.898 / 0.950 | 0.83 (n=510) | 0.18 |
| GSE169749 | mouse / colon | 4,395 (1,872) | 0.804 / 0.899 / 0.950 | 0.69 (n=443) | 0.05 |

**Max deviation across all datasets × levels × modes (global + locally adaptive): 0.5 pp.** The framework's statistical guarantees do not depend on the PAS caller, the species, or the tissue. Full numbers: `scripts/multicaller_validation/report.md`.

### Domain Recovery (Unsupervised)

MOB (mouse olfactory bulb): **ARI = 0.60, NMI = 0.68** (unsupervised Leiden, no labels needed — vs spvAPA which requires supervised labels). Fair comparison: identical Leiden pipeline, only the APA source differs (mean-imputed vs GP-imputed).

---

## Datasets

### Visium (8 GSE, 32 samples)

| GSE | Tissue | Species | Samples | PAS Sites |
|-----|--------|---------|---------|-----------|
| GSE237183 | Glioma | Human | 18 | 20k–76k |
| GSE183456 | Kidney | Human | 1 | 53,572 |
| GSE179572 | Brain metastasis | Human | 1 | 26,478 |
| GSE220442 | AD brain (PFC) | Human | 6 (3 ctrl + 3 AD) | 27k–51k |
| GSE338525 | Liver (normal) | Human | 2 | 49k–55k |
| GSE206391 | Skin (psoriasis) | Human | 2 | 229–276 |
| GSE169749 | Colon (DSS) | **Mouse** | 1 | 40,795 |
| GSE263303 | Brain (Nf1) | **Mouse** | 1 | 36,586 |

### Stereo-seq (3 GSE, subcellular)

| GSE | Tissue | Species | Samples | PAS × Spots |
|-----|--------|---------|---------|-------------|
| GSE263789 | AD brain | Mouse | 5 (AD/WT/3mo) | 21,455 PAS × 20.7M DNB |
| GSE293464 | Retinal organoids (RA± × 16/26 wk) | **Human** | 4 | 18k–32k PAS × ~14–16k bin200 spots |
| GSE333693 | Thymus | **Rat** | 1 | 23,138 PAS × 13,917 bin200 spots |

Processed via SAW 8.2.2 → BAM retag adapter (CB=Cx_Cy, UB padded, fixed-length RX) → scAPAtrap → bin200 aggregation. Real spatial APA examples: mouse **Cdk8**, **Apoe** (Alzheimer's APOE), **Gnb1l**. The GSE293464/GSE333693 expansion (2026-09) carries deposited `barcodeToPos.h5` masks (GEO supplementary) — no STOmics retention-window dependency. Conformal coverage on all 5 new samples (11,259,201 test points, 30 sample×level×mode combinations): 80/90/95% all within **±0.07 pp** of nominal (e.g. D4: 0.7999/0.8999/0.9499). Cross-sample transfer (leave-one-out, 110 pairs): mean decay <0.4 pp at all levels, 95% level max deviation 6.75 pp. See `pipeline_output/stereo_expansion_downstream/conformal_expansion_final5.csv` and `pipeline_output/conformal_transfer/transfer_summary.json`.

### MOB (external APA matrix)

ST11 mouse olfactory bulb, 260 spots, 5 annotated layers — a published movAPA-prepared APA/RUD matrix (non-scAPAtrap input), used for domain recovery and cross-caller triangulation.

**3'-bias validation**: Mouse 53.4%, Human 47.0% within 500bp of TES (polyA-capture signature confirmed).

---

## Installation

### Requirements

```bash
# Python 3.10+
conda create -n spagapa python=3.10
conda activate spagapa
pip install numpy scipy scikit-learn pandas matplotlib pysam
pip install statsmodels

# For Stereo-seq processing: SAW 8.2.2, samtools, umi_tools, featureCounts
# For PAS calling: scAPAtrap (R) or Sierra (R); both supported
```

### From source

```bash
git clone <repo>
cd spaGAPA/spaGAPA
pip install -e .
```

### Environment (S91 server)

```bash
export OPENBLAS_NUM_THREADS=8     # CRITICAL: 64/32 causes segfault on heavy linalg
export TMPDIR=/s3/mengzijun/tmp
```

See `CLAUDE.md` / `docs/REPRODUCIBILITY_MANIFEST.md` for the per-host (S90/S91/S97/S98) environment table.

## Documentation

Full usage manual (installation → data preparation → CLI/API → presets →
uncertainty → Stereo-seq guide → troubleshooting) builds with Sphinx and is
hosted-ready for Read the Docs (`.readthedocs.yaml`):

```bash
cd docs && make html   # or connect the repo on readthedocs.org
```

---

## Quick Start

### Basic APA Analysis Pipeline

```python
from spagapa import SpaGAPA

pipeline = SpaGAPA(analysis_preset="highres_fast")
results = pipeline.run(
    apa_matrix="data/processed/gse183456_gsm6047774_scapatrap/apa_matrix.csv",
    coordinates="data/processed/gse183456_gsm6047774_scapatrap/coordinates.csv",
)
pipeline.save_results("out_run/")

print(f"Domains: {results['domains']['n_domains']}")
print(f"Imputed: {results['imputed_values'].shape}")
print(f"Uncertainty: {results['uncertainty'].shape}")
```

### Conformal Uncertainty Calibration

```python
from spagapa.imputation import SparseGPImputer, CalibratedInterval
from spagapa.imputation.calibration import ConformalCalibrator

# Fit GP
imputer = SparseGPImputer(n_inducing=200, length_scale=100.0, noise_level=0.1)
batch = imputer.fit_batch(coords, apa_matrix, mask=observed_mask)
predictions, uncertainty = batch.impute(return_uncertainty=True)

# Conformal calibration (split held-out into cal/test)
calibrator = ConformalCalibrator(alpha=0.1)  # 90% intervals
calibrator.fit(calibration_errors)  # |y_true - y_pred| on calibration set
interval = calibrator.predict(predictions, uncertainty)  # CalibratedInterval(lower, upper)

print(f"Coverage: {interval.coverage(test_truth)}")  # ≈ 0.90
```

### Differential APA Analysis

```python
from spagapa.analysis import DifferentialAPAAnalyzer

# IMPORTANT: for valid statistical inference, the replicate unit must be the
# biological sample (donor/animal), not the individual spot. Pooling spots
# from multiple samples inflates n (pseudoreplication). Per-gene per-sample
# aggregation should precede group comparison when n_samples >= 3.
analyzer = DifferentialAPAAnalyzer(method="t-test")
results = analyzer.test_differential_apa(
    apa_matrix=sample_level_matrix,  # gene × n_samples (NOT gene × all_spots)
    group1_indices=control_sample_indices,
    group2_indices=ad_sample_indices,
    gene_names=gene_names,
)
results = analyzer.adjust_pvalues(results, method="fdr_bh")
```

### Spatially Variable APA (SVAPA)

```python
from spagapa.analysis.svapa import svapa, morans_i

# Moran's I per gene
score, pval, n = morans_i(apa_values, coords, k=8)
svapa_results = svapa(apa_matrix, coords, k=8, n_perm=200, fdr="fdr_bh")
```

### APA Batch Correction

```python
from spagapa.analysis.bias_correction import quantile_normalize, linear_batch_correction

# Quantile normalization across batches
corrected = quantile_normalize(apa_matrix, group_labels=batch_labels)

# Linear batch removal (preserves biology)
corrected = linear_batch_correction(
    apa_matrix, batch_labels=batch_labels, preserve_labels=condition_labels
)
```

---

## Architecture

### Package Structure

```
spagapa/
├── core/                   # APA dataset, site, coordinate management
│   ├── apa_dataset.py      # APADataset: load, validate, align
│   ├── apa_site.py         # APASite: peak annotation, gene assignment
│   └── ...
├── calling/                # PAS quality + spatial validation
│   ├── spatial_validator.py
│   └── quality_filter.py
├── imputation/             # GP + uncertainty (CORE INNOVATION)
│   ├── sparse_gp.py        # SparseGPImputer + SparseGPImputerBatch
│   ├── calibration.py      # ConformalCalibrator + CalibratedInterval
│   ├── gp_imputer.py       # Dense GP (legacy)
│   └── expression_gp.py    # Expression-informed GP
├── bioml/                  # BioML graph + domain detection
│   ├── multiview_graph.py  # MultiViewGraphBuilder (spatial+APA+expr)
│   ├── factorization.py    # GraphRegularizedAPAFactorizer (chunked)
│   ├── highres.py          # HighresBioMLConfig + highres_bioml_recover()
│   └── domain.py           # BioMLDomainDetector (Leiden/spectral)
├── analysis/               # Downstream analysis
│   ├── differential.py     # DifferentialAPAAnalyzer
│   ├── svapa.py            # SVAPA (Moran's I + Geary's C)
│   ├── bias_correction.py  # QN + linear batch correction
│   ├── spatial_pattern.py  # Spatial pattern detection
│   ├── trajectory.py       # Spatial trajectory analysis
│   └── gp_trend_detector.py
├── benchmark/              # Benchmark suite
│   ├── evaluator.py        # Metric computation (RMSE, Pearson, coverage)
│   ├── simulator.py        # Synthetic APA data generator
│   └── real_data_benchmark.py
├── visualization/          # Plotting
└── pipeline.py             # SpaGAPA pipeline orchestrator
```

### Key Modules

#### `spagapa/imputation/sparse_gp.py`
- `SparseGPImputer`: inducing-point GP with configurable n_inducing, length_scale, noise_level
- `SparseGPImputerBatch`: batch fit across all genes with K_nm precomputation + inducing-point reuse
- `fit_batch(coords, values, mask)`: O(nm²) per gene, K_nm computed once for all genes
- `impute(return_uncertainty=True)`: returns predictions + posterior std

#### `spagapa/imputation/calibration.py`
- `ConformalCalibrator`: split-conformal prediction with global and locally-adaptive modes
- `CalibratedInterval`: lower/upper bounds with guaranteed marginal coverage
- `evaluate_coverage(predictions, truths, intervals)`: empirical coverage computation
- Uses finite-sample-corrected quantile: ceil((n+1)(1-α))/n

#### `spagapa/analysis/bias_correction.py`
- `quantile_normalize(apa, group_labels)`: rank-based marginal alignment across batches
- `linear_batch_correction(apa, batch_labels, preserve_labels)`: limma-style removeBatchEffect
- Auto-detects collinear designs (warns when batch ⊄ biology)

#### `spagapa/bioml/highres.py`
- `highres_bioml_recover()`: full high-resolution pipeline (GP + BioML + domain)
- Auto-enables: sparse GP (n_spots > 5k), chunked factorizer (n_spots > 20k), Leiden domain (>5k)
- kNN-Laplacian fix: factorizer uses spatial_laplacian (not fused graph) for O(seconds) solve

---

## Manuscript Figures (BIB)

All figures are **vector PDF + 300-DPI PNG**, colorblind-friendly Okabe-Ito palette, DejaVu Sans, sized at BIB print width (~178 mm) with no downscaling. Generation scripts are version-controlled; rendered outputs live under `pipeline_output/` (gitignored).

### Main figures 1–7 — `scripts/main_figures/` (index: `figure_index.md`)

| # | Title | Key content |
|---|-------|-------------|
| 1 | Framework overview | gaps · input pipeline · sparse GP · split-conformal |
| 2 | Spatial vs mean benchmark | 5 methods × 2 datasets, honest mean-baseline framing |
| 3 | Conformal marginal coverage | 11 samples, 523k points, 80/90/95% |
| 4 | Uncertainty + risk-coverage | 4 noise models, subgroup coverage, −23% RMSE triage |
| 5 | Domain recovery | MOB layers, fair mean-vs-GP comparison, ARI/NMI |
| 6 | Scalability | 1k–100k runtime/memory/completion/Pareto |
| 7 | Stereo-seq pilot | 20.7M DNB, QC, Cdk8/Apoe/Gnb1l maps, binning robustness |

### Supplementary figures S1–S16 — `scripts/supplementary_figures/` (legends: `supp_figure_legends.md`)

Dataset overview, inducing-point & masking sensitivity, mean stratification, benchmark parameter table, full conformal coverage, leakage/LOOCV audit, per-gene uncertainty correlation, spatial-block conformal, domain-recovery details, pseudoreplication, Stereo-seq QC/binning, batch correction, AD-vs-WT effect sizes (descriptive), coverage deviation forest, MAE risk-coverage.

---

## Reproducibility

### Benchmark & Validation Scripts

| Script | Purpose |
|--------|---------|
| `scripts/benchmark_stapaminer_headtohead.py` | 5-method imputation head-to-head (GP/stAPAminer/spvAPA/KNN/mean) |
| `scripts/benchmark_runtime_scalability.py` | Runtime/memory at 1k–100k spots |
| `scripts/benchmark_mean_transparent_report.py` | Transparent mean-baseline comparison |
| `scripts/calibrate_uncertainty_all_datasets.py` | Multi-dataset conformal coverage validation |
| `scripts/conditional_coverage_validation.py` | Subgroup + conditional coverage |
| `scripts/result12_riskcoverage_blockconformal.py` | Risk-coverage + spatial-block conformal |
| `scripts/audit_uncertainty_leakage.py` (+ 3 more `audit*.py`) | Uncertainty audits (leakage/LOOCV/per-gene/width) |
| `scripts/analyze_gse220442_diff_apa_unified.py` | Unified-peak differential APA (3v3, donor-level) |
| `scripts/mob_domain_recovery.py` | MOB unsupervised domain ARI/NMI |
| `scripts/run_bias_correction.py` | APA batch correction demo |
| `scripts/multicaller_validation/` | Sierra × spaGAPA caller-robustness chain (one-command `run_dataset.sh`) |
| `scripts/main_figures/` · `scripts/supplementary_figures/` | All figure generation scripts |

### Processing Logs & Manifests
- `logs/20260710_数据处理记录.md` – `logs/20260726_数据处理记录.md`
- `docs/REPRODUCIBILITY_MANIFEST.md` – Full environment + command manifest
- `docs/apa_evidence_source_manifest.md` – APA-evidence decision rules per source
- `docs/real_data_benchmark_protocol.md` – Real-data benchmark track definitions
- `docs/TABLE1_tool_comparison.md` – Manuscript Table 1

### Tests
```bash
OPENBLAS_NUM_THREADS=8 python -m pytest tests/ -q
# 439 tests collected; 1 pre-existing collection error (tests/benchmark/test_performance.py)
```

---

## Innovation Summary

### What spaGAPA Does That No Competitor Can

| Capability | stAPAminer | spvAPA | metaAPA | **spaGAPA** |
|-----------|-----------|--------|---------|------------|
| Imputation method | KNN (heuristic) | WNN (heuristic) | N/A | **Sparse GP (probabilistic)** |
| Uncertainty quantification | ❌ | ❌ | caller-level only | **✅ Conformal (guaranteed coverage)** |
| Caller-agnostic guarantees | ❌ | ❌ | integrates callers | **✅ validated (scAPAtrap + Sierra)** |
| Batch correction | ❌ | ❌ | ❌ | **✅ QN + linear** |
| Stereo-seq (subcellular) | ❌ (Visium only) | ❌ (Visium only) | N/A | **✅ 21k PAS × 20.7M DNB** |
| Scalability | ❌ OOM @ 42k | ❌ FAIL @ 42k | N/A | **✅ 100k spots** |
| Supervised analysis | ❌ | ✅ (sPLS-DA, needs labels) | ❌ | ❌ (unsupervised) |
| Validation data | MOB × 3 | 9 datasets (sc + ST) | 4 datasets | **9 GSE × 32 samples + MOB + Stereo-seq** |

### Contributions to the Field

1. **Statistical rigor**: First spatial APA tool with conformal-calibrated uncertainty (distribution-free coverage guarantee)
2. **Caller robustness**: Coverage guarantees empirically invariant to the PAS caller (scAPAtrap/Sierra), species, and tissue
3. **Scalability**: Sparse GP framework extends spatial APA to subcellular-resolution platforms (competitors are fundamentally O(n²))
4. **Cross-study integration**: APA-specific batch correction (existing tools are expression-only)
5. **High-resolution validation**: First complete APA pipeline on Stereo-seq subcellular data
6. **Competitive benchmarking**: Head-to-head comparison on shared data with transparent baseline reporting (GP > stAPAminer > spvAPA on spatial fidelity + speed; mean baseline reported honestly)

---

## Limitations (Honest)

1. **GP doesn't beat per-gene mean on entry-wise RMSE** — the gene-level distal-usage index is bimodal; mean predicts the dominant mode. GP's value is in spatial fidelity (0.42 vs 0.00), uncertainty, and scalability, not raw RMSE. This pattern holds on every caller input tested.
2. **Raw GP uncertainty correlation is modest within genes** (per-gene median r ≈ 0.10–0.19; pooled r ≈ 0.5 partly driven by cross-gene heterogeneity). Residual-based heteroscedastic estimation (method D) clears pooled r ≥ 0.3 and enables risk-coverage triage; conformal calibration provides distribution-free marginal coverage regardless.
3. **Marginal ≠ conditional coverage** — high-expression bins undercover (~0.82); the guarantee is marginal, not per-subgroup (shown explicitly in Fig 4E / Supp Fig S9).
4. **Stereo-seq raw FASTQ + mask availability is a field-wide bottleneck, but improving** — STOmics chip masks are retained on OSS only for a limited window after each experiment and expire for published datasets (documented with official reply: [STOmics/SAW#268](https://github.com/STOmics/SAW/issues/268); GSE269906 human AD brain is mask-blocked for this reason). Mitigation demonstrated (2026-09): a systematic GEO survey (75 series → 28 with DNBSEQ raw data) identified datasets whose authors deposited `barcodeToPos.h5` directly in GEO supplementary; via that route spaGAPA now runs on human retinal-organoid (GSE293464, 4 samples) and rat thymus (GSE333693) Stereo-seq data with unchanged conformal guarantees.
5. **sAPA-RegNet perturbation model not validated** — cis-regression go/no-go found the simple stability model unsupported at spot level. The regulatory annotation remains descriptive.
6. **No supervised analysis** — spaGAPA is unsupervised; spvAPA offers supervised sPLS-DA.
7. **Batch-correction module still under evaluation** — Harmony comparison pending; listed for completeness.

---

## Citation

```
spaGAPA: A statistical framework for spatial alternative polyadenylation analysis
with calibrated uncertainty quantification.
```

## License

MIT

## Contact

For questions about the pipeline, datasets, or collaboration:
- Project: spaGAPA (Spatial Gaussian-Process and Graph-Aware APA Analyzer)
- Target journal: Briefings in Bioinformatics (BIB)

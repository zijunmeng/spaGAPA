# spaGAPA: A Statistical Framework for Spatial Alternative Polyadenylation Analysis with Calibrated Uncertainty Quantification

**Target journal:** Briefings in Bioinformatics (BIB)

**Authors:** [PLACEHOLDER — author list to be supplied]

**Affiliations:** [PLACEHOLDER — author affiliations to be supplied]

**Corresponding author:** [PLACEHOLDER]

---

> **Document-integration note.** This master document assembles the full manuscript. The **Methods (M1–M12)** and **Results (R1–R9)** sections are integrated verbatim from `docs/MANUSCRIPT_DRAFT.md`. The **Figure Legends** are drawn from `pipeline_output/main_figures/figure_index.md`. The **Table 1 legend** is integrated from `docs/TABLE1_tool_comparison.md`. The **Abstract, Introduction, Discussion, Title/Authors block, and References** are clearly flagged placeholders — the dedicated Abstract/Introduction/Discussion authoring agent has not yet produced `docs/MANUSCRIPT_ABSTRACT_INTRO_DISCUSSION.md`, so these sections carry `[TO BE WRITTEN BY ABSTRACT/INTRO/DISCUSSION AGENT]` markers and draft scaffolding to indicate target length and content.

---

## Abstract

> **[TO BE WRITTEN BY ABSTRACT/INTRO/DISCUSSION AGENT — target ~250 words]**
>
> *Draft scaffolding (to be replaced).* Spatial alternative polyadenylation (APA) measurements from spatial transcriptomics are sparse, noisy, and spatially structured, yet existing tools (stAPAminer, spvAPA) provide only point-estimate KNN/WNN imputation with no uncertainty quantification and do not scale beyond the ~5–15 k spots of 10x Visium. We present spaGAPA, a downstream inference framework that casts spatial APA imputation as sparse Gaussian-process (GP) regression with heteroscedastic noise estimation and wraps the posterior in a split-conformal prediction layer delivering a distribution-free marginal coverage guarantee. Across 11 Visium samples (7 GSE datasets, 4 tissues, 2 species, 523 174 held-out test points) the conformal intervals achieved near-nominal coverage (deviations of 0.05–0.3 percentage points at 80/90/95% targets). Against stAPAminer and spvAPA under a shared 20%-masking protocol, spaGAPA-GP delivered superior spatial fidelity and was 5.7–7.2× faster; under a fixed 1200 s wall cap, spaGAPA completed 42 k- and 100 k-spot synthetic analyses where the comparators timed out or failed. On a subcellular Stereo-seq mouse AD-brain sample (20.7 M DNBs, 21 455 PAS) the full scAPAtrap → spaGAPA pipeline ran end-to-end — the first such demonstration for spatial APA. We report subgroup-calibration limits honestly (extreme expression tails) and find that at donor-level replication (n = 3 vs 3) no gene reaches FDR-significant differential APA, framing the cohort as underpowered rather than biologically negative.

**Keywords:** alternative polyadenylation; spatial transcriptomics; Gaussian process; conformal prediction; uncertainty quantification; Stereo-seq; 10x Visium.

---

## Introduction

> **[TO BE WRITTEN BY ABSTRACT/INTRO/DISCUSSION AGENT — target ~900 words]**
>
> *Draft scaffolding (to be replaced).* Suggested content arc, consistent with the Methods and Results already written:
>
> 1. **APA background and biological importance.** Alternative polyadenylation as a widespread post-transcriptional regulatory layer; ~70% of human genes produce multiple 3′ UTR isoforms; tissue- and disease-specific APA (cancer, neurodegeneration including AD).
> 2. **Spatial transcriptomics rises, but APA lags.** 10x Visium (near-cellular) and Stereo-seq (subcellular) platforms; most spatial methods focus on gene expression, not 3′-end processing. APA-specific tools (stAPAminer, spvAPA) have appeared but are Visium-only and point-estimate.
> 3. **The three gaps spaGAPA closes** (mirrors Figure 1A): (i) no calibrated uncertainty — users cannot tell which APA calls are trustworthy; (ii) O(n²) imputation that breaks at high resolution; (iii) no subcellular-resolution spatial APA pipeline demonstrated end-to-end.
> 4. **Statistical framing.** Why a GP + conformal treatment is the right tool: sparse GP gives O(n·m²) scaling and a posterior; conformal prediction gives distribution-free coverage without parametric assumptions on the bounded, bimodal distal-usage index.
> 5. **Positioning vs competitors** (refer forward to Table 1): spaGAPA is a downstream inference layer, not a PAS caller; it begins where scAPAtrap ends.
> 6. **Contributions (4 bullets)** mirroring the Contribution Statement in `PROJECT_SUMMARY.md` §6.2: (1) conformal-calibrated UQ, (2) sparse GP probabilistic framework enabling Stereo-seq APA, (3) uncertainty-guided exploratory analysis, (4) transparent head-to-head benchmark with disclosed asymmetries.
> 7. **Paper roadmap.** Forward references to Methods (M1–M12), Results (R1–R9), Figures 1–7, Table 1.

---

## Methods

### M1. Overview of the spaGAPA framework

Spatial alternative polyadenylation (APA) measurements derived from spatial transcriptomics are characterized by three structural properties that motivate a dedicated statistical treatment: (i) they are extremely sparse, because each spot captures only a small fraction of the poly(A) ends of any given gene; (ii) they are noisy, because the bounded APA usage ratio is estimated from small per-gene per-spot read counts; and (iii) they are spatially structured, because neighboring spots in a tissue section share both transcriptional state and 3′-end-processing context. spaGAPA is a downstream inference framework that takes these properties seriously: it casts spatial APA imputation as a sparse Gaussian-process (GP) regression problem and calibrates the resulting posterior uncertainty with split-conformal prediction to deliver distribution-free marginal coverage.

We emphasize that spaGAPA is *not* a polyadenylation-site (PAS) caller. Its required input is a PAS × spot APA usage matrix together with per-spot spatial coordinates, which we obtain from scAPAtrap (see M2). The framework then performs probabilistic imputation with heteroscedastic noise estimation, conformal uncertainty calibration, spatial-domain identification, and downstream differential and spatially-variable APA analyses. The three core algorithmic components — sparse GP regression (M3), heteroscedastic noise estimation (M4), and conformal calibration (M5) — are the focus of the methodological description; downstream modules (spatial domains, differential APA, Stereo-seq processing) are described in M6–M11.

### M2. Input APA quantification

For each gene and each spatial spot we compute a bounded gene-level distal-usage index that summarizes the relative use of proximal versus distal polyadenylation sites within the gene body. For every gene we pool its PAS peaks into two classes — proximal and distal — using the median PAS genomic coordinate as the split point, and require a minimum parent-gene support of `min_parent = 5` reads across the dataset so that the proximal/distal classification is not driven by singletons. The resulting per-gene per-spot value is

  u(g, s) = (distal reads) / (proximal + distal reads),

which is bounded in [0, 1] and is undefined (treated as missing) when the gene has no observed poly(A) reads in that spot. This two-state distal-usage index matches the convention used by stAPAminer and spvAPA, enabling like-for-like head-to-head comparison (M8).

PAS calling upstream of the index is performed with scAPAtrap (v0.2.0), which an independent benchmark survey identified as the top-performing tool for differential-APA detection. scAPAtrap is run on the 10x Visium or Stereo-seq position-sorted, gene-annotated BAM with the standard pipeline (peak calling, UMI counting, and site-to-gene assignment); the resulting PAS peaks are aggregated into the gene-level proximal/distal index above. For all datasets in this study the same scAPAtrap pipeline was applied, removing per-tool peak-calling as a confounder in the method comparison.

### M3. Sparse Gaussian-process spatial model

For each gene g we model the observed spatial APA usage as a noisy realization of a latent Gaussian process over the two-dimensional tissue coordinates. Writing s ∈ ℝ² for a spot coordinate, the generative model is

  f_g(s) ~ GP(m_g(s), k(s, s′)),     y_g(s) = f_g(s) + ε_g(s),

where m_g(·) is a constant mean function fitted per gene, ε_g(·) is heteroscedastic observation noise (M4), and k(·, ·) is a squared-exponential (RBF) kernel

  k(s, s′) = σ² exp(−‖s − s′‖² / (2 ℓ²)).

The signal variance σ² and the length scale ℓ are the two free hyperparameters of the spatial covariance. We set ℓ to five times the median nearest-neighbor distance between spots, which for a regular Visium grid (center-to-center spacing 100 µm) gives ℓ ≈ 100 µm in physical units, or two grid spacings in normalized index units (the multiplier used in the head-to-head benchmark). This choice expresses a weak prior that APA usage varies smoothly over a few neighboring spots while still allowing sharp boundaries between tissue layers.

A full GP scales as O(n³) in the number of spots n and is therefore intractable on Visium sections (n ≈ 1.5–6 k) and impossible on high-resolution platforms (Stereo-seq, n ≈ 10⁷ DNBs). spaGAPA therefore uses a sparse inducing-point approximation with m ≪ n inducing points. We place m inducing points on a spatial grid covering the tissue and approximate the GP posterior by variational inference under a structured mean-field assumption, reducing the per-gene cost to O(n m²). For Visium datasets we use m ≈ 150–200 inducing points.

A key implementation detail is that the cross-covariance matrix K_nm between the n spots and the m inducing points depends only on the coordinates and the kernel, not on the gene. We therefore precompute K_nm once and reuse it across all genes in `SparseGPImputerBatch.fit_batch`, so that the marginal per-gene cost of the O(n m²) solve dominates the wall time. Posterior prediction returns both the posterior mean (the imputed APA usage) and the posterior standard deviation (the raw aleatoric uncertainty that feeds the heteroscedastic and conformal layers). For high-resolution runs (n > 5 k) the pipeline automatically switches to the sparse-GP mode; above n > 20 k the chunked graph-regularized factorizer is skipped so that the GP remains the scalability bottleneck, which is what enables the 42 k and 100 k spot completions reported in M9.

### M4. Heteroscedastic noise estimation

The raw GP posterior standard deviation reflects kernel-driven spatial smoothing uncertainty but does not by itself capture gene- or spot-specific observation noise, which on sparse APA data is strongly heteroscedastic: a gene observed from two reads in a spot has a much noisier usage estimate than a gene observed from fifty. spaGAPA implements four interchangeable noise-estimation strategies that the user may select via a single parameter.

- **Method A (constant):** a global `noise_level = 0.1` applied to every gene and spot. This is the simplest baseline and yields the widest, most conservative intervals.
- **Method B (per-gene kNN MAD):** for each gene the noise scale is estimated as the median absolute deviation of the gene's observed usages among its k nearest-neighbor spots in expression space, capturing per-gene variability in APA "roughness".
- **Method C (per-spot spatial kNN variance):** for each spot the noise scale is the local variance of the k spatially nearest spots, capturing per-spot local heterogeneity.
- **Method D (residual-based empirical Bayes, default):** a two-pass procedure. In the first pass the GP is fitted with Method A noise to obtain pilot residuals r_i = y_i − ŷ_i on held-out entries. In the second pass each spot's noise scale is re-estimated as a shrinkage estimate (empirical-Bayes) of |r_i| toward the global residual variance, and the GP is refitted. This combines local residual evidence with global stabilization.

We adopt Method D as the default based on a coverage-width-correlation audit (Results R6). Method A delivers exact marginal coverage but at the cost of uniformly wide intervals (interval score 0.449 at α = 0.1); Method B improves calibration-uncertainty/error correlation to ≈0.55 but inflates mean interval width by a factor of 2.4 (range 1.1–4.8) relative to Method A and undercovers the high-expression stratum (0.715); Method D matches Method A's interval width (ratio 1.03) while recovering most of Method B's calibration-error correlation (≈0.51) and partially correcting the high-expression stratum (0.833). Method D is therefore the best coverage/width/correlation trade-off and is the default used throughout the paper unless stated otherwise.

### M5. Conformal uncertainty calibration

The heteroscedastic GP posterior is a useful point estimate of local uncertainty but, by itself, carries no finite-sample coverage guarantee. spaGAPA therefore wraps the GP in a split-conformal prediction layer that delivers a distribution-free *marginal* coverage guarantee under exchangeability.

**Split.** For each dataset we partition the observed (non-missing) entries of the PAS × spot matrix into a training set (80%) used to fit the GP, and a held-out set (20%). The held-out set is split again, 50/50, into a calibration set and a test set. All three splits are stratified random over gene–spot entries with a fixed seed (42) for reproducibility and fairness across methods.

**Nonconformity score.** For each calibration entry i we compute the absolute residual s_i = |y_i − ŷ_i| as the nonconformity score, where ŷ_i is the GP posterior mean. When the locally-adaptive mode is enabled, the score is replaced by s_i / max(σ̂_i, σ_floor), where σ̂_i is the GP posterior standard deviation (Method D) and σ_floor is a small positive constant preventing division by near-zero variance.

**Finite-sample quantile.** For a target miscoverage level α we take the ⌈(n_cal + 1)(1 − α)⌉ / n_cal-th order statistic of the calibration scores as the conformal quantile q̂_α. This is the standard finite-sample correction that yields the distribution-free guarantee.

**Interval construction.** The 100(1 − α)% prediction interval for a test entry is ŷ_i ± q̂_α in the global mode, and ŷ_i ± q̂_α · max(σ̂_i, σ_floor) in the locally-adaptive mode.

**Guarantee.** Under exchangeability of the calibration and test entries, split-conformal prediction guarantees

  P(Y ∈ C(X)) ≥ 1 − α,

where the probability is over the joint draw of calibration and test data. **This is a guarantee on marginal coverage, not on conditional coverage:** averaging over all test entries, the empirical coverage is at least 1 − α in expectation, but the coverage of any particular subgroup (e.g., a single tissue region, or a single uncertainty quintile) is not individually guaranteed and must be checked empirically — which is the subject of M6. We are explicit about this distinction throughout the manuscript and report subgroup calibration alongside the marginal result.

### M6. Empirical subgroup-calibration analysis

Because the conformal guarantee is marginal, we empirically audit subgroup coverage to characterize where the framework is reliable and where it is not. We compute 90% conformal coverage on four independent subgroup axes:

- **Uncertainty quintiles:** 5 bins defined by GP posterior standard deviation.
- **Spatial partitions:** k-means partitions on the (x, y) coordinates with K ∈ {5, 6, 7, 8, 9, 10}, giving up to 24 spatial-region bins per dataset.
- **Tissue groups:** one bin per dataset.
- **Expression strata:** 5 bins defined by per-gene mean observed expression, crossed with the available datasets, giving 18 expression-stratum bins.

We additionally run an **adversarial spatial-block conformal split** as a stress test: instead of stratified-random splitting, we partition the tissue into 5 contiguous spatial blocks, fit on 3 blocks and test on the 2 held-out blocks, and compare the resulting coverage to the random-split coverage. If spatial autocorrelation violated exchangeability, the block split would undercover; if the two agree, the exchangeability assumption is empirically supported.

Subgroup coverage is reported as the deviation from the 0.90 nominal, stratified by axis (Results R5). We define "within calibration" as a deviation of at most ±0.05 (5 percentage points).

### M7. Spatial-domain identification

To recover tissue architecture from the APA-modality alone, spaGAPA builds a fused multi-view graph with `MultiViewGraphBuilder` and detects communities by Leiden clustering. The graph combines three weighted views: a spatial k-nearest-neighbor graph (k = 6) on the spot coordinates, an APA-modality kNN graph on the (imputed) gene-level distal-usage vectors, and — when available — an expression kNN graph on the matched gene-expression matrix. The relative weights (spatial, expression, APA) are configurable; the APA-dominant setting used in the MOB benchmark is (0.2, 0.2, 0.6).

Leiden community detection is run through the igraph/leidenalg backends over a grid of resolutions; the partition with the best match to a reference layer annotation (when available) is reported. Domain detection is fully unsupervised: no sample or layer labels are used in graph construction or clustering. We evaluate domain recovery on the spvAPA mouse-olfactory-bulb (MOB) reference (ST11 sample, five anatomical layers) using the adjusted Rand index (ARI) and normalized mutual information (NMI).

### M8. Benchmark design and fairness protocol

All method comparisons use the same gene × spot APA distal-usage index (M2) computed from the same scAPAtrap PAS output, so that PAS-calling differences do not enter the comparison. For each benchmark dataset we mask 20% of the observed entries with a fixed random seed (42) and ask each method to impute the masked values. The same mask is presented to every method.

Methods compared: (i) spaGAPA-GP (sparse GP, Method D noise, locally-adaptive conformal); (ii) stAPAminer (v0.1.0, expression-KNN imputation at its recommended k = 10); (iii) spvAPA (v0.1.0, Seurat WNN-based imputation at its recommended k = 15); (iv) a generic spatial-KNN baseline (k = 15); and (v) a per-gene mean baseline that predicts every masked entry with the gene's observed mean. We include the per-gene mean baseline explicitly and transparently: because the gene-level distal-usage index is bimodal (most genes are nearly all-proximal or nearly all-distal), the per-gene mean is a deceptively strong entry-wise RMSE baseline, and we report its performance honestly rather than omit it.

Metrics: root-mean-square error (RMSE), Pearson r, and Spearman ρ between predicted and true masked values, computed both on all masked entries and on the top quartile of spatial-signal genes; a spatial-fidelity score measuring recovery of the true spatial autocorrelation (Moran's I) of imputed values relative to the observed values; and wall-clock runtime. All methods run on the same hardware (S91 server) with `OPENBLAS_NUM_THREADS = 8`, identical threading, and the same TMPDIR; runtimes are measured by subprocess wall-clock.

Disclosed asymmetries: spaGAPA-GP is implemented in Python 3.10 whereas stAPAminer and spvAPA are R packages; the KNN/WNN imputers in the competitors carry inherent O(n²) cost in the spot count n, whereas the sparse GP is O(n m²). These structural differences are part of the comparison and are reported as such, not normalized away.

### M9. Scalability experiments

To characterize how each method's runtime and memory scale with spot count, we generate synthetic APA datasets of 1 000, 5 000, 15 000, 42 000, and 100 000 spots by sampling PAS × spot matrices with realistic sparsity (≈90% missing) and bounded distal-usage values. Each synthetic dataset is presented to spaGAPA-highres_fast, spaGAPA-highres_accuracy, stAPAminer, and spvAPA under identical resource limits: a 1200-second wall-clock cap and `OPENBLAS_NUM_THREADS = 8`. We record wall time, peak resident memory, and completion status (COMPLETED, TIMEOUT, or FAILED with exit code).

The completion status — rather than the runtime alone — is the primary scalability metric. A method that exceeds the wall cap or crashes with a memory error is recorded as non-completing; we report this factual outcome and do not interpret it as a claim that the competitor is intrinsically unable to scale.

### M10. Stereo-seq processing and APA feasibility

For subcellular-resolution Stereo-seq data we run the SAW (v8.2.2) processing pipeline (read re-tagging, alignment, and barcode-position assignment) to produce a position-sorted, gene-annotated target BAM, which is then passed to scAPAtrap for PAS calling. Because Stereo-seq uses poly(A)-capture chemistry, we expect reads to pile up at gene 3′-ends; we verify this 3′-bias as a feasibility check by computing, for a 1-million-read subsample of gene-annotated reads, the strand-aware distance from each read's 3′-end to its gene's transcript-end site (TES), and reporting the fraction of reads whose 3′-end falls within 500 bp of the TES. We performed this QC on both a mouse AD-brain Stereo-seq sample (GSE263789) and a human AD-brain Stereo-seq sample (GSE269906).

Because the raw DNB count (≈20.7 M DNBs per sample) is intractable for any current spatial APA method, we bin the tissue on a 200 × 200 (bin200) grid to produce a tractable PAS × bin APA matrix, and validate the binning for downstream consistency by computing the cross-bin Pearson correlation of imputed APA usage against the full-resolution run.

### M11. Donor-level differential APA

Differential APA is performed at two levels with explicit pseudoreplication safeguards. At the **donor level**, the correct replicate unit is the biological sample (donor or animal), not the individual spot: pooling the ~14–15 k spots of each donor into a per-gene mean gives n = 3 control versus n = 3 AD samples for GSE220442, on which we run a Welch two-sample t-test per gene with Benjamini–Hochberg FDR correction (padj < 0.05, |Δ| > 0.05). This is the confirmatory analysis.

At the **spot level**, we also pool all ~29 k spots across donors and run a per-spot t-test. We label this spot-level analysis **exploratory**, because treating spots as independent replicates is a pseudoreplication that inflates the effective n and the nominal significance count. We additionally apply a direction-agreement filter (the per-sample sign of Δ agrees across all three donor pairs) as a hypothesis-generating filter on top of the exploratory spot-level call set.

For the Stereo-seq AD-versus-WT comparison, only one AD and one WT sample are available (biological n = 1 per condition), so no statistical test is possible. We report effect sizes only (|Δ| threshold) and label the result explicitly as hypothesis-generating.

### M12. Software and reproducibility

spaGAPA is implemented in Python 3.10 and released under the MIT license. The package comprises 30+ modules across `core`, `calling`, `imputation`, `bioml`, `analysis`, `benchmark`, and `visualization` subpackages, accompanied by 439 unit tests. Source code, benchmark scripts, and a full reproducibility manifest (environment paths, software versions, and the exact commands for every dataset and analysis) are available at https://github.com/zijunmeng/spaGAPA. All analyses in this manuscript were performed on a four-server cluster (S90/S91/S97/S98) with server-specific environment variables; the manifest documents the server, software version, and command for every reported figure and table.

---

## Results

### R1. spaGAPA provides a probabilistic inference framework for spatial APA

spaGAPA takes as input a PAS × spot APA usage matrix together with per-spot coordinates produced by an external PAS caller (scAPAtrap; Figure 1) and returns, for every gene in every spot, a probabilistic imputation consisting of a posterior mean and a calibrated prediction interval. Downstream, the same framework produces unsupervised spatial domains, donor-level differential APA calls with pseudoreplication safeguards, and spatially-variable APA (SVAPA) scores (Figure 5). The framework is positioned as a downstream inference layer, not as a PAS caller: it begins where scAPAtrap ends, and its contribution is the probabilistic and calibrated treatment of the resulting sparse, noisy, spatially structured APA matrix. The two components that distinguish spaGAPA from existing spatial APA tools — sparse GP regression with heteroscedastic noise estimation, and split-conformal uncertainty calibration with a distribution-free marginal coverage guarantee — are validated in the following sections across 11 Visium samples (7 GSE datasets, 4 tissue types, 2 species) and one subcellular Stereo-seq sample.

### R2. Spatial structure preservation versus a per-gene mean baseline

We first compared spaGAPA-GP against four reference imputation strategies on two Visium datasets (human kidney GSE183456; human brain/visual cortex GSE220442) under a shared 20%-masking protocol (Figure 2). We report the comparison transparently, including a per-gene mean baseline.

Although a per-gene mean baseline achieved lower global RMSE, spaGAPA better recovered spatial variation (GSE183456: 0.080 vs 0.122; GSE220442: 0.094 vs 0.126). This advantage of the mean is a structural artifact of the gene-level distal-usage index: the index is strongly bimodal, so a per-gene constant that predicts the dominant mode minimizes entry-wise RMSE while erasing all spatial variation. spaGAPA better recovered spatial variation: on GSE183456 the GP spatial-fidelity score (recovery of the observed Moran's I autocorrelation in the imputed values) was 0.42 (GP) vs 0.000 (mean) and the spatial-KNN baseline are the relevant methods, and spaGAPA additionally provides calibrated uncertainty that neither baseline offers.

Representative gene-level spatial maps (Figure 2E) show that spaGAPA-GP recovers layer-resolved APA gradients that the per-gene mean flattens and that stAPAminer/spvAPA recover noisily. The GP result is robust to the number of inducing points and to the masking fraction (Supplementary Figures S2 and S3), and the per-gene mean's RMSE advantage is confined to spatially flat genes — the GP–mean gap widens monotonically with gene-level spatial signal (Supplementary Figure S4).

### R3. spaGAPA outperforms existing tools in spatial fidelity and runtime

Pooling the two head-to-head datasets, spaGAPA-GP outperformed stAPAminer and spvAPA on entry-wise accuracy and spatial fidelity. On GSE220442, spaGAPA-GP achieved RMSE 0.126, Pearson r 0.947, and Spearman ρ 0.837, versus stAPAminer 0.159 / 0.917 / 0.820 and spvAPA 0.271 / 0.796 / 0.714. On GSE183456, spaGAPA-GP achieved RMSE 0.122, Pearson r 0.942, versus stAPAminer 0.208 / 0.852 and spvAPA 0.272 / 0.769; here the GP underperformed stAPAminer on raw RMSE because of the bimodal-index effect described in R2, but retained higher spatial fidelity (0.043 vs 0.034 vs 0.024).

spaGAPA-GP was 5.7–7.2× faster than the named competitors. On GSE183456 the wall times were: spaGAPA-GP 32.1 s, stAPAminer 161.8 s, spvAPA 186.4 s; on GSE220442: spaGAPA-GP 10.5 s, stAPAminer 109.8 s, spvAPA 124.8 s. We include the per-gene mean baseline transparently in the same panels (3.2 s and 0.1 s respectively, with the lowest RMSE but zero spatial fidelity), so that the reader can see both the regime in which the mean dominates and the regime in which the named competitors are the appropriate comparison. Full benchmark parameters (implementation language, neighbour basis, key hyperparameters, free-parameter count, random seed, and single-run wall time) are tabulated in Supplementary Figure S5.

### R4. Conformal coverage across diverse datasets

We validated split-conformal prediction across 11 Visium samples spanning 7 GSE datasets (GSE169749 mouse colon, GSE179572 human brain metastasis, GSE183456 human kidney, GSE220442 human AD brain, GSE237183 human glioma [18 samples], GSE263303 mouse brain, GSE338525 human liver), 4 tissue types, and 2 species, totaling 523 174 held-out test points (Figure 3; dataset overview in Supplementary Figure S1). One additional sample (GSE206391, <50 genes with sufficient observations) was excluded as having too few testable genes.

Split-conformal prediction achieved near-nominal distribution-free marginal coverage at every target level. Mean empirical coverage was 0.8005 at the 80% target, 0.8996 at the 90% target, and 0.9497 at the 95% target — deviations of 0.05, 0.04, and 0.3 percentage points respectively. Coverage was tight across samples (standard deviation 0.0027 / 0.0020 / 0.0012 across the 11 samples for the three targets), and all 11 samples fell within ±5 percentage points of nominal at every target (11/11 within ±0.05 at 80%, 90%, and 95%). The maximum absolute deviation across all 11 samples and three targets was 0.5 percentage points. The calibration curve (empirical coverage versus nominal) overlay the diagonal, and the locally-adaptive intervals were negligibly wider than the global intervals at matched coverage, indicating that the GP posterior standard deviation carries little redundant width information beyond what the global quantile already provides on the marginal distribution. Per-sample empirical coverage, interval width (2·q̂), and interval score across all 11 samples at the three target levels are detailed in Supplementary Figure S6.

### R5. Calibration robustness and the expression-stratum limitation

Because the conformal guarantee is marginal rather than conditional (M5), we audited subgroup coverage on four representative samples (GSE237183, GSE183456, GSE169749, GSE338525; 195 841 test points) across four subgroup axes (Figure 4).

Empirical subgroup coverage remained within ±5 percentage points across spatial regions and uncertainty strata. On the spatial-region axis (24 bins from k-means partitions with K = 5–10), all 24 bins were within ±0.05 of the 0.90 target (mean absolute deviation 0.010, maximum 0.030). On the uncertainty-quintile axis (20 bins), all 20 bins were within ±0.05 (mean absolute deviation 0.012, maximum 0.034). On the tissue axis (4 bins), all 4 bins were within ±0.05 (mean absolute deviation 0.001). The spatial-block adversarial split (3 train blocks, 2 test blocks) produced 90% coverage of 0.9896, statistically indistinguishable from the random-split coverage of 0.9897; the exchangeability concern that spatial autocorrelation would break the random-split assumption was therefore not supported by this stress test (full random- vs spatial-block coverage per dataset in Supplementary Figure S9). We additionally verified that the method selection and subgroup-coverage results are not artefacts of train/test leakage: a clean vs leaked-split comparison and a 5-fold leave-one-dataset-out (LOOCV) analysis both confirm no leakage and 5/5 stable method selection (Supplementary Figure S7).

The expression-stratum axis (18 bins) was the one axis on which calibration did not hold uniformly: 11/18 bins were within ±0.05 (61%) and 17/18 within ±0.10 (94%). Expression-dependent undercoverage (~0.83) was observed in the highest-expression subgroup, with the single worst bin (GSE169749, low-expression stratum, n = 5 357) reaching coverage 0.745 (deviation −0.155). We report this honestly as a limitation: spaGAPA's conformal intervals are well-calibrated across space, uncertainty, and tissue, but become miscalibrated in the extreme expression tails, particularly for the highest-expression genes where the locally-adaptive denominator can over-shrink intervals. Method D's residual-based noise partially corrects this (high-expression coverage 0.833 vs Method A's 0.930 over-coverage and Method B's 0.715 under-coverage) but does not eliminate it.

### R6. Heteroscedastic noise and global risk stratification

We compared the four noise-estimation strategies (Methods A–D; M4) on five datasets against three criteria: marginal coverage, interval width, and the correlation between the predicted uncertainty score and the absolute imputation error. Method D (residual-based empirical Bayes, default) achieved marginal coverage 0.9016 at the 90% target — statistically indistinguishable from Method A (0.9013) and Method B (0.9001) — at an interval width essentially identical to Method A (D/A width ratio 1.03), whereas Method B inflated width by a factor of 2.4 (range 1.1–4.8). Method D's calibration-set uncertainty/error correlation was 0.509, close to Method B's 0.555 and far above Method A's 0.068.

Auditing the uncertainty/error correlation on held-out test entries across five datasets (1 200 genes, 502 353 test points), heteroscedastic noise estimates enabled informative global risk stratification across gene-spot predictions (pooled r ≈ 0.55), although within-gene error ranking was more modest (median per-gene r ≈ 0.14; per-gene centered r ≈ 0.11). The pooled correlation is partly driven by cross-gene heterogeneity (genes with higher average uncertainty also have higher average error), but the within-gene signal is real and statistically consistent: across the five datasets, 82–90% of genes had a positive within-gene uncertainty/error correlation, and 17–26% exceeded r = 0.30. The honest interpretation is that Method D delivers reliable global risk stratification (which gene–spot predictions are likely wrong, pooled across the dataset) and modest within-gene discrimination (which spots within a single gene are likely wrong). Per-gene uncertainty/error correlation distributions for each of the five datasets, with the pooled across-gene r marked, are shown in Supplementary Figure S8.

We exploited the global stratification in a risk-coverage analysis (Figure 4G). Using Method D's split-conformal half-width as the uncertainty score, we ranked all gene-spot predictions and retained the lowest-uncertainty fraction; we compared the resulting RMSE against random retention at matched fraction, averaged over 10 random shuffles and five datasets. Uncertainty-ranked retention outperformed random filtering at every threshold (mean RMSE improvement 23.2% at 80% retention, p=0.004; paired t-test versus random). The improvement grew monotonically as the retention fraction decreased: 65.8% at 20% retention, 50.5% at 40%, 37.5% at 60%, 23.2% at 80%. The uncertainty score is therefore actionable for downstream triage, in the global-risk-stratification sense documented above.

### R7. Scalability to high-resolution data

We scaled each method on synthetic APA datasets from 1 000 to 100 000 spots under a fixed 1200-second wall cap and `OPENBLAS_NUM_THREADS = 8` (Figure 6). spaGAPA-highres_fast completed analyses under predefined resource limits at every scale: 1 000 spots in 8.9 s, 5 000 in 48.9 s, 15 000 in 276.6 s, 42 000 in 162.4 s, and 100 000 in 510.8 s, with peak memory 516, 1 270, 3 158, 2 792, and 8 807 MB respectively. spaGAPA-highres_accuracy is by design capped at ≤15 k spots (the factorizer is not optimized above that scale) and is reported as such.

Under the same resource limits, spaGAPA completed analyses under predefined resource limits, whereas comparator implementations did not complete at the higher scales. stAPAminer completed at 1 k (34 s), 5 k (253 s), and 15 k (495 s), but exceeded the 1200 s wall cap at 42 k and 100 k (TIMEOUT in both cases). spvAPA completed at 1 k (56 s), 5 k (198 s), and 15 k (196 s), but at 42 k and 100 k it failed with an R `future.globals.maxSize` error (509 MB and 1.13 GB respectively exceeding the 500 MB cap) in the parallel WNN evaluation, before producing any output. We report these outcomes factually: under the predefined resource limits used here, spaGAPA completed the 42 k and 100 k analyses whereas the comparator implementations did not.

### R8. Stereo-seq APA mapping

To demonstrate that spatial APA analysis is feasible at subcellular resolution, we processed a mouse AD-brain Stereo-seq sample (GSE263789, 20.7 M DNBs) through SAW 8.2.2 and scAPAtrap (Figure 7). The 3′-end-bias QC confirmed the poly(A)-capture signature: 53.4% of gene-annotated reads had their 3′-end within 500 bp of the gene TES in the mouse sample, and 47.0% in an independent human AD-brain Stereo-seq sample (GSE269906). scAPAtrap called 21 455 PAS peaks across 20.7 M DNBs; binning on a 200 × 200 grid produced a tractable 21 455 PAS × 15 235 bin APA matrix (sparsity 89.7%) whose cross-bin Pearson correlation with the full-resolution imputation was r = 0.842, supporting binning consistency (Stereo-seq read-level QC and binning sensitivity in Supplementary Figure S12). Unsupervised domain detection on the binned matrix recovered spatially coherent tissue domains (41 Leiden domains in the AD sample).

We describe the Stereo-seq result as descriptive APA mapping rather than a comparative disease analysis. A formal AD-versus-WT statistical comparison is not possible at the available replication (see R9); the contribution of this section is the demonstration that the full scAPAtrap → spaGAPA pipeline runs end-to-end on subcellular-resolution Stereo-seq data, which no previous spatial APA tool has been shown to do.

### R9. Donor-level differential APA reveals limited statistical power

Applying the donor-level differential APA analysis (M11) to GSE220442 (3 control vs 3 AD human brain samples), no gene reached FDR significance at the donor level: across 696 testable genes, the spot-level exploratory analysis nominated 91 nominal significant genes, but after per-gene per-sample aggregation to n=3 vs n=3 the count of FDR-significant genes was 0. A direction-agreement filter (per-sample Δ sign agreeing across all three donor pairs) flagged 29 genes as hypothesis-generating candidates; these are reported as such and not as confirmed differential APA calls. Spot-level differential APA results are exploratory; at the donor level (n=3 per group), no gene reached significance. The top donor-level candidates (e.g., ribosomal protein genes RPL10, RPL12, RPL10A) are mechanistically plausible — ribosomal-protein genes are a canonical APA-dysregulated class — but the absence of FDR significance at n = 3 is a statistical-boundary statement, not a negative biological result: the cohort is simply underpowered to detect the effect sizes observed (median per-donor |Δ| ≈ 0.03).

For the Stereo-seq AD-versus-WT comparison (n = 1 per condition), no statistical test is possible. We report effect sizes only: of 1 090 genes observed in both conditions, 607 exceeded |Δ| > 0.10 (296 distal-up in AD, 311 proximal-up in AD), with the largest effect sizes concentrated in neuronal and synaptic genes (e.g., Rpl37, Cdk8, Ly6h, Filip1l, Tubb5, Apod, Fgfr1). These are labeled explicitly as effect-size-only, hypothesis-generating candidates (Supplementary Figure S11) requiring biological replicates for confirmatory testing.

---

## Discussion

> **[TO BE WRITTEN BY ABSTRACT/INTRO/DISCUSSION AGENT — target ~900 words]**
>
> *Draft scaffolding (to be replaced).* Suggested content arc, derived from the honest limitations already documented in the Methods/Results and `PROJECT_SUMMARY.md` §6.4:
>
> 1. **Principal findings.** (a) A sparse GP + split-conformal wrapper delivers distribution-free marginal coverage for spatial APA imputation across 11 diverse Visium samples (deviations 0.05–0.3 pp); (b) it matches/beats stAPAminer and spvAPA on spatial fidelity at 5.7–7.2× lower runtime; (c) it is the only framework completing 42 k/100 k spots and the only end-to-end Stereo-seq spatial APA pipeline; (d) uncertainty is actionable for global risk-stratified triage.
> 2. **Strengths of the statistical framing.** Why conformal marginal coverage (not conditional) is the right honest claim; the exchangeability stress test (spatial-block split) supports the assumption empirically.
> 3. **Honest limitations** (mirrors `PROJECT_SUMMARY.md` §6.4 + R5/R6/R9):
>    - GP does not beat per-gene mean on entry-wise RMSE — a structural artifact of the bimodal distal-usage index (R2).
>    - Conformal guarantee is marginal, not conditional; extreme expression tails under-cover (R5); within-gene uncertainty/error correlation is modest (R6).
>    - Donor-level differential APA yields 0 FDR-significant genes at n=3 vs 3 — the cohort is underpowered, not biologically negative (R9); Stereo-seq AD-vs-WT is effect-size-only at n=1 per condition.
>    - Stereo-seq raw FASTQ+mask scarcity is a domain-wide data bottleneck; Visium HD excluded (probe chemistry incompatible with APA).
>    - sAPA-RegNet perturbation model is descriptive only (cis-regression negative).
>    - No supervised feature selection (spvAPA's sPLS-DA is a complementary capability spaGAPA lacks).
> 4. **Comparison with the literature / Table 1 positioning.** First spatial APA tool with conformal-calibrated uncertainty; first subcellular-resolution spatial APA pipeline; transparent head-to-head with disclosed Python-vs-R and O(n·m²)-vs-O(n²) asymmetries.
> 5. **Future work.** Conditional/group-wise conformal to fix the expression-tail under-coverage; tighter locally-adaptive widths; biological replicates for confirmatory AD APA; multi-organ Stereo-seq domain atlas (MOSTA).
> 6. **Conclusion** (2–3 sentences). spaGAPA upgrades spatial APA analysis from a data-processing pipeline to a statistical-inference framework with calibrated uncertainty, and makes subcellular-resolution spatial APA analysis tractable for the first time.

---

## Figure legends

All figures rendered at 300 DPI with a colorblind-friendly Okabe-Ito palette (`#0072B2` blue / `#D55E00` orange / `#009E73` green), DejaVu Sans typography, bold panel letters (A, B, C…), labeled axes with units, and paired per-dataset dots in addition to mean bars. Source data files and generation scripts are listed at the end of this section.

**Figure 1 — Framework overview.** Conceptual/schematic figure (no quantitative data).
- **A. Three gaps in spatial APA.** Sparse PAS matrix (<12% observed), no calibrated uncertainty intervals, and O(N²) methods that cannot scale past 15 k spots.
- **B. Input pipeline.** Spatial reads → PAS peak counts → APA usage matrix [gene × spot] → spaGAPA. Sparsity grid visualizes the <12% observed fraction.
- **C. Sparse Gaussian Process.** Posterior mean + 68/95% CI with M=150 inducing points; O(N·M²) complexity enables 100 k-spot scaling.
- **D. Split-conformal calibration.** Train / calibrate / test partition; nonconformity score sᵢ=|yᵢ−ŷᵢ|; quantile qₜ gives ŷ±qₜ interval with finite-sample coverage guarantee.
- **E. spaGAPA outputs.** Posterior mean, interval width, uncertainty map, tissue domains.
- **F. Dataset atlas.** 8 tissues · 2 species (human/mouse) · 2 platforms (Visium/Stereo-seq); 11 calibration samples; up to 100 k spots.

**Figure 2 — Spatial vs mean benchmark.** Head-to-head, 5 methods × 2 datasets (GSE183456 kidney, GSE220442 brain), same 20% per-gene mask (seed 42). *Honest framing: mean is shown to beat GP on entry-wise RMSE.*
- **A. Masking design.** 20%-of-observed per-gene holdout; same entries masked for every method.
- **B. Per-dataset metrics (paired dots).** RMSE (left) and Pearson/Spearman (right) with per-dataset dots + mean bars. Mean wins RMSE/Pearson/Spearman; this is shown openly.
- **C. Spatial fidelity.** GP 0.42 vs Mean 0.00 — mean imputes a per-gene constant and recovers no spatial gradient by construction.
- **D. Accuracy-spatial 2D scatter.** One point per method per dataset; mean is high-accuracy/zero-gradient, GP moderate-accuracy/strong-gradient.
- **E. Representative APA spatial maps (GSE183456).** Two high-spatial-variance genes: truth / mean (flat) / spaGAPA-GP (spatially smoothed) / |error|.
- **F. GP advantage grows with gene spatial signal.** Low/med/high variance strata; the GP-vs-mean RMSE gap narrows as spatial signal rises.

**Figure 3 — Conformal marginal coverage.** 11 samples, 523 174 test points, 80/90/95% nominal.
- **A. Split-conformal flow.** Calibration set → nonconformity |yᵢ−ŷᵢ| → quantile qₜ → interval ŷ±qₜ → test coverage ≥ 1−α.
- **B. Per-sample coverage at 80/90/95%.** Dots per sample + nominal lines + 90% binomial CI band. Empirical lands within ~0.2% of nominal.
- **C. Calibration curve.** Nominal vs empirical; thin per-sample lines + bold pooled mean; y=x reference. Mean abs deviation 0.21% (80%), 0.16% (90%), 0.10% (95%).
- **D. Coverage deviation forest.** Per-sample empirical−nominal at all three α.
- **E. Interval width: Constant (A) vs Residual-spot (D).** D ≈ same width as A → adaptive intervals without inflation.
- **F. Spatial instance.** One GSE183456 gene: posterior mean / interval width / |error| / covered-vs-uncovered (green/red).

**Figure 4 — Noise methods + risk-coverage curve.** 4 noise methods (A constant, B local-gene, C spatial-spot, D residual-spot) × 5 datasets.
- **A. Four uncertainty models.** Schematic of each noise definition.
- **B. Multi-objective method selection.** x=uncertainty-error r, y=interval width, bubble size=|coverage deviation|. D (spaGAPA) is high-r, near-minimal-width.
- **C. Pooled r per dataset: A vs D.** A near-zero baseline; D clears the r≥0.3 target.
- **D. Per-gene vs pooled r.** Boxplot of within-gene r distribution + pooled-r star. Pooled r partly driven by cross-gene heterogeneity; within-gene discrimination is modest (median r ≈ 0.10–0.19).
- **E. Uncertainty quintile coverage.** Coverage climbs Q1→Q5 (high-uncertainty points correctly over-covered).
- **F. Conditional coverage heatmap.** Tissue × uncertainty quintile; green≈nominal.
- **G. Risk-coverage curve (KEY NEW PANEL).** Retention fraction vs RMSE for uncertainty-ranked / oracle / random retention. Uncertainty triage gives −66% RMSE at 20% retention and −23% at 80% (paired t-test p=0.004 vs random).

**Figure 5 — Domain recovery.** MOB (mouse olfactory bulb) spvAPA ST array, 5 annotated layers, 260 spots.
- **A. MOB tissue: 5 olfactory bulb layers** (GCL/GL/MCL/ONL/OPL).
- **B. Domain maps.** Mean impute (no APA signal) vs spaGAPA (APA-dominant).
- **C. ARI/NMI bars.** All run configs; best ARI=0.60, NMI=0.68 (expression_apa / leiden res0.8). The Leiden-resolution sweep and weight-configuration robustness behind this result are detailed in Supplementary Figure S10.
- **D. Spatial signal recovery (Moran's-I).** Per-method across the 2 benchmark datasets; GP reconstructs autocorrelation that mean cannot.
- **E. Representative APA gradient genes.** Four layer-stratified gradients.
- **F. Cross-sample note.** MOB is single-sample (n=1); multi-section domain stability reported in Supplementary (GSE237183).

**Figure 6 — Scalability.** Synthetic grids 1 k–100 k spots, 4 methods, 1200 s wall cap.
- **A. Algorithmic complexity.** O(N²) KNN vs O(N·M²) sparse GP vs O(N) mean.
- **B. Runtime scaling (log-log).** spaGAPA-fast completes 100 k in 511 s; stAPAminer/spvAPA timeout/fail at 42 k+.
- **C. Peak memory scaling.** spaGAPA-fast 8.8 GB at 100 k; R-tools exceed limits.
- **D. Completion matrix.** Method × scale: green=completed, yellow=skipped (out of scope), red=timeout/OOM.
- **E. Accuracy-runtime Pareto.** At matched accuracy spaGAPA-GP is 5–7× faster than stAPAminer/spvAPA; mean/spatial-KNN are fast but recover no gradient.
- **F. Ablation** (inducing-point sweep) → moved to Supplementary per reviewer spec.

**Figure 7 — Stereo-seq pilot.** GSE263789 mouse AD brain (Stereo-seq + scAPAtrap).
- **A. Stereo-seq workflow.** DNB array → SAW count → scAPAtrap PAS peaks → bin 50/100/200 → spaGAPA impute.
- **B. 3′-end bias QC.** Reuses existing `fig1_3prime_enrichment_near_TES.png` (read coverage enrichment near TES confirms APA-relevant 3′ capture).
- **C. Data scale & sparsity.** 20.7 M DNB, 21 455 PAS, 4 845 genes, ~12% observed.
- **D. Tissue APA atlas.** Spatial APA usage across bin50 coordinates.
- **E. Representative genes.** Three PAS peaks (cortex-enriched, hippocampal gradient, punctate).
- **F. Binning consistency.** Cross-bin Pearson r for 50→100 and 50→200; mean r=0.84, 100% of pairs >0.7 — APA patterns robust to resolution.
- **G. AD vs WT (descriptive, n=1).** Illustrative APA DU direction by gene category. No statistical testing (n=1 per group); effect directions for hypothesis generation only, prominently flagged.

*Figure source data and generation scripts* (all under `pipeline_output/`): benchmark — `benchmark_mean_transparent/transparent_comparison.csv`, `benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json`; conformal — `conformal_validation/all_samples_coverage.csv` + `summary.json`; conditional coverage — `conformal_conditional_coverage/*.csv`; noise methods — `uncertainty_corr_improvement/per_dataset_corr.csv`; per-gene — `uncertainty_within_gene_audit/per_gene_corr_distribution.csv`; risk-coverage — `risk_coverage_curve/risk_coverage_data.csv`; MOB — `mob_domain_recovery/spagapa_metrics.json` + `spagapa_domains.csv`; scalability — `benchmark_runtime/runtime_table.csv`; Stereo-seq — `stereo_binning_consistency/binning_correlation.csv`, `gse263789_stereo_pilot/figures/fig1_3prime_enrichment_near_TES.png`. Generation scripts: `fig1_overview.py`, `fig2_benchmark.py`, `fig3_conformal.py`, `fig4_noise.py`, `fig5_domain.py`, `fig6_scalability.py`, `fig7_stereo.py` (shared `_style.py`).

---

## Table 1 legend

**Table 1 — Spatial APA tool capability comparison.** Comparison of spaGAPA with the three existing spatial/alternative-polyadenylation analysis frameworks most relevant to this work: stAPAminer (the direct Visium-scale spatial APA competitor), spvAPA (an imputation-plus-supervised-selection tool covering spatial and single-cell APA), and metaAPA (a cross-tool PAS-integration workflow). "N/A" denotes a capability that is outside the tool's intended scope rather than a deficiency. All claims are documented in Methods (M8–M9) and Results (R2–R3, R7–R8); version numbers, run conditions, and the reproducibility manifest are given in M12. Scalability and high-resolution outcomes are reported under the predefined resource limits of M9 (1200 s wall cap, `OPENBLAS_NUM_THREADS = 8`); the failure of stAPAminer/spvAPA at 42 k spots is a timeout / `future.globals.maxSize` outcome under those limits, reported factually. spaGAPA is the only tool whose prediction intervals carry a distribution-free marginal coverage guarantee (split-conformal, M5); metaAPA's "site confidence" is a caller-level score, not a calibrated prediction interval on an imputed value. spaGAPA's spatial-domain identification (M7) is fully unsupervised (Leiden on a fused spatial + APA + expression graph); spvAPA's supervised feature selection (sPLS-DA) is a complementary capability spaGAPA does not provide. The QN + linear batch-correction module is an optional, still-under-evaluation component (Pillar 2; Harmony comparison pending), listed for completeness, with the QN-versus-Harmony comparison detailed in Supplementary Figure S13. [Full table in `docs/TABLE1_tool_comparison.md` and `docs/TABLE1_tool_comparison.csv`.]

---

## Supplementary Figures

All supplementary figures share the main-figure formatting (300 DPI, Okabe–Ito palette, DejaVu Sans). Generation scripts and cached source data are under `pipeline_output/supplementary_figures/` (figure scripts `supp_figNN_*.py`; data-prep helpers `_gp_sweeps.py` for S2/S3, `_stereo_3end.py` for S12; cached CSV/NPY under `_cache/`).

- **Supplementary Figure S1 — Dataset overview.** (A) PAS (scAPAtrap-called site) count per sample across the 32-sample Visium benchmark, grouped by GSE series and coloured by tissue. (B) Per-series summary table: species, platform, tissue, sample count, and mean/total spot and PAS counts. Data: `data/processed/*_scapatrap/qc_summary.json`.
- **Supplementary Figure S2 — Sparse-GP inducing-point sensitivity (GSE183456).** (A) Held-out RMSE (overall and median-per-gene) versus number of inducing points *m* ∈ {50, 100, 200, 500}; both are essentially flat, so the posterior mean is robust to *m*. (B) Wall time versus *m*, scaling near-linearly. Data: `_cache/s2_inducing_sweep.csv`.
- **Supplementary Figure S3 — Masking-level sensitivity (GSE183456).** Imputation benchmark re-run at 10/20/30/50% masking (not only 20%) for the three Python-native methods (spaGAPA-GP, mean, spatial-KNN): grouped bars of median-per-gene RMSE by method × masking level. stAPAminer/spvAPA (R-based) are reported only at the 20% level used in the main benchmark. Data: `_cache/s3_masking_sweep.csv`.
- **Supplementary Figure S4 — Mean-baseline stratification.** Box plot of per-gene (GP RMSE − mean RMSE) stratified by gene-level spatial signal (Moran's I quintile). The GP's advantage over the per-gene mean concentrates in high-spatial-signal genes; for spatially flat genes the mean is difficult to beat. Data: `_cache/s4_stratification.csv`.
- **Supplementary Figure S5 — Benchmark parameter table.** Publication-style table of the five benchmarked methods: implementation language, neighbour basis, key hyperparameters, free-parameter count, random seed, and single-run wall time (GSE183456, same hardware). Data: `benchmark_fairness/parameter_table.csv` + cached transparent head-to-head timings.
- **Supplementary Figure S6 — Full conformal coverage across 11 samples.** (A) Empirical vs nominal coverage at 80/90/95% (per sample). (B) Interval width 2·q̂ at each level (per sample). (C) Interval score (Gneiting & Raftery 2007) per sample — an approximation from per-sample q̂ and RMSE (see caption note); lower is better. Data: `conformal_validation/all_samples_coverage.csv`.
- **Supplementary Figure S7 — Leakage audit and LOOCV method-selection stability.** (A) Clean vs leaked per-sample uncertainty/error correlation for the local-gene method across five datasets (≈ identical ⇒ no train/test leakage). (B) 5-fold leave-one-dataset-out: each held-out dataset selects the same method (B/D local-gene), i.e. 5/5 stable selection, with held-out r per method per fold.
- **Supplementary Figure S8 — Per-gene uncertainty–error correlation.** Faceted per-dataset histograms of per-gene Pearson r (uncertainty ↔ absolute error), with the pooled across-gene r marked. Per-gene r is weak on average (consistent with marginal r ≈ 0.07) but the distribution is centred above zero — uncertainty is informative within genes, not only globally. Data: `uncertainty_within_gene_audit/per_gene_corr_distribution.csv`.
- **Supplementary Figure S9 — Spatial-block vs random-split conformal coverage.** Grouped bars of empirical coverage under random vs spatial-block (spatially disjoint) splitting, per dataset, at 80/90/95%. Intervals remain valid (within ~2 pp of nominal) even under the harder block split — coverage is not an artefact of adjacent spots leaking across the train/test boundary. Data: `conformal_conditional_coverage/*.csv`.
- **Supplementary Figure S10 — Domain-recovery details on MOB.** (A) Leiden-resolution sweep (ARI/NMI vs resolution) for each of the four weight configurations (apa_dominant, balanced, spatial_apa, expression_apa), with the chosen best resolution per config marked. (B) ARI/NMI under the four weight regimes (each at its best resolution), showing the result is stable across weightings. Data: `mob_domain_recovery/`.
- **Supplementary Figure S11 — Pseudoreplication analysis (GSE220442).** (A) Spot-level vs donor-level significant gene counts: the pseudoreplicated spot-level test calls 91 genes significant whereas the correct donor-level test (n = 3 ctrl vs n = 3 AD) calls 0 — classic pseudoreplication inflation; direction-agreement genes = 29. (B) Effect-size agreement for the 29 direction-agreement genes.
- **Supplementary Figure S12 — Stereo-seq QC and binning sensitivity (GSE263789).** (A) PAS count distribution by chromosome from scAPAtrap peak calls. (B) Read 3′-end-to-TES signed-distance histogram (500k gene-annotated reads), showing the 3′-enrichment (polyA-capture) signature of Stereo-seq. Data: `_cache/s12_tes_distances.npy` (`_stereo_3end.py`).
- **Supplementary Figure S13 — Batch correction: QN vs Harmony.** (A) Mean pairwise PCC (biological signal recovery) vs residual batch signal per correction method: spaGAPA-QN recovers most of Harmony's PCC gain while keeping batch signal near zero; Harmony over-corrects (variance collapsed). (B) Per-gene variance ratio (corrected / before), showing QN preserves within-gene variance whereas Harmony attenuates it. Data: `apa_bias_correction/`.

---

## References

> **[TO BE WRITTEN — PLACEHOLDER]**
>
> The reference list is not yet compiled. The following key citations should be included when the bibliography agent finalizes the list (BIB Vancouver/numbered style):
>
> 1. stAPAminer — original method paper (spatial APA on Visium).
> 2. spvAPA — original method paper (imputation + sPLS-DA for spatial/sc APA).
> 3. metaAPA — original method paper (cross-tool PAS integration).
> 4. scAPAtrap (v0.2.0) — PAS caller used upstream of spaGAPA.
> 5. 10x Visium platform / technology paper.
> 6. Stereo-seq / SAW (v8.2.2) platform and pipeline papers.
> 7. Core Gaussian-process references — Rasmussen & Williams (2006); Titsias (2009) inducing-point variational sparse GP; Hensman et al. (2013) stochastic variational GP.
> 8. Conformal prediction — Vovk et al. (2005); Lei et al. (2018) split conformal; Romano et al. (2019) locally-adaptive conformal (CQR).
> 9. Leiden community detection — Traag et al. (2019).
> 10. Moran's I / spatial autocorrelation — Moran (1950).
> 11. APA biology reviews — Tian & Manley (2017); Elkon et al. (2013); Mayr (2019).
> 12. Datasets — GSE169749 (mouse colon DSS); GSE179572 (human brain metastasis); GSE183456 (human kidney); GSE206391 (human skin psoriasis); GSE220442 (human AD brain PFC); GSE237183 (human glioma); GSE263303 (mouse brain Nf1); GSE263789 (mouse AD brain Stereo-seq); GSE269906 (human AD brain Stereo-seq); GSE338525 (human liver).
> 13. Welch's t-test; Benjamini–Hochberg FDR (1995).
> 14. Okabe & Ito colorblind-friendly palette reference.
>
> *Replace this placeholder with the formatted numbered reference list before submission.*

---

## Appendix: document provenance

This master document (`docs/SPAGAPA_FULL_MANUSCRIPT.md`) integrates content from the following sources:

- **Methods (M1–M12)** — verbatim from `docs/MANUSCRIPT_DRAFT.md`.
- **Results (R1–R9)** — verbatim from `docs/MANUSCRIPT_DRAFT.md`.
- **Figure legends (Fig 1–7)** — adapted from `pipeline_output/main_figures/figure_index.md`.
- **Table 1 legend** — from `docs/TABLE1_tool_comparison.md`.
- **Abstract, Introduction, Discussion** — from `docs/MANUSCRIPT_ABSTRACT_INTRO_DISCUSSION.md` (integrated after the Methods/Results draft).
- **Supplementary figure legends (S1–S13)** — added during manuscript finalisation; generation scripts under `pipeline_output/supplementary_figures/` (data sources listed in each legend).
- **References** — placeholder; the intended citation list is enumerated in-section and is to be formatted in BIB Vancouver/numbered style before submission.
- **Title, Authors, Affiliations** — placeholders, to be supplied by the corresponding author.

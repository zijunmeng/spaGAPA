# spaGAPA Project Summary

**Project name**: spaGAPA
**Full name**: spatial Gaussian-process and graph-aware APA analyzer
**Language**: Python
**Main environment**: conda env `spagapa`
**Current manuscript target**: Briefings in Bioinformatics (BIB)
**Last updated**: 2026-07-13

---

## 1. Executive Summary

spaGAPA is a Python toolkit for alternative polyadenylation (APA) analysis in spatial transcriptomics data. The project was initially designed around spatially aware APA validation and Gaussian process (GP) imputation. It has now evolved into a broader statistical and machine-learning framework for sparse spatial APA signals.

The current BIB-oriented thesis is:

> Spatial APA signals are sparse, spatially structured, uncertain, and highly dependent on tissue context. A strong spatial APA framework should not only call APA sites, but also validate spatial support, recover missing APA usage with calibrated uncertainty, separate value recovery from biological-domain discovery, and provide reproducible benchmarks on real spatial transcriptomics datasets.

spaGAPA addresses this through four connected layers:

1. Real-data APA input construction from Space Ranger BAM and scAPAtrap-compatible APA/PAS evidence.
2. Spatial APA validation, quantification, and GP-based imputation with uncertainty.
3. CPU-friendly BioML graph learning for biological-domain recovery.
4. Benchmark and visualization infrastructure for simulation, MOB/layer-like data, high-resolution pseudo-bin data, and real datasets.

The project is currently in the phase of converting algorithmic prototypes into manuscript-grade evidence, especially real-data benchmarks and high-resolution validation.

---

## 2. Scientific Motivation

APA changes the 3' end of transcripts and can alter transcript stability, localization, translation, and regulatory interactions. Existing single-cell and spatial transcriptomics APA tools are useful, but several problems remain under-addressed:

- APA measurements are sparse and dropout-prone.
- Spatial transcriptomics adds physical neighborhood structure that should be modeled explicitly.
- KNN-style imputation often lacks calibrated uncertainty.
- Biological domain recovery and APA value recovery are related but not identical tasks.
- High-resolution spatial transcriptomics introduces stronger sparsity and scalability constraints.
- Published tools often lack comprehensive, reproducible benchmark pipelines.

spaGAPA is designed around these gaps.

---

## 3. Core Contributions

### 3.1 Spatial APA data model and workflow

spaGAPA provides dedicated data structures and workflow utilities for spatial APA analysis:

- APA site representation.
- APA count/usage matrix handling.
- Spatial coordinate handling.
- Result writing.
- Integration with AnnData-like analysis ecosystems.
- CLI and pipeline presets.

This makes spaGAPA a package rather than a one-off benchmark script.

### 3.2 Spatial validation before downstream modeling

APA sites are not treated as independent molecular events detached from tissue geometry. spaGAPA includes spatial validation and filtering modules that evaluate:

- read support
- spot support
- spatial support
- neighbor consistency
- spatial autocorrelation

This layer is important because noisy APA calls can otherwise dominate imputation and downstream SVAPA detection.

### 3.3 GP-based APA value recovery with uncertainty

spaGAPA implements Gaussian process-based imputation for APA usage:

- spatial kernels
- sparse/fast variants
- batch imputation
- posterior standard deviation
- uncertainty-aware downstream analysis

The main advantage of GP over simple KNN is not just lower error in favorable settings. The key statistical advantage is that GP provides uncertainty estimates, which can be used for:

- confidence filtering
- uncertainty-aware SVAPA ranking
- interpretation of poorly supported spatial regions
- calibration analysis

### 3.4 Systematic exploration of expression-informed GP

The project explicitly tested whether expression-informed GP should become the central method. Explored variants include:

- additive expression/spatial GP
- radial GP
- product-kernel GP
- adaptive expression weighting
- layer-local expression kernel

Benchmark conclusion:

- Expression-aware GP can improve some domain metrics in selected settings.
- It did not robustly dominate spatial GP or graph/KNN baselines.
- It sometimes risks worse calibration or runtime.

Strategic decision:

- Expression-informed GP remains an experimental/ablation component.
- It is not the main BIB method.
- Biological-domain recovery is handled by a separate BioML graph route.

This is important because it prevents the method from being overclaimed.

### 3.5 Decoupled BioML route for biological-domain consistency

The major algorithmic evolution is the BioML route:

- CPU-only.
- No GPU requirement.
- No deep learning dependency.
- Multi-view graph construction from:
  - spatial coordinates
  - expression matrix
  - APA usage matrix
- Graph-regularized factor/domain recovery.

The conceptual shift is:

- GP is responsible for APA value recovery and uncertainty.
- BioML graph learning is responsible for biological-domain consistency.

This decoupling is currently the most defensible route toward comprehensive performance across RMSE, layer/domain ARI/NMI, and runtime.

### 3.6 High-resolution spatial transcriptomics support

spaGAPA includes a high-resolution route designed for sparse pseudo-bin or high-resolution spatial data:

- `highres_accuracy`
- `highres_fast`
- `auto` preset dispatch
- pseudo-bin benchmark suites
- multiscale and graph-based domain handling

The current high-resolution candidate default is:

- `gp_blend = 0.1`
- graph weight `spatial = 0.1`
- graph weight `expression = 0.7`
- graph weight `apa = 0.2`

The goal is to make spaGAPA useful not only for classic Visium-scale data, but also for emerging platforms such as Stereo-seq and Visium HD where appropriate APA evidence exists.

Important technical caveat:

- Visium HD FFPE/probe-based data are not automatically suitable for primary APA/PAS calling.
- Stereo-seq and poly(A)-compatible high-resolution data are more promising for true high-resolution APA discovery.

---

## 4. Package Architecture

Current main modules:

| Module | Role |
|--------|------|
| `spagapa.core` | APA dataset and APA site data structures |
| `spagapa.io` | BAM, coordinate, matrix, and result I/O |
| `spagapa.spatial` | spatial graph construction |
| `spagapa.calling` | spatial validation and quality filtering |
| `spagapa.imputation` | GP, sparse GP, fast GP, expression GP prototypes |
| `spagapa.quantification` | RUD, PDUI, WUL, PAI, QC |
| `spagapa.analysis` | differential APA, domain detection, SVAPA |
| `spagapa.visualization` | spatial/statistical/QC plots |
| `spagapa.benchmark` | simulation and real benchmark infrastructure |
| `spagapa.bioml` | CPU-friendly multi-view graph/domain route |
| `spagapa.pipeline` | end-to-end workflow |
| `spagapa.cli` | command-line interface |
| `spagapa.presets` | analysis preset definitions |

User-facing presets:

| Preset | Intended use |
|--------|--------------|
| `auto` | default; dispatch based on matrix size/sparsity |
| `standard` | conventional low/moderate-resolution spatial transcriptomics |
| `highres_accuracy` | high-resolution route prioritizing domain stability and RMSE |
| `highres_fast` | high-resolution route prioritizing speed |

---

## 5. Competitor and Reference Tools

The main reference packages are:

- `scAPAtrap`
- `stAPAminer`
- `metaAPA`

### 5.1 scAPAtrap

Role:

- Strong APA/PAS calling front end from BAM.
- Used by spaGAPA as one major route for real APA evidence generation.

Limitations relative to spaGAPA:

- Not designed as a spatial APA analysis framework.
- Does not provide spaGAPA's spatial validation, GP uncertainty, BioML domain route, or high-resolution benchmark framework.

### 5.2 stAPAminer

Role:

- Primary spatial APA competitor.
- Uses expression-neighborhood/KNN-like logic for APA imputation and spatial analysis.

Strength:

- Expression-informed neighborhood structure can produce strong biological-domain consistency.
- Runtime can be faster than GP-heavy methods.

Limitations:

- No GP posterior uncertainty.
- KNN imputation has limited calibration.
- Less explicit separation between APA value recovery and biological-domain discovery.
- Benchmark and high-resolution scalability are less systematic.

spaGAPA response:

- Do not try to beat stAPAminer-like KNN on every domain metric using GP alone.
- Use GP where GP is strongest.
- Use BioML graph route for biological consistency.

### 5.3 metaAPA

Role:

- Workflow-oriented APA/meta-analysis framework.

Relationship to spaGAPA:

- More complementary than directly competitive.
- spaGAPA focuses on spatial APA validation, imputation, uncertainty, domains, and high-resolution spatial benchmarking.

---

## 6. Benchmark Status

### 6.1 Simulation and MOB/layer-like benchmarks

Completed:

- Spatial dropout masks.
- Random masks.
- Spatial block masks.
- Ring/sector masks.
- Layer-aware masks.
- Multi-seed benchmark.
- Runtime and memory tracking.
- Biological consistency metrics:
  - ARI
  - NMI
- Uncertainty metrics:
  - interval coverage
  - uncertainty-error correlation
  - uncertainty-filtered accuracy

Main conclusion:

- GP is strong for APA value recovery and uncertainty.
- Biological-domain recovery requires graph/expression-aware structure.
- Decoupling value recovery and domain recovery is necessary.

### 6.2 Expression GP sweeps

Completed:

- spatial GP
- additive expression GP
- product kernel GP
- adaptive GP
- layer-local GP

Main conclusion:

- Expression GP is not robust enough to be the central method.
- These experiments are still useful as ablation evidence.

### 6.3 High-resolution benchmark

Completed:

- pseudo-bin generation
- 2x/4x/8x tests
- multi-seed smoke suites
- `highres_accuracy`
- `highres_fast`
- BioML decoupled route

Representative result pattern:

- `highres_fast` substantially improves runtime.
- `highres_accuracy` preserves better value recovery.
- BioML graph route improves or stabilizes biological-domain metrics in selected high-resolution settings.

### 6.4 Real-data benchmark

Completed:

- GSE179572/GSM5420751 real Space Ranger + scAPAtrap dataset.
- GSE183456/GSM6047774 real Space Ranger + scAPAtrap dataset.
- GSE179572 high-resolution smoke.
- GSE179572 external validation smoke.

In progress:

- GSE237183 multi-sample Space Ranger processing.
- GSE220442 SRA split and Space Ranger route.
- Additional brain/high-resolution candidates.

---

## 7. Real Dataset Status

### 7.1 GSE179572 / GSM5420751

Status:

- Space Ranger completed.
- BAM barcode/UMI tags verified.
- scAPAtrap completed.
- Expression matrix exported.
- Processed dataset available.

Use:

- First real-data benchmark.
- High-resolution smoke validation.
- External validation smoke.
- Marker weak-label exploration.

Current limitation:

- No expert pathology ROI label yet.

### 7.2 GSE183456 / GSM6047774

Status:

- Space Ranger completed.
- scAPAtrap completed.
- Processed dataset available.

Use:

- Strong candidate for biological-domain benchmark.
- Potentially better than GSE179572 for domain label construction.

### 7.3 GSE237183

Status:

- FASTQ complete.
- Dry-runs passed.
- Space Ranger formal counts submitted.
- Most completed samples pass; one observed failure is image/fiducial alignment-related.

Use:

- Brain/Visium candidate for multi-real-data benchmark.

### 7.4 GSE220442

Status:

- Non-standard SRA technical read structure decoded.
- Correct read mapping found.
- Pilot Space Ranger dry-run completed.
- Full suite processing started.

Use:

- Important brain dataset candidate.
- Demonstrates spaGAPA workflow can handle complicated SRA formats.

### 7.5 High-resolution candidate search

Candidate documents:

- `docs/highres_st_geo_candidates.md`
- `docs/stereo_visiumhd_geo_candidates.md`

Priority:

- Stereo-seq first for true APA raw-read evidence.
- Visium HD for high-resolution expression/domain validation, with APA chemistry caveat.

---

## 8. Current Strengths

spaGAPA currently has several credible strengths for a BIB-level package/method paper:

1. Clear biological problem: spatial APA remodeling.
2. Statistical angle: sparse, uncertain, spatially structured APA usage.
3. Methodological components:
   - spatial validation
   - GP imputation
   - posterior uncertainty
   - BioML graph domain route
4. CPU-only design:
   - no GPU dependency
   - no deep learning requirement
5. Reproducible benchmark framework.
6. Real-data processing path from Space Ranger BAM to APA matrix.
7. Active high-resolution strategy.
8. Honest ablation story: expression-informed GP was tested and not overclaimed.

---

## 9. Current Weaknesses

The project is not yet manuscript-complete. Remaining weaknesses are:

1. Biological gold-standard labels are still limited.
2. GP runtime can be slower than KNN baselines.
3. Domain ARI/NMI is not always superior unless the BioML route is used.
4. Real-data benchmarks are still being expanded.
5. High-resolution APA evidence requires careful dataset selection.
6. Some real datasets have difficult image alignment or SRA read structures.
7. Documentation still needs to be consolidated around current presets and workflows.

These are solvable, but they must be handled before submission.

---

## 10. Publication Strategy

### 10.1 BIB positioning

The strongest BIB framing is:

- spaGAPA is a reproducible spatial APA analysis framework.
- It is not only an implementation of one model.
- It systematically evaluates APA recovery, uncertainty, biological consistency, and scalability.
- It is CPU-friendly and accessible to ordinary bioinformatics labs.
- It includes real spatial transcriptomics workflows and high-resolution readiness.

### 10.2 What not to claim

Avoid claims such as:

- "GP alone fully solves biological domain recovery."
- "spaGAPA beats every competitor on every metric."
- "Visium HD automatically supports APA calling."
- "Weak labels are pathology gold standards."

### 10.3 What to claim if supported by final benchmark

Defensible claims:

- spaGAPA improves or matches APA value recovery in real spatial data.
- spaGAPA provides calibrated uncertainty unavailable in KNN-only tools.
- BioML graph route improves biological-domain consistency in sparse/high-resolution settings.
- High-resolution presets provide an accuracy/speed tradeoff without GPU dependency.
- The package provides a reproducible benchmark and real-data pipeline for spatial APA analysis.

---

## 11. Roadmap

### Short term

1. Finish GSE220442 Space Ranger processing.
2. Decide repair/skip strategy for failed GSE237183 sample.
3. Run scAPAtrap on additional GSE237183/GSE220442 BAMs.
4. Build formal processed datasets.
5. Run multi-real-data benchmark.
6. Produce benchmark visualizations.

### Medium term

1. Build marker-defined and ROI/domain labels.
2. Add one high-resolution raw-read dataset, preferably Stereo-seq.
3. Add robust runtime/memory scaling figure.
4. Consolidate CLI tutorials.
5. Prepare manuscript figure panels.

### Submission-ready target

Minimum evidence for BIB:

- Two or more real spatial APA datasets.
- One convincing biological-domain validation.
- One high-resolution or pseudo-high-resolution validation.
- Clear uncertainty calibration figure.
- Runtime and scalability figure.
- Reproducible code/data manifest.

---

## 12. Final Project Direction

The current best direction is:

> Build spaGAPA as a CPU-friendly, uncertainty-aware, benchmark-driven spatial APA framework where GP handles APA recovery and uncertainty, while BioML graph learning handles biological-domain consistency.

This direction is technically coherent, defensible to reviewers, and aligned with the user's constraint of avoiding deep learning/GPU dependence.

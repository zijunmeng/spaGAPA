# spaGAPA Development Progress

**Project**: spaGAPA - spatial transcriptomics APA analysis toolkit
**Current target**: Briefings in Bioinformatics (BIB)
**Last updated**: 2026-07-15
**Current status**: BIB-oriented algorithm refinement + multi-real-dataset benchmark expansion (18+1 real Visium APA datasets; Stereo-seq APA proven end-to-end)

---

## 1. Current Position

spaGAPA has moved from an initial package prototype into a BIB-oriented method project. The current scientific claim is no longer only "a package for spatial APA", but:

> Spatial APA analysis requires dedicated statistical treatment of sparse, spatially structured, and uncertain APA usage measurements. spaGAPA provides a reproducible framework that integrates spatial validation, Gaussian process imputation, uncertainty propagation, CPU-friendly graph machine learning, and downstream SVAPA/differential analysis for robust interpretation of spatially organized APA remodeling in real tissues.

The active main line is:

1. Use scAPAtrap or compatible APA callers to generate real APA/PAS evidence from Space Ranger BAM.
2. Use spaGAPA to validate, impute, quantify, benchmark, and interpret spatial APA signals.
3. Use Gaussian process models primarily for APA value recovery and uncertainty.
4. Use BioML multi-view graph models for biological domain consistency, especially in high-resolution or pseudo-bin data.
5. Validate on multiple real spatial transcriptomics datasets instead of relying on one simulated benchmark.

---

## 2. Repository and Package Status

### 2.1 Core package

Implemented and actively used modules:

- `spagapa.core`
  - `APADataset`
  - `APASite`
  - `APASiteCollection`
- `spagapa.io`
  - coordinate readers
  - BAM/spatial file readers
  - result writers
  - scAPAtrap wrapper/integration utilities
- `spagapa.spatial`
  - KNN graph
  - radius graph
  - Delaunay graph
  - spatial weights
- `spagapa.calling`
  - spatial validation
  - quality filtering
- `spagapa.imputation`
  - standard GP
  - sparse/fast GP variants
  - expression-informed GP prototypes
  - product/adaptive/layer-local GP experiments
- `spagapa.quantification`
  - RUD
  - PDUI
  - WUL
  - PAI
  - QC metrics
- `spagapa.analysis`
  - spatial domains
  - differential APA
  - SVAPA/spatial pattern analysis
  - GP trend detection
- `spagapa.visualization`
  - spatial plots
  - statistical plots
  - QC plots
  - benchmark plots
- `spagapa.benchmark`
  - simulation benchmark
  - MOB/stAPAminer-style benchmark
  - high-resolution pseudo-bin benchmark
  - real-data benchmark runners
- `spagapa.bioml`
  - CPU-only multi-view graph construction
  - graph-regularized domain route
  - high-resolution domain recovery route
- `spagapa.pipeline`
  - end-to-end pipeline
  - preset dispatch
  - BioML/highres route integration
- `spagapa.cli`
  - command-line entry
  - analysis presets

### 2.2 Presets and pipeline routes

Current user-facing analysis presets:

- `auto`
  - Inspects input size/sparsity and dispatches to a suitable route.
- `standard`
  - Spatial APA analysis for conventional Visium-scale data.
- `highres_accuracy`
  - High-resolution/pseudo-bin route prioritizing biological-domain stability.
- `highres_fast`
  - Faster high-resolution route with small accuracy tradeoff.

Current preferred high-resolution candidate parameters:

- `gp_blend = 0.1`
- graph weight `spatial = 0.1`
- graph weight `expression = 0.7`
- graph weight `apa = 0.2`

This configuration was selected because it improved real-data high-resolution smoke RMSE while preserving the decoupled biological-domain graph route.

---

## 3. Milestone History

### Milestone 0: Initial package foundation - completed

Completed:

- Python package structure.
- Core APA data classes.
- I/O layer.
- Spatial graph utilities.
- Spatial validation.
- GP imputation.
- APA quantification.
- Differential APA and SVAPA analysis.
- Visualization and benchmark framework.
- Initial documentation and examples.

Outcome:

- spaGAPA became a functional APA analysis package rather than a collection of scripts.

### Milestone 1: BIB-ready foundation audit - completed

Main issues addressed:

- Matrix orientation consistency.
- Reader return type consistency.
- Sparse GP parameter handling.
- Quantification placeholder cleanup.
- End-to-end toy tests.
- Pipeline audit for basic BIB credibility.

Outcome:

- The package can run a controlled toy workflow and expose predictable input/output behavior.

### Milestone 2: Real APA data path - completed for first datasets

Completed:

- Identified real Visium datasets with accessible raw data.
- Ran Space Ranger for GSE179572/GSM5420751.
- Checked BAM barcode/UMI tags.
- Ran scAPAtrap on Space Ranger BAM.
- Exported expression matrix from Space Ranger `filtered_feature_bc_matrix.h5`.
- Built processed real-data directories with:
  - `apa_matrix.csv`
  - `apa_sites.csv`
  - `coordinates.csv`
  - `expression_matrix.csv`
  - `qc_summary.json`

Key processed datasets:

- `data/processed/gse179572_gsm5420751_scapatrap/`
- `data/processed/gse183456_gsm6047774_scapatrap/`

Outcome:

- spaGAPA is no longer only benchmarked on simulated or tutorial data. It can consume real Space Ranger + scAPAtrap APA matrices.

### Milestone 3: Formal benchmark v2 - completed

Completed benchmark components:

- More genes.
- Multiple random seeds.
- Random mask.
- Spatial block mask.
- Low-coverage mask.
- Ring/sector mask.
- Layer-aware mask.
- Runtime/memory tracking.
- Layer ARI/NMI for biological consistency.
- Uncertainty calibration:
  - interval coverage
  - uncertainty-error correlation
  - accuracy after high-uncertainty filtering

Key lessons:

- Spatial GP is strong for RMSE and uncertainty-aware APA recovery.
- Pure expression KNN can remain strong for layer/domain clustering.
- A single model should not be forced to optimize both APA value recovery and domain discovery.

Outcome:

- This motivated the current decoupled route: GP for APA values, BioML graph for biological domains.

### Milestone 4: Expression-informed GP exploration - completed and deprioritized

Explored variants:

- additive spatial + expression GP
- radial GP
- product kernel GP
- adaptive/gated expression weighting
- layer-local expression kernel

Main conclusion:

- Expression-informed GP variants improved some biological-consistency metrics in limited settings but did not robustly dominate spatial GP.
- Product/adaptive/layer-local GP did not provide a clean enough improvement to become the main BIB story.
- The main risk was losing GP calibration or runtime while still not beating graph/KNN baselines on domain recovery.

Decision:

- Keep expression-informed GP as an experimental route.
- Do not make it the core BIB method.
- Move biological consistency into a dedicated CPU-friendly BioML graph route.

### Milestone 5: BioML and high-resolution route - completed as main candidate

Implemented:

- `spagapa.bioml` module.
- Multi-view graph builder.
- Graph-regularized APA factor/domain route.
- Decoupled high-resolution strategy:
  - APA value recovery handled separately from biological-domain recovery.
  - Domain labels are derived from graph structure, not by forcing GP to learn every biological manifold.

Representative benchmark evidence:

- High-resolution BioML suite showed stronger layer/domain recovery than raw APA and expression KNN in selected pseudo-bin settings.
- `highres_fast` gave substantial runtime improvement with limited RMSE cost.
- Decoupled route is currently the most promising high-resolution BIB narrative.

Outcome:

- spaGAPA has a credible CPU-only route for high-resolution spatial transcriptomics, avoiding deep learning and GPU dependency.

### Milestone 6: Pipeline/CLI integration - mostly completed

Completed:

- Main pipeline supports analysis presets.
- CLI exposes preset selection.
- BioML/highres route connected to pipeline-level execution.
- Smoke and small suite tests added for high-resolution route.
- Real-data highres smoke benchmark completed on GSE179572.

Remaining:

- Further hardening of CLI documentation.
- More examples showing exact command-line use for real Space Ranger/scAPAtrap outputs.
- More stable benchmark summary tables for manuscript use.

### Milestone 7: Real dataset expansion - completed for core BIB set

Current real datasets:

#### GSE179572 / GSM5420751

Status:

- Space Ranger completed.
- scAPAtrap completed.
- Processed dataset ready.
- Expression matrix exported.
- APA matrix and sites available.
- Marker-defined weak labels explored.

Current QC:

- `apa_ready: true`
- `expression_ready: true`
- `n_spots: 4992`
- `n_called_sites: 26478`
- `n_gene_annotated_sites: 22860`
- `n_apa_usage_sites: 14101`

Use:

- First real-data benchmark.
- High-resolution smoke validation.
- External smoke validation.
- Candidate dataset for ROI/pathology weak/gold labels.

#### GSE183456 / GSM6047774

Status:

- Space Ranger completed successfully.
- scAPAtrap completed.
- Processed dataset fully ready (apa + expression + coordinates all in place).

Current QC:

- `apa_ready: true`
- `expression_ready: true` (38,606 genes)
- `n_spots: 3010`
- `n_called_sites: 53572`
- `n_gene_annotated_sites: 41747`
- `n_apa_usage_sites: 34118`

Tissue confirmation:

- Tissue = kidney, confirmed by three-way cross-check (GEO metadata, KPMP human kidney atlas, local SOFT).

Use:

- Fully ready for spaGAPA real-data benchmark (apa + expression + coords).
- Strong candidate for biological consistency/domain recovery.
- Better chance than GSE179572 for interpretable tissue/domain labels.

#### GSE237183

Status:

- 76/76 FASTQ files downloaded.
- 19/19 Space Ranger dry-runs passed.
- Formal Space Ranger count complete (submitted with `JOBS=2`, `THREADS=16`, `MEM_GB=96`).
- Final count summary:
  - 18 samples PASS (complete `outs/`, usable for scAPAtrap)
  - 1 sample FAIL: GSM7596588 (GBM ZH881 infiltrating)
- scAPAtrap batch 18/18 complete. Pilot 2 samples (GSM7596590 GBM bulk + GSM7596601 IDHm astrocytoma) plus batch of 16 (3-way tmux parallel, reads-balanced greedy bins) all finished.
- 18 Visium APA datasets produced in `data/processed/gse237183_gsm*_scapatrap/` (each contains `apa_matrix.csv`, `apa_sites.csv.gz`, `coordinates.csv`, `qc_summary.json`).
- GSM7596588 skipped (fiducial failure, see note below); ZH881 still covered by 4 passing sections (GSM7596589, GSM7596597, GSM7596598, GSM7596599).

Important note:

- GSM7596588 failed at `ALIGN_FIDUCIALS` ("too many points rotated out of the image space") - a deterministic fiducial/tissue-image alignment failure, not a FASTQ/reference/resource issue. Re-run or `--resume` would fail identically; Loupe manual alignment is the only fix.
- Decision: skip GSM7596588 in the first formal analysis. ZH881 is still covered by 4 passing sections (GSM7596589, GSM7596597, GSM7596598, GSM7596599).
- BAM tag audit on the 18 PASS samples (first 200 records each): CB 190-200/200, UB 199-200/200, all `samtools quickcheck` OK. CB/UB are effectively complete, so all 18 are scAPAtrap-ready. GX/GN presence is low/variable as expected (only gene-annotated reads carry these tags; not required by scAPAtrap). Audit written to `spaGAPA/pipeline_output/gse237183_spaceranger_counts/_logs/gse237183_bam_tag_audit_summary.tsv`.

scAPAtrap spot/site validation (sample checks):

- GSM7596590: 54,131 sites x 4,992 spots; value range [0,1]; spot count matches coordinates.
- GSM7596601: 26,216 sites x 3,311 spots; value range [0,1]; spot count matches coordinates.

Use:

- Important brain/Visium candidate for multi-real-data BIB benchmark (18 usable Visium sections, now all APA-ready).

#### GSE220442

Status:

- Original ENA FASTQ appeared as single FASTQ and was not directly Space Ranger-ready.
- Investigation showed SRA has four technical reads.
- Correct mapping:
  - read 3 = barcode/UMI, use as R1
  - read 4 = cDNA, use as R2
  - read 1/read 2 are index reads
- Pilot sample `GSM6801755` completed split + dry-run.
- Full suite submitted in tmux:
  - session: `spagapa_gse220442_sra_split_all`
  - log: `logs/20260713_gse220442_sra_split_suite.log`

Use:

- Brain Visium candidate.
- Important for showing spaGAPA can handle non-standard SRA read structures.

### Milestone 8: High-resolution dataset scouting - audited + Stereo-seq APA pilot proven end-to-end

#### 8.1 High-resolution candidate audit (GEO suppl + SRA raw FASTQ via ENA)

The top candidates were audited at both the GEO supplementary level (what is actually inside each `RAW.tar`) and the SRA raw-FASTQ level (true sizes from ENA `filereport`). Candidates have moved from paper lists to evidence-backed entries.

Audited candidates (GEO suppl content + SRA raw FASTQ total size):

| GSE | Platform / tissue | GEO suppl content | SRA raw FASTQ | SRA project |
|-----|-------------------|-------------------|---------------|-------------|
| GSE269906 | Stereo-seq / human PFC AD | processed `*.cellbin.gef` + `*.tissue.gef` (no raw FASTQ/BAM) | 557 GB (12 runs) | PRJNA1120963 |
| GSE263789 | Stereo-seq / AD mouse | `barcodeToPos.h5` + `bin1.max.gem.txt.gz` + h5ad | 511 GB (32 runs) | PRJNA1099247 |
| GSE311383 | Visium HD / human liver | complete multi-resolution Space Ranger output (2/8/16um) | 100 GB (28 runs) | PRJNA1369387 |
| GSE268519 | Slide-seqV2 / mouse hippocampus + cerebellum | series-level RCTD `.rds.gz` only (per-sample = NONE) | 267 GB (16 runs) | PRJNA1117495 |
| **Total** | | | **~1.44 TB** | |

Demoted candidate:

- GSE290724 (Visium HD GBM): `RAW.tar` only 252 MB; GSM-level contains only `barcodes/features/matrix.mtx.gz` (processed matrix, no raw, no spatial files, no GSM-level SRA link). Demoted from the high-resolution list; Visium HD expression/domain candidate switched to GSE311383.

Verdict (download priority):

- GSE263789 (mouse Stereo-seq AD) = APA raw-read pilot first choice (smallest SRR28637909 is only 7 GB; carries its own `barcodeToPos.h5` + h5ad; strong AD biology).
- GSE269906 (human Stereo-seq AD/normal) = upgrade target after the mouse pipeline is proven (557 GB; smallest run 180 GB).
- GSE311383 (human liver Visium HD) = expression/domain only; raw may not preserve 3' (FFPE/probe), not an APA benchmark.
- GSE268519 (Slide-seqV2) = one puck to test Slide-seqV2 APA feasibility (smallest 1.7 GB).

Stereo-seq read structure (verified on SRR28637909):

- read1 = 35 bp barcode/UMI; read2 = 100 bp cDNA; layout = PAIRED; instrument = DNBSEQ-T7.
- The spatial barcode is encoded as a chip coordinate in the read name (e.g. `E100043442L1C001R001...`), resolved to physical x/y via the supplementary `barcodeToPos.h5` - not a 10x-style `CB` read.
- read2 (100 bp) is cDNA derived from polyA in-situ capture -> 3'-derived, carrying 3'/polyA evidence suitable for PAS/3'-end peak calling.

#### 8.2 Stereo-seq APA pilot (GSE263789 / SRR28637909) - PROVEN end-to-end

The full Stereo-seq APA chain was built and run successfully on the smallest GSE263789 run (SRR28637909, mouse AD brain):

```text
SRA raw (SRR28637909) -> fasterq-dump (read1 35bp barcode/UMI + read2 100bp cDNA)
  -> SAW 8.2.2 count (bcSTAR align to GRCm38 + spatial-coordinate annotation + GEF)
  -> retag (CB=Cx_Cy, UB=UR padded to len 5, GX=GI, GN=GS) + chr-prefix via samtools reheader
  -> scAPAtrap (mouse chr1-19+XY)
  -> 8,659 PAS sites x 20,677,477 spots spatial APA matrix
```

Installed infrastructure:

- SAW 8.2.2 (BGI Stereo-seq Analysis Workflow) at `/s1/SHARE/01_software/saw-8.2.2` (self-contained conda envs with Python/R/STAR/samtools; no docker or network needed).
- Mouse SAW reference at `/s1/SHARE/01_software/SAW_refs/Mus_musculus_index` (GRCm38 / Ensembl93).

SAW count results (pilot GSM8199179 / SRR28637909, ~4h, 13/15 steps successful):

- ALIGNMENT: sorted BAM (7.2 GB).
- ANNOTATION: dedup + gene-annotated target BAM (5.9 GB).
- EXPRESSION_MATRIX: `tissue.gef` spatial expression matrix (2.7 GB).
- REPORT: `gef` + `report.html`.
- visualization: cosmetic failure only (`invalid serial` SN-format parse), does not affect BAM/GEF.

BAM tag structure (APA inputs present, tag names differ from scAPAtrap expectations):

- `Cx`/`Cy` = spatial coordinate (x, y) - the Stereo-seq "cell barcode".
- `UR` = UMI.
- `GS`/`GI`/`GE` = gene annotation (~79% of reads carry it).
- Chromosomes named `1,2,...,X,Y` (Ensembl, no `chr` prefix); 80.4M deduplicated reads.

Retag step (adapts SAW BAM tags to scAPAtrap's expected 10x tags):

- `CB` = `Cx_Cy`, `UB` = `UR` padded to fixed length 5, `GX` = `GI`, `GN` = `GS`.
- `samtools reheader` adds the `chr` prefix to the BAM header (header-only edit, seconds) so scAPAtrap's `chr1`-style peaks match.
- UMI length was non-uniform (97.8% length 5, 2.2% length 1-4) which broke umi_tools' fixed-length assertion; padding UB to length 5 fixed it.

APA feasibility evidence (50k-read sample of the SAW annotated BAM, read 3' end to gene TES distance):

| Distance from gene TES | Fraction |
|------------------------|----------|
| 0-500 bp (flush against gene 3' end) | 33.8% |
| 500 bp - 2 kb | 21.4% |
| 2-10 kb | 9.4% |
| >10 kb | 34.8% |

- 55% of gene-annotated reads fall within 2 kb of the gene 3' end; 33.8% flush within 500 bp. This is the characteristic polyA-capture 3'-enrichment signature (random-fragmentation RNA-seq does not pile up at 3' ends) -> the read 3'-end enrichment marks polyA cleavage sites (PAS), which is exactly the signal scAPAtrap relies on for PAS calling.
- polyA-in-read 1.27% (normal; polyA tails are often soft-clipped); 66.4% of reads carry gene annotation.
- Conclusion: Stereo-seq cDNA is 3'-enriched; APA calling is feasible.

Final spatial APA matrix:

- n_sites = 8,659 PAS sites (polyA-tail confirmed, `min.count`/`min.cells=10` filters passed).
- n_spots = 20,677,477 (Stereo-seq subcellular DNB resolution).
- n_nonzero = 23,019,011 (sparse).
- Outputs in `pipeline_output/gse263789_stereo_pilot/scapatrap_raw/`: `scAPAtrapData.rda` (149M), `peaks_meta.csv.gz` (8,659 peaks), `apa_site_counts.csv.gz` (sparse long format, 127M), `peaks.saf.reduced`, `counts.tsv.gz.reduced` (880M), `scapatrap_qc.json`.
- Dense matrix export OOMs (exit 137) at 20.7M spots; sparse long-format export (peak x spot x count triples) was used instead.

Outcome:

- This proves Stereo-seq can do APA at subcellular resolution, providing real-data evidence for the spaGAPA high-resolution APA story.
- Next steps: bin the subcellular APA matrix into the spaGAPA downstream pipeline; scale from the pilot to the full GSE263789 study and then to human GSE269906.

Important caveat:

- Visium HD FFPE/probe-based data are valuable for high-resolution expression/domain validation, but may not provide primary APA/PAS evidence. For true APA calling, Stereo-seq or poly(A)-compatible high-resolution protocols are more important.

---

## 4. Representative Benchmark Evidence

### 4.1 Real GSE179572 high-resolution smoke

File:

- `benchmark_results/real/gse179572_scapatrap_highres_smoke/highres_default_validation_summary.csv`

Observed result:

- Current default RMSE: approximately `0.05215`
- Candidate highres default RMSE: approximately `0.05019`

Interpretation:

- The candidate highres parameters slightly improve real-data RMSE.
- This supports changing the highres default toward `gp_blend=0.1`, `spatial=0.1`, `expression=0.7`, `apa=0.2`.

### 4.2 Highres pipeline smoke suite

File:

- `benchmark_results/real/pipeline_highres_smoke_suite_v1/pipeline_highres_smoke_overall_summary.csv`

Representative result:

- `pipeline_highres_accuracy`
  - RMSE approximately `0.08277`
  - layer ARI approximately `0.3174`
  - runtime approximately `2.16s`
- `pipeline_highres_fast`
  - RMSE approximately `0.08668`
  - layer ARI approximately `0.3174`
  - runtime approximately `0.20s`

Interpretation:

- `highres_fast` is much faster while preserving domain recovery in this small suite.
- This supports exposing both high-accuracy and fast presets.

### 4.3 BioML highres suite

File:

- `benchmark_results/real/highres_bioml_suite_v1/highres_bioml_suite_overall_summary.csv`

Representative result:

- highres BioML improved layer/domain recovery over raw APA and expression KNN in selected pseudo-bin benchmarks.
- It became the strongest candidate for high-resolution biological consistency.

Interpretation:

- BioML should be framed as the domain-consistency route, not as a replacement for GP-based APA value recovery.

### 4.4 Product/adaptive expression GP suite

File:

- `benchmark_results/real/product_gp_suite_v1/product_sweep_summary.csv`

Interpretation:

- Product/adaptive expression GP did not provide enough robust improvement to justify becoming the main method.
- It remains useful as an ablation and negative result showing why the decoupled BioML route was necessary.

---

## 5. Known Problems and Risks

### 5.1 Biological labels remain the largest evaluation bottleneck

For BIB, RMSE alone is not enough. We need biological consistency evidence:

- pathology ROI labels
- tumor/non-tumor/stroma/necrosis labels
- layer/domain labels
- marker-defined weak labels
- spatially coherent domain recovery

Current status:

- GSE179572 has weak-label potential but no pathology gold standard yet.
- GSE183456 is fully ready (apa + expression + coords, 38,606 genes x 3,010 spots, kidney confirmed) and is the strongest current biological-domain candidate.
- GSE237183 scAPAtrap batch is complete (18/18), giving 18 brain/Visium sections for domain labels.
- GSE220442 SRA split suite is still running.
- Stereo-seq pilot (GSE263789) produced a subcellular APA matrix (8,659 PAS x 20.7M spots) but no pathology/domain gold standard yet.

### 5.2 GP runtime is still a weakness

Observed issue:

- Standard GP can have strong RMSE but slower runtime.
- Fast GP reduces runtime but needs careful validation to avoid accuracy loss.

Current strategy:

- Keep standard GP for accuracy route.
- Use fast/sparse GP for scalable route.
- Do not claim runtime dominance until larger real-data results support it.

### 5.3 Domain ARI is not automatically solved by GP

Repeated benchmark finding:

- GP can recover APA values well.
- Expression KNN or graph methods may recover biological domains better.

Decision:

- Use decoupled route:
  - GP for APA values and uncertainty.
  - BioML graph for biological domains.

This is a methodological strength if presented clearly.

### 5.4 High-resolution APA data are harder than high-resolution expression data

Visium HD and many FFPE/probe-based platforms are not automatically suitable for primary APA/PAS calling.

For high-resolution APA evidence, priority should be:

- Stereo-seq or other poly(A)-compatible raw reads.
- BAM/FASTQ retaining barcode/UMI/cDNA structure.
- Datasets where 3' end or poly(A) evidence can be extracted.

Update: the Stereo-seq APA pilot (GSE263789 / SRR28637909) has now demonstrated this end-to-end - SAW count + retag + scAPAtrap produced an 8,659 PAS x 20.7M-spot subcellular APA matrix, with 55% of gene-annotated reads within 2 kb of the gene 3' end (33.8% within 500 bp) confirming the polyA-capture 3'-enrichment signature. This validates Stereo-seq as a viable high-resolution APA data source.

### 5.5 Real-data processing is still fragile

Observed issues:

- GEO/SRA read structures are inconsistent.
- Space Ranger image/fiducial alignment can fail.
- Some datasets have expression matrices but no usable APA evidence.
- Some datasets need manual read remapping or custom FASTQ reconstruction.
- Stereo-seq data require a SAW pipeline and a retag/chr-prefix/UMI-padding adapter (now built and working).

Mitigation:

- Maintain per-dataset processing notes.
- Keep reusable scripts for each abnormal data format.
- Document all failure modes for reproducibility.
- The Stereo-seq SAW-to-scAPAtrap adapter is reusable for other Stereo-seq datasets (e.g. GSE269906).

---

## 6. Immediate Next Steps

### 6.1 Finish current real-data processing

Priority:

1. Monitor and complete GSE220442 SRA split suite.
2. GSE237183 scAPAtrap is complete (18/18); the failed GSM7596588 remains skipped (Loupe manual alignment is the only fix and is not blocking).
3. Standardize the 18 GSE237183 processed outputs in `data/processed/gse237183_gsm*_scapatrap/`.
4. Bin the Stereo-seq subcellular APA matrix (GSE263789, 20.7M spots) into the spaGAPA downstream pipeline.
5. Scale the Stereo-seq pilot to the full GSE263789 study and then to human GSE269906.

### 6.2 Build biological labels without waiting for pathology gold standard

Parallel routes:

1. Marker-defined weak labels.
2. Spatial expression domain labels.
3. Tumor/stroma/immune marker labels for cancer datasets.
4. ROI labels if pathology support becomes available.

Important:

- Weak labels must be clearly described as weak or marker-derived, not pathology gold standard.

### 6.3 Formal multi-real-data benchmark

Required datasets:

- at least one completed cancer Visium dataset: GSE179572 (done)
- one additional completed Visium dataset: GSE183456 (fully ready) or GSE237183 (18/18 complete) - both available
- one brain/layer or brain-region dataset: GSE237183 (18 sections, APA-ready) and/or GSE220442 (still processing)
- one high-resolution candidate: Stereo-seq pilot (GSE263789) proven end-to-end; full-study + human GSE269906 pending

Required metrics:

- RMSE / MAE
- Pearson / Spearman
- layer/domain ARI and NMI
- SVAPA ranking stability
- uncertainty calibration
- runtime and memory

### 6.4 High-resolution main story

Next technical work:

- High-resolution candidate datasets are now audited (GEO suppl + SRA sizes, ~1.44 TB total).
- Stereo-seq APA pilot proven end-to-end on GSE263789 (8,659 PAS x 20.7M spots); bin the subcellular matrix into the downstream pipeline.
- Prioritize Stereo-seq raw-read datasets for true APA evidence (GSE263789 full study, then GSE269906 human brain).
- Keep Visium HD (GSE311383) as high-resolution expression/domain validation only (FFPE/probe chemistry not APA-suitable).
- Add formal highres benchmark figures with:
  - 2x / 4x / 8x pseudo-bins
  - runtime scaling
  - domain stability
  - uncertainty filtering

### 6.5 Manuscript-facing figures

Target figure set:

1. Method overview.
2. Simulation/MOB benchmark.
3. Real Visium APA benchmark.
4. Biological-domain consistency benchmark.
5. High-resolution pseudo-bin benchmark.
6. Uncertainty calibration and downstream filtering.
7. Runtime/scalability.
8. Real tissue biological case study.

---

## 7. BIB Readiness Checklist

Current state:

- [x] Functional package.
- [x] CLI and pipeline presets.
- [x] Real Space Ranger + scAPAtrap path.
- [x] GP imputation and uncertainty.
- [x] BioML graph route.
- [x] High-resolution route prototype.
- [x] Formal benchmark scripts.
- [x] First real-data benchmarks.
- [x] Candidate high-resolution dataset search.
- [x] At least two completed real APA benchmark datasets (GSE179572 + GSE183456 fully ready; GSE237183 18/18 complete).
- [x] At least one strong biological-domain validation dataset (GSE183456 kidney, fully ready; GSE237183 brain 18 sections).
- [x] At least one high-resolution or pseudo-high-resolution APA validation (Stereo-seq pilot GSE263789 proven end-to-end, 8,659 PAS x 20.7M spots).
- [ ] Robust runtime/memory scaling figure.
- [ ] Final documentation and tutorials.
- [ ] Manuscript-quality figure set.
- [ ] Clean reproducibility manifest.

---

## 8. Practical Commands

Activate development environment:

```bash
conda activate spagapa
```

Run tests:

```bash
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
pytest
```

Run CLI with automatic preset:

```bash
spagapa run \
  --apa-matrix data/processed/<dataset>/apa_matrix.csv \
  --coordinates data/processed/<dataset>/coordinates.csv \
  --expression-matrix data/processed/<dataset>/expression_matrix.csv \
  --analysis-preset auto \
  --output pipeline_output/<dataset>_spagapa
```

Use high-resolution accuracy route:

```bash
spagapa run \
  --apa-matrix data/processed/<dataset>/apa_matrix.csv \
  --coordinates data/processed/<dataset>/coordinates.csv \
  --expression-matrix data/processed/<dataset>/expression_matrix.csv \
  --analysis-preset highres_accuracy \
  --output pipeline_output/<dataset>_highres_accuracy
```

Use high-resolution fast route:

```bash
spagapa run \
  --apa-matrix data/processed/<dataset>/apa_matrix.csv \
  --coordinates data/processed/<dataset>/coordinates.csv \
  --expression-matrix data/processed/<dataset>/expression_matrix.csv \
  --analysis-preset highres_fast \
  --output pipeline_output/<dataset>_highres_fast
```

---

## 9. Current Strategic Decision

The project should not attempt to force a single GP model to win every metric. The stronger BIB-level framing is:

- spaGAPA uses GP where GP is statistically appropriate: APA value recovery and calibrated uncertainty.
- spaGAPA uses BioML graph learning where graph methods are appropriate: biological domain consistency in sparse spatial data.
- The package is CPU-friendly and avoids deep learning/GPU dependency.
- The benchmark explicitly tests where each component is useful, instead of hiding negative results.

This makes the method more defensible to reviewers and gives a clearer route to comprehensive improvement.

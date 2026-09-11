# spaGAPA Project Summary

**Project name**: spaGAPA
**Full name**: spatial Gaussian-process and graph-aware APA analyzer
**Language**: Python (R used only for external callers/competitors)
**Main environment**: conda env `spagapa` (S91 primary host)
**Current manuscript target**: Briefings in Bioinformatics (BIB)
**Last updated**: 2026-09-09
**Current status**: submission preparation complete — manuscript figure set (7 main + 16 supplementary, print-width vector PDF with submission-ready legends), Table 1, multi-caller robustness validation (Sierra × spaGAPA, 3 datasets / 2 species / 3 tissues), and honest-narrative review fixes are all in place; remaining work is the submission itself and any post-review response.

---

## 1. Executive Summary

spaGAPA is a Python toolkit for alternative polyadenylation (APA) analysis in spatial transcriptomics data. The project evolved from a package prototype into a statistical and machine-learning framework for sparse spatial APA signals, and has now passed its publication freeze: the public repository carries the package, the full benchmark/validation suite, and all figure-generation scripts, while the manuscript text itself is maintained outside the repo (its drafting history is preserved in git).

The BIB-oriented thesis is:

> Spatial APA signals are sparse, spatially structured, uncertain, and highly dependent on tissue context. A strong spatial APA framework should not only call APA sites, but also validate spatial support, recover missing APA usage with calibrated uncertainty, separate value recovery from biological-domain discovery, and provide reproducible benchmarks on real spatial transcriptomics datasets.

spaGAPA addresses this through five connected layers:

1. Real-data APA input construction from Space Ranger BAM via scAPAtrap or Sierra (caller-agnostic by design, validated by caller swap).
2. Spatial APA validation, quantification, and sparse-GP imputation with conformal-calibrated uncertainty.
3. CPU-friendly BioML graph learning for biological-domain recovery (decoupled from value recovery).
4. Uncertainty-aware downstream analysis: risk-coverage triage, SVAPA, donor-level differential APA, optional APA batch correction.
5. Benchmark and visualization infrastructure spanning simulation, MOB/layer data, 9 real GSE datasets, and subcellular Stereo-seq.

---

## 2. Scientific Motivation

APA changes the 3' end of transcripts and can alter transcript stability, localization, translation, and regulatory interactions. Existing single-cell and spatial transcriptomics APA tools leave several problems under-addressed:

- APA measurements are sparse and dropout-prone; high-resolution platforms make this worse.
- Spatial transcriptomics adds physical neighborhood structure that should be modeled explicitly.
- KNN-style imputation lacks calibrated uncertainty; no spatial APA tool offers prediction intervals with coverage guarantees.
- Biological domain recovery and APA value recovery are related but not identical tasks; forcing one model to win both is fragile.
- Published tools often lack comprehensive, reproducible benchmark pipelines, and conclusions may be entangled with the choice of PAS caller.

spaGAPA is designed around these gaps, and its claims are deliberately scoped to what the benchmarks support.

---

## 3. Core Contributions

### 3.1 Spatial APA data model and workflow

- APA site representation, count/usage matrix handling, spatial coordinates, result writing.
- CLI and pipeline presets (`auto`, `standard`, `highres_accuracy`, `highres_fast`).
- Integration with AnnData-like ecosystems; this makes spaGAPA a package rather than a one-off benchmark script.

### 3.2 Spatial validation before downstream modeling

APA sites are not treated as independent molecular events detached from tissue geometry. The validation layer evaluates read support, spot support, spatial support, neighbor consistency, and spatial autocorrelation — noisy APA calls would otherwise dominate imputation and SVAPA detection.

### 3.3 Sparse-GP value recovery with conformal uncertainty (core innovation)

- Inducing-point sparse GP: O(n·m²), batch fit across genes with K_nm precomputation.
- Split-conformal prediction (global + locally adaptive) with finite-sample-corrected quantiles.
- Validated empirically: 11 samples, 523,174 test points — mean |coverage deviation| 0.21 / 0.16 / 0.10 pp at 80/90/95%; max deviation ≤ 0.5 pp; holds under spatial-block splits.
- Four uncertainty models (constant, local-gene, spatial-spot, residual-spot) benchmarked on 5 datasets; residual-spot (method D) clears pooled uncertainty–error r ≥ 0.3 without inflating width.
- Risk-coverage triage: filtering by uncertainty removes 23% of RMSE at 80% retention (paired across 5 datasets, p = 0.004).

### 3.4 Systematic exploration of expression-informed GP (negative result, kept as ablation)

Additive / radial / product-kernel / adaptive / layer-local expression GPs were all benchmarked. They improved some domain metrics in selected settings but did not robustly dominate spatial GP, and risked worse calibration or runtime. Decision: kept as ablation evidence, not the main method. This prevents overclaiming and motivated the decoupled BioML route.

### 3.5 Decoupled BioML route for biological-domain consistency

- CPU-only, no GPU, no deep learning dependency.
- Multi-view graph (spatial coordinates + expression + APA usage) → graph-regularized factorization → Leiden domains.
- Conceptual split: GP handles APA value recovery and uncertainty; BioML handles biological-domain consistency.
- MOB validation: unsupervised ARI = 0.60, NMI = 0.68 against 5 annotated layers, with a fair-comparison design (identical Leiden pipeline; only the APA source differs between mean-imputed and GP-imputed).

### 3.6 Caller-agnostic design with multi-caller validation (new)

Sierra 0.99.27 was run on the same Space Ranger BAMs of three Visium datasets (GSE183456 human kidney, GSE220442 human AD-brain PFC, GSE169749 mouse colon), converted to spaGAPA input format (usage matrix + shared coordinates, ≥2 sites per gene, min_parent_count = 5) with **no re-tuning**, and pushed through the identical inference:

- Conformal coverage stays at nominal on every dataset: 0.801/0.900/0.950 (kidney), 0.796/0.898/0.950 (brain), 0.804/0.899/0.950 (colon); max deviation 0.5 pp across all datasets × levels × modes.
- Gene-level distal-usage concordance across callers: median per-gene r = 0.51 / 0.83 / 0.69.
- PAS overlap itself is caller-dependent (Jaccard 0.05–0.18 at ±50 bp; 81–88% of Sierra peaks within 500 bp of a scAPAtrap site), consistent with the independent-caller benchmark literature.
- The "per-gene mean is a strong RMSE baseline" pattern holds on both caller inputs in every dataset.
- Companion evidence: the MOB analysis runs on a movAPA-prepared (non-scAPAtrap) APA matrix.

Full numbers and a Methods-ready paragraph: `scripts/multicaller_validation/report.md`; one-command driver `scripts/multicaller_validation/run_dataset.sh` (hostname-aware).

### 3.7 High-resolution / subcellular Stereo-seq support

End-to-end chain proven on GSE263789 (mouse AD brain): SRA FASTQ → SAW 8.2.2 count (bcSTAR + spatial annotation) → BAM retag adapter (CB=Cx_Cy, UB padded to length 5, chr-prefix reheader) → scAPAtrap → 21,455 PAS × 20.7M DNB subcellular matrix → bin200 aggregation (15,235 spots, ~10.3% observed) → spaGAPA downstream.

- 3'-enrichment signature confirmed: 53.4% of reads within 500 bp of TES (73.5% within 2 kb) — polyA-capture, not random fragmentation.
- Real spatial APA maps recovered for Cdk8, Apoe (Alzheimer's APOE), Gnb1l.
- Binning robustness quantified (cross-bin Pearson r 50→100/50→200 + Moran's I vs bin size).
- High-resolution presets: `highres_accuracy` / `highres_fast` with candidate defaults gp_blend=0.1, spatial=0.1, expression=0.7, apa=0.2.
- Caveat maintained: Visium HD FFPE/probe data are not APA-suitable (GSE311383 blocked on chemistry, kept as expression/domain-only candidate).

### 3.8 APA batch correction (optional module)

Quantile normalization + limma-style linear removal with collinearity auto-detection. Evaluated against baselines (QN vs Harmony comparison pending); deliberately labeled optional/under-evaluation rather than a headline claim.

---

## 4. Package Architecture

| Module | Role |
|--------|------|
| `spagapa.core` | APA dataset and APA site data structures |
| `spagapa.io` | BAM, coordinate, matrix, and result I/O; scAPAtrap integration |
| `spagapa.spatial` | KNN / radius / Delaunay spatial graphs |
| `spagapa.calling` | spatial validation and quality filtering |
| `spagapa.imputation` | dense/sparse GP, calibration, expression-GP prototypes |
| `spagapa.quantification` | RUD, PDUI, WUL, PAI, QC |
| `spagapa.analysis` | differential APA, domains, SVAPA, bias correction |
| `spagapa.visualization` | spatial/statistical/QC/benchmark plots |
| `spagapa.benchmark` | simulation + real benchmark infrastructure |
| `spagapa.bioml` | CPU-friendly multi-view graph/domain route |
| `spagapa.pipeline` | end-to-end workflow + preset dispatch |
| `spagapa.cli` / `spagapa.presets` | CLI and preset definitions |

Presets: `auto` (size/sparsity dispatch), `standard` (Visium-scale), `highres_accuracy`, `highres_fast`.

Tests: 439 collected (1 pre-existing collection error in `tests/benchmark/test_performance.py`).

---

## 5. Competitor and Reference Tools

- **scAPAtrap / Sierra** — PAS calling front ends. spaGAPA consumes either; the multi-caller validation shows its guarantees are caller-independent. Sierra requires a splice-junction BED (extracted from the BAM by pysam, ≥25 supporting reads) and a chr-consistent GTF.
- **stAPAminer** — primary spatial APA competitor (KNN-style imputation). Strong domain consistency via expression neighborhoods; no posterior uncertainty, limited calibration, less systematic scalability. spaGAPA does not chase it on every domain metric — GP where GP is strongest, BioML graph for biological consistency.
- **spvAPA** — imputation + supervised (sPLS-DA) selection; requires labels. Complementary supervised capability spaGAPA deliberately does not provide.
- **metaAPA** — cross-tool PAS integration workflow; complementary rather than competitive.
- **movAPA** — source of the MOB published APA/RUD matrix used for domain recovery.

Manuscript Table 1 (`docs/TABLE1_tool_comparison.md`) compares stAPAminer / spvAPA / metaAPA / spaGAPA across 10 capability rows, with scope notes distinguishing "not validated" from "out of scope".

---

## 6. Validation and Benchmark Status

All headline numbers below are the frozen, cross-checked versions used in the figures (see `scripts/main_figures/figure_index.md` for panel-level provenance).

### 6.1 Head-to-head imputation benchmark

5 methods (spaGAPA-GP, stAPAminer, spvAPA, spatial-KNN, mean) × 2 real datasets (GSE183456 kidney, GSE220442 brain), identical 20% per-gene masks (seed 42), parameter-fairness table published (Supp Fig S5).

- GP: RMSE 0.122, Pearson 0.942, 19 s (GSE183456) — beats stAPAminer (0.208) and spvAPA (0.272) on accuracy and is 5.7–7.2× faster.
- Mean: RMSE 0.080 — wins pointwise RMSE openly (bimodal index; mean predicts the dominant mode).
- Spatial fidelity: GP 0.42 vs mean 0.00 (a per-gene constant recovers no gradient by construction).
- ΔRMSE stratification: the GP−mean gap does not shrink with spatial signal (Moran's-I quintiles; Supp Fig S4) — the honest framing is "spatial reconstruction + calibrated uncertainty + scalable inference", not raw-RMSE dominance.

### 6.2 Conformal coverage validation

11 samples / 523,174 test points / 7 GSE / 4 tissues / 2 species: mean |deviation| 0.21 / 0.16 / 0.10 pp at 80/90/95%, max ≤ 0.5 pp, all 11 samples within ±5% (Fig 3; full table Supp Fig S6; per-sample deviation forest Supp Fig S15).

### 6.3 Uncertainty audits (reviewer-driven)

- **Leakage audit + LOOCV method-selection stability** (Supp Fig S7): calibration-set hygiene verified; method choice stable under leave-one-dataset-out.
- **Per-gene vs pooled correlation** (Supp Fig S8): within-gene median r ≈ 0.10–0.19; pooled r ≈ 0.5 partly cross-gene — reported transparently.
- **Conditional/subgroup coverage** (Fig 4E, Supp Fig S9): spatial-block splits stay at nominal; high-expression bins undercover (~0.82) — shown as an honest marginal-vs-conditional limitation.
- **Coverage–width tradeoff** (audit 4): interval width tracked, not just coverage.

### 6.4 Risk-coverage triage

Uncertainty filtering removes 23% of RMSE at 80% retention, paired across 5 datasets, p = 0.004 (Fig 4F; MAE variant Supp Fig S16).

### 6.5 MOB domain recovery (unsupervised)

ARI = 0.60, NMI = 0.68 vs 5 annotated layers. Fair comparison: same Leiden pipeline, only the APA source differs (mean- vs GP-imputed); mean-imputation also recovers substantial structure on MOB — reported honestly. Global-best config (expression_apa) moved to Supp Fig S10; confusion matrix replaces a prior Moran's-I panel whose caption had contradicted its data.

### 6.6 Scalability

Synthetic grids 1k–100k spots, 1200 s wall cap: spaGAPA-fast completes 100k in 511 s / 8.8 GB peak; stAPAminer and spvAPA time out or fail at 42k. Empirical power-law slopes reported (spaGAPA-fast 0.84, accuracy 1.01, stAPAminer 1.01, spvAPA 0.49) with the explicit note that deployed-implementation slopes can differ from textbook complexity. Spatial-fidelity–runtime Pareto: spatial-KNN highest fidelity (0.91) but non-probabilistic and non-scaling; spaGAPA-GP = moderate fidelity + uncertainty + scalability.

### 6.7 Stereo-seq pilot (GSE263789)

See §3.7. Figure 7 shows workflow, native-redraw APA QC, single-definition scale/sparsity (20.7M DNB · 21,455 PAS · 15,235 bin200 spots · ~10.3% observed), real Cdk8/Apoe/Gnb1l maps, and binning robustness. An earlier 41-domain + raw-σ uncertainty panel was removed (bin-grid × inducing-point geometry artifacts; over-segmented domains) — documented rather than hidden. AD-vs-WT effect sizes are descriptive only (n=1 per condition; Supp Fig S14, no p-values).

### 6.8 Multi-caller robustness

See §3.6. Answers the anticipated reviewer question "why scAPAtrap / how do you know its results are correct?" empirically.

### 6.9 Statistical-practice fixes (review-driven)

- Pseudoreplication eliminated: differential APA aggregated to donor/sample level (GSE220442 3v3; Supp Fig S11); GSE263789 AD-vs-WT explicitly labeled n=1 descriptive.
- Mean baseline included transparently in every benchmark where it wins pointwise RMSE.
- Publication freeze (2026-07-28) unified numbers across README/summary/manuscript; three data discrepancies and stale claims were fixed in that pass.

---

## 7. Real Dataset Inventory

### 7.1 Visium (8 GSE, 32 samples, scAPAtrap-processed, `data/processed/`)

| GSE | Tissue | Species | Samples | Notes |
|-----|--------|---------|---------|-------|
| GSE237183 | Glioma | Human | 18 | 18/18 PASS (GSM7596588 skipped: fiducial failure, Loupe-only fix; ZH881 covered by 4 other sections) |
| GSE183456 | Kidney | Human | 1 | fully ready (APA + expression + coords); kidney confirmed 3-way |
| GSE179572 | Brain metastasis | Human | 1 | first real-data benchmark; weak-label exploration |
| GSE220442 | AD brain (PFC) | Human | 6 (3+3) | non-standard SRA 4-read structure decoded (read3=R1 barcode, read4=R2 cDNA); donor-level diff-APA |
| GSE338525 | Liver | Human | 2 | normal |
| GSE206391 | Skin (psoriasis) | Human | 2 | very low PAS counts (229–276) |
| GSE169749 | Colon (DSS) | Mouse | 1 | multi-caller validation set |
| GSE263303 | Brain (Nf1) | Mouse | 1 | |

### 7.2 Stereo-seq (subcellular)

| GSE | Tissue | Species | Status |
|-----|--------|---------|--------|
| GSE263789 | AD brain | Mouse | end-to-end proven; 5 samples (AD/WT/3mo); outputs in `pipeline_output/gse263789_stereo_pilot/` |
| GSE269906 | AD brain | Human | feasibility only (STAR BAM 61 GB + 3'-bias check); chip mask not obtainable |

### 7.3 MOB (external APA matrix)

`data/processed/mob_st11` — ST11 mouse olfactory bulb, 260 spots, 5 layers, published movAPA-prepared APA/RUD matrix; domain-recovery anchor and non-scAPAtrap triangulation.

### 7.4 Blocked / supportive candidates

- GSE311383 (Visium HD liver): probe/FFPE chemistry not APA-suitable; expression/domain-only candidate.
- 10x V1_Human_Brain_Section_1, GSE153859: expression-only; excluded from APA claims per `docs/apa_evidence_source_manifest.md` decision rules.

3'-bias validation across the corpus: Mouse 53.4%, Human 47.0% within 500 bp of TES.

---

## 8. Manuscript and Figure Status

### 8.1 Manuscript

- Methods + Results (5,009 words), Abstract + Intro + Discussion (3,267 words), and the integrated full manuscript were drafted in-repo (git: `4fffc68`, `5bc3a5a`, `028b6dc`) with Table 1.
- The 2026-07-29 public-repo prune (`45588a2`) removed the manuscript drafts (and 20 other redundant docs) to keep the public repo code+figures+manifests; the text is recoverable from git history and maintained outside the repo.
- Known post-freeze figure-side revisions (honest-narrative rebuild, print-width pass, S14–S16 moves) are reflected in `figure_index.md` and the supplementary legends; the manuscript text should be re-checked against these before submission.

### 8.2 Figures (submission-ready)

- **Main 1–7** — `scripts/main_figures/fig{1..7}_*.py` + `figure_index.md`. Vector PDF + 300-DPI PNG, Okabe-Ito, DejaVu Sans, ~178 mm print width, bold panel letters, panels show data only (prose in captions).
- **Supplementary S1–S16** — `scripts/supplementary_figures/supp_fig{01..16}_*.py` + `supp_figure_legends.md` (submission-ready legends, added `ec275a5`).
- Rebuild history: per-review rebuild with real data only (`c12aa42`), fairness/provenance P0 fixes (`d353fa1`), supplementary print-width unification (`591adcf`).
- Reproduce: `OPENBLAS_NUM_THREADS=8 TMPDIR=/s3/mengzijun/tmp ~/anaconda3/envs/spagapa/bin/python figN_*.py` (S91; per-host env in CLAUDE.md).
- All figure source data live under `pipeline_output/` (gitignored, on-disk); data-source map at the end of `figure_index.md`.

### 8.3 Reproducibility artifacts

- `docs/REPRODUCIBILITY_MANIFEST.md` — environments, versions, per-dataset commands (2026-07-25).
- `docs/apa_evidence_source_manifest.md`, `docs/real_data_benchmark_protocol.md` — evidence rules and benchmark tracks.
- Processing logs `logs/2026*_数据处理记录.md`; runlogs/ for long suites.

---

## 9. Current Strengths

1. Clear biological problem (spatial APA remodeling) with a statistical angle (sparse, uncertain, spatially structured measurements).
2. Only spatial APA tool with conformal-calibrated uncertainty — and the guarantee is empirically caller-, species-, and tissue-invariant.
3. Scalable to 100k spots and to subcellular Stereo-seq; CPU-only, no GPU/deep-learning dependency.
4. Decoupled GP (values + uncertainty) / BioML graph (domains) architecture, with the negative result (expression-GP) kept as honest ablation.
5. Broad, real-data-only benchmark corpus: 9 GSE × 32 Visium samples + MOB + Stereo-seq, plus a 3-dataset caller-swap validation.
6. Review-driven statistical hygiene: pseudoreplication fixed, mean baseline reported transparently, conditional-coverage limits shown, n=1 analyses kept descriptive.
7. Fully version-controlled figure generation with print-width vector output and submission-ready legends.

---

## 10. Current Weaknesses and Honest Limitations

1. Per-gene mean wins entry-wise RMSE on the bimodal distal-usage index (holds on both callers); spaGAPA's claim is spatial fidelity, uncertainty, and scalability.
2. Within-gene uncertainty–error correlation is modest (median r ≈ 0.10–0.19); pooled r ≈ 0.5 is partly cross-gene ranking.
3. Marginal coverage ≠ conditional coverage: high-expression bins undercut (~0.82).
4. PAS overlap between callers is low (Jaccard 0.05–0.18) — caller-dependent site catalogs are a field-wide reality; spaGAPA mitigates this at the usage/guarantee level, not by unifying catalogs.
5. Stereo-seq raw FASTQ + mask availability is a field-wide bottleneck (1/66 GEO datasets); human Stereo-seq (GSE269906) is feasibility-only.
6. sAPA-RegNet perturbation model unvalidated (descriptive annotation only).
7. No supervised mode (spvAPA has sPLS-DA); batch-correction module still pending Harmony comparison.
8. Biological gold-standard labels remain limited (MOB layers are the main annotated anchor; pathology ROI labels never materialized).

---

## 11. Roadmap

### Short term (submission)

1. Final manuscript pass: re-check text numbers against the rebuilt figures (esp. conformal deviations 0.21/0.16/0.10 pp, multi-caller section, S14–S16 references).
2. Assemble submission package: main 1–7 + supp S1–S16 PDFs, legends, Table 1, cover letter.
3. Decide Zenodo/GitHub release tagging for the code + reproducibility manifest.

### Medium term (post-submission / revision-ready)

1. Batch-correction Harmony comparison (closes limitation 7).
2. Full GSE263789 study (beyond pilot) if compute allows; GSE269906 human upgrade remains mask-blocked.
3. CLI tutorial consolidation (`docs/user_guide.md`, `docs/tutorial_basic.ipynb` refresh).
4. Reviewer-response experiments from the BIB trail.

### Not planned

- Supervised selection module; deep-learning/GPU components; expression-GP promotion (all deliberately out of scope).

---

## 12. Publication Strategy

### 12.1 BIB positioning

- spaGAPA is a reproducible spatial APA analysis framework, not a single-model implementation.
- It systematically evaluates APA recovery, uncertainty, biological consistency, and scalability — with transparent baselines and negative results.
- CPU-friendly and accessible to ordinary bioinformatics labs; real-data workflows from Space Ranger/SAW BAM to APA matrix.

### 12.2 What not to claim

- "GP alone fully solves biological domain recovery."
- "spaGAPA beats every competitor on every metric."
- "Visium HD automatically supports APA calling."
- "Weak labels are pathology gold standards."
- "Conformal coverage is conditional/subgroup-wise."

### 12.3 Defensible claims (all benchmark-backed)

- Calibrated uncertainty unavailable in KNN-only tools, with distribution-free marginal coverage validated on 523k real test points and invariant to PAS caller, species, and tissue.
- GP > stAPAminer/spvAPA on spatial fidelity and runtime (5.7–7.2×), with the mean baseline reported honestly.
- BioML graph route delivers unsupervised domain recovery (MOB ARI 0.60) without labels.
- First complete APA pipeline on subcellular Stereo-seq (21,455 PAS × 20.7M DNB).
- Scales to 100k spots where both competitors fail at 42k under identical limits.

---

## 13. Final Project Direction

> Build spaGAPA as a CPU-friendly, uncertainty-aware, benchmark-driven spatial APA framework where GP handles APA recovery and uncertainty, while BioML graph learning handles biological-domain consistency — with every claim scoped to what real-data benchmarks support.

This direction is technically coherent, defensible to reviewers, aligned with the no-deep-learning/no-GPU constraint, and has now been carried through to a submission-ready state.

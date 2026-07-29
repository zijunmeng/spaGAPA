# Stereo-seq and Visium HD GEO candidates for spaGAPA

Date: 2026-07-13

This document records a platform-focused GEO search for Stereo-seq and 10x Visium HD datasets. The goal is to identify high-resolution spatial transcriptomics datasets that can support the spaGAPA BIB benchmark story.

## Practical Distinction

### Stereo-seq

Stereo-seq is sequencing-based high-resolution spatial transcriptomics. Many datasets include raw sequencing files or platform-specific processed files such as GEF/H5/H5AD and barcode-to-position files.

For spaGAPA:

- Better candidate for real high-resolution APA/PAS exploration.
- More likely to preserve transcriptome-wide read information.
- Preprocessing is more difficult than 10x Visium.
- Need to validate barcode/UMI/read structure before claiming APA compatibility.

### Visium HD

Visium HD is widely used and important for package adoption, but most datasets are FFPE/probe-style spatial gene expression.

For spaGAPA:

- Strong candidate for high-resolution expression/domain validation.
- Good for testing scalability, coordinate handling, binning, graph construction, and uncertainty-aware downstream.
- Weak as a primary APA/PAS raw-read benchmark unless the chemistry and raw reads can be shown to preserve 3' end evidence.

Conclusion:

```text
Stereo-seq should be prioritized for high-resolution APA benchmark.
Visium HD should be prioritized for high-resolution user-facing support and expression/domain validation.
```

## Stereo-seq Candidates

### Tier S1: Highest priority for spaGAPA

#### GSE269906 - Human prefrontal cortex AD/normal Stereo-seq

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269906

Why it matters:

- Platform: Stereo-seq.
- Tissue: human prefrontal cortex BA10.
- Design: six AD and six normal samples.
- Disease relevance: Alzheimer's disease.
- Strong biological labels: AD vs normal, human brain, cortical structure.
- GEO supplementary: `GSE269906_RAW.tar`.

Potential role in spaGAPA:

- Best candidate for a high-resolution human brain disease validation route.
- Could test whether spaGAPA recovers spatial APA/domain patterns in AD-relevant cortex.

Main risks:

- Need to inspect RAW.tar contents.
- GEO text did not expose SRA links in the quick series view; raw files may be supplementary rather than SRA.
- Stereo-seq preprocessing will require format-specific handling.

Recommended next step:

```text
Download/file-list only: inspect GSE269906_RAW.tar contents before downloading full data.
Prioritize if it contains FASTQ/BAM or GEF/H5 with coordinates.
```

#### GSE263789 - AD mouse model Stereo-seq spatial transcriptomics

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263789

Why it matters:

- Platform: Stereo-seq.
- Tissue: mouse brain.
- Disease model: amyloid plaque niche in AD model.
- GEO supplementary:
  - `GSE263789_RAW.tar`
  - `GSE263789_StereoSeq_mouse_combined.h5ad`
- Strong biological context: plaque niche, microglia-astrocyte interaction, AD mouse model.

Potential role in spaGAPA:

- Excellent high-resolution biological consistency dataset.
- H5AD may allow fast expression/domain ingestion before raw-read APA work.

Main risks:

- Mixed CosMx and Stereo-seq project; must isolate Stereo-seq subset.
- Need raw read or processed read-level evidence for APA.

Recommended next step:

```text
Inspect H5AD metadata first.
If coordinates and labels are clean, use as high-resolution expression/domain benchmark.
Then inspect RAW.tar for APA-compatible raw reads.
```

#### GSE256319 - Nova-ST / Stereo-seq mouse brain

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256319

Why it matters:

- Platform: Nova-ST derived from Stereo-seq.
- Tissue: mouse coronal brain section.
- GEO supplementary: `GSE256319_RAW.tar`.
- Previously identified as a strong high-resolution brain candidate.

Potential role in spaGAPA:

- Technical high-resolution benchmark.
- Good for multiscale binning and graph/coordinate handling.

Main risks:

- Nova-ST may have custom barcode/array structure.
- Need to inspect RAW.tar for GEF/H5/barcodeToPos/raw reads.

Recommended next step:

```text
Inspect RAW.tar file list; if it includes counts.gef.h5/barcodeToPos.h5 or FASTQ, build an ingestion pilot.
```

### Tier S2: Strong but less direct

#### GSE299386 - Comparative cortical architecture, Stereo-seq

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299386

Why it matters:

- Platform: Stereo-seq.
- Tissue: opossum and mouse primary visual cortex.
- Biological value: cortical architecture and evolutionary comparison.
- GEO supplementary: `GSE299386_RAW.tar`.

Use case:

- Good cortical architecture / layer-like biological consistency benchmark.

Risks:

- Cross-species design adds complexity.
- Not necessarily easiest first raw-read pilot.

#### GSE333479 - Unified Stereo-seq atlas of ten mouse organs

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333479

Why it matters:

- Platform: Stereo-seq.
- Tissue: ten mouse organs including brain, kidney, lung, thymus, intestine, skin, spleen, ovary, testis, uterus.
- GEO supplementary: `GSE333479_RAW.tar`.
- Very useful for broad user-facing high-resolution support.

Use case:

- Good package robustness dataset.
- Good for showing spaGAPA handles multiple Stereo-seq tissues.

Risks:

- Published/public in 2026; likely large.
- Atlas scope may be too broad for first APA benchmark.

#### GSE297068 - Olfactory receptor spatial code, nose and brain

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297068

Why it matters:

- Includes spatial data related to olfactory epithelium/brain mapping.
- GEO supplementary includes `GSE297068_RAW.tar` and processed CSV files.

Use case:

- Potential biological structure benchmark.

Risks:

- Mixed methods including MERFISH and scRNA; must isolate the sequencing-based spatial component.

### Tier S3: Useful for tumor/high-res user support

#### GSE328481 - Stereo-XCR-seq lung adenocarcinoma

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328481

Why it matters:

- Platform: Stereo-XCR-seq.
- Tissue: fresh-frozen LUAD patient tissues.
- GEO supplementary: `GSE328481_RAW.tar`.
- Tumor immune niche / TCR-BCR enrichment.

Use case:

- High-resolution tumor spatial immune benchmark.

Risks:

- Raw reads are reported as available in GSA, not necessarily GEO/SRA.
- XCR enrichment complicates APA interpretation.

#### GSE317755 - Axolotl brain Stereo-seq

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317755

Why it matters:

- Platform: Stereo-seq.
- Tissue: axolotl brain.
- Interesting regeneration biology.

Use case:

- Not first-line for spaGAPA due non-model organism annotation complexity.

## Visium HD Candidates

### Tier V1: Best match for spaGAPA high-resolution validation

#### GSE290724 - Visium HD glioblastoma infiltration

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290724

Why it matters:

- Platform: 10x Genomics Visium HD.
- Tissue: GBM patient tissue.
- Design: peritumoral region, tumor edge, tumor core.
- GEO supplementary: `GSE290724_RAW.tar`.
- Strong biological/domain labels for tumor infiltration.

Potential role in spaGAPA:

- Best Visium HD candidate for high-resolution domain validation.
- Strong fit for biological consistency, graph-domain detection, and uncertainty-aware domain stability.

APA caveat:

- Visium HD is likely not a primary APA/PAS raw-read benchmark.
- Use for high-resolution expression/domain route unless raw chemistry proves APA-compatible.

#### GSE311383 - Human liver Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311383

Why it matters:

- Platform: 10x Genomics Visium HD.
- Tissue: healthy human liver.
- GEO supplementary: `GSE311383_RAW.tar`.
- Strong anatomical structure in liver.

Potential role:

- High-resolution tissue architecture / coordinate scaling benchmark.

#### GSE318916 - High-grade T1 bladder cancer Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318916

Why it matters:

- Platform: Visium HD.
- Tissue: FFPE high-grade T1 bladder cancer.
- GEO supplementary: `GSE318916_RAW.tar`.
- Tumor heterogeneity and immune microenvironment.

Potential role:

- Tumor domain and high-resolution package usability benchmark.

### Tier V2: Useful secondary Visium HD candidates

#### GSE301973 - EGFR-mutant NSCLC Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301973

Why it matters:

- Platform: Visium HD.
- Tissue: NSCLC before/after osimertinib.
- GEO supplementary: `GSE301973_RAW.tar`.

Use case:

- Tumor microenvironment and treatment comparison.

#### GSE280315 - Colorectal cancer Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280315

Why it matters:

- Platform: Visium HD.
- Tissue: human CRC and adjacent mucosa.
- GEO supplementary: `GSE280315_RAW.tar`.

Use case:

- High-resolution tumor/normal domain validation.

#### GSE314260 - Human trunk embryoid model Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314260

Why it matters:

- Platform: Visium HD.
- Model: patterned human trunk embryoid.
- GEO supplementary includes `GSE314260_RAW.tar` and RDS objects.

Use case:

- Spatial patterning/domain validation.

Limitation:

- Organoid model, less ideal than tissue benchmark for BIB biological validation.

#### GSE293042 - Human/mouse adrenal gland Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293042

Why it matters:

- Platform: Visium HD.
- Tissue: adrenal gland.
- GEO supplementary: `GSE293042_RAW.tar`.

Use case:

- Tissue organization and species comparison.

#### GSE260926 - Developing human knee Visium HD

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260926

Why it matters:

- Platform: Visium HD.
- Tissue: developing human knee.
- GEO text reports raw FASTQ files and Space Ranger output.
- GEO supplementary: `GSE260926_RAW.tar`.

Use case:

- Strong for user-facing Visium HD input support.

Limitation:

- Not brain/tumor; APA suitability uncertain.

## Recommended Strategy

### For high-resolution APA evidence

Prioritize Stereo-seq:

```text
1. GSE269906 - human PFC AD/normal Stereo-seq
2. GSE263789 - AD mouse model Stereo-seq
3. GSE256319 - Nova-ST/Stereo-seq mouse brain
4. GSE299386 - mouse/opossum visual cortex Stereo-seq
```

Why:

```text
These are sequencing-based high-resolution ST datasets and are more plausible for APA/PAS raw-read analysis.
```

### For high-resolution package support and biological consistency

Prioritize Visium HD:

```text
1. GSE290724 - GBM infiltration, Visium HD
2. GSE311383 - human liver, Visium HD
3. GSE318916 - bladder cancer, Visium HD
4. GSE280315 - colorectal cancer, Visium HD
```

Why:

```text
These datasets are likely more user-facing and easier to support as high-resolution input formats.
They are excellent for domain detection, graph construction, binning, runtime, memory, and uncertainty-aware downstream validation.
They should not be overclaimed as APA raw-read benchmarks unless read chemistry supports APA calling.
```

## Immediate Next Step

Run metadata/file-list audits for the top two from each platform:

```text
Stereo-seq:
  GSE269906
  GSE263789

Visium HD:
  GSE290724
  GSE311383
```

For each:

```text
1. Download or inspect GEO supplementary file list.
2. Identify whether RAW.tar contains FASTQ, BAM, GEF/H5/H5AD, Parquet, Space Ranger output, or only matrices.
3. Record file sizes.
4. Decide whether the dataset is:
   A. APA raw-read candidate;
   B. high-resolution expression/domain validation candidate;
   C. unsuitable or too costly.
```


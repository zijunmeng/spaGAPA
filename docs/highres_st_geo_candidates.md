# High-resolution ST GEO candidates for spaGAPA

Date: 2026-07-13

This document records a first-pass GEO search for high-resolution spatial transcriptomics datasets that may support the spaGAPA high-resolution benchmark narrative. The filtering criterion is not just spatial resolution; the dataset must be plausibly useful for APA-aware analysis.

## Selection Criteria

Primary requirements:

1. Sequencing-based spatial transcriptomics, not purely image-based targeted FISH.
2. Raw FASTQ/SRA or barcode/UMI-preserving BAM available.
3. Spot/bead/bin/cell coordinates available.
4. Poly(A), cDNA, or whole-transcriptome capture that can plausibly preserve 3' end evidence.
5. Tissue structure labels or interpretable anatomy available.
6. CPU-feasible preprocessing path.

Exclusion or low-priority criteria:

1. MERFISH/Xenium/CosMx-only data: high spatial resolution but not APA-calling friendly.
2. Probe-based FFPE assays: useful for expression/domain benchmarking, but may not preserve true APA/PAS evidence.
3. Datasets with processed matrices only and no usable raw reads.

## Priority A: Strong Candidates

### 1. GSE268519 - Slide-seqV2 mouse hippocampus/cerebellum

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268519

Why it is strong:

- Platform: Slide-seqV2.
- Tissue: mouse hippocampus and cerebellum.
- Raw data: SRA available.
- Library source/selection: transcriptomic single-cell / cDNA.
- Biological labels: hippocampus and cerebellum have strong anatomical structure.
- Good fit for spaGAPA: high-resolution bead-level sequencing data, brain tissue, suitable for multiscale binning.

Risks:

- Requires Slide-seq-specific preprocessing.
- Supplementary files are RDS/RCTD objects rather than a simple Space Ranger output.
- Need to confirm FASTQ read structure and barcode/UMI extraction.

Recommended first action:

Download SRA metadata for one hippocampus puck, then inspect FASTQ read structure and processed RDS contents.

### 2. GSE130682 - HDST mouse olfactory bulb

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130682

Why it is strong:

- Platform: High-Density Spatial Transcriptomics (HDST).
- Tissue: adult mouse main olfactory bulb.
- Resolution: reported 2 um.
- Raw data: SRA available.
- Processed files: supplementary TAR of TIFF/TSV.
- Biological labels: olfactory bulb has clear layered anatomy.

Risks:

- Older/custom HDST pipeline.
- Need to reconstruct or parse bead barcode to coordinate mapping.
- May require more custom preprocessing than 10x/Slide-seq.

Recommended first action:

Download supplementary TAR first, inspect TSV coordinate/count files, then decide whether raw SRA is needed.

### 3. GSE256319 - Nova-ST / Stereo-seq mouse brain

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256319

Why it is strong:

- Platform: Nova-ST / Stereo-seq.
- Tissue: mouse brain.
- Raw data: SRA available.
- Processed files: RAW TAR includes H5/TAR; sample-level files include barcodeToPos.h5 and counts.gef.h5.
- Library: polyA RNA reported at sample level.
- Strong high-resolution narrative value.

Risks:

- BGI/Stereo-seq format requires SAW/GEF/H5 parsing.
- Large files: RAW TAR approximately 12 GB; sample files can be multiple GB.
- APA calling from raw reads may require custom handling of Stereo-seq barcodes.

Recommended first action:

Start from processed counts.gef.h5 + barcodeToPos.h5 to validate coordinate/expression ingestion before raw FASTQ.

### 4. GSE251926 - Open-ST 3D high-resolution ST

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251926

Why it is strong:

- Platform: Open-ST, sequencing-based high-resolution ST.
- Tissues: adult mouse hippocampus, mouse head, primary HNSCC, healthy/metastatic lymph node.
- Raw data: SRA available.
- Processed files: H5AD files with coordinates, staining image, segmentation mask, and cell type annotation.
- Strong for tumor/domain and 3D spatial narrative.

Risks:

- Very large: RAW TAR approximately 35 GB; metastatic 3D H5AD approximately 2.2 GB.
- Custom Open-ST/spacemake pipeline.
- More complex than Slide-seqV2 for immediate APA calling.

Recommended first action:

Use adult mouse hippocampus or one metastatic lymph node section H5AD as expression/domain ingestion test; defer full raw FASTQ until format is understood.

## Priority B: Useful But Not First Choice

### 5. GSE169706 - Seq-Scope mouse liver/colon

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169706

Why it is useful:

- Platform: Seq-Scope.
- Resolution: submicrometer-scale spatial barcoding.
- Raw data: SRA available.
- Processed files: supplementary TAR of MTX/TSV.
- Strong technical high-resolution benchmark.

Limitations:

- Tissue is liver/colon, not brain.
- Biological domain labels may be weaker than brain layer/region datasets.
- Seq-Scope workflow has two sequencing steps and custom barcode handling.

Recommended use:

Use as a second-phase technical scalability dataset, not the first high-resolution biological benchmark.

### 6. GSE299816 - 10x Visium HD mouse tibia

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299816

Why it is useful:

- Platform: 10x Visium HD.
- Raw data: SRA available.
- Processed files: RAW TAR includes CLOUPE/H5/HTML/JPG/JSON/PARQUET/PNG/TIFF.
- Easiest among high-resolution platforms for 10x-style parsing.

Limitations:

- Tissue is mouse tibia, not brain/tumor.
- If probe-based FFPE, APA/PAS evidence may not be reliable.
- More suitable for high-resolution expression/domain validation than APA calling.

Recommended use:

Use only if we need a 10x-format high-resolution smoke test.

## Not Suitable As APA Benchmarks

### GSE246717 - whole mouse brain MERFISH atlas

URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246717

Reason:

- Very valuable high-resolution brain atlas.
- However, MERFISH is image-based/targeted and not suitable for APA site calling.
- Could be used as anatomical/cell-type reference, not as a spaGAPA APA benchmark.

## Recommended Order

1. GSE268519 Slide-seqV2 hippocampus/cerebellum.
2. GSE130682 HDST olfactory bulb.
3. GSE256319 Nova-ST/Stereo-seq mouse brain.
4. GSE251926 Open-ST adult mouse hippocampus or metastatic lymph node.
5. GSE169706 Seq-Scope as a technical stress test.
6. GSE299816 Visium HD as a 10x-format smoke test only.

## Immediate Next Step

The best next action is to run a metadata-only audit for:

```text
GSE268519
GSE130682
GSE256319
GSE251926
```

For each dataset, collect:

```text
GSE
GSM
SRX/SRR
organism
tissue
library layout
FASTQ URLs
read lengths
supplementary file list
coordinate file type
processed matrix format
expected preprocessing route
APA suitability
```

Then select one dataset for a small raw-read pilot.

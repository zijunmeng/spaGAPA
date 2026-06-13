# APA Evidence Source Manifest

This manifest records candidate real-data sources for spaGAPA's BIB-oriented
APA benchmark. The goal is to separate true APA/PAS evidence from expression-only
spatial data, and to make the next download/calling step reproducible.

Last updated: 2026-06-13

## 1. Decision Rules

A dataset can enter the true APA benchmark only when it has at least one of:

- published APA usage matrix;
- published PAS/poly(A) site table with site-by-spot or site-by-cell counts;
- public BAM suitable for scAPAtrap, Sierra, polyApipe, SCAPE/metaAPA, or a
  similar APA caller;
- public FASTQ that can be aligned into such a BAM.

Expression matrices, spatial coordinates, molecule-info files, or long-read
isoform matrices are useful supporting evidence, but they are not sufficient
for `apa_matrix.csv` unless a real APA/PAS calling step is performed.

## 2. Source Status Summary

| Source | Tissue | Evidence Found | APA Benchmark Status | Decision |
| --- | --- | --- | --- | --- |
| `stapaminer_mob` | mouse olfactory bulb | prepared APA/RUD matrix, layers, site-level partial evidence | ready | keep as current real-data anchor |
| 10x `V1_Human_Brain_Section_1` | human brain | expression H5, spatial files, molecule info | expression-only candidate | do not count until BAM/FASTQ/PAS is found |
| `GSE153859` CBS1/CBS2 | mouse hippocampal brain | GEO expression/spatial files, SRA Illumina FASTQ source, Nanopore isoform matrices | APA-callable candidate | use SRA Illumina FASTQ/BAM route |

## 3. 10x `V1_Human_Brain_Section_1`

Local candidate directory:

```text
data/processed/visium_human_brain_section1_candidate/
```

Checked canonical 10x CDN paths:

| File | URL Status | Notes |
| --- | --- | --- |
| `V1_Human_Brain_Section_1_filtered_feature_bc_matrix.h5` | available | expression matrix |
| `V1_Human_Brain_Section_1_spatial.tar.gz` | available | spatial images and coordinates |
| `V1_Human_Brain_Section_1_raw_feature_bc_matrix.h5` | available | raw expression matrix |
| `V1_Human_Brain_Section_1_molecule_info.h5` | available, 129,941,082 bytes | UMI molecule support, not sufficient for APA site calling |
| `V1_Human_Brain_Section_1_possorted_genome_bam.bam` | not publicly available at canonical path | cannot run APA callers from this source |
| `V1_Human_Brain_Section_1_possorted_genome_bam.bam.bai` | not publicly available at canonical path | no index |

Conclusion:

- Keep this dataset as an expression/spatial ingestion candidate.
- Do not use it in true APA benchmark claims.
- Revisit only if a public BAM/FASTQ/PAS table is found.

## 4. `GSE153859`: Mouse Brain/Hippocampus Candidate

Primary source:

- GEO: `GSE153859`
- SRA study: `SRP270322`
- BioProject: `PRJNA644362`
- Paper title on GEO: "The spatial landscape of gene expression isoforms in
  tissue sections"

Local files downloaded:

```text
data/raw/gse153859/GSE153859_RAW.tar
data/raw/gse153859/GSE153859_runinfo.csv
```

Downloaded sizes:

```text
GSE153859_RAW.tar       176,271,360 bytes
GSE153859_runinfo.csv       8,214 bytes
```

`GSE153859_RAW.tar` contents:

```text
GSM4656179_CBS2-illumina.tar.gz
GSM4656180_CBS1-illumina.tar.gz
GSM4656181_MOB-illumina.tar.gz
GSM4656182_CBS2-SiT.tar.gz
GSM4656183_CBS1-SiT.tar.gz
GSM4656184_MOB-SiT.tar.gz
```

Nested Illumina contents:

```text
CBS1/CBS2/MOB illumina:
  filtered_feature_bc_matrix.h5
  spatial/scalefactors_json.json
  spatial/tissue_hires_image.png
  spatial/tissue_lowres_image.png
  spatial/tissue_positions_list.csv
```

Nested Nanopore/SiT contents:

```text
CBS1/CBS2/MOB SiT:
  *_genematrix.txt
  *_isomatrix.txt
  *_juncmatrix.txt
  *_snpmatrix.csv
```

Interpretation:

- The GEO RAW tar gives expression/spatial scaffolds and long-read isoform
  support.
- It does not provide a direct PAS/poly(A) call table.
- It does not provide a BAM.
- The SRA Illumina runs are therefore the real APA entry point.

## 5. `GSE153859` Illumina Runs

| GEO Sample | Biological Sample | SRA Run | Platform | Layout | SRA Size | APA Route |
| --- | --- | --- | --- | --- | --- | --- |
| `GSM4656179` | CBS2 Illumina | `SRR12157782` | NextSeq 500 | paired | 7,482 MB | FASTQ -> BAM -> APA caller |
| `GSM4656180` | CBS1 Illumina | `SRR12157783` | NextSeq 500 | paired | 7,742 MB | FASTQ -> BAM -> APA caller |
| `GSM4656181` | MOB Illumina | `SRR12157784` | NextSeq 500 | paired | 8,426 MB | FASTQ -> BAM -> APA caller |

Recommended next target:

```text
GSM4656180 / CBS1 / SRR12157783
```

Rationale:

- It is a brain/hippocampal Visium sample, giving a stronger second-tissue
  benchmark than another MOB-only replicate.
- metaAPA explicitly used CBS1/CBS2 as spatial transcriptomics samples with
  BAM input, making them directly relevant competitors/precedents.
- Once APA calls are generated, CBS2 can become replicate validation.

## 6. Reproducible Commands

Download metadata and RAW tar:

```bash
mkdir -p data/raw/gse153859
wget -O data/raw/gse153859/GSE153859_runinfo.csv \
  'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP270322'
wget -c -O data/raw/gse153859/GSE153859_RAW.tar \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153859/suppl/GSE153859_RAW.tar'
```

Inspect nested file contents:

```bash
tar -tf data/raw/gse153859/GSE153859_RAW.tar
tar -xOf data/raw/gse153859/GSE153859_RAW.tar GSM4656180_CBS1-illumina.tar.gz | tar -tzf -
tar -xOf data/raw/gse153859/GSE153859_RAW.tar GSM4656183_CBS1-SiT.tar.gz | tar -tzf -
```

Planned SRA download and FASTQ conversion for CBS1:

```bash
mkdir -p data/raw/gse153859/sra data/raw/gse153859/fastq/CBS1
prefetch SRR12157783 --output-directory data/raw/gse153859/sra
fasterq-dump data/raw/gse153859/sra/SRR12157783/SRR12157783.sra \
  --split-files \
  --threads 8 \
  --temp /tmp \
  --outdir data/raw/gse153859/fastq/CBS1
```

Before running the SRA download, check available storage. The SRA file is about
7.7 GB and the expanded FASTQ plus downstream BAM may require tens of GB.

## 7. Planned APA Calling Route

Preferred route:

1. Download `SRR12157783` FASTQ.
2. Align and barcode/UMI-preserve the reads into a BAM compatible with APA
   callers.
3. Run one or more APA callers:
   - scAPAtrap, matching stAPAminer as closely as possible;
   - Sierra or polyApipe, for metaAPA-style caller comparison;
   - optional metaAPA integration if multiple caller outputs are available.
4. Export spaGAPA-ready files:

```text
data/processed/gse153859_cbs1_apa/
  apa_matrix.csv
  apa_sites.csv
  apa_site_counts.csv
  coordinates.csv
  expression_matrix.csv
  metadata.csv
  qc_summary.json
```

Acceptance criteria:

- `apa_matrix.csv` is derived from real PAS/poly(A) evidence, not expression
  smoothing.
- Coordinates align exactly to APA matrix columns.
- At least one site-level table is retained for audit.
- Biological region labels are added before using layer/domain metrics.

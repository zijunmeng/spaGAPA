# Real-Data Benchmark Protocol

This document defines the real-data benchmark standard for spaGAPA's
BIB-targeted validation. Its purpose is to keep real-data evidence reproducible,
auditable, and clearly separated from simulation or expression-derived proxies.

## 1. Scope

spaGAPA real-data benchmarks must evaluate spatial APA measurements, not only
spatial gene expression. A dataset can enter different tracks depending on the
available evidence:

| Track | Required evidence | Counts toward BIB real-data benchmark |
| --- | --- | --- |
| Prepared spatial APA | APA usage matrix + coordinates | Yes |
| External biological validation | Prepared APA + biological labels | Yes |
| High-resolution validation | Prepared APA + expression + labels | Yes |
| Site-level APA readiness | poly(A) sites + site counts | Supportive |
| Expression-only candidate | expression + coordinates only | No |

Expression-only Visium, Stereo-seq, Slide-seq, MERFISH, or Xenium datasets are
useful for ingestion and coordinate handling tests, but they must not be
reported as true APA benchmarks unless real APA/PAS evidence is added.

## 2. Dataset Directory Layout

Each prepared dataset lives under:

```text
data/processed/<dataset_id>/
```

Use lowercase, underscore-separated identifiers, for example:

```text
data/processed/stapaminer_mob/
data/processed/visium_human_brain_section1_candidate/
```

### 2.1 Required For Prepared APA

```text
apa_matrix.csv
coordinates.csv
```

`apa_matrix.csv`:

- rows: genes, APA events, or APA sites depending on the analysis level.
- columns: spots/cells/bins.
- values: APA usage such as RUD, PDUI, WUL, distal usage, or normalized site
  usage.
- missing values: `NaN` for unobserved or insufficient-coverage events.

`coordinates.csv`:

- one row per spot/cell/bin.
- required columns:
  - `spot_id`
  - `x`
  - `y`

The column names in `apa_matrix.csv` must match `coordinates.csv:spot_id`.

### 2.2 Required For Biological Validation

```text
metadata.csv
```

Recommended columns:

- `spot_id`
- `layer`, `region`, `domain`, or equivalent biological label.
- `dataset`
- `source`
- `sample_id` or `replicate` when available.

The default validation label column is `layer`.

### 2.3 Required For Expression-Aware Baselines

```text
expression_matrix.csv
```

Format:

- rows: genes.
- columns: spots/cells/bins.
- values: raw counts, normalized expression, or log-normalized expression.

This file enables stAPAminer-like expression-KNN and spaGAPA BioML/highres
routes. Expression alone is not an APA measurement.

### 2.4 Recommended For stAPAminer Alignment

```text
stapaminer_rud_raw.csv
stapaminer_rud_imputed.csv
```

These files are optional but valuable when reproducing stAPAminer exactly.

### 2.5 Recommended For metaAPA / Site-Level Validation

```text
apa_sites.csv(.gz)
apa_site_counts.csv(.gz)
```

`apa_sites.csv` should contain site-level annotations such as:

- `site_id`
- `gene`
- `chrom`
- `start`
- `end`
- `strand`
- `source_caller`

`apa_site_counts.csv` should contain site-by-spot or site-by-cell counts.

For full metaAPA-style validation, also record:

- caller outputs from Sierra/polyApipe/SCAPE or equivalent.
- genome FASTA/GTF versions.
- PAS motif or sequence-context evidence.
- optional long-read support.

## 3. Provenance

Each prepared dataset should include:

```text
qc_summary.json
```

Recommended fields:

```json
{
  "dataset": "dataset_id",
  "source": "publication or repository",
  "platform": "Visium / Stereo-seq / Slide-seq / ...",
  "species": "mouse / human / ...",
  "tissue": "brain / MOB / tumor / ...",
  "apa_source": "published_APA_index / PAS_counts / BAM_calling / expression_only_candidate",
  "apa_ready": true,
  "expression_ready": true,
  "biological_labels_ready": true,
  "site_level_ready": false,
  "n_spots": 0,
  "n_genes": 0,
  "notes": []
}
```

## 4. Inclusion Criteria

A dataset counts as a BIB real-data APA benchmark only if:

1. `apa_matrix.csv` exists and is derived from real APA/PAS evidence.
2. `coordinates.csv` exists and aligns to APA columns.
3. The source, preprocessing commands, and filtering thresholds are recorded.
4. The benchmark uses the same genes/spots/masks across methods.

For biological-consistency claims, the dataset must also contain:

1. `metadata.csv`.
2. A biological label column such as `layer` or `region`.

For high-resolution claims, the dataset should contain:

1. enough spots/cells/bins or a defensible pseudo-bin stress test;
2. `expression_matrix.csv`;
3. biological labels for validation.

## 5. Exclusion Criteria

Do not count a dataset as true APA benchmark if:

1. APA values are simulated from expression.
2. APA values are random, smoothed, or generated only for software smoke tests.
3. only gene-expression matrices are available.
4. spot coordinates cannot be aligned to APA columns.
5. the source of APA/PAS calls is unknown.

Expression-only datasets may be stored as candidates using:

```text
qc_summary.json: apa_ready = false
qc_summary.json: apa_source = "expression_only_candidate"
```

## 6. Benchmark Tracks

### 6.1 Readiness

Command:

```bash
conda run -n spagapa python scripts/run_real_data_suite.py \
  --processed-root data/processed \
  --output-dir benchmark_results/real/real_data_suite_v1
```

Outputs:

```text
real_data_readiness.csv
real_data_readiness_summary.json
decision_summary.json
```

### 6.2 External Biological Validation

Command:

```bash
conda run -n spagapa python scripts/run_real_data_suite.py \
  --processed-root data/processed \
  --output-dir benchmark_results/real/real_data_suite_v1 \
  --run-external
```

Core methods:

- raw
- stAPAminer original imputed, when available
- stAPAminer-like expression KNN
- spaGAPA GP
- spaGAPA BioML/highres

### 6.3 High-Resolution Default Validation

Command:

```bash
conda run -n spagapa python scripts/run_real_data_suite.py \
  --processed-root data/processed \
  --output-dir benchmark_results/real/real_data_suite_v1 \
  --run-highres-default
```

This compares current/candidate high-resolution defaults over pseudo-bin
density, dropout stress, and multiple random seeds.

## 7. Current BIB Gate

spaGAPA should not claim BIB-scale real-data validation until the suite reaches:

```text
n_external_validation_ready >= 3
n_highres_validation_ready >= 1
n_site_level_partial_ready >= 2
```

Current status is tracked by:

```text
benchmark_results/real/real_data_suite_v1/real_data_readiness_summary.json
```

## 8. Second Dataset Strategy

Priority order:

1. brain/layer dataset with published APA/PAS usage or accessible BAM/FASTQ for
   APA calling;
2. tumor or tissue-region dataset with biological region labels;
3. high-resolution dataset with 3-prime/poly(A)-informative reads.

For a brain Visium dataset with expression and coordinates but no APA evidence:

1. standardize it as an expression-only candidate;
2. keep `apa_ready=false`;
3. do not include it in APA benchmark claims;
4. later add APA by running scAPAtrap/metaAPA/Sierra/polyApipe on BAM/FASTQ or
   by importing published APA/PAS calls.

## 9. Current APA Evidence Manifest

Detailed evidence checks are tracked in:

```text
docs/apa_evidence_source_manifest.md
```

Current conclusions:

1. The 10x `V1_Human_Brain_Section_1` sample has expression/spatial files and
   `molecule_info.h5`, but no publicly accessible BAM at the canonical 10x CDN
   path. It remains an expression-only candidate.
2. `GSE153859` was audited as a possible second brain/layer candidate:
   - GEO `GSE153859`
   - SRA `SRP270322`
   - BioProject `PRJNA644362`
   - CBS1 Illumina `SRR12157783`
   - CBS2 Illumina `SRR12157782`
3. `GSE153859_RAW.tar` has been downloaded locally and contains Visium
   expression/spatial files plus Nanopore isoform matrices, but no direct
   PAS/APA call table.
4. CBS1 Illumina `SRR12157783` has been downloaded through the public SRA S3
   object route and converted to FASTQ:
   `data/raw/gse153859/fastq/CBS1/SRR12157783_1.fastq.gz`.
5. The conversion produced one FASTQ file rather than a conventional paired
   FASTQ pair. The next true-APA step is therefore a read-structure audit:
   inspect whether spatial barcode and UMI information are represented in read
   names, sequence segments, or another SRA-derived layout.
6. Initial audit found plain SRA headers and uniform 91 bp reads in a 1,000-read
   sample, with no obvious barcode/UMI field. RunInfo also reports
   `spots_with_mates=0`. This may be an SRA-lite/cDNA-only object rather than a
   full Visium R1/R2 submission.
7. If original R1/R2 FASTQ or BAM can be located, generate BAM/PAS calls using
   scAPAtrap, Sierra, polyApipe, or metaAPA-style integration.
8. Since the public SRA object appears to expose only a cDNA-like read, keep
   `GSE153859` as expression/spatial biological support unless original R1/R2
   FASTQ or BAM is found outside the current public SRA route.

Do not count `GSE153859` as APA-ready until the resulting `apa_matrix.csv` and
site-level audit files are produced.

### `GSE179572` Human Brain Metastasis Visium Replacement Candidate

`GSE179572` is now the preferred second real-data APA calling candidate.

Reasons:

1. It is fresh frozen human brain metastasis Visium, not FFPE probe-based
   spatial data.
2. GEO SOFT records six SRA-linked samples with matched expression/spatial
   supplementary files.
3. SRA XML confirms Visium-compatible public FASTQ evidence:
   - `R1` / barcode-UMI read: 28 bp;
   - `R2` / cDNA read: 90 bp;
   - some lanes also include two 10 bp index reads.
4. The smallest recommended validation run is:
   `GSM5420751 / SRX11362761 / SRR15052395`, approximately 4.8 GB.

Minimal validation command:

```bash
mkdir -p data/raw/gse179572/sra data/raw/gse179572/fastq/GSM5420751
aria2c \
  --dir=data/raw/gse179572/sra \
  --out=SRR15052395.sra \
  --check-certificate=false \
  --max-connection-per-server=8 \
  --split=8 \
  --min-split-size=20M \
  --continue=true \
  'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15052395/SRR15052395'
fasterq-dump \
  -O data/raw/gse179572/fastq/GSM5420751 \
  -t /tmp/spagapa_fasterq_srr15052395 \
  -e 8 \
  -p \
  --split-files \
  --include-technical \
  data/raw/gse179572/sra/SRR15052395.sra
```

Expected decision:

- This validation passed on 2026-06-14.
- `fasterq-dump` produced 4 FASTQ files from `SRR15052395`:
  - `_1.fastq`: 28 bp barcode/UMI read;
  - `_2.fastq`: 90 bp cDNA read;
  - `_3.fastq`: 10 bp index read;
  - `_4.fastq`: 10 bp index read.
- Proceed to Space Ranger/Cell Ranger BAM generation and
  scAPAtrap/Sierra/polyApipe APA calling.
- `GSE179572` now replaces `GSE153859` as the active second true-APA dataset
  route.

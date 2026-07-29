# APA Evidence Source Manifest

This manifest records candidate real-data sources for spaGAPA's BIB-oriented
APA benchmark. The goal is to separate true APA/PAS evidence from expression-only
spatial data, and to make the next download/calling step reproducible.

Last updated: 2026-06-14

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
| `GSE153859` CBS1/CBS2 | mouse hippocampal brain | GEO expression/spatial files, SRA Illumina cDNA-only-like FASTQ source, Nanopore isoform matrices | expression/spatial support | public SRA object appears to lack Visium barcode/UMI mate |
| `GSE179572` | human brain metastasis | fresh frozen Visium expression/spatial files plus public SRA R1/R2 FASTQ | APA-callable candidate | prioritize as second real-data APA source |

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

Completed SRA download and FASTQ conversion for CBS1:

```bash
mkdir -p data/raw/gse153859/sra_aria2 data/raw/gse153859/fastq/CBS1
printf '%s\n' \
  'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12157783/SRR12157783' \
  > data/raw/gse153859/sra_aria2/SRR12157783.links.txt

aria2c \
  --input-file=data/raw/gse153859/sra_aria2/SRR12157783.links.txt \
  --dir=data/raw/gse153859/sra_aria2 \
  --check-certificate=false \
  --max-concurrent-downloads=1 \
  --max-connection-per-server=8 \
  --split=8 \
  --min-split-size=20M \
  --continue=true \
  --auto-file-renaming=false

mv data/raw/gse153859/sra_aria2/SRR12157783 \
  data/raw/gse153859/sra_aria2/SRR12157783.sra

fasterq-dump \
  -O /tmp/spagapa_fasterq_srr12157783/out \
  -t /tmp/spagapa_fasterq_srr12157783/tmp \
  -e 8 \
  -p \
  --split-files \
  --include-technical \
  data/raw/gse153859/sra_aria2/SRR12157783.sra

pigz -p 8 /tmp/spagapa_fasterq_srr12157783/out/SRR12157783_1.fastq
mv /tmp/spagapa_fasterq_srr12157783/out/SRR12157783_1.fastq.gz \
  data/raw/gse153859/fastq/CBS1/SRR12157783_1.fastq.gz
```

Implementation notes:

- The S3/`aria2c` route follows the local reference workflow in
  `01_ref_codes/GSE171306_src`.
- `SRR12157783.sra` was downloaded and matched the expected S3 object size:
  `8,118,305,016` bytes.
- A sandboxed `fasterq-dump` run failed because this SRA needs RefSeq
  dependencies such as `NC_000072.6`, `NC_000074.6`, and `NC_005089.1`.
  Re-running `fasterq-dump` with network access resolved the dependency lookup.
- `fasterq-dump` reported `217,143,116` spots/reads read and written.
- The resulting compressed FASTQ is:
  `data/raw/gse153859/fastq/CBS1/SRR12157783_1.fastq.gz`
  (`7,096,935,823` bytes).
- `pigz -t` completed without error.
- Initial read-structure audit:
  - `GSE153859_runinfo.csv` reports `LibraryLayout=PAIRED`, but also
    `spots_with_mates=0` and `avgLength=91`.
  - `fasterq-dump` produced one FASTQ file.
  - A 1,000-read gzip sample found all reads are length 91.
  - Example headers are plain SRA headers such as
    `@SRR12157783.1 1 length=91`.
  - No obvious spatial barcode or UMI field was observed in read headers.
- This suggests the public SRA object may only expose the cDNA read needed for
  expression/spatial counting, not the barcode/UMI mate needed to reconstruct a
  molecule-level Visium BAM. The next decision point is to locate original
  R1/R2 FASTQ or BAM for CBS1/CBS2, or to treat this dataset as
  expression/spatial support rather than true APA evidence.

## 7. `GSE179572`: Human Brain Metastasis Visium Candidate

Primary source:

- GEO: `GSE179572`
- BioProject: `PRJNA744217`
- SRA study: `SRP327191`
- Tissue: surgically resected human brain metastases
- Platform: fresh frozen 10x Visium Spatial Gene Expression
- Samples: 6 Visium samples with expression H5/matrix and spatial folders

Why this dataset is higher priority than `GSE153859`:

- It is fresh frozen Visium rather than FFPE probe-based Visium.
- GEO SOFT records direct SRA relations for all 6 samples.
- SRA XML confirms public original FASTQ evidence with Visium-compatible read
  structure.
- At least one lane per sample has R1/R2 FASTQ; several lanes also expose
  I1/I2 index reads.
- The read structure includes the expected 28 bp barcode/UMI read and a 90 bp
  cDNA read, making Space Ranger/Cell Ranger BAM generation feasible.

Sample-to-SRA mapping:

| GEO Sample | Patient | Primary Tumor | SRA Experiment | Runs | Read Structure | SRA Size Range |
| --- | --- | --- | --- | --- | --- | --- |
| `GSM5420749` | patient 15 | lung adenocarcinoma | `SRX11362758` | `SRR15052390`, `SRR15052391` | 28 + 90 + 10 + 10 bp | 7.0-7.4 GB |
| `GSM5420750` | patient 16 | melanoma | `SRX11362760` | `SRR15052392`, `SRR15052393` | 28 + 90 bp | 5.5-5.6 GB |
| `GSM5420751` | patient 19 | lung adenocarcinoma | `SRX11362761` | `SRR15052394`, `SRR15052395` | 28 + 90 + 10 + 10 bp | 4.8-5.1 GB |
| `GSM5420752` | patient 24 | renal cell carcinoma | `SRX11362762` | `SRR15052396`, `SRR15052397` | 28 + 90 bp | 8.5-8.7 GB |
| `GSM5420753` | patient 26 | breast carcinoma | `SRX11362764` | `SRR15052398`, `SRR15052399` | 28 + 90 bp | 16.6-17.4 GB |
| `GSM5420754` | patient 27 | lung adenocarcinoma | `SRX11362765` | `SRR15052401`, `SRR15052402` | 28 + 90 bp | 7.6-8.0 GB |

Recommended minimal validation target:

```text
GSM5420751 / SRX11362761 / SRR15052395
```

Rationale:

- It is the smallest SRA run found in this series, about 4.8 GB.
- It has a full Visium-like 4-read structure: 28 bp, 90 bp, 10 bp, 10 bp.
- It belongs to a brain metastasis sample and has matched GEO expression and
  spatial files.

Planned minimal validation route:

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

Expected validation output:

```text
SRR15052395_1.fastq  # 28 bp barcode/UMI read
SRR15052395_2.fastq  # 90 bp cDNA read
SRR15052395_3.fastq  # 10 bp index read, if emitted by fasterq-dump
SRR15052395_4.fastq  # 10 bp index read, if emitted by fasterq-dump
```

Validation completed on 2026-06-14:

```text
Downloaded SRA:
  data/raw/gse179572/sra/SRR15052395.sra
  size on disk: 4.5 GB

fasterq-dump summary:
  spots read:    113,316,797
  reads read:    453,267,188
  reads written: 453,267,188

FASTQ output:
  SRR15052395_1.fastq  28 bp  barcode/UMI read
  SRR15052395_2.fastq  90 bp  cDNA read
  SRR15052395_3.fastq  10 bp  index read
  SRR15052395_4.fastq  10 bp  index read
```

Example first-read lengths:

```text
SRR15052395_1.fastq  28
SRR15052395_2.fastq  90
SRR15052395_3.fastq  10
SRR15052395_4.fastq  10
```

Decision:

- `GSE179572` is now the active second real-data APA calling route.
- `GSE153859` remains useful as expression/spatial support, but should not be
  counted as true APA evidence unless original R1/R2 FASTQ or BAM is found.
- Next step: prepare Space Ranger/Cell Ranger-compatible FASTQ names for
  `GSM5420751`, run alignment with barcode/UMI preservation, then run
  scAPAtrap/Sierra/polyApipe to produce `apa_matrix.csv` and site-level audit
  tables.

Space Ranger and scAPAtrap progress as of 2026-06-16:

- Space Ranger completed successfully for `GSM5420751 / SRR15052395`.
- Barcode/UMI-preserving BAM:

```text
pipeline_output/gse179572_GSM5420751_sr/outs/possorted_genome_bam.bam
```

- BAM tag audit confirmed 10x-style barcode and UMI tags:
  - `CB` / `CR`: corrected/raw spot barcode;
  - `UB` / `UR`: corrected/raw UMI;
  - `GX` / `GN`: gene annotation tags.
- The first scAPAtrap attempt entered `findUniqueMap` and reached
  `dedupByPos`, but it was launched in a foreground validation session and was
  interrupted before `scAPAtrapData.rda` was produced.
- scAPAtrap was relaunched with a staging BAM path so all large intermediate
  files are isolated from the Space Ranger output directory. The staging
  directory and BAM stem intentionally avoid strings such as `_bam`, because
  scAPAtrap builds intermediate filenames with `gsub(".bam", "", input)`:

```text
pipeline_output/gse179572_GSM5420751_scapatrap/stage/spaceranger_input.bam
```

- The first staging attempt used `work_bam/possorted_genome_bam.bam` and failed
  at `separateBamBystrand`, because scAPAtrap's regex-style `gsub(".bam", "")`
  rewrote the output path to a non-existent `work/` directory. The current
  `stage/spaceranger_input.bam` path fixes this.

- Launcher scripts:

```text
scripts/run_scapatrap_spaceranger.py
scripts/run_scapatrap_spaceranger.sh
```

- Current target outputs after successful scAPAtrap completion:

```text
data/processed/gse179572_gsm5420751_scapatrap/
  apa_matrix.csv
  apa_sites.csv.gz
  apa_site_counts.csv.gz
  coordinates.csv
  metadata.csv
  qc_summary.json
```

- Monitor commands:

```bash
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
tail -f pipeline_output/gse179572_GSM5420751_scapatrap/logs/run_scapatrap_spaceranger.nohup.log
tail -f pipeline_output/gse179572_GSM5420751_scapatrap/logs/scapatrap_internal.log
pgrep -af 'run_scapatrap_spaceranger|run_scapatrap.R|scAPAtrap|Rscript'
```

## 8. Planned APA Calling Route

Preferred route:

1. Search for original CBS1/CBS2 R1/R2 FASTQ or BAM outside the SRA-lite/public
   cDNA object.
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

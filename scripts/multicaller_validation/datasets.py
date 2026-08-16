#!/usr/bin/env python
"""Dataset registry for the multi-caller (Sierra vs scAPAtrap) robustness
validation. Single source of truth for paths shared by all pipeline steps.

Three datasets: human kidney (original), human AD-brain PFC, mouse colon —
2 species x 3 tissues. Each entry:
  key         output dir name under pipeline_output/multicaller_validation/
  label       human-readable (GEO accession, tissue, platform, n_spots)
  species / tissue
  bam         Space Ranger possorted_genome_bam.bam (CB/UB tagged)
  gtf_gz      reference GTF (the exact Space Ranger reference GTF the
              scAPAtrap baseline was annotated with; chr-prefixed contig
              names matching the BAM header)
  baseline    data/processed/<scapatrap baseline> (usage matrix +
              coordinates + gene-annotated apa_sites for caller-2 gene ids)
  outdir      working/output dir (junctions, sierra peaks/counts, converted
              spaGAPA tables, conformal runs, metric JSONs)

Chromosome-naming note (GSE169749): the Space Ranger BAM uses UCSC-style
chr-prefixed contigs (chr1..chrY, chrM; un-prefixed scaffolds). The Ensembl
GRCm39.111 GTF uses un-prefixed contigs (1..19, X, Y, MT) and would yield
empty peak annotation. We therefore use the Space Ranger reference's own
GTF (refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz, GENCODE M33 = Ensembl
110, chr-prefixed) — the identical file the scAPAtrap baseline was
annotated with (scripts/run_scapatrap_gse169749_gsm5213483.sh), so no
renaming is needed and gene ids match by construction.
"""

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"

DATASETS = {
    "gse183456_gsm6047774": {
        "label": "GSE183456 / GSM6047774",
        "species": "human",
        "tissue": "kidney",
        "detail": "human kidney Visium, 3010 spots",
        "bam": f"{ROOT}/pipeline_output/gse183456_GSM6047774_sr/outs/possorted_genome_bam.bam",
        "gtf_gz": "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz",
        "baseline": f"{ROOT}/data/processed/gse183456_gsm6047774_scapatrap",
        "outdir": f"{ROOT}/pipeline_output/multicaller_validation/sierra",
    },
    "gse220442_gsm6801751": {
        "label": "GSE220442 / GSM6801751",
        "species": "human",
        "tissue": "brain (AD PFC, control)",
        "detail": "human AD-brain prefrontal cortex Visium (control), 4179 spots",
        "bam": f"{ROOT}/pipeline_output/gse220442_GSM6801751_sr/outs/possorted_genome_bam.bam",
        "gtf_gz": "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz",
        "baseline": f"{ROOT}/data/processed/gse220442_gsm6801751_scapatrap",
        "outdir": f"{ROOT}/pipeline_output/multicaller_validation/gse220442_gsm6801751",
    },
    "gse169749_gsm5213483": {
        "label": "GSE169749 / GSM5213483",
        "species": "mouse",
        "tissue": "colon (DSS d0)",
        "detail": "mouse colon Visium (DSS injury, day 0), 2715 spots",
        "bam": f"{ROOT}/pipeline_output/gse169749/gsm5213483_d0_sr/outs/possorted_genome_bam.bam",
        "gtf_gz": "/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz",
        "baseline": f"{ROOT}/data/processed/gse169749_gsm5213483_scapatrap",
        "outdir": f"{ROOT}/pipeline_output/multicaller_validation/gse169749_gsm5213483",
    },
}


def get(key):
    if key not in DATASETS:
        raise SystemExit(f"unknown dataset '{key}'; choose from {sorted(DATASETS)}")
    return DATASETS[key]

#!/usr/bin/env bash
# GSE169749 (mouse colon, DSS injury) Visium -> spaceranger count -> scAPAtrap
# Run on S91. Both samples (day0 = SRR14083626, day14 = SRR14083627).
#
# FASTQ conversion: fasterq-dump with 10x-aware splitting (R1=barcode+UMI, R2=CDNA).
# GSE169749 is standard 10x Visium: R1 28bp (bc1+bc2+umi = 16+12), R2 90bp (cDNA).
set -euo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/home/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/home/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s3/mengzijun/tmp"
        PYTHON="/home/mengzijun/anaconda3/envs/spagapa/bin/python"
        SRATK="/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin"
        SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
        ;;
    *)
        echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac

export OPENBLAS_NUM_THREADS=8
PROJECT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
RAW="$PROJECT/data/raw/gse169749"
OUT="$PROJECT/pipeline_output/gse169749"
LOGS="$PROJECT/logs"
MMREF="/s1/SHARE/00_ref_genecode/refdata-gex-mm10-2020-A"

mkdir -p "$OUT" "$LOGS" "$TMPDIR"

# --- fasterq-dump conversion per sample (10x layout: --split-files) ---
declare -A SAMPLES=(
    [SRR14083626]="gsm5213483_d0"
    [SRR14083627]="gsm5213484_d14"
)

for srr in "${!SAMPLES[@]}"; do
    label="${SAMPLES[$srr]}"
    fastq_dir="$RAW/fastq_${label}"
    if [ -f "$fastq_dir/${srr}_2.fastq.gz" ]; then
        echo "[skip] FASTQ already exists for $srr ($label)"
        continue
    fi
    echo "[fasterq] Converting $srr ($label) -> $fastq_dir"
    mkdir -p "$fastq_dir"
    cd "$fastq_dir"
    # 10x Visium: fasterq-dump --split-files produces _1 (R1) and _2 (R2).
    # --include-technical keeps R1 (barcode). -e 8 threads, -p progress.
    "$SRATK/fasterq-dump" "$RAW/${srr}.sra" \
        --split-files --include-technical \
        --threads 8 --progress \
        --temp "$TMPDIR" \
        --outdir "$fastq_dir" 2>&1 | tee "$LOGS/fasterq_${srr}.log"
    # gzip the outputs
    pigz -p 8 "$fastq_dir"/${srr}_*.fastq 2>/dev/null || gzip "$fastq_dir"/${srr}_*.fastq
    cd "$PROJECT"
    echo "[fasterq] done $srr -> $(ls -la "$fastq_dir")"
done

echo "=== FASTQ conversion complete ==="

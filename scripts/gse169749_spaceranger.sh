#!/usr/bin/env bash
# Run spaceranger count for both GSE169749 mouse colon Visium samples.
# Assumes fasterq-dump output has been renamed into bcl2fastq convention.
set -euo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91) export TMPDIR="/s3/mengzijun/tmp" ;;
    *) echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
export OPENBLAS_NUM_THREADS=8

PROJECT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
RAW="$PROJECT/data/raw/gse169749"
OUT="$PROJECT/pipeline_output/gse169749"
LOGS="$PROJECT/logs"
MMREF="/s1/SHARE/00_ref_genecode/refdata-gex-mm10-2020-A"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"

mkdir -p "$OUT" "$LOGS"

# Per-sample FASTQ dirs (one sample per dir, with renamed _S1_L001_ files).
declare -A FASTQDIRS=(
    [gsm5213483_d0]="$RAW/fastq_gsm5213483_d0_renamed"
    [gsm5213484_d14]="$RAW/fastq_gsm5213484_d14_renamed"
)

for sample in "${!FASTQDIRS[@]}"; do
    fqdir="${FASTQDIRS[$sample]}"
    if [ -d "$OUT/${sample}_sr/outs" ]; then
        echo "[skip] spaceranger output exists: $OUT/${sample}_sr/outs"
        continue
    fi
    echo "[spaceranger] $sample <- $fqdir"
    # Build --sample (SampleName) and --lanes auto. spaceranger count expects:
    #   fastq_dir/SampleName_S1_L001_R1_001.fastq.gz (+R2)
    "$SPACERANGER" count \
        --id="${sample}_sr" \
        --transcriptome="$MMREF" \
        --fastqs="$fqdir" \
        --sample="$sample" \
        --create-bam=true \
        --localcores=8 \
        --localmem=32 \
        --maxjobs=4 \
        --jobmode=local
    # move output if it landed in cwd
    if [ -d "${sample}_sr" ] && [ ! -d "$OUT/${sample}_sr" ]; then
        mv "${sample}_sr" "$OUT/"
    fi
done

echo "=== spaceranger complete ==="
ls -d "$OUT"/*/outs 2>/dev/null

#!/usr/bin/env bash
# Run mouse-adapted scAPAtrap on both GSE169749 spaceranger outputs.
set -euo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91)
        export TMPDIR="/s3/mengzijun/tmp"
        PYTHON="/home/mengzijun/anaconda3/envs/spagapa/bin/python"
        ;;
    *) echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
export OPENBLAS_NUM_THREADS=8

PROJECT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SR="$PROJECT/pipeline_output/gse169749"
LAUNCHER="$PROJECT/scripts/run_scapatrap_spaceranger_mouse.py"
MMGTF="/s1/SHARE/00_ref_genecode/refdata-gex-mm10-2020-A/genes/genes.gtf"

# Process each sample
declare -A SAMPLES=(
    [gsm5213483_d0]="GSM5213483 Colon DSS day0"
    [gsm5213484_d14]="GSM5213484 Colon DSS day14"
)

for sample in "${!SAMPLES[@]}"; do
    desc="${SAMPLES[$sample]}"
    bam="$SR/${sample}_sr/outs/possorted_genome_bam.bam"
    spatial="$SR/${sample}_sr/outs/spatial"
    outroot="$SR/${sample}_scapatrap"
    processed="$PROJECT/data/processed/gse169749_${sample}_scapatrap"
    dsname="gse169749_${sample}_scapatrap"

    if [ ! -f "$bam" ]; then
        echo "[ERROR] BAM not found for $sample: $bam" >&2
        continue
    fi

    echo "[scapatrap] $sample -> $outroot"
    $PYTHON "$LAUNCHER" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --gtf "$MMGTF" \
        --output-root "$outroot" \
        --processed-dir "$processed" \
        --dataset-name "$dsname" \
        --source-label "GSE169749 $desc (mouse Visium) Space Ranger BAM + scAPAtrap" \
        --species mouse \
        --tissue "mouse colon" \
        --threads 12 \
        --readlength 90 \
        --cov-cutoff 10 \
        --min-cells 10 \
        --min-count 10 \
        --min-parent-count 5 \
        --tails-search peaks
done

echo "=== scAPAtrap complete ==="

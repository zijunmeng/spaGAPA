#!/usr/bin/env bash
# Rename fasterq-dump outputs to spaceranger's expected bcl2fastq naming and gzip.
#   SRR_<12>_1.fastq  ->  <sample>_S1_L001_R1_001.fastq.gz
#   SRR_<12>_2.fastq  ->  <sample>_S1_L001_R2_001.fastq.gz
set -euo pipefail

RAW="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/data/raw/gse169749"

rename_sample() {
    local srr="$1" sample="$2"
    local srcdir="$RAW/fastq_${sample}"
    local dstdir="$RAW/fastq_${sample}_renamed"
    mkdir -p "$dstdir"
    if [ -f "$dstdir/${sample}_S1_L001_R2_001.fastq.gz" ]; then
        echo "[skip] $sample already renamed+gzipped"
        return
    fi
    echo "[rename] $srr -> $sample"
    for r in 1 2; do
        src="$srcdir/${srr}_${r}.fastq"
        dst="$dstdir/${sample}_S1_L001_R${r}_001.fastq"
        if [ ! -f "$src" ]; then echo "[warn] missing $src"; continue; fi
        if [ -f "$dst.gz" ]; then continue; fi
        echo "  gzipping $src -> $dst.gz"
        pigz -p 6 -c "$src" > "$dst.gz"
    done
    echo "[done] $sample -> $(ls -la "$dstdir")"
}

rename_sample SRR14083626 gsm5213483_d0
rename_sample SRR14083627 gsm5213484_d14

echo "=== rename+gzip complete ==="
ls -la "$RAW"/*_renamed/ 2>/dev/null

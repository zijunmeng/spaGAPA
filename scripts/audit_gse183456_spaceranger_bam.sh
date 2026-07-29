#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
BAM="${BAM:-$ROOT/pipeline_output/gse183456_GSM6047774_sr/outs/possorted_genome_bam.bam}"
LOG_DIR="$ROOT/pipeline_output/gse183456_GSM6047774_scapatrap/logs"
LOG="$LOG_DIR/bam_tag_audit.log"

mkdir -p "$LOG_DIR"

if [[ ! -s "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi

set +u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate samtools
set -u

{
  echo "[$(date '+%F %T')] BAM tag audit"
  echo "BAM: $BAM"
  samtools quickcheck -v "$BAM" || true
  echo
  echo "Header preview:"
  set +o pipefail
  samtools view -H "$BAM" | head -n 30 || true
  echo
  echo "First alignments with key tags:"
  samtools view "$BAM" | head -n 200 | awk '
    BEGIN { cb=0; ub=0; gx=0; gn=0; n=0 }
    {
      n++;
      has_cb=0; has_ub=0; has_gx=0; has_gn=0;
      for (i=12; i<=NF; i++) {
        if ($i ~ /^CB:Z:/) has_cb=1;
        if ($i ~ /^UB:Z:/) has_ub=1;
        if ($i ~ /^GX:Z:/) has_gx=1;
        if ($i ~ /^GN:Z:/) has_gn=1;
      }
      cb+=has_cb; ub+=has_ub; gx+=has_gx; gn+=has_gn;
    }
    END {
      printf("records_checked=%d\nCB=%d\nUB=%d\nGX=%d\nGN=%d\n", n, cb, ub, gx, gn)
    }' || true
  set -o pipefail
} > "$LOG" 2>&1

cat "$LOG"

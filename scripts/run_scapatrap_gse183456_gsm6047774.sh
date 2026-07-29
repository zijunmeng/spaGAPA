#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
BAM="${BAM:-$ROOT/pipeline_output/gse183456_GSM6047774_sr/outs/possorted_genome_bam.bam}"
SPATIAL_DIR="${SPATIAL_DIR:-$ROOT/pipeline_output/gse183456_GSM6047774_sr/outs/spatial}"
OUT_ROOT="$ROOT/pipeline_output/gse183456_GSM6047774_scapatrap"
PROCESSED_DIR="$ROOT/data/processed/gse183456_gsm6047774_scapatrap"
LOG_DIR="$OUT_ROOT/logs"
OUT_LOG="$LOG_DIR/run_scapatrap_gse183456.nohup.log"
PID_FILE="$OUT_ROOT/scapatrap_launcher.pid"

mkdir -p "$LOG_DIR"

if [[ ! -s "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi
if [[ ! -d "$SPATIAL_DIR" ]]; then
  echo "Spatial directory not found: $SPATIAL_DIR" >&2
  exit 1
fi

set +u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate spagapa
set -u

setsid nohup python "$ROOT/scripts/run_scapatrap_spaceranger.py" \
  --bam "$BAM" \
  --spatial-dir "$SPATIAL_DIR" \
  --output-root "$OUT_ROOT" \
  --processed-dir "$PROCESSED_DIR" \
  --dataset-name "gse183456_gsm6047774_scapatrap" \
  --source-label "GSE183456 GSM6047774 Space Ranger BAM + scAPAtrap" \
  --species "human" \
  --tissue "kidney" \
  --threads "${THREADS:-12}" \
  --force \
  >> "$OUT_LOG" 2>&1 < /dev/null &

echo $! > "$PID_FILE"
echo "scAPAtrap started: $(cat "$PID_FILE")"
echo "Monitor with:"
echo "  tail -f $OUT_LOG"
echo "  tail -f $OUT_ROOT/logs/scapatrap_internal.log"

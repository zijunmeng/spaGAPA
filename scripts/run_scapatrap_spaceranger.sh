#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
PID_FILE="$ROOT/pipeline_output/gse179572_GSM5420751_scapatrap/scapatrap_launcher.pid"
LOG_DIR="$ROOT/pipeline_output/gse179572_GSM5420751_scapatrap/logs"
OUT_LOG="$LOG_DIR/run_scapatrap_spaceranger.nohup.log"

mkdir -p "$LOG_DIR"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate spagapa
python "$ROOT/scripts/run_scapatrap_spaceranger.py" >> "$OUT_LOG" 2>&1

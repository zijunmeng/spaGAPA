#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/s1/SHARE/mengzijun/01_project/26_spaGAPA}"
CONDA_ENV="${CONDA_ENV:-samtools}"
MODE="${1:-all}"
MAX_FASTQ_DOWNLOADS="${MAX_FASTQ_DOWNLOADS:-4}"
MAX_SUPP_DOWNLOADS="${MAX_SUPP_DOWNLOADS:-2}"

ACCS=(GSE220442 GSE237183)

log() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

generate_manifest() {
  local acc="$1"
  local out_dir="$ROOT/data/raw/${acc,,}"
  mkdir -p "$out_dir/logs" "$out_dir/geo_supplementary" "$out_dir/fastq_ena"

  log "Generating GEO/SRA/ENA manifest for $acc"
  conda run -n "$CONDA_ENV" python \
    "$ROOT/spaGAPA/scripts/build_visium_candidate_manifest.py" \
    "$acc" \
    "$out_dir"
}

download_supplementary() {
  local acc="$1"
  local out_dir="$ROOT/data/raw/${acc,,}"
  local log_dir="$out_dir/logs"
  local url_file="$log_dir/geo_supplementary_urls.txt"
  if [[ ! -s "$url_file" ]]; then
    log "No GEO supplementary URLs for $acc"
    return 0
  fi
  log "Downloading GEO supplementary files for $acc"
  aria2c \
    --input-file="$url_file" \
    --dir="$out_dir/geo_supplementary" \
    --check-certificate=false \
    --max-concurrent-downloads="$MAX_SUPP_DOWNLOADS" \
    --max-connection-per-server=2 \
    --split=2 \
    --min-split-size=20M \
    --continue=true \
    --file-allocation=none \
    --auto-file-renaming=false \
    --allow-overwrite=false \
    --max-tries=0 \
    --retry-wait=30 \
    --console-log-level=warn \
    --download-result=full \
    --summary-interval=120 \
    --log="$log_dir/geo_supplementary_aria2c.log"
}

download_fastq() {
  local acc="$1"
  local out_dir="$ROOT/data/raw/${acc,,}"
  local log_dir="$out_dir/logs"
  local url_file="$log_dir/ena_fastq_urls.txt"
  if [[ ! -s "$url_file" ]]; then
    log "No ENA FASTQ URLs for $acc"
    return 0
  fi
  log "Downloading ENA FASTQ files for $acc"
  aria2c \
    --input-file="$url_file" \
    --dir="$out_dir/fastq_ena" \
    --check-certificate=false \
    --max-concurrent-downloads="$MAX_FASTQ_DOWNLOADS" \
    --max-connection-per-server=2 \
    --split=2 \
    --min-split-size=20M \
    --continue=true \
    --file-allocation=none \
    --auto-file-renaming=false \
    --allow-overwrite=false \
    --max-tries=0 \
    --retry-wait=30 \
    --console-log-level=warn \
    --download-result=full \
    --summary-interval=120 \
    --log="$log_dir/ena_fastq_aria2c.log"
}

need_cmd aria2c
need_cmd conda

case "$MODE" in
  all|metadata|supp|fastq) ;;
  *)
    echo "Usage: $0 [all|metadata|supp|fastq]" >&2
    exit 2
    ;;
esac

for acc in "${ACCS[@]}"; do
  if [[ "$MODE" == "all" || "$MODE" == "metadata" ]]; then
    generate_manifest "$acc"
  fi
done

for acc in "${ACCS[@]}"; do
  if [[ "$MODE" == "all" || "$MODE" == "supp" ]]; then
    download_supplementary "$acc"
  fi
done

for acc in "${ACCS[@]}"; do
  if [[ "$MODE" == "all" || "$MODE" == "fastq" ]]; then
    download_fastq "$acc"
  fi
done

log "Finished mode=$MODE for ${ACCS[*]}"

#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
OUT_DIR="$ROOT/data/raw/gse183456/GSM6047774"
FASTQ_RAW_DIR="$OUT_DIR/fastq_ena"
FASTQ_SR_DIR="$OUT_DIR/fastq_spaceranger_ena"

mkdir -p "$FASTQ_SR_DIR"

for fq in \
  "$FASTQ_RAW_DIR/SRR18794263_1.fastq.gz" \
  "$FASTQ_RAW_DIR/SRR18794263_2.fastq.gz" \
  "$FASTQ_RAW_DIR/SRR18794264_1.fastq.gz" \
  "$FASTQ_RAW_DIR/SRR18794264_2.fastq.gz"
do
  if [[ ! -s "$fq" ]]; then
    echo "Missing or empty FASTQ: $fq" >&2
    exit 1
  fi
  gzip -t "$fq"
done

rm -f "$FASTQ_SR_DIR"/GSM6047774_S1_L*_R*_001.fastq.gz
ln -sf "$FASTQ_RAW_DIR/SRR18794263_1.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_L001_R1_001.fastq.gz"
ln -sf "$FASTQ_RAW_DIR/SRR18794263_2.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_L001_R2_001.fastq.gz"
ln -sf "$FASTQ_RAW_DIR/SRR18794264_1.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_L002_R1_001.fastq.gz"
ln -sf "$FASTQ_RAW_DIR/SRR18794264_2.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_L002_R2_001.fastq.gz"

echo "Space Ranger ENA FASTQ links ready: $FASTQ_SR_DIR"

#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
OUT_DIR="$ROOT/data/raw/gse183456/GSM6047774"
FASTQ_RAW_DIR="$OUT_DIR/fastq_ena"
FASTQ_SR_DIR="$OUT_DIR/fastq_spaceranger_ena"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$FASTQ_RAW_DIR" "$FASTQ_SR_DIR" "$LOG_DIR"

cat > "$LOG_DIR/ena_fastq_urls.txt" <<'URLS'
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR187/063/SRR18794263/SRR18794263_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR187/063/SRR18794263/SRR18794263_2.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR187/064/SRR18794264/SRR18794264_1.fastq.gz
https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR187/064/SRR18794264/SRR18794264_2.fastq.gz
URLS

aria2c \
  --input-file="$LOG_DIR/ena_fastq_urls.txt" \
  --dir="$FASTQ_RAW_DIR" \
  --check-certificate=false \
  --max-concurrent-downloads=4 \
  --max-connection-per-server=2 \
  --split=2 \
  --min-split-size=20M \
  --continue=true \
  --file-allocation=none \
  --auto-file-renaming=false \
  --allow-overwrite=false \
  --max-tries=5 \
  --retry-wait=30 \
  --console-log-level=warn \
  --download-result=full \
  --summary-interval=60 \
  --log="$LOG_DIR/ena_aria2c.log"

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

echo "ENA FASTQs ready: $FASTQ_SR_DIR"

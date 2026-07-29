#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
OUT_DIR="$ROOT/data/raw/gse183456/GSM6047774"
PREFETCH_DIR="$OUT_DIR/prefetch_sra"
FASTQ_RAW_DIR="$OUT_DIR/fastq_raw"
FASTQ_SR_DIR="$OUT_DIR/fastq_spaceranger"
LOG_DIR="$OUT_DIR/logs"
TMP_DIR="$OUT_DIR/tmp_fasterq_prefetch"
THREADS="${THREADS:-8}"

mkdir -p "$PREFETCH_DIR" "$FASTQ_RAW_DIR" "$FASTQ_SR_DIR" "$LOG_DIR" "$TMP_DIR"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate samtools

for run in SRR18794263 SRR18794264; do
  echo "[$(date '+%F %T')] prefetch $run"
  prefetch \
    --transport http \
    --resume yes \
    --max-size 100G \
    --output-directory "$PREFETCH_DIR" \
    "$run"
done

for run in SRR18794263 SRR18794264; do
  sra_path="$PREFETCH_DIR/$run/$run.sra"
  if [[ ! -f "$sra_path" ]]; then
    echo "Missing $sra_path after prefetch" >&2
    exit 1
  fi
  echo "[$(date '+%F %T')] vdb-validate $run"
  vdb-validate "$sra_path"

  if compgen -G "$FASTQ_RAW_DIR/${run}_*.fastq.gz" > /dev/null; then
    echo "[$(date '+%F %T')] FASTQ already exists for $run; skipping fasterq-dump"
    continue
  fi
  echo "[$(date '+%F %T')] fasterq-dump $run"
  fasterq-dump \
    --split-files \
    --include-technical \
    --threads "$THREADS" \
    --temp "$TMP_DIR" \
    --outdir "$FASTQ_RAW_DIR" \
    "$sra_path"
  echo "[$(date '+%F %T')] pigz $run"
  pigz -p "$THREADS" "$FASTQ_RAW_DIR/${run}"_*.fastq
done

echo "[$(date '+%F %T')] Creating Space Ranger FASTQ symlinks"
rm -f "$FASTQ_SR_DIR"/GSM6047774_S1_L*_R*_001.fastq.gz
rm -f "$FASTQ_SR_DIR"/GSM6047774_S1_L*_I*_001.fastq.gz

lane=1
for run in SRR18794263 SRR18794264; do
  lane_label=$(printf "L%03d" "$lane")
  ln -sf "$FASTQ_RAW_DIR/${run}_1.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_R1_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_2.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_R2_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_3.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_I1_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_4.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_I2_001.fastq.gz"
  lane=$((lane + 1))
done

rm -rf "$TMP_DIR"
echo "[$(date '+%F %T')] prefetch FASTQ conversion complete"
echo "Space Ranger FASTQs: $FASTQ_SR_DIR"

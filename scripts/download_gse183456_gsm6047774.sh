#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
OUT_DIR="$ROOT/data/raw/gse183456/GSM6047774"
SRA_DIR="$OUT_DIR/sra"
FASTQ_RAW_DIR="$OUT_DIR/fastq_raw"
FASTQ_SR_DIR="$OUT_DIR/fastq_spaceranger"
LOG_DIR="$OUT_DIR/logs"
THREADS="${THREADS:-8}"

mkdir -p "$SRA_DIR" "$FASTQ_RAW_DIR" "$FASTQ_SR_DIR" "$LOG_DIR"

cat > "$LOG_DIR/srr_urls.txt" <<'URLS'
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18794263/SRR18794263
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18794264/SRR18794264
URLS

cat > "$LOG_DIR/download_and_fastq.worker.sh" <<'WORKER'
#!/usr/bin/env bash
set -euo pipefail

source ~/anaconda3/etc/profile.d/conda.sh
conda activate samtools

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
OUT_DIR="$ROOT/data/raw/gse183456/GSM6047774"
SRA_DIR="$OUT_DIR/sra"
FASTQ_RAW_DIR="$OUT_DIR/fastq_raw"
FASTQ_SR_DIR="$OUT_DIR/fastq_spaceranger"
LOG_DIR="$OUT_DIR/logs"
THREADS="${THREADS:-8}"

echo "[$(date '+%F %T')] Downloading SRR18794263/SRR18794264"
aria2c \
  --input-file="$LOG_DIR/srr_urls.txt" \
  --dir="$SRA_DIR" \
  --check-certificate=false \
  --max-concurrent-downloads=2 \
  --max-connection-per-server=2 \
  --split=2 \
  --min-split-size=20M \
  --continue=true \
  --auto-file-renaming=false \
  --allow-overwrite=true \
  --max-tries=5 \
  --retry-wait=30 \
  --console-log-level=warn \
  --download-result=full \
  --summary-interval=60 \
  --log="$LOG_DIR/aria2c.log"

for run in SRR18794263 SRR18794264; do
  if [[ -f "$SRA_DIR/$run" && ! -f "$SRA_DIR/$run.sra" ]]; then
    mv "$SRA_DIR/$run" "$SRA_DIR/$run.sra"
  fi
  if [[ ! -f "$SRA_DIR/$run.sra" ]]; then
    echo "Missing $SRA_DIR/$run.sra" >&2
    exit 1
  fi
done

echo "[$(date '+%F %T')] Converting SRA to FASTQ"
TMP_DIR="$OUT_DIR/tmp_fasterq"
mkdir -p "$TMP_DIR"
for run in SRR18794263 SRR18794264; do
  if compgen -G "$FASTQ_RAW_DIR/${run}_*.fastq.gz" > /dev/null; then
    echo "[$(date '+%F %T')] FASTQ already exists for $run; skipping fasterq-dump"
    continue
  fi
  fasterq-dump \
    --split-files \
    --include-technical \
    --threads "$THREADS" \
    --temp "$TMP_DIR" \
    --outdir "$FASTQ_RAW_DIR" \
    "$SRA_DIR/$run.sra"
  pigz -p "$THREADS" "$FASTQ_RAW_DIR/${run}"_*.fastq
done
rm -rf "$TMP_DIR"

echo "[$(date '+%F %T')] Creating Space Ranger FASTQ names"
rm -f "$FASTQ_SR_DIR"/GSM6047774_S*_L001_R*.fastq.gz
lane=1
for run in SRR18794263 SRR18794264; do
  lane_label=$(printf "L%03d" "$lane")
  ln -sf "$FASTQ_RAW_DIR/${run}_1.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_R1_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_2.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_R2_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_3.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_I1_001.fastq.gz"
  ln -sf "$FASTQ_RAW_DIR/${run}_4.fastq.gz" "$FASTQ_SR_DIR/GSM6047774_S1_${lane_label}_I2_001.fastq.gz"
  lane=$((lane + 1))
done

echo "[$(date '+%F %T')] Done"
echo "Space Ranger FASTQs: $FASTQ_SR_DIR"
WORKER

chmod +x "$LOG_DIR/download_and_fastq.worker.sh"
nohup "$LOG_DIR/download_and_fastq.worker.sh" > "$LOG_DIR/download_and_fastq.nohup.log" 2>&1 &
echo $! > "$LOG_DIR/download_and_fastq.pid"

echo "Download/conversion started: $(cat "$LOG_DIR/download_and_fastq.pid")"
echo "Monitor with:"
echo "  tail -f $LOG_DIR/download_and_fastq.nohup.log"
echo "  tail -f $LOG_DIR/aria2c.log"

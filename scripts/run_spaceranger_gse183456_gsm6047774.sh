#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
PACKAGE_ROOT="$ROOT/spaGAPA"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A"
FASTQ_DIR="${FASTQ_DIR:-$ROOT/data/raw/gse183456/GSM6047774/fastq_spaceranger_ena}"
IMAGE="$ROOT/data/raw/gse183456/GSM6047774/extracted/spatial/detected_tissue_image.jpg"
OUT_ROOT="$PACKAGE_ROOT/pipeline_output"
ID="gse183456_GSM6047774_sr"
LOG_DIR="$OUT_ROOT/${ID}_logs"
DRY_ROOT="$OUT_ROOT/${ID}_dryrun"
THREADS="${THREADS:-24}"
MEM_GB="${MEM_GB:-128}"

mkdir -p "$LOG_DIR"
mkdir -p "$OUT_ROOT"

if [[ ! -x "$SPACERANGER" ]]; then
  echo "Space Ranger not found or not executable: $SPACERANGER" >&2
  exit 1
fi
if [[ ! -d "$TRANSCRIPTOME" ]]; then
  echo "Transcriptome reference not found: $TRANSCRIPTOME" >&2
  exit 1
fi
if [[ ! -f "$IMAGE" ]]; then
  echo "Image not found: $IMAGE" >&2
  exit 1
fi
if ! compgen -G "$FASTQ_DIR/GSM6047774_S1_L001_R1_001.fastq.gz" > /dev/null; then
  echo "Space Ranger FASTQs not ready in $FASTQ_DIR" >&2
  echo "Run/monitor scripts/download_ena_gse183456_gsm6047774_fastq.sh first." >&2
  exit 1
fi

COMMON_ARGS=(
  count
  --id "$ID"
  --description "GSE183456 GSM6047774 IU-F59 human kidney Visium"
  --transcriptome "$TRANSCRIPTOME"
  --fastqs "$FASTQ_DIR"
  --sample "GSM6047774"
  --image "$IMAGE"
  --unknown-slide "visium-1"
  --create-bam true
  --localcores "$THREADS"
  --localmem "$MEM_GB"
  --disable-cell-annotation
  --disable-ui
)

echo "[$(date '+%F %T')] Running Space Ranger dry-run"
rm -rf "$DRY_ROOT"
mkdir -p "$DRY_ROOT"
(
  cd "$DRY_ROOT"
  "$SPACERANGER" "${COMMON_ARGS[@]}" --dry > "$LOG_DIR/spaceranger_dry_run.log" 2>&1
)

echo "[$(date '+%F %T')] Dry-run passed; launching Space Ranger"
if [[ -e "$OUT_ROOT/$ID" ]]; then
  echo "Output pipestance already exists: $OUT_ROOT/$ID" >&2
  echo "Move/remove it before launching a new Space Ranger count." >&2
  exit 1
fi
(
  cd "$OUT_ROOT"
  setsid nohup "$SPACERANGER" "${COMMON_ARGS[@]}" > "$LOG_DIR/spaceranger_count.nohup.log" 2>&1 < /dev/null &
  echo $! > "$LOG_DIR/spaceranger_count.pid"
)

echo "Space Ranger started: $(cat "$LOG_DIR/spaceranger_count.pid")"
echo "Monitor with:"
echo "  tail -f $LOG_DIR/spaceranger_count.nohup.log"
echo "  tail -f $OUT_ROOT/$ID/_log"

#!/usr/bin/env bash
# =============================================================================
# Run spaceranger count for GSE169749 GSM5213483 (mouse colon, day 0 Visium).
#
# Inputs (already on disk):
#   - FASTQ: data/raw/gse169749/fastq_gsm5213483_d0_renamed/
#       gsm5213483_d0_S1_L001_R1_001.fastq.gz   (28bp barcode+UMI)
#       gsm5213483_d0_S1_L001_R2_001.fastq.gz   (120bp cDNA, 3' GE)
#   - HE image: data/raw/gse169749/images/GSM5213483_d0_he_8bit.tif
#
# Slide serial number V19S23-097 / capture area A1 (from GSM5213483 GEO suppl
# filename GSM5213483_V19S23-097_A1_S1_tissue_positions_list.csv.gz). The 10x
# slide serial is no longer in the spaceranger 4.1 bundled registry, so we use
# --unknown-slide visium-1 (the standard 11mm Visium layout). This is identical
# to the layout V19S23-097 used and produces a valid tissue_positions.csv.
#
# Output: pipeline_output/gse169749/gsm5213483_d0_sr/outs/
#   - possorted_genome_bam.bam (CB/UB tagged)  [scAPAtrap input]
#   - spatial/tissue_positions.csv             [scAPAtrap input]
#
# Hostname-aware (CLAUDE.md S90/S91/S97/S98).
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md) ---
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) export TMPDIR="/s2/mengzijun/tmp" ;;
    S91) export TMPDIR="/s3/mengzijun/tmp" ;;
    S97) export TMPDIR="/s972/mengzijun/tmp" ;;
    S98) export TMPDIR="/s982/mengzijun/tmp" ;;
    *)  echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac
export TMPDIR
export OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

# --- 2. Paths ---
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
FASTQ_DIR="$ROOT/data/raw/gse169749/fastq_gsm5213483_d0_renamed"
SAMPLE_NAME="gsm5213483_d0"
HE_IMAGE="$ROOT/data/raw/gse169749/images/GSM5213483_d0_he_8bit.tif"
OUT_ROOT="$ROOT/pipeline_output/gse169749"
OUT_ID="gsm5213483_d0_sr"

# GRCm39 (GENCODE M33 / 10x 2024-A) - newer reference, matches the GTF used by scAPAtrap.
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"

LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SAMPLE_LOG="$LOG_ROOT/${DATE_TAG}_gsm5213483_spaceranger_count.log"
mkdir -p "$LOG_ROOT" "$OUT_ROOT"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$SAMPLE_LOG"; }

# --- 3. Preflight ---
for f in "$FASTQ_DIR/${SAMPLE_NAME}_S1_L001_R1_001.fastq.gz" \
         "$FASTQ_DIR/${SAMPLE_NAME}_S1_L001_R2_001.fastq.gz"; do
    if [[ ! -s "$f" ]]; then
        log "ERROR missing FASTQ: $f"
        exit 2
    fi
done
if [[ ! -s "$HE_IMAGE" ]]; then
    log "ERROR missing HE image: $HE_IMAGE"
    exit 2
fi
if [[ ! -x "$SPACERANGER" ]]; then
    log "ERROR spaceranger not executable: $SPACERANGER"
    exit 2
fi

# --- 4. Skip if already complete ---
if [[ -f "$OUT_ROOT/$OUT_ID/outs/web_summary.html" || \
      -f "$OUT_ROOT/$OUT_ID/outs/count_summary.html" ]]; then
    log "SKIP spaceranger: output already complete at $OUT_ROOT/$OUT_ID/outs/"
    exit 0
fi

# --- 5. Clean any partial prior run (the previous attempt was terminated mid-align) ---
if [[ -d "$OUT_ROOT/$OUT_ID" ]]; then
    log "WARN removing prior partial output $OUT_ROOT/$OUT_ID"
    rm -rf "$OUT_ROOT/$OUT_ID"
fi

log "=== spaceranger count GSE169749 GSM5213483 (mouse colon d0, V19S23-097 A1) ==="
log "host=$HOSTNAME transcriptome=$TRANSCRIPTOME fastq=$FASTQ_DIR"
log "spaceranger binary: $SPACERANGER"

# --- 6. Run spaceranger ---
# --unknown-slide visium-1 = standard 11mm Visium capture-area layout (matches
# V19S23-097 A1). --create-bam true is REQUIRED so we get CB/UB-tagged reads
# for scAPAtrap. localmem=48 leaves ample headroom on the 1TB box; localcores=8
# (other jobs share the box).
(
    cd "$OUT_ROOT" && \
    "$SPACERANGER" count \
        --id "$OUT_ID" \
        --description "GSE169749 GSM5213483 d0 mouse colon Visium" \
        --transcriptome "$TRANSCRIPTOME" \
        --fastqs "$FASTQ_DIR" \
        --sample "$SAMPLE_NAME" \
        --unknown-slide visium-1 \
        --image "$HE_IMAGE" \
        --create-bam true \
        --localcores 8 \
        --localmem 48 \
        --disable-ui
) > "$SAMPLE_LOG" 2>&1
rc=$?

if [[ $rc -ne 0 ]]; then
    log "FAIL spaceranger rc=$rc; tail of log:"
    tail -n 50 "$SAMPLE_LOG" | sed 's/^/      /' | tee -a "$SAMPLE_LOG.err"
    exit $rc
fi

log "OK spaceranger -> $OUT_ROOT/$OUT_ID/outs"
log "--- outs/ listing ---"
ls -la "$OUT_ROOT/$OUT_ID/outs/" 2>/dev/null | tee -a "$SAMPLE_LOG"
log "--- spatial/ listing ---"
ls -la "$OUT_ROOT/$OUT_ID/outs/spatial/" 2>/dev/null | tee -a "$SAMPLE_LOG"
exit 0

#!/usr/bin/env bash
# =============================================================================
# Run scAPAtrap (via run_scapatrap_spaceranger.py --species mouse) on the
# GSE169749 GSM5213483 (mouse colon, day 0) Space Ranger BAM.
#
# Prereq: scripts/gse169749_gsm5213483_spaceranger.sh must have completed and
#         produced:
#           pipeline_output/gse169749/gsm5213483_d0_sr/outs/possorted_genome_bam.bam
#           pipeline_output/gse169749/gsm5213483_d0_sr/outs/spatial/tissue_positions.csv
#
# Uses GRCm39 mouse GTF (matches the spaceranger transcriptome).
# Hostname-aware (CLAUDE.md S90/S91/S97/S98).
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md) ---
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) export TMPDIR="/s2/mengzijun/tmp"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    S91) export TMPDIR="/s3/mengzijun/tmp"; PY="$HOME/anaconda3/envs/spagapa/bin/python" ;;
    S97) export TMPDIR="/s972/mengzijun/tmp"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    S98) export TMPDIR="/s982/mengzijun/tmp"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    *)  echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac
export TMPDIR
export OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

# --- 2. Paths ---
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SR_OUT="$ROOT/pipeline_output/gse169749/gsm5213483_d0_sr/outs"
BAM="$SR_OUT/possorted_genome_bam.bam"
SPATIAL="$SR_OUT/spatial"
GTF="/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"

OUT_ROOT="$ROOT/pipeline_output/gse169749_gsm5213483_scapatrap"
PROC_DIR="$ROOT/data/processed/gse169749_gsm5213483_scapatrap"
DATASET="gse169749_gsm5213483_scapatrap"

LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
LOG="$LOG_ROOT/${DATE_TAG}_gsm5213483_scapatrap.log"
mkdir -p "$LOG_ROOT"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

# --- 3. Preflight ---
if [[ ! -s "$BAM" ]]; then
    log "ERROR BAM not found: $BAM"
    log "Did scripts/gse169749_gsm5213483_spaceranger.sh complete?"
    exit 2
fi
if [[ ! -d "$SPATIAL" ]]; then
    log "ERROR spatial dir not found: $SPATIAL"
    exit 2
fi
if [[ ! -s "$GTF" ]]; then
    log "ERROR mouse GTF not found: $GTF"
    exit 2
fi

log "=== scAPAtrap GSE169749 GSM5213483 (mouse colon d0) ==="
log "BAM=$BAM"
log "SPATIAL=$SPATIAL"
log "GTF=$GTF (mouse GRCm39)"
log "OUT_ROOT=$OUT_ROOT"
log "PROC_DIR=$PROC_DIR"

# --- 4. Run launcher ---
# R2 is 120bp cDNA. scAPAtrap readlength param controls coverage window sizing;
# 90 is the conservative default used for 90bp Visium reads and works for 120bp.
"$PY" "$ROOT/scripts/run_scapatrap_spaceranger.py" \
    --bam "$BAM" \
    --spatial-dir "$SPATIAL" \
    --gtf "$GTF" \
    --output-root "$OUT_ROOT" \
    --processed-dir "$PROC_DIR" \
    --dataset-name "$DATASET" \
    --source-label "GSE169749 GSM5213483 d0 mouse colon (Space Ranger BAM + scAPAtrap)" \
    --tissue "mouse colon" \
    --species mouse \
    --readlength 90 \
    --threads 12 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}

if [[ $rc -ne 0 ]]; then
    log "FAIL scAPAtrap rc=$rc"
    exit $rc
fi

log "OK scAPAtrap complete"
log "--- qc_summary.json ---"
if [[ -f "$PROC_DIR/qc_summary.json" ]]; then
    cat "$PROC_DIR/qc_summary.json" | tee -a "$LOG"
fi
log "--- processed files ---"
ls -lh "$PROC_DIR/" 2>/dev/null | tee -a "$LOG"
exit 0

#!/usr/bin/env bash
# =============================================================================
# MASTER ORCHESTRATOR: Process both MOUSE Visium datasets for APA analysis.
#
#   Dataset 1: GSE169749  GSM5213483  (mouse colon, day 0)
#              - FASTQ already on disk (renamed to 10x bcl2fastq convention)
#              - pipeline: spaceranger count (mouse GRCm39) -> scAPAtrap
#
#   Dataset 2: GSE263303  GSM8189356 + GSM8189359  (mouse brain Nf1+/-)
#              - FASTQ must be downloaded via aria2c S3 then fasterq-dump
#              - pipeline: download -> fasterq-dump -> spaceranger count -> scAPAtrap
#
# Each phase is idempotent (skips completed work). Run inside tmux.
#
# Usage:
#   ./run_scapatrap_mouse_batch.sh                       # all phases, both datasets
#   ./run_scapatrap_mouse_batch.sh gse169749             # dataset 1 only
#   ./run_scapatrap_mouse_batch.sh gse263303             # dataset 2 only
#   ./run_scapatrap_mouse_batch.sh gse169749 scapatrap   # dataset 1, scapatrap phase only
#
# Phases: download | spaceranger | scapatrap  (gse169749 has no download phase)
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md mandatory) ---
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

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
MASTER_LOG="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs/${DATE_TAG}_mouse_batch_master.log"
mkdir -p "$(dirname "$MASTER_LOG")"

mlog() { printf '[%s] [MASTER] %s\n' "$(date '+%F %T')" "$*" | tee -a "${MASTER_LOG}"; }

DATASET="${1:-all}"
PHASE="${2:-all}"

mlog "=== MOUSE APA BATCH START (host=$HOSTNAME, dataset=$DATASET, phase=$PHASE) ==="

rc=0

# ---------------------------------------------------------------------------
# Dataset 1: GSE169749 GSM5213483 (mouse colon d0)
# ---------------------------------------------------------------------------
run_gse169749() {
    local phase="$1"
    if [[ "$phase" == "all" || "$phase" == "spaceranger" ]]; then
        mlog "--- GSE169749 GSM5213483: spaceranger count (mouse GRCm39) ---"
        if bash "$ROOT/scripts/gse169749_gsm5213483_spaceranger.sh"; then
            mlog "GSE169749 spaceranger OK"
        else
            mlog "GSE169749 spaceranger FAIL (rc=$?)"; rc=1
        fi
    fi
    if [[ "$phase" == "all" || "$phase" == "scapatrap" ]]; then
        mlog "--- GSE169749 GSM5213483: scAPAtrap (mouse) ---"
        if bash "$ROOT/scripts/run_scapatrap_gse169749_gsm5213483.sh"; then
            mlog "GSE169749 scAPAtrap OK"
        else
            mlog "GSE169749 scAPAtrap FAIL (rc=$?)"; rc=1
        fi
    fi
}

# ---------------------------------------------------------------------------
# Dataset 2: GSE263303 (mouse brain Nf1+/-) - delegates to existing master script
# ---------------------------------------------------------------------------
run_gse263303() {
    local phase="$1"
    case "$phase" in
        download)
            mlog "--- GSE263303: download (aria2c S3) ---"
            bash "$ROOT/scripts/download_gse263303_sra.sh" || { mlog "GSE263303 download FAIL"; rc=1; }
            ;;
        spaceranger|fasterq)
            mlog "--- GSE263303: fasterq-dump + spaceranger count (mouse GRCm39) ---"
            bash "$ROOT/scripts/run_gse263303_fasterq_spaceranger.sh" || { mlog "GSE263303 fasterq+SR FAIL"; rc=1; }
            ;;
        scapatrap)
            mlog "--- GSE263303: scAPAtrap (mouse) ---"
            bash "$ROOT/scripts/run_scapatrap_gse263303.sh" || { mlog "GSE263303 scAPAtrap FAIL"; rc=1; }
            ;;
        all)
            bash "$ROOT/scripts/run_gse263303_full_pipeline.sh" || { mlog "GSE263303 full pipeline FAIL"; rc=1; }
            ;;
        *) mlog "unknown phase $phase"; rc=1 ;;
    esac
}

case "$DATASET" in
    gse169749) run_gse169749 "$PHASE" ;;
    gse263303) run_gse263303 "$PHASE" ;;
    all)
        run_gse169749 "$PHASE"
        run_gse263303 "$PHASE"
        ;;
    *) echo "Usage: $0 [gse169749|gse263303|all] [download|spaceranger|scapatrap|all]"; exit 2 ;;
esac

mlog "=== MOUSE APA BATCH DONE (rc=$rc) ==="

# --- summary ---
mlog "=== SUMMARY ==="
for qc in \
    "$ROOT/data/processed/gse169749_gsm5213483_scapatrap/qc_summary.json" \
    "$ROOT/data/processed/gse263303_GSM8189356_scapatrap/qc_summary.json" \
    "$ROOT/data/processed/gse263303_GSM8189359_scapatrap/qc_summary.json"; do
    if [[ -f "$qc" ]]; then
        "$HOME/anaconda3/envs/spagapa/bin/python" -c "
import json
q=json.load(open('$qc'))
print('  $qc:')
print('    species=', q.get('species'), 'tissue=', q.get('tissue'))
print('    n_called_sites=', q.get('n_called_sites'), 'n_gene_annotated=', q.get('n_gene_annotated_sites'), 'n_apa_usage_sites=', q.get('n_apa_usage_sites'), 'n_spots=', q.get('n_spots'))
" 2>/dev/null | tee -a "$MASTER_LOG"
    else
        mlog "  (no qc_summary.json yet: $qc)"
    fi
done
exit $rc

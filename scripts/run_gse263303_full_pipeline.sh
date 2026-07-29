#!/usr/bin/env bash
# =============================================================================
# MASTER ORCHESTRATOR: GSE263303 (mouse brain, Nf1+/- Visium) -> APA analysis.
#
# Chains: download (aria2c S3) -> fasterq-dump -> spaceranger count -> scAPAtrap
#
# This driver runs the three phase scripts in sequence. Each phase is idempotent
# (skips completed work via done markers / output checks). Run inside tmux.
#
# Usage:
#   ./run_gse263303_full_pipeline.sh              # all phases
#   ./run_gse263303_full_pipeline.sh download     # phase only
#   ./run_gse263303_full_pipeline.sh fasterq      # phase only
#   ./run_gse263303_full_pipeline.sh scapatrap    # phase only
#
# Selected samples (2 of 4, spanning both Nf1+/- animals):
#   GSM8189356  K73-6-FMFC  SRR28566759 + SRR28566760   (animal K73)
#   GSM8189359  K75-2-FMFC  SRR28566753 + SRR28566754   (animal K75)
# Other samples NOT processed: GSM8189357 (K73-8), GSM8189358 (K75-1).
# =============================================================================
set -uo pipefail

# --- hostname detection (CLAUDE.md) ---
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
MASTER_LOG="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs/${DATE_TAG}_gse263303_master.log"
mkdir -p "$(dirname "$MASTER_LOG")"

mlog() { printf '[%s] [MASTER] %s\n' "$(date '+%F %T')" "$*" | tee -a "${MASTER_LOG}"; }

PHASE="${1:-all}"

mlog "=== GSE263303 master pipeline START (host=$HOSTNAME, phase=$PHASE) ==="
mlog "samples: GSM8189356 (K73-6), GSM8189359 (K75-2)  [mouse brain Nf1+/-]"

rc=0

# --- Phase 1: download ---
if [[ "$PHASE" == "all" || "$PHASE" == "download" ]]; then
    mlog "--- PHASE 1: download (aria2c S3) ---"
    if bash "$ROOT/scripts/download_gse263303_sra.sh"; then
        mlog "phase 1 OK"
    else
        mlog "phase 1 FAIL (rc=$?)"; rc=1
    fi
fi

# --- Phase 2: fasterq-dump + spaceranger ---
if [[ "$PHASE" == "all" || "$PHASE" == "fasterq" ]]; then
    mlog "--- PHASE 2: fasterq-dump + spaceranger count (mouse GRCm39) ---"
    if bash "$ROOT/scripts/run_gse263303_fasterq_spaceranger.sh"; then
        mlog "phase 2 OK"
    else
        mlog "phase 2 FAIL (rc=$?)"; rc=1
    fi
fi

# --- Phase 3: scAPAtrap (mouse) ---
if [[ "$PHASE" == "all" || "$PHASE" == "scapatrap" ]]; then
    mlog "--- PHASE 3: scAPAtrap (mouse-adapted launcher) ---"
    if bash "$ROOT/scripts/run_scapatrap_gse263303.sh"; then
        mlog "phase 3 OK"
    else
        mlog "phase 3 FAIL (rc=$?)"; rc=1
    fi
fi

mlog "=== MASTER DONE (rc=$rc) ==="

# --- summary ---
mlog "=== SUMMARY ==="
mlog "SRA downloads:"
ls -lh "$ROOT"/data/raw/gse263303/ 2>/dev/null | grep -E "SRR" | tee -a "$MASTER_LOG"
mlog "spaceranger outputs:"
ls -d "$ROOT"/pipeline_output/gse263303_*_sr 2>/dev/null | tee -a "$MASTER_LOG"
mlog "scapatrap processed dirs:"
ls -d "$ROOT"/data/processed/gse263303_*_scapatrap 2>/dev/null | tee -a "$MASTER_LOG"
for gsm in GSM8189356 GSM8189359; do
    qc="$ROOT/data/processed/gse263303_${gsm}_scapatrap/qc_summary.json"
    if [[ -f "$qc" ]]; then
        python3 -c "
import json
q=json.load(open('$qc'))
print('  $gsm: n_called_sites=', q['n_called_sites'], 'n_gene_annotated=', q['n_gene_annotated_sites'], 'n_apa_usage_sites=', q['n_apa_usage_sites'], 'n_spots=', q['n_spots'])
" 2>/dev/null | tee -a "$MASTER_LOG"
    fi
done
exit $rc

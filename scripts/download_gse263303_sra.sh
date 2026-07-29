#!/usr/bin/env bash
# =============================================================================
# Download GSE263303 (mouse brain, Nf1+/- model Visium) FASTQ via aria2c S3.
#
# Selected 2 representative samples (spanning both animals K73 & K75):
#   GSM8189356 (K73-6-FMFC) : SRR28566759 + SRR28566760
#   GSM8189359 (K75-2-FMFC) : SRR28566753 + SRR28566754
# Each GSM has 2 technical-split SRA runs (NovaSeq); both must be merged.
#
# Hostname-aware (CLAUDE.md S90/S91/S97/S98).
# aria2c with proxy DISABLED (proven method, ~60-70 MiB/s vs ~250 KiB/s through socks).
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md mandatory) ---
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) TMPDIR="/s2/mengzijun/tmp" ;;
    S91) TMPDIR="/s3/mengzijun/tmp" ;;
    S97) TMPDIR="/s972/mengzijun/tmp" ;;
    S98) TMPDIR="/s982/mengzijun/tmp" ;;
    *)  echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac
export TMPDIR
mkdir -p "$TMPDIR"

# --- 2. Paths ---
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT="$ROOT/data/raw/gse263303"
LOG_DIR="$OUT/logs"
SUITE_LOG="$LOG_DIR/download.log"
mkdir -p "$OUT" "$LOG_DIR"

# proxy DISABLED (system socks proxy throttles aria2c)
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

# GSM -> array of SRR run accessions (each GSM = 2 technical splits)
declare -A SRRS
SRRS[GSM8189356]="SRR28566759 SRR28566760"   # K73-6-FMFC
SRRS[GSM8189359]="SRR28566753 SRR28566754"   # K75-2-FMFC

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }

log "=== GSE263303 FASTQ download (mouse brain Nf1+/- Visium) ==="
log "host=$HOSTNAME  out=$OUT  tmp=$TMPDIR"
log "samples: GSM8189356 (K73-6), GSM8189359 (K75-2)"

# --- 3. Download each SRR from NCBI SRA S3 (Open Data Program) ---
for gsm in GSM8189356 GSM8189359; do
    log "--- ${gsm} ---"
    for SRR in ${SRRS[$gsm]}; do
        if [[ -f "$OUT/${SRR}" && ! -f "$OUT/${SRR}.aria2" ]]; then
            sz=$(du -h "$OUT/${SRR}" 2>/dev/null | cut -f1)
            log "SKIP ${SRR}: already complete (${sz})"
            continue
        fi
        log "START ${SRR} (https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR})"
        t0=$(date +%s)
        aria2c -c -x4 -s4 -d "$OUT" \
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}" \
            2>&1 | tail -3 | tee -a "${SUITE_LOG}"
        rc=${PIPESTATUS[0]}
        t1=$(date +%s)
        sz=$(du -h "$OUT/${SRR}" 2>/dev/null | cut -f1)
        log "END   ${SRR} (rc=${rc}, elapsed=$(( (t1-t0)/60 ))min, size=${sz})"
    done
done

log "=== DOWNLOAD COMPLETE ==="
ls -lh "$OUT"/*.sra "$OUT"/SRR28566* 2>/dev/null | tee -a "${SUITE_LOG}"
log "total: $(du -sh "$OUT" | cut -f1)"

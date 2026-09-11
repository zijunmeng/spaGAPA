#!/usr/bin/env bash
# ENA FASTQ downloader for GSE293464 (32 runs). ENA scales to ~6 MB/s per file
# with 16 connections (vs ODP's decayed <1 MB/s); files are ready-made gzip
# FASTQ with official MD5s (skips fasterq-dump entirely).
#
# Layout matches the fasterq phase output convention:
#   <base>/fasterq/<SRR>/<SRR>_{1,2}.fastq.gz   (merge phase works unchanged)
#
# Usage: bash download_gse293464_ena.sh [n_workers]
set -uo pipefail
HOSTNAME=$(hostname -s); [ "$HOSTNAME" = "S91" ] || { echo "not S91" >&2; exit 1; }
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

NW=${1:-6}
BASE=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse293464
MAN=/tmp/ena_manifest.tsv
LOGDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/logs
MLOG="$LOGDIR/20260909_stereo_ena_download.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$MLOG"; }

worker() {
    local wid=$1
    tail -n +2 "$MAN" | awk -v w="$wid" -v n="$NW" 'NR % n == w' | while IFS=$'\t' read -r srr fq1 fq2 b1 b2 m1 m2; do
        [ -z "$srr" ] && continue
        local d="$BASE/fasterq/$srr"; mkdir -p "$d"
        local ok=1
        for spec in "1:$fq1:$b1:$m1" "2:$fq2:$b2:$m2"; do
            IFS=':' read -r idx url bytes md5 <<< "$spec"
            local f="$d/${srr}_${idx}.fastq.gz"
            if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then
                continue
            fi
            local tries=0
            while :; do
                tries=$((tries+1))
                timeout 600 aria2c -c -x16 -s16 -k20M --check-certificate=false \
                    --allow-overwrite=true --auto-file-renaming=false \
                    --console-log-level=error --summary-interval=0 \
                    --max-tries=2 --retry-wait=5 --timeout=60 --connect-timeout=20 \
                    -d "$d" -o "${srr}_${idx}.fastq.gz" "https://$url" > /dev/null 2>&1
                if [ -s "$f" ] && [ "$(stat -c%s "$f")" = "$bytes" ] \
                   && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then
                    log "w$wid OK $srr _${idx}.fastq.gz ($bytes bytes, $tries tries)"; break
                fi
                [ $tries -ge 15 ] && { log "w$wid GIVEUP $srr _${idx} (got $(stat -c%s "$f" 2>/dev/null || echo 0)/$bytes)"; ok=0; break; }
            done
        done
        [ $ok -eq 1 ] && touch "$d/.ena_complete"
    done
    log "w$wid queue finished"
}

log "=== ENA start: $NW workers, $(($(wc -l < "$MAN") - 1)) runs x 2 files ==="
pids=()
for w in $(seq 0 $((NW-1))); do worker $w & pids+=($!); done
for p in "${pids[@]}"; do wait "$p"; done
log "=== ENA ALL DONE: $(find "$BASE/fasterq" -name .ena_complete | wc -l)/32 runs complete ==="

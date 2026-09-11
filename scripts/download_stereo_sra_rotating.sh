#!/usr/bin/env bash
# Rotating-parallel SRA downloader for the stereo expansion (2026-09-09).
#
# Empirically (S91, 2026-09-09): fresh ODP transfers start at ~20-70 MB/s but
# DECAY to ~1 MB/s after several minutes (per-transfer, not per-IP: a fresh
# parallel transfer runs at full speed while an old one crawls). Fix: N workers
# in parallel, each restarting its aria2c every ROTATE s (timeout) — aria2c -c
# resumes from the control file with fresh connections.
#
# Usage: bash download_stereo_sra_rotating.sh [n_workers] [rotate_seconds]
set -uo pipefail
HOSTNAME=$(hostname -s)
[ "$HOSTNAME" = "S91" ] || { echo "ERROR: not S91" >&2; exit 1; }
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

NW=${1:-4}; ROT=${2:-300}
BASE=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw
LOGDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/logs
MLOG="$LOGDIR/20260909_stereo_sra_rotating.log"
log() { echo "[$(date '+%F %T')] $*" >> "$MLOG"; }

ALL=$( (for i in $(seq 32936115 32936146); do echo "SRR$i $BASE/gse293464/sra"; done; echo "SRR38888148 $BASE/gse333693/sra") )

expected_bytes() { curl -sI "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" --max-time 20 | grep -i content-length | tr -dc '0-9'; }

worker() {
    local wid=$1
    local items
    items=$(echo "$ALL" | awk -v w="$wid" -v n="$NW" 'NR % n == w {print}')
    while read -r srr dest; do
        [ -z "$srr" ] && continue
        local f="$dest/$srr" exp
        exp=$(expected_bytes "$srr")
        if [[ -n "$exp" && -s "$f" && ! -f "$f.aria2" && $(stat -c%s "$f" 2>/dev/null || echo 0) -eq "$exp" ]]; then
            log "w$wid skip $srr (complete)"; continue
        fi
        local tries=0
        while :; do
            tries=$((tries+1))
            timeout "$ROT" aria2c -c -x4 -s4 --check-certificate=false \
                --allow-overwrite=false --auto-file-renaming=false \
                --console-log-level=error --summary-interval=0 \
                --max-tries=2 --retry-wait=5 --timeout=60 --connect-timeout=20 \
                -d "$dest" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$srr/$srr" \
                > /dev/null 2>&1
            local sz; sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
            if [[ -n "$exp" && "$sz" = "$exp" && ! -f "$f.aria2" ]]; then
                log "w$wid DONE $srr ($sz bytes, $tries rotations)"; break
            fi
            [ $tries -ge 60 ] && { log "w$wid GIVEUP $srr after $tries rotations (size=$sz/exp=$exp)"; break; }
        done
    done <<< "$items"
    log "w$wid finished its queue"
}

log "=== start: $NW workers, rotate=${ROT}s, $(echo "$ALL" | wc -l) files ==="
pids=()
for w in $(seq 0 $((NW-1))); do worker $w & pids+=($!); done
for p in "${pids[@]}"; do wait "$p"; done
log "=== ALL DONE: gse293464 $(ls "$BASE/gse293464/sra" | grep -c '^SRR')/32 gse333693 $(ls "$BASE/gse333693/sra" | grep -c '^SRR')/1 ==="
log "sizes: $(du -sh "$BASE/gse293464/sra" "$BASE/gse333693/sra" 2>/dev/null | tr '\n' ' ')"

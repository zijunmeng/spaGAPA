#!/usr/bin/env bash
# Sequential single-file SRA downloads for the stereo expansion (2026-09-09).
# Replaces the -j3 batch mode that ODP throttled to ~1 MiB/s; the proven
# single-file pattern (download_gse263303_sra.sh) sustains 19-70 MB/s.
# Idempotent: skips complete files (no .aria2 + nonzero size), aria2c -c resumes.
# Includes stall guard: kills+retries an aria2c whose bytes-written stop growing.
set -uo pipefail
HOSTNAME=$(hostname -s)
[ "$HOSTNAME" = "S91" ] || { echo "ERROR: not S91" >&2; exit 1; }
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY all_proxy ALL_PROXY

BASE=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw
LOGDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/logs
MLOG="$LOGDIR/20260909_stereo_sra_sequential.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$MLOG"; }

declare -A DEST
for srr in $(seq 32936115 32936146); do DEST[SRR$srr]="$BASE/gse293464/sra"; done
DEST[SRR38888148]="$BASE/gse333693/sra"

expected_bytes() { curl -sI "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" --max-time 15 | grep -i content-length | tr -dc '0-9'; }

for srr in $(seq 32936115 32936146) 38888148; do
    srr="SRR$srr"; d="${DEST[$srr]}"; f="$d/$srr"
    exp=$(expected_bytes "$srr")
    if [[ -n "$exp" && -s "$f" && ! -f "$f.aria2" && $(stat -c%s "$f") -eq "$exp" ]]; then
        log "skip $srr (complete, $exp bytes)"; continue
    fi
    ok=0
    for attempt in 1 2 3 4 5; do
        log "dl $srr attempt $attempt (target $exp bytes)"
        aria2c -c -x4 -s4 --check-certificate=false --allow-overwrite=false \
            --auto-file-renaming=false --console-log-level=warn --summary-interval=120 \
            --max-tries=3 --retry-wait=10 --timeout=120 --connect-timeout=30 \
            --log="$LOGDIR/20260909_aria2_sra_single.log" \
            -d "$d" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$srr/$srr" \
            >> "$MLOG" 2>&1 &
        apid=$!
        # stall guard: check bytes-written every 45 s; kill if no growth twice
        last=-1; stalls=0
        while kill -0 $apid 2>/dev/null; do
            sleep 45
            cur=$(stat -c%s "$f" 2>/dev/null || echo 0)
            # aria2 preallocates; use .aria2 control mtime + proc io instead
            io=$(grep write_bytes /proc/$apid/io 2>/dev/null | head -1 | tr -dc '0-9' || echo 0)
            cur=${io:-0}
            if [ "$cur" = "$last" ]; then stalls=$((stalls+1)); else stalls=0; fi
            last=$cur
            [ $stalls -ge 2 ] && { log "  stall detected on $srr (io=$cur), killing"; kill $apid 2>/dev/null; break; }
        done
        wait $apid 2>/dev/null; rc=$?
        sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
        if [[ -n "$exp" && "$sz" = "$exp" && ! -f "$f.aria2" ]]; then
            log "  $srr OK ($sz bytes, attempt $attempt)"; ok=1; break
        fi
        log "  $srr incomplete (size=$sz, rc=$rc), retrying"
        sleep 10
    done
    [ $ok -eq 1 ] || log "  $srr FAILED after 5 attempts"
done
log "ALL DONE: gse293464 $(ls "$BASE/gse293464/sra" | grep -c '^SRR')/32, gse333693 $(ls "$BASE/gse333693/sra" | grep -c '^SRR')/1"
log "total: $(du -sh "$BASE/gse293464/sra" "$BASE/gse333693/sra" | tr '\n' ' ')"

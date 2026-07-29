#!/usr/bin/env bash
# =============================================================================
# Normalize R1 read length to a fixed value (default 28bp) for spaceranger.
#
# Background: GSE220442 SRA archives store R1 with VARIABLE length
# (read_length="variable" in sra-stat). The 16bp barcode + 12bp UMI = 28bp
# nominal, but NCBI trimmed low-quality 3' bases per-spot during archiving,
# producing reads of 18-28bp. spaceranger errors with TXRNGR10009 ("mixture
# of different R1 lengths"). --r1-length N cannot fix this when reads are
# shorter than N (hard-trim only).
#
# Fix: 3'-pad each R1 read with N to TARGET_LEN (default 28). Barcode (5',
# positions 1-16) is preserved exactly; only 3' UMI positions are padded,
# which spaceranger tolerates. Quality strings are padded with '#' (phred 2).
# Reads longer than TARGET_LEN are 3'-trimmed (rare; preserves 5' barcode).
#
# Handles both uncompressed (.fastq) and gzipped (.fastq.gz) inputs.
# Writes gzipped output (uniform, space-efficient).
# =============================================================================
set -euo pipefail

# --- Host detection (CLAUDE.md) ---
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90)  export TMPDIR="/s2/mengzijun/tmp"  ;;
    S91)  export TMPDIR="/s3/mengzijun/tmp"  ;;
    S97)  export TMPDIR="/s972/mengzijun/tmp" ;;
    S98)  export TMPDIR="/s982/mengzijun/tmp" ;;
    *) echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
DATA_ROOT="${ROOT}/data/raw/gse220442/spaceranger_inputs"
LOG_ROOT="${ROOT}/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse220442_normalize_r1.log"
mkdir -p "$LOG_ROOT"

TARGET_LEN="${TARGET_LEN:-28}"
SAMPLES=(
  GSM6801751 GSM6801752 GSM6801753
  GSM6801754 GSM6801755 GSM6801756
)

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }

normalize_one_file() {
    # $1 = input fastq (.fastq or .fastq.gz)
    # $2 = output fastq.gz path
    local in="$1" out="$2"
    local n="${TARGET_LEN}"
    # Build padding strings once
    local pad_seq pad_qual
    pad_seq=$(printf 'N%.0s' $(seq 1 "$n"))
    pad_qual=$(printf '#%.0s' $(seq 1 "$n"))
    # awk reads 4 lines at a time; pads seq (line 2) and qual (line 4) to n.
    # Handles both cases: shorter (3'-pad N/#) and longer (3'-trim to n).
    local reader
    if [[ "$in" == *.gz ]]; then reader="gunzip -c $in"; else reader="cat $in"; fi
    eval "$reader" | awk -v n="$n" -v padS="$pad_seq" -v padQ="$pad_qual" '
        BEGIN { padS = substr(padS, 1, n); padQ = substr(padQ, 1, n) }
        {
            hdr = $0
            getline seq
            getline plus
            getline qual
            L = length(seq)
            if (L < n) {
                seq  = seq  substr(padS, 1, n - L)
                qual = qual substr(padQ, 1, n - L)
            } else if (L > n) {
                seq  = substr(seq,  1, n)
                qual = substr(qual, 1, n)
            }
            print hdr
            print seq
            print plus
            print qual
        }
    ' | gzip -c > "$out"
}

normalize_sample() {
    local gsm="$1"
    local in_dir="${DATA_ROOT}/${gsm}/fastqs_sra_split"
    local out_dir="${DATA_ROOT}/${gsm}/fastqs_sra_split_r1norm28"
    rm -rf "$out_dir"; mkdir -p "$out_dir"
    log "=== ${gsm}: normalize R1 -> ${TARGET_LEN}bp (out: ${out_dir}) ==="
    # R1: normalize. R2: pass-through (symlink if gz, gzip if uncompressed).
    local f base stem r1_in r1_out r2_in r2_out
    for f in "$in_dir"/*; do
        base=$(basename "$f")
        # Compute a clean stem (strip trailing .fastq.gz or .fastq)
        if [[ "$base" == *.fastq.gz ]]; then
            stem="${base%.fastq.gz}"
        elif [[ "$base" == *.fastq ]]; then
            stem="${base%.fastq}"
        else
            stem="$base"
        fi
        case "$base" in
            *_R1_001.fastq|*_R1_001.fastq.gz)
                r1_in="$f"
                r1_out="${out_dir}/${stem}.fastq.gz"
                log "  normalize R1: $base -> $(basename "$r1_out")"
                t0=$(date +%s)
                normalize_one_file "$r1_in" "$r1_out"
                t1=$(date +%s)
                log "    done (elapsed=$((t1-t0))s, $(stat -L -c '%s' "$r1_out" | numfmt --to=iec))"
                ;;
            *_R2_001.fastq|*_R2_001.fastq.gz)
                r2_in="$f"
                r2_out="${out_dir}/${stem}.fastq.gz"
                # If input is already gzipped, symlink; else gzip-copy.
                if [[ "$r2_in" == *.gz ]]; then
                    ln -sf "$r2_in" "$r2_out"
                    log "  symlink R2: $base -> $(basename "$r2_out")"
                else
                    log "  gzip R2: $base -> $(basename "$r2_out")"
                    t0=$(date +%s)
                    gzip -c "$r2_in" > "$r2_out"
                    t1=$(date +%s)
                    log "    done (elapsed=$((t1-t0))s)"
                fi
                ;;
        esac
    done
    # Verify all 4 files exist + R1 uniform length
    local nf
    nf=$(find "$out_dir" -maxdepth 1 -type f -o -type l | wc -l)
    if (( nf != 4 )); then
        log "  ERROR(${gsm}): expected 4 output files, found ${nf}"
        return 2
    fi
    # Quick length check on first R1
    local r1
    r1=$(ls "$out_dir"/*_R1_001.fastq.gz | head -1)
    local lens
    lens=$(zcat "$r1" | head -4000 | awk 'NR%4==2{print length($0)}' | sort -u | tr '\n' ',')
    log "  verify R1 lengths (first 1000 reads): ${lens}"
}

run_all() {
    log "HOSTNAME=${HOSTNAME}  TMPDIR=${TMPDIR}  TARGET_LEN=${TARGET_LEN}"
    log "samples=(${SAMPLES[*]})  PARALLEL=${PARALLEL:-1}"
    echo "----------------------------------------" | tee -a "${SUITE_LOG}"
    local par="${PARALLEL:-1}"
    if (( par <= 1 )); then
        for gsm in "${SAMPLES[@]}"; do
            normalize_sample "$gsm"
        done
    else
        # Parallel: launch up to $par samples at once
        declare -A pid2gsm
        running=0; i=0
        while (( i < ${#SAMPLES[@]} )) || (( running > 0 )); do
            while (( i < ${#SAMPLES[@]} )) && (( running < par )); do
                local gsm="${SAMPLES[$i]}"
                i=$((i+1))
                ( normalize_sample "$gsm" ) &
                pid2gsm[$!]=$gsm
                running=$((running+1))
                log "launched ${gsm} (pid=$!)"
            done
            if (( running > 0 )); then
                wait -n
                for pid in "${!pid2gsm[@]}"; do
                    if ! kill -0 "$pid" 2>/dev/null; then
                        log "finished ${pid2gsm[$pid]} (pid=$pid)"
                        unset 'pid2gsm[$pid]'
                        running=$((running-1))
                    fi
                done
                (( running > 0 )) && sleep 2
            fi
        done
    fi
    log "ALL NORMALIZE DONE"
}

usage() {
    cat <<EOF
Usage:
  $0 run-all      Normalize R1 for all 6 samples.
  $0 run-one GSM  Normalize a single sample.
Env: TARGET_LEN=${TARGET_LEN}  PARALLEL=${PARALLEL:-1}
EOF
}

MODE="${1:-help}"
case "$MODE" in
    run-all)   run_all ;;
    run-one)   normalize_sample "${2:-}" ;;
    tmux-submit)
        SESSION="${SESSION:-spagapa_gse220442_normr1}"
        if tmux has-session -t "${SESSION}" >/dev/null 2>&1; then
            echo "tmux session already exists: ${SESSION}" >&2; exit 1
        fi
        tmux new-session -d -s "${SESSION}" -c "${ROOT}" \
            "PARALLEL=${PARALLEL:-3} TARGET_LEN=${TARGET_LEN} \
             ${0} run-all > ${SUITE_LOG}.tmux 2>&1"
        echo "Submitted tmux session: ${SESSION}"
        echo "Suite log: ${SUITE_LOG} (and ${SUITE_LOG}.tmux)"
        ;;
    help|-h|--help) usage ;;
    *) echo "Unknown mode: $MODE" >&2; usage; exit 2 ;;
esac

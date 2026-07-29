#!/usr/bin/env bash
# =============================================================================
# Formal Space Ranger count for GSE220442 (human AD brain Visium).
# Reads SRA-split FASTQs (read3->R1 barcode/UMI, read4->R2 cDNA).
# Follows CLAUDE.md host-detection rules (S90/S91/S97/S98).
# Runs N samples concurrently with per-sample cores/mem; times every command.
# =============================================================================
set -uo pipefail

# -----------------------------------------------------------------------------
# 1. Hostname detection (CLAUDE.md mandatory)
# -----------------------------------------------------------------------------
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s2/mengzijun/tmp"
        ;;
    S91)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s3/mengzijun/tmp"
        ;;
    S97)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s972/mengzijun/tmp"
        ;;
    S98)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s982/mengzijun/tmp"
        ;;
    *)
        echo "ERROR: unknown host $HOSTNAME" >&2
        exit 1
        ;;
esac
mkdir -p "$TMPDIR"

# -----------------------------------------------------------------------------
# 2. Paths + parameters
# -----------------------------------------------------------------------------
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A"
DATA_ROOT="${ROOT}/data/raw/gse220442"
OUT_ROOT="${ROOT}/spaGAPA/pipeline_output"
LOG_ROOT="${ROOT}/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse220442_spaceranger_count.log"

# Per-sample resources
THREADS="${THREADS:-16}"
MEM_GB="${MEM_GB:-64}"
# Concurrency (samples run in parallel)
CONCURRENCY="${CONCURRENCY:-2}"

mkdir -p "${OUT_ROOT}" "${LOG_ROOT}" "$TMPDIR"

SAMPLES=(
  GSM6801751
  GSM6801752
  GSM6801753
  GSM6801754
  GSM6801755
  GSM6801756
)

declare -A TITLE
declare -A VISIUM_ID
TITLE[GSM6801751]="control 1";   VISIUM_ID[GSM6801751]="1-1"
TITLE[GSM6801752]="control 2";   VISIUM_ID[GSM6801752]="2-5"
TITLE[GSM6801753]="control 3";   VISIUM_ID[GSM6801753]="18-64"
TITLE[GSM6801754]="AD 1";        VISIUM_ID[GSM6801754]="2-3"
TITLE[GSM6801755]="AD 2";        VISIUM_ID[GSM6801755]="2-8"
TITLE[GSM6801756]="AD 3";        VISIUM_ID[GSM6801756]="T4857"

# -----------------------------------------------------------------------------
# 3. Timing helper
# -----------------------------------------------------------------------------
SECONDS=0
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }
timed() {
    # Run a command, print elapsed. $1=label, rest=cmd.
    local label="$1"; shift
    local t0 t1 rc
    t0=$(date +%s.%N)
    log "START ${label}: $*"
    "$@"
    rc=$?
    t1=$(date +%s.%N)
    log "END   ${label} (rc=${rc}, elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f min",(b-a)/60}'))"
    return $rc
}

# -----------------------------------------------------------------------------
# 4. Per-sample count
# -----------------------------------------------------------------------------
run_count() {
    local gsm="$1"
    local title="${TITLE[$gsm]}"
    local visium="${VISIUM_ID[$gsm]}"
    # FASTQ_SUBDIR defaults to the R1-normalized dir (28bp) which fixes the
    # TXRNGR10009 variable-R1-length error from SRA's per-spot quality trimming.
    # Override with FASTQ_SUBDIR=fastqs_sra_split to use the raw split dir.
    local fastq_sub="${FASTQ_SUBDIR:-fastqs_sra_split_r1norm28}"
    local fastq_dir="${DATA_ROOT}/spaceranger_inputs/${gsm}/${fastq_sub}"
    local image="${DATA_ROOT}/spaceranger_inputs/${gsm}/detected_tissue_image.jpg"
    local out_id="gse220442_${gsm}_sr"
    local sample_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_spaceranger_count.log"

    log "=== ${gsm} (${title}, visium=${visium}) ==="
    log "fastq_dir=${fastq_dir}"
    log "image=${image}"
    log "out_id=${out_id} (output dir: ${OUT_ROOT}/${out_id})"
    log "log=${sample_log}"

    # Preflight
    if [[ ! -d "$fastq_dir" ]]; then
        log "ERROR(${gsm}): fastq dir missing -> ${fastq_dir}; skipping"
        return 2
    fi
    if [[ ! -f "$image" ]]; then
        log "ERROR(${gsm}): image missing -> ${image}; skipping"
        return 2
    fi
    local nfiles
    nfiles=$(find "$fastq_dir" -maxdepth 1 -type f -o -type l | wc -l)
    if (( nfiles < 4 )); then
        log "ERROR(${gsm}): expected >=4 fastq files/symlinks in ${fastq_dir}, found ${nfiles}; skipping"
        return 2
    fi
    # Skip if already complete
    if [[ -f "${OUT_ROOT}/${out_id}/outs/count_summary.html" || -f "${OUT_ROOT}/${out_id}/outs/web_summary.html" ]]; then
        log "SKIP(${gsm}): output already complete at ${OUT_ROOT}/${out_id}/outs"
        return 0
    fi
    # Clean prior partial
    if [[ -d "${OUT_ROOT}/${out_id}" ]]; then
        log "WARN(${gsm}): removing prior partial output ${OUT_ROOT}/${out_id}"
        rm -rf "${OUT_ROOT}/${out_id}"
    fi

    # Run from OUT_ROOT so --id lands inside it.
    local t0 t1 rc
    t0=$(date +%s.%N)
    log "START spaceranger count ${gsm} (cores=${THREADS}, mem=${MEM_GB}GB)"
    (
        cd "${OUT_ROOT}" && \
        "${SPACERANGER}" count \
            --id "${out_id}" \
            --description "GSE220442 ${gsm} ${title} (visium slide ${visium})" \
            --transcriptome "${TRANSCRIPTOME}" \
            --fastqs "${fastq_dir}" \
            --sample "${gsm}" \
            --image "${image}" \
            --unknown-slide visium-1 \
            --create-bam true \
            --localcores "${THREADS}" \
            --localmem "${MEM_GB}" \
            --disable-cell-annotation \
            --disable-ui
    ) > "${sample_log}" 2>&1
    rc=$?
    t1=$(date +%s.%N)
    log "END   spaceranger count ${gsm} (rc=${rc}, elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f min",(b-a)/60}'))"

    if (( rc != 0 )); then
        log "FAIL(${gsm}): spaceranger rc=${rc}; tail of log:"
        tail -n 30 "${sample_log}" 2>/dev/null | sed 's/^/      /' | tee -a "${SUITE_LOG}"
    else
        log "OK(${gsm}): output at ${OUT_ROOT}/${out_id}"
    fi
    return $rc
}

# -----------------------------------------------------------------------------
# 5. Concurrency runner (max $CONCURRENCY parallel samples)
# -----------------------------------------------------------------------------
run_all() {
    log "HOSTNAME=${HOSTNAME}  TMPDIR=${TMPDIR}"
    log "CONCURRENCY=${CONCURRENCY}  THREADS/sample=${THREADS}  MEM_GB/sample=${MEM_GB}"
    log "samples=(${SAMPLES[*]})"
    log "suite log: ${SUITE_LOG}"
    echo "----------------------------------------" | tee -a "${SUITE_LOG}"

    declare -A pid2gsm
    running=0
    i=0
    while (( i < ${#SAMPLES[@]} )) || (( running > 0 )); do
        # Launch while under concurrency and samples remain
        while (( i < ${#SAMPLES[@]} )) && (( running < CONCURRENCY )); do
            local gsm="${SAMPLES[$i]}"
            i=$((i + 1))
            (
                run_count "$gsm"
            ) &
            pid2gsm[$!]=$gsm
            running=$((running + 1))
            log "launched ${gsm} (pid=$!)"
        done
        # Wait for any one to finish
        if (( running > 0 )); then
            wait -n
            # find which pid exited by scanning
            local still=""
            for pid in "${!pid2gsm[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    still+=" $pid"
                else
                    log "finished ${pid2gsm[$pid]} (pid=$pid)"
                    unset 'pid2gsm[$pid]'
                    running=$((running - 1))
                fi
            done
            (( running > 0 )) && sleep 2
        fi
    done
    log "ALL DONE. Total elapsed: ${SECONDS}s ($(awk -v s="$SECONDS" 'BEGIN{printf "%.1f min",s/60}'))"
}

status() {
    log "=== STATUS ==="
    for gsm in "${SAMPLES[@]}"; do
        local out_id="gse220442_${gsm}_sr"
        local out="${OUT_ROOT}/${out_id}"
        local state="MISSING"
        if [[ -f "${out}/outs/web_summary.html" ]]; then
            state="COMPLETE"
        elif [[ -f "${out}/_log" || -d "${out}" ]]; then
            state="PARTIAL/RUNNING"
        fi
        log "${gsm}\t${state}\t${out}"
    done
}

usage() {
    cat <<EOF
Usage:
  $0 run-all      Launch all 6 samples (concurrency=${CONCURRENCY}).
  $0 run-one GSM  Launch a single sample (foreground).
  $0 status       Print per-sample output status.
  $0 tmux-submit  Launch run-all inside persistent tmux session
                  'spagapa_gse220442_count'.

Env overrides:
  CONCURRENCY=${CONCURRENCY}  THREADS=${THREADS}  MEM_GB=${MEM_GB}  DATE_TAG=...
EOF
}

MODE="${1:-help}"
case "$MODE" in
    run-all)   run_all ;;
    run-one)   run_count "${2:-}" ;;
    status)    status ;;
    tmux-submit)
        SESSION="spagapa_gse220442_count"
        if tmux has-session -t "${SESSION}" >/dev/null 2>&1; then
            echo "tmux session already exists: ${SESSION}" >&2
            exit 1
        fi
        tmux new-session -d -s "${SESSION}" -c "${ROOT}" \
            "CONCURRENCY=${CONCURRENCY} THREADS=${THREADS} MEM_GB=${MEM_GB} \
             ${0} run-all > ${SUITE_LOG}.tmux 2>&1"
        echo "Submitted tmux session: ${SESSION}"
        echo "Suite log: ${SUITE_LOG} (and ${SUITE_LOG}.tmux)"
        echo "Attach:    tmux attach -t ${SESSION}"
        ;;
    help|-h|--help) usage ;;
    *) echo "Unknown mode: $MODE" >&2; usage; exit 2 ;;
esac

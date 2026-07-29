#!/usr/bin/env bash
# =============================================================================
# Space Ranger count for GSE206391 (human skin, psoriasis lesional) Visium.
# 2 chosen samples (smallest psoriasis-lesional FASTQ):
#   GSM6252925 (Patient 4 lesional psoriasis, rep 1) -- SRR19737306-9
#   GSM6252926 (Patient 4 lesional psoriasis, rep 2) -- SRR19737310-3
#
# ENA read layout (standard Visium, NO SRA-split needed):
#   _1.fastq.gz = 28 bp  -> Visium R1 (16 bp barcode + 12 bp UMI)
#   _2.fastq.gz = 120 bp -> Visium R2 (cDNA)
# Each GSM has 4 SRR runs that we assign to lanes L001-L004.
#
# Follows CLAUDE.md host-detection rules (S90/S91/S97/S98).
# OPENBLAS_NUM_THREADS=8 (per task). GRCh38 reference.
# Do NOT git commit.
# =============================================================================
set -uo pipefail

# -----------------------------------------------------------------------------
# 1. Hostname detection (CLAUDE.md mandatory)
# -----------------------------------------------------------------------------
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/home/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/home/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s3/mengzijun/tmp"
        ;;
    *)
        echo "ERROR: this script is tuned for S91 (host=$HOSTNAME)" >&2
        exit 1
        ;;
esac
export OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

# -----------------------------------------------------------------------------
# 2. Paths + parameters
# -----------------------------------------------------------------------------
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A"
DATA_ROOT="${ROOT}/data/raw/gse206391"
FASTQ_SRC="${DATA_ROOT}/fastq_ena"
SUPP_SRC="${DATA_ROOT}/geo_supplementary"
OUT_ROOT="${ROOT}/spaGAPA/pipeline_output"
LOG_ROOT="${ROOT}/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse206391_spaceranger_count.log"

THREADS="${THREADS:-16}"
MEM_GB="${MEM_GB:-48}"
CONCURRENCY="${CONCURRENCY:-2}"

mkdir -p "${OUT_ROOT}" "${LOG_ROOT}" "$TMPDIR"

# Sample -> ordered SRR list (becomes L001..L004)
SAMPLES=(GSM6252925 GSM6252926)
declare -A SRUN_LIST
SRUN_LIST[GSM6252925]="SRR19737306 SRR19737307 SRR19737308 SRR19737309"
SRUN_LIST[GSM6252926]="SRR19737310 SRR19737311 SRR19737312 SRR19737313"
declare -A TITLE
TITLE[GSM6252925]="P4 lesional psoriasis rep1"
TITLE[GSM6252926]="P4 lesional psoriasis rep2"

SECONDS=0
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }

# -----------------------------------------------------------------------------
# 3. Prepare spaceranger input dir: symlink ENA FASTQ into 10x naming.
#    Pattern: <sample>_S1_L00X_R{1,2}_001.fastq.gz
# -----------------------------------------------------------------------------
prepare_links() {
    local gsm="$1"
    local in_dir="${DATA_ROOT}/spaceranger_inputs/${gsm}/fastqs"
    rm -f "${in_dir}"/*.fastq.gz
    mkdir -p "${in_dir}"
    local lane=1
    for srr in ${SRUN_LIST[$gsm]}; do
        local r1="${FASTQ_SRC}/${srr}_1.fastq.gz"
        local r2="${FASTQ_SRC}/${srr}_2.fastq.gz"
        if [[ ! -s "$r1" || ! -s "$r2" ]]; then
            log "ERROR(${gsm}): missing FASTQ ${r1} or ${r2}"
            return 1
        fi
        local L=$(printf '%03d' "$lane")
        ln -sf "$r1" "${in_dir}/${gsm}_S1_L${L}_R1_001.fastq.gz"
        ln -sf "$r2" "${in_dir}/${gsm}_S1_L${L}_R2_001.fastq.gz"
        lane=$((lane + 1))
    done
    echo "$in_dir"
}

# -----------------------------------------------------------------------------
# 4. Extract tissue image for spaceranger --image
# -----------------------------------------------------------------------------
prepare_image() {
    local gsm="$1"
    local out_dir="${DATA_ROOT}/spaceranger_inputs/${gsm}"
    local img="${out_dir}/detected_tissue_image.jpg"
    mkdir -p "$out_dir"
    if [[ -s "$img" ]]; then echo "$img"; return 0; fi
    # Locate the downloaded .jpg.gz for this GSM
    local gz
    gz=$(ls "${SUPP_SRC}"/${gsm}_*.jpg.gz 2>/dev/null | head -1)
    if [[ -z "$gz" ]]; then
        gz=$(ls "${SUPP_SRC}"/*${gsm}*.jpg.gz 2>/dev/null | head -1)
    fi
    if [[ -z "$gz" || ! -s "$gz" ]]; then
        log "ERROR(${gsm}): no tissue image jpg.gz found in ${SUPP_SRC}"
        return 1
    fi
    log "Extracting image: $gz -> $img"
    gunzip -c "$gz" > "$img"
    echo "$img"
}

# -----------------------------------------------------------------------------
# 5. Per-sample count
# -----------------------------------------------------------------------------
run_count() {
    local gsm="$1"
    local title="${TITLE[$gsm]}"
    local fastq_dir out_id sample_log img
    fastq_dir=$(prepare_links "$gsm") || return 2
    img=$(prepare_image "$gsm") || return 2
    out_id="gse206391_${gsm}_sr"
    sample_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_spaceranger_count.log"

    log "=== ${gsm} (${title}) ==="
    log "fastq_dir=${fastq_dir}  image=${img}"
    log "out_id=${out_id} (output: ${OUT_ROOT}/${out_id})  log=${sample_log}"

    if [[ -f "${OUT_ROOT}/${out_id}/outs/web_summary.html" ]]; then
        log "SKIP(${gsm}): already complete"
        return 0
    fi
    if [[ -d "${OUT_ROOT}/${out_id}" ]]; then
        log "WARN(${gsm}): removing prior partial output"
        rm -rf "${OUT_ROOT}/${out_id}"
    fi

    local t0 t1 rc
    t0=$(date +%s.%N)
    log "START spaceranger count ${gsm} (cores=${THREADS}, mem=${MEM_GB}GB)"
    (
        cd "${OUT_ROOT}" && \
        "${SPACERANGER}" count \
            --id "${out_id}" \
            --description "GSE206391 ${gsm} ${title}" \
            --transcriptome "${TRANSCRIPTOME}" \
            --fastqs "${fastq_dir}" \
            --sample "${gsm}" \
            --image "${img}" \
            --unknown-slide visium-1 \
            --create-bam true \
            --localcores "${THREADS}" \
            --localmem "${MEM_GB}" \
            --disable-cell-annotation \
            --disable-ui
    ) > "${sample_log}" 2>&1
    rc=$?
    t1=$(date +%s.%N)
    log "END spaceranger count ${gsm} (rc=${rc}, elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f min",(b-a)/60}'))"
    if (( rc != 0 )); then
        log "FAIL(${gsm}): spaceranger rc=${rc}; tail of log:"
        tail -n 30 "${sample_log}" 2>/dev/null | sed 's/^/      /' | tee -a "${SUITE_LOG}"
    else
        log "OK(${gsm}): ${OUT_ROOT}/${out_id}"
    fi
    return $rc
}

# -----------------------------------------------------------------------------
# 6. Concurrency runner
# -----------------------------------------------------------------------------
run_all() {
    log "HOSTNAME=${HOSTNAME}  TMPDIR=${TMPDIR}  OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
    log "CONCURRENCY=${CONCURRENCY}  THREADS/sample=${THREADS}  MEM_GB/sample=${MEM_GB}"
    log "TRANSCRIPTOME=${TRANSCRIPTOME}"
    log "samples=(${SAMPLES[*]})"
    echo "----------------------------------------" | tee -a "${SUITE_LOG}"

    declare -A pid2gsm
    running=0; i=0
    while (( i < ${#SAMPLES[@]} )) || (( running > 0 )); do
        while (( i < ${#SAMPLES[@]} )) && (( running < CONCURRENCY )); do
            local gsm="${SAMPLES[$i]}"
            i=$((i + 1))
            ( run_count "$gsm" ) &
            pid2gsm[$!]=$gsm
            running=$((running + 1))
            log "launched ${gsm} (pid=$!)"
        done
        if (( running > 0 )); then
            wait -n
            for pid in "${!pid2gsm[@]}"; do
                if ! kill -0 "$pid" 2>/dev/null; then
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
        local out_id="gse206391_${gsm}_sr"
        local out="${OUT_ROOT}/${out_id}"
        local state="MISSING"
        if [[ -f "${out}/outs/web_summary.html" ]]; then state="COMPLETE"
        elif [[ -f "${out}/_log" || -d "${out}" ]]; then state="PARTIAL/RUNNING"; fi
        log "${gsm}\t${state}\t${out}"
    done
}

usage() {
    cat <<EOF
Usage:
  $0 run-all      Launch both samples (concurrency=${CONCURRENCY}).
  $0 run-one GSM  Single sample (foreground).
  $0 status       Print per-sample output status.
  $0 prepare GSM  Build fastqs/ + extract image (no spaceranger).

Env: CONCURRENCY=${CONCURRENCY}  THREADS=${THREADS}  MEM_GB=${MEM_GB}
EOF
}

MODE="${1:-help}"
case "$MODE" in
    run-all)  run_all ;;
    run-one)  run_count "${2:-}" ;;
    status)   status ;;
    prepare)  gsm="${2:-}"; prepare_links "$gsm"; prepare_image "$gsm" ;;
    help|-h|--help) usage ;;
    *) echo "Unknown mode: $MODE" >&2; usage; exit 2 ;;
esac

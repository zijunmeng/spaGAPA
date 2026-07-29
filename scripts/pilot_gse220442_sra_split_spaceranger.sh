#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
DATA_ROOT="${ROOT}/data/raw/gse220442"
PACKAGE_ROOT="${ROOT}/spaGAPA"
LOG_ROOT="${ROOT}/logs"
SPACERANGER="${SPACERANGER:-/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger}"
TRANSCRIPTOME="${TRANSCRIPTOME:-/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A}"
PREFETCH="${PREFETCH:-prefetch}"
FASTERQ_DUMP="${FASTERQ_DUMP:-fasterq-dump}"
VDB_VALIDATE="${VDB_VALIDATE:-vdb-validate}"

MODE="${1:-help}"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"

# Pilot sample: smallest GSE220442 sample by downloaded ENA FASTQ size.
# The inferred Visium directory is based on sample order and metadata patientID.
GSM="${GSM:-GSM6801755}"
VISIUM_ID="${VISIUM_ID:-2-8}"
RUNS="${RUNS:-SRR22561784 SRR22561785}"
THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-48}"
PREFETCH_MAX_SIZE="${PREFETCH_MAX_SIZE:-100G}"

SRA_DIR="${DATA_ROOT}/sra"
SPLIT_DIR="${DATA_ROOT}/fastq_sra_split/${GSM}"
INPUT_DIR="${DATA_ROOT}/spaceranger_inputs/${GSM}"
FASTQ_LINK_DIR="${INPUT_DIR}/fastqs_sra_split"
IMAGE_OUT="${INPUT_DIR}/detected_tissue_image.jpg"
SUPP_TAR="${DATA_ROOT}/geo_supplementary/GSE220442_counts_and_images.tar.gz"
OUT_ROOT="${PACKAGE_ROOT}/pipeline_output/gse220442_sra_split_pilot"
WORK_DIR="${OUT_ROOT}/${GSM}"
DRYRUN_LOG="${WORK_DIR}/logs/spaceranger_dry_run.log"
AUDIT_TSV="${WORK_DIR}/logs/sra_split_fastq_audit.tsv"
NOHUP_LOG="${LOG_ROOT}/${DATE_TAG}_gse220442_sra_split_pilot.nohup.log"
PID_FILE="${LOG_ROOT}/${DATE_TAG}_gse220442_sra_split_pilot.pid"

mkdir -p "${LOG_ROOT}" "${SRA_DIR}" "${SPLIT_DIR}" "${FASTQ_LINK_DIR}" "${WORK_DIR}/logs"

usage() {
  cat <<EOF
Usage:
  $0 check
  $0 foreground
  $0 nohup
  $0 monitor

Purpose:
  Pilot whether GSE220442 SRA can recover Space Ranger-compatible R1/R2.

Default pilot:
  GSM=${GSM}
  RUNS="${RUNS}"
  VISIUM_ID=${VISIUM_ID}

Key behavior:
  1. prefetch each SRR into ${SRA_DIR}
  2. vdb-validate the .sra files
  3. fasterq-dump --split-files --include-technical
  4. keep uncompressed split FASTQ to avoid heavy recompression
  5. link read 3 as Space Ranger R1 and read 4 as Space Ranger R2
  6. run Space Ranger --dry only

Override example:
  GSM=GSM6801754 VISIUM_ID=2-3 RUNS="SRR22561786 SRR22561787" $0 nohup

Monitor:
  $0 monitor
  tail -f ${NOHUP_LOG}
  tail -f ${DRYRUN_LOG}

Outputs:
  Split FASTQ: ${SPLIT_DIR}
  Space Ranger input links: ${FASTQ_LINK_DIR}
  Dry-run log: ${DRYRUN_LOG}
  Audit TSV: ${AUDIT_TSV}
EOF
}

check_tools() {
  command -v "${PREFETCH}" >/dev/null
  command -v "${FASTERQ_DUMP}" >/dev/null
  command -v "${VDB_VALIDATE}" >/dev/null
  [[ -x "${SPACERANGER}" ]]
  [[ -d "${TRANSCRIPTOME}" ]]
  [[ -f "${SUPP_TAR}" ]]
}

print_plan() {
  echo "GSM=${GSM}"
  echo "RUNS=${RUNS}"
  echo "VISIUM_ID=${VISIUM_ID}"
  echo "THREADS=${THREADS}"
  echo "MEM_GB=${MEM_GB}"
  echo "PREFETCH_MAX_SIZE=${PREFETCH_MAX_SIZE}"
  echo "SRA_DIR=${SRA_DIR}"
  echo "SPLIT_DIR=${SPLIT_DIR}"
  echo "FASTQ_LINK_DIR=${FASTQ_LINK_DIR}"
  echo "IMAGE_OUT=${IMAGE_OUT}"
  echo "OUT_ROOT=${OUT_ROOT}"
  echo "DRYRUN_LOG=${DRYRUN_LOG}"
  echo "AUDIT_TSV=${AUDIT_TSV}"
}

prefetch_and_split() {
  # shellcheck disable=SC2206
  local run_array=(${RUNS})
  for run in "${run_array[@]}"; do
    local sra_path="${SRA_DIR}/${run}/${run}.sra"
    local lock_path="${sra_path}.lock"
    if [[ -e "${lock_path}" ]]; then
      echo "[$(date '+%F %T')] moving stale lock ${lock_path}"
      mv "${lock_path}" "${lock_path}.stale.$(date +%s)"
    fi
    if [[ ! -s "${sra_path}" ]]; then
      echo "[$(date '+%F %T')] prefetch ${run}"
      "${PREFETCH}" --max-size "${PREFETCH_MAX_SIZE}" "${run}" -O "${SRA_DIR}"
    else
      echo "[$(date '+%F %T')] existing SRA ${sra_path}"
    fi

    echo "[$(date '+%F %T')] vdb-validate ${run}"
    "${VDB_VALIDATE}" "${sra_path}"

    local read1="${SPLIT_DIR}/${run}_1.fastq.gz"
    local read4="${SPLIT_DIR}/${run}_4.fastq.gz"
    local read1_raw="${SPLIT_DIR}/${run}_1.fastq"
    local read4_raw="${SPLIT_DIR}/${run}_4.fastq"
    if [[ ! -s "${read1}" && ! -s "${read1_raw}" || ! -s "${read4}" && ! -s "${read4_raw}" ]]; then
      if [[ -s "${read1_raw}" && -s "${read4_raw}" ]]; then
        echo "[$(date '+%F %T')] existing uncompressed split FASTQ for ${run}; keep as-is"
      else
        echo "[$(date '+%F %T')] fasterq-dump --split-files --include-technical ${run}"
        "${FASTERQ_DUMP}" \
          --split-files \
          --include-technical \
          --threads "${THREADS}" \
          --outdir "${SPLIT_DIR}" \
          --temp "${SPLIT_DIR}/tmp_${run}" \
          "${sra_path}"
      fi
    else
      echo "[$(date '+%F %T')] existing split FASTQ for ${run}"
    fi
  done
}

first_header() {
  python - "$1" <<'PY'
import gzip
import sys

opener = gzip.open if sys.argv[1].endswith(".gz") else open
with opener(sys.argv[1], "rt", encoding="utf-8", errors="replace") as handle:
    print(handle.readline().rstrip("\n"))
PY
}

first_length() {
  python - "$1" <<'PY'
import gzip
import sys

opener = gzip.open if sys.argv[1].endswith(".gz") else open
with opener(sys.argv[1], "rt", encoding="utf-8", errors="replace") as handle:
    handle.readline()
    print(len(handle.readline().rstrip("\n")))
PY
}

audit_split_fastq() {
  echo -e "run\tread_index\tfile\tfirst_seq_len\tfirst_header" > "${AUDIT_TSV}"
  # shellcheck disable=SC2206
  local run_array=(${RUNS})
  for run in "${run_array[@]}"; do
    for read_idx in 1 2 3 4; do
      local fq len header
      fq="$(select_split_fastq "${run}" "${read_idx}")" || continue
      len="$(first_length "${fq}")"
      header="$(first_header "${fq}")"
      echo -e "${run}\t${read_idx}\t${fq}\t${len}\t${header}" >> "${AUDIT_TSV}"
    done
  done
  cat "${AUDIT_TSV}"
}

select_split_fastq() {
  local run="$1"
  local read_index="$2"
  local raw="${SPLIT_DIR}/${run}_${read_index}.fastq"
  local gz="${SPLIT_DIR}/${run}_${read_index}.fastq.gz"
  if [[ -s "${raw}" ]]; then
    echo "${raw}"
  elif [[ -s "${gz}" ]]; then
    echo "${gz}"
  else
    return 1
  fi
}

extract_image() {
  if [[ -s "${IMAGE_OUT}" ]]; then
    echo "[$(date '+%F %T')] existing image ${IMAGE_OUT}"
    return
  fi

  local gz_member="counts_and_images/${VISIUM_ID}/spatial/detected_tissue_image.jpg.gz"
  local plain_member="counts_and_images/${VISIUM_ID}/spatial/detected_tissue_image.jpg"
  if tar -tzf "${SUPP_TAR}" "${gz_member}" >/dev/null 2>&1; then
    tar -xOzf "${SUPP_TAR}" "${gz_member}" | gzip -cd > "${IMAGE_OUT}"
  elif tar -tzf "${SUPP_TAR}" "${plain_member}" >/dev/null 2>&1; then
    tar -xOzf "${SUPP_TAR}" "${plain_member}" > "${IMAGE_OUT}"
  else
    echo "No detected_tissue_image found for VISIUM_ID=${VISIUM_ID}" >&2
    exit 1
  fi
}

prepare_spaceranger_links() {
  rm -f "${FASTQ_LINK_DIR}"/*.fastq.gz
  # shellcheck disable=SC2206
  local run_array=(${RUNS})
  local lane=1
  for run in "${run_array[@]}"; do
    local read1 read4 ext1 ext2
    read1="$(select_split_fastq "${run}" 3)" || {
      echo "Missing ${SPLIT_DIR}/${run}_3.fastq(.gz)" >&2
      exit 1
    }
    read4="$(select_split_fastq "${run}" 4)" || {
      echo "Missing ${SPLIT_DIR}/${run}_4.fastq(.gz)" >&2
      exit 1
    }
    ext1="fastq"
    ext2="fastq"
    [[ "${read1}" == *.gz ]] && ext1="fastq.gz"
    [[ "${read4}" == *.gz ]] && ext2="fastq.gz"
    ln -s "${read1}" "${FASTQ_LINK_DIR}/${GSM}_S1_L$(printf '%03d' "${lane}")_R1_001.${ext1}"
    ln -s "${read4}" "${FASTQ_LINK_DIR}/${GSM}_S1_L$(printf '%03d' "${lane}")_R2_001.${ext2}"
    lane=$((lane + 1))
  done
}

run_spaceranger_dry() {
  local dry_id="gse220442_${GSM}_sra_split_dry"
  rm -rf "${WORK_DIR}/work/${dry_id}"
  mkdir -p "${WORK_DIR}/work"
  {
    echo "[$(date '+%F %T')] Running Space Ranger dry-run"
    echo "GSM=${GSM}"
    echo "RUNS=${RUNS}"
    echo "VISIUM_ID=${VISIUM_ID}"
    echo "FASTQ_LINK_DIR=${FASTQ_LINK_DIR}"
    echo "IMAGE_OUT=${IMAGE_OUT}"
    "${SPACERANGER}" count \
      --id "${dry_id}" \
      --description "GSE220442 ${GSM} SRA split dry-run" \
      --transcriptome "${TRANSCRIPTOME}" \
      --fastqs "${FASTQ_LINK_DIR}" \
      --sample "${GSM}" \
      --image "${IMAGE_OUT}" \
      --unknown-slide visium-1 \
      --create-bam true \
      --localcores "${THREADS}" \
      --localmem "${MEM_GB}" \
      --disable-cell-annotation \
      --disable-ui \
      --dry
  } > "${DRYRUN_LOG}" 2>&1
}

run_all() {
  check_tools
  print_plan
  prefetch_and_split
  audit_split_fastq
  extract_image
  prepare_spaceranger_links
  run_spaceranger_dry
  echo "[$(date '+%F %T')] Done"
  echo "Dry-run log: ${DRYRUN_LOG}"
  echo "Audit TSV: ${AUDIT_TSV}"
}

case "${MODE}" in
  check)
    check_tools
    print_plan
    ;;

  foreground)
    run_all
    ;;

  nohup|start)
    if [[ -s "${PID_FILE}" ]]; then
      old_pid="$(cat "${PID_FILE}")"
      if ps -p "${old_pid}" >/dev/null 2>&1; then
        echo "Existing pilot is running: PID=${old_pid}"
        echo "Use: $0 monitor"
        exit 1
      fi
    fi
    nohup "${SCRIPT_PATH}" foreground > "${NOHUP_LOG}" 2>&1 &
    pid="$!"
    echo "${pid}" > "${PID_FILE}"
    echo "Submitted GSE220442 SRA split pilot PID: ${pid}"
    echo "PID file: ${PID_FILE}"
    echo "Log: ${NOHUP_LOG}"
    echo "Monitor: ${SCRIPT_PATH} monitor"
    ;;

  monitor)
    if [[ -s "${PID_FILE}" ]]; then
      pid="$(cat "${PID_FILE}")"
      ps -p "${pid}" -o pid,ppid,stat,etime,%mem,%cpu,cmd || true
    else
      echo "No PID file: ${PID_FILE}"
    fi
    echo
    echo "NOHUP tail:"
    [[ -f "${NOHUP_LOG}" ]] && tail -n 80 "${NOHUP_LOG}" || true
    echo
    echo "Split FASTQ audit:"
    [[ -f "${AUDIT_TSV}" ]] && cat "${AUDIT_TSV}" || true
    echo
    echo "Space Ranger dry-run tail:"
    [[ -f "${DRYRUN_LOG}" ]] && tail -n 80 "${DRYRUN_LOG}" || true
    ;;

  help|-h|--help)
    usage
    ;;

  *)
    echo "Unknown mode: ${MODE}" >&2
    usage
    exit 2
    ;;
esac

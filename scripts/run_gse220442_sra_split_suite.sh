#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
PILOT="${ROOT}/spaGAPA/scripts/pilot_gse220442_sra_split_spaceranger.sh"
DATA_ROOT="${ROOT}/data/raw/gse220442"
OUT_ROOT="${ROOT}/spaGAPA/pipeline_output/gse220442_sra_split_pilot"
LOG_ROOT="${ROOT}/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse220442_sra_split_suite.log"
SESSION="${SESSION:-spagapa_gse220442_sra_split_all}"
MODE="${1:-help}"
THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-48}"
PREFETCH_MAX_SIZE="${PREFETCH_MAX_SIZE:-100G}"

mkdir -p "${LOG_ROOT}"

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
declare -A RUNS

TITLE[GSM6801751]="control 1, biol rep 1"
TITLE[GSM6801752]="control 2, biol rep 2"
TITLE[GSM6801753]="control 3, biol rep 3"
TITLE[GSM6801754]="AD 1, biol rep 1"
TITLE[GSM6801755]="AD 2, biol rep 2"
TITLE[GSM6801756]="AD 3, biol rep 3"

VISIUM_ID[GSM6801751]="1-1"
VISIUM_ID[GSM6801752]="2-5"
VISIUM_ID[GSM6801753]="18-64"
VISIUM_ID[GSM6801754]="2-3"
VISIUM_ID[GSM6801755]="2-8"
VISIUM_ID[GSM6801756]="T4857"

RUNS[GSM6801751]="SRR22561792 SRR22561793"
RUNS[GSM6801752]="SRR22561790 SRR22561791"
RUNS[GSM6801753]="SRR22561788 SRR22561789"
RUNS[GSM6801754]="SRR22561786 SRR22561787"
RUNS[GSM6801755]="SRR22561784 SRR22561785"
RUNS[GSM6801756]="SRR22561782 SRR22561783"

usage() {
  cat <<EOF
Usage:
  $0 check
  $0 run-one <GSM>
  $0 run-all
  $0 tmux-submit
  $0 monitor

Purpose:
  Reconstruct Space Ranger-compatible GSE220442 FASTQ from SRA using:
    fasterq-dump --split-files --include-technical
  and link technical read 3 as Space Ranger R1, read 4 as Space Ranger R2.

Default resources per sample:
  THREADS=${THREADS}
  MEM_GB=${MEM_GB}
  PREFETCH_MAX_SIZE=${PREFETCH_MAX_SIZE}

Recommended background full suite:
  cd ${ROOT}
  THREADS=8 MEM_GB=48 ${0} tmux-submit

Monitor:
  ${0} monitor
  tail -f ${SUITE_LOG}

Pilot output root:
  ${OUT_ROOT}/<GSM>/
EOF
}

sample_exists() {
  local gsm="$1"
  [[ -n "${RUNS[${gsm}]:-}" ]]
}

status_for_sample() {
  local gsm="$1"
  local split_dir="${DATA_ROOT}/fastq_sra_split/${gsm}"
  local link_dir="${DATA_ROOT}/spaceranger_inputs/${gsm}/fastqs_sra_split"
  local dry_log="${OUT_ROOT}/${gsm}/logs/spaceranger_dry_run.log"
  local read3_ok=0
  local read4_ok=0
  local link_count=0
  # shellcheck disable=SC2206
  local run_array=(${RUNS[${gsm}]})
  for run in "${run_array[@]}"; do
    [[ -s "${split_dir}/${run}_3.fastq" || -s "${split_dir}/${run}_3.fastq.gz" ]] && read3_ok=$((read3_ok + 1))
    [[ -s "${split_dir}/${run}_4.fastq" || -s "${split_dir}/${run}_4.fastq.gz" ]] && read4_ok=$((read4_ok + 1))
  done
  if [[ -d "${link_dir}" ]]; then
    link_count="$(find "${link_dir}" -maxdepth 1 -type l | wc -l)"
  fi
  printf "%s\t%s\t%s\t%s\tread3=%s/%s\tread4=%s/%s\tlinks=%s\tdry_log=%s\n" \
    "${gsm}" "${TITLE[${gsm}]}" "${VISIUM_ID[${gsm}]}" "${RUNS[${gsm}]}" \
    "${read3_ok}" "${#run_array[@]}" "${read4_ok}" "${#run_array[@]}" \
    "${link_count}" "$([[ -s "${dry_log}" ]] && echo yes || echo no)"
}

check_suite() {
  echo -e "gsm\ttitle\tvisium_id\truns\tread3_status\tread4_status\tlinks\tdry_log"
  for gsm in "${SAMPLES[@]}"; do
    status_for_sample "${gsm}"
  done
}

run_one() {
  local gsm="$1"
  if ! sample_exists "${gsm}"; then
    echo "Unknown GSM: ${gsm}" >&2
    exit 2
  fi
  echo "[$(date '+%F %T')] Start ${gsm} ${TITLE[${gsm}]}"
  GSM="${gsm}" \
  VISIUM_ID="${VISIUM_ID[${gsm}]}" \
  RUNS="${RUNS[${gsm}]}" \
  THREADS="${THREADS}" \
  MEM_GB="${MEM_GB}" \
  PREFETCH_MAX_SIZE="${PREFETCH_MAX_SIZE}" \
    "${PILOT}" foreground
  echo "[$(date '+%F %T')] Finish ${gsm}"
}

run_all() {
  check_suite
  for gsm in "${SAMPLES[@]}"; do
    run_one "${gsm}"
    check_suite
  done
}

case "${MODE}" in
  check)
    check_suite
    ;;

  run-one)
    run_one "${2:-}"
    ;;

  run-all)
    run_all
    ;;

  tmux-submit)
    if tmux has-session -t "${SESSION}" >/dev/null 2>&1; then
      echo "tmux session already exists: ${SESSION}" >&2
      exit 1
    fi
    tmux new-session -d -s "${SESSION}" -c "${ROOT}" \
      "THREADS=${THREADS} MEM_GB=${MEM_GB} PREFETCH_MAX_SIZE=${PREFETCH_MAX_SIZE} ${0} run-all > ${SUITE_LOG} 2>&1"
    echo "Submitted tmux session: ${SESSION}"
    echo "Suite log: ${SUITE_LOG}"
    ;;

  monitor)
    tmux list-sessions 2>/dev/null | grep "${SESSION}" || true
    echo
    check_suite
    echo
    [[ -f "${SUITE_LOG}" ]] && tail -n 120 "${SUITE_LOG}" || true
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

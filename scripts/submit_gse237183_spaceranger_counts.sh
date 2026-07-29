#!/usr/bin/env bash
set -euo pipefail

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA"
RUNNER="${ROOT}/spaGAPA/scripts/run_gse237183_spaceranger_counts.py"
OUT_ROOT="${ROOT}/spaGAPA/pipeline_output/gse237183_spaceranger_counts"
LOG_ROOT="${ROOT}/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
NOHUP_LOG="${LOG_ROOT}/${DATE_TAG}_gse237183_spaceranger_counts.nohup.log"
PID_FILE="${LOG_ROOT}/${DATE_TAG}_gse237183_spaceranger_counts.pid"

MODE="${1:-help}"
JOBS="${JOBS:-1}"
THREADS="${THREADS:-24}"
MEM_GB="${MEM_GB:-128}"
RESUME_EXISTING="${RESUME_EXISTING:-1}"
FORCE="${FORCE:-0}"
PYTHON_BIN="${PYTHON_BIN:-python}"

mkdir -p "${LOG_ROOT}" "${OUT_ROOT}"

build_cmd() {
  CMD=(
    "${PYTHON_BIN}" "${RUNNER}"
    --jobs "${JOBS}"
    --threads "${THREADS}"
    --mem-gb "${MEM_GB}"
  )

  if [[ "${RESUME_EXISTING}" == "1" ]]; then
    CMD+=(--resume-existing)
  fi
  if [[ "${FORCE}" == "1" ]]; then
    CMD+=(--force)
  fi
  if [[ -n "${SAMPLES:-}" ]]; then
    # shellcheck disable=SC2206
    SAMPLE_ARRAY=(${SAMPLES})
    CMD+=(--samples "${SAMPLE_ARRAY[@]}")
  fi
}

print_usage() {
  cat <<EOF
Usage:
  $0 check
  $0 nohup
  $0 foreground
  $0 monitor

Default resources:
  JOBS=${JOBS}
  THREADS=${THREADS}
  MEM_GB=${MEM_GB}
  RESUME_EXISTING=${RESUME_EXISTING}
  FORCE=${FORCE}

Recommended full run:
  cd ${ROOT}
  JOBS=1 THREADS=24 MEM_GB=128 ${0} nohup

Optional two-sample parallel run, only if CPU/I/O/memory are sufficient:
  cd ${ROOT}
  JOBS=2 THREADS=16 MEM_GB=96 ${0} nohup

Optional subset:
  SAMPLES="GSM7596587 GSM7596588" ${0} nohup

Monitor:
  ${0} monitor
  tail -f ${NOHUP_LOG}
  tail -f ${OUT_ROOT}/count_summary.tsv

Output:
  Space Ranger pipestances: ${OUT_ROOT}/gse237183_<GSM>_sr
  Per-sample logs:          ${OUT_ROOT}/_logs/<GSM>/spaceranger_count.log
  Suite summary:            ${OUT_ROOT}/count_summary.tsv
EOF
}

check_existing_job() {
  if [[ -s "${PID_FILE}" ]]; then
    PID="$(cat "${PID_FILE}")"
    if ps -p "${PID}" >/dev/null 2>&1; then
      echo "Existing job appears to be running: PID=${PID}"
      echo "PID file: ${PID_FILE}"
      return 0
    fi
  fi
  return 1
}

case "${MODE}" in
  check)
    "${PYTHON_BIN}" "${RUNNER}" --check-only
    ;;

  nohup|start)
    build_cmd
    if check_existing_job; then
      echo "Refusing to submit a second job. Use monitor, or remove stale PID file if needed."
      exit 1
    fi
    {
      echo "[$(date '+%F %T')] Starting GSE237183 Space Ranger count suite"
      echo "ROOT=${ROOT}"
      echo "OUT_ROOT=${OUT_ROOT}"
      echo "JOBS=${JOBS}"
      echo "THREADS=${THREADS}"
      echo "MEM_GB=${MEM_GB}"
      echo "RESUME_EXISTING=${RESUME_EXISTING}"
      echo "FORCE=${FORCE}"
      echo "COMMAND=${CMD[*]}"
    } > "${NOHUP_LOG}"
    nohup "${CMD[@]}" >> "${NOHUP_LOG}" 2>&1 &
    PID="$!"
    echo "${PID}" > "${PID_FILE}"
    echo "Submitted PID: ${PID}"
    echo "PID file: ${PID_FILE}"
    echo "Main log: ${NOHUP_LOG}"
    echo "Summary: ${OUT_ROOT}/count_summary.tsv"
    ;;

  foreground)
    build_cmd
    echo "Running foreground command:"
    printf ' %q' "${CMD[@]}"
    echo
    exec "${CMD[@]}"
    ;;

  monitor)
    if [[ -s "${PID_FILE}" ]]; then
      PID="$(cat "${PID_FILE}")"
      echo "PID file: ${PID_FILE}"
      ps -p "${PID}" -o pid,ppid,stat,etime,%mem,%cpu,cmd || true
    else
      echo "No PID file found: ${PID_FILE}"
    fi
    echo
    echo "Completed BAM count:"
    find "${OUT_ROOT}" -path "*/outs/possorted_genome_bam.bam" -type f 2>/dev/null | wc -l
    echo
    if [[ -f "${OUT_ROOT}/count_summary.tsv" ]]; then
      echo "Current summary:"
      tail -n 30 "${OUT_ROOT}/count_summary.tsv"
      echo
    fi
    if [[ -f "${NOHUP_LOG}" ]]; then
      echo "Main log tail:"
      tail -n 80 "${NOHUP_LOG}"
    else
      echo "No main log found: ${NOHUP_LOG}"
    fi
    ;;

  help|-h|--help)
    print_usage
    ;;

  *)
    echo "Unknown mode: ${MODE}" >&2
    print_usage
    exit 2
    ;;
esac

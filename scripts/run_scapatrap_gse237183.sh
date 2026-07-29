#!/usr/bin/env bash
# Run scAPAtrap on GSE237183 samples sequentially.
#
# Modes:
#   ./run_scapatrap_gse237183.sh pilot          # GSM7596590 then GSM7596601 (2 representative samples)
#   ./run_scapatrap_gse237183.sh batch          # all 18 valid samples
#   ./run_scapatrap_gse237183.sh gsm GSM7596590 # one specific sample
#   ./run_scapatrap_gse237183.sh dryrun         # dry-validate pilot samples (no R run)
#
# Designed to run inside a persistent tmux session. Logs to:
#   spaGAPA/pipeline_output/gse237183_scapatrap_pilot/logs/  (pilot)
#   spaGAPA/pipeline_output/gse237183_scapatrap_batch/logs/  (batch)
set -euo pipefail

# ---------------------------------------------------------------------------
# Server detection per CLAUDE.md (S90/S91/S97/S98)
# ---------------------------------------------------------------------------
HOSTNAME_SHORT=$(hostname -s | tr 'a-z' 'A-Z')
case "$HOSTNAME_SHORT" in
    S90)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s2/mengzijun/tmp"
        PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
        ;;
    S91)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="$HOME/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="$HOME/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s3/mengzijun/tmp"
        PYTHON_BIN="$HOME/anaconda3/envs/spagapa/bin/python"
        ;;
    S97)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s972/mengzijun/tmp"
        PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
        ;;
    S98)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
        export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
        export TMPDIR="/s982/mengzijun/tmp"
        PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
        ;;
    *)
        echo "错误：未知服务器 $HOSTNAME_SHORT，请检查环境配置 (expected S90/S91/S97/S98)." >&2
        exit 1
        ;;
esac
export TMPDIR
mkdir -p "$TMPDIR"

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DRIVER="$ROOT/scripts/run_scapatrap_gse237183.py"

MODE="${1:-pilot}"
THREADS="${THREADS:-8}"

# ---------------------------------------------------------------------------
# Pick log dir + sample set based on mode
# ---------------------------------------------------------------------------
case "$MODE" in
    pilot)
        LOG_BASE="$ROOT/pipeline_output/gse237183_scapatrap_pilot/logs"
        GSMS="GSM7596590 GSM7596601"
        DRIVER_ARGS="--gsm GSM7596590 GSM7596601"
        ;;
    batch)
        LOG_BASE="$ROOT/pipeline_output/gse237183_scapatrap_batch/logs"
        GSMS="all 18 valid samples"
        DRIVER_ARGS="--all"
        ;;
    dryrun)
        LOG_BASE="$ROOT/pipeline_output/gse237183_scapatrap_pilot/logs"
        GSMS="GSM7596590 GSM7596601 (dry-run validation)"
        DRIVER_ARGS="--gsm GSM7596590 GSM7596601 --dry-run"
        ;;
    gsm)
        LOG_BASE="$ROOT/pipeline_output/gse237183_scapatrap_pilot/logs"
        GSM="${2:?usage: $0 gsm GSM7596590}"
        GSMS="$GSM"
        DRIVER_ARGS="--gsm $GSM"
        ;;
    *)
        echo "Unknown mode: $MODE" >&2
        echo "Usage: $0 {pilot|batch|dryrun|gsm <GSM>}" >&2
        exit 1
        ;;
esac

mkdir -p "$LOG_BASE"
TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="$LOG_BASE/scapatrap_${MODE}_${TS}.log"

echo "########################################################" | tee "$OUT_LOG"
echo "# scAPAtrap GSE237183 run" | tee -a "$OUT_LOG"
echo "# host    : $HOSTNAME_SHORT" | tee -a "$OUT_LOG"
echo "# mode    : $MODE" | tee -a "$OUT_LOG"
echo "# samples : $GSMS" | tee -a "$OUT_LOG"
echo "# threads : $THREADS" | tee -a "$OUT_LOG"
echo "# python  : $PYTHON_BIN" | tee -a "$OUT_LOG"
echo "# Rscript : $RSCRIPT" | tee -a "$OUT_LOG"
echo "# R_LIBS  : $R_LIBS" | tee -a "$OUT_LOG"
echo "# TMPDIR  : $TMPDIR" | tee -a "$OUT_LOG"
echo "# log     : $OUT_LOG" | tee -a "$OUT_LOG"
echo "# started : $(date)" | tee -a "$OUT_LOG"
echo "########################################################" | tee -a "$OUT_LOG"

# All subsequent output (python + R) appended to $OUT_LOG.
# stdout/stderr of python's subprocess (the R run) flows through.
exec "$PYTHON_BIN" "$DRIVER" $DRIVER_ARGS --threads "$THREADS" --force >> "$OUT_LOG" 2>&1

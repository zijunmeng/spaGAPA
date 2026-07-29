#!/usr/bin/env bash
# Run scAPAtrap on all 6 GSE220442 Visium (human AD brain) Space Ranger BAMs.
# 3 control (GSM6801751-3) + 3 AD (GSM6801754-6).
# Hostname-aware env (CLAUDE.md), per-sample timing, concurrency-limited.
set -uo pipefail

# ---- hostname-aware environment (CLAUDE.md) ----
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s2/mengzijun/tmp" ;;
    S91) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/home/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/home/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s3/mengzijun/tmp" ;;
    S97) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s972/mengzijun/tmp" ;;
    S98) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; PY="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"; export TMPDIR="/s982/mengzijun/tmp" ;;
    *) echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

B=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOGROOT="$B/pipeline_output/gse220442_scapatrap_batch_logs"
mkdir -p "$LOGROOT"
SUMMARY="$LOGROOT/batch_summary.tsv"
PROGRESS="$LOGROOT/batch_progress.log"
echo -e "sample\tcondition\texit\twall_s" > "$SUMMARY"

# sample -> condition mapping
run_sample () {
    local s="$1" cond="$2"
    local t0 t1 rc
    local gsm_lc=$(echo "$s" | tr 'A-Z' 'a-z')
    local bam="$B/pipeline_output/gse220442_${s}_sr/outs/possorted_genome_bam.bam"
    local spatial="$B/pipeline_output/gse220442_${s}_sr/outs/spatial"
    local outroot="$B/pipeline_output/gse220442_${s}_scapatrap"
    local processed="$B/data/processed/gse220442_${gsm_lc}_scapatrap"
    local slabel="GSE220442 ${s} Visium human AD brain PFC (${cond}) + scAPAtrap"
    local tissue="human AD brain PFC (${cond})"
    local slog="$LOGROOT/${s}.log"

    echo "[$(date '+%F %T')] START $s ($cond) -> $slog" | tee -a "$PROGRESS"
    t0=$(date +%s)
    "$PY" "$B/scripts/run_scapatrap_spaceranger.py" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --output-root "$outroot" \
        --processed-dir "$processed" \
        --dataset-name "gse220442_${gsm_lc}_scapatrap" \
        --source-label "$slabel" \
        --tissue "$tissue" \
        > "$slog" 2>&1
    rc=$?
    t1=$(date +%s)
    local wall=$((t1 - t0))
    echo -e "${s}\t${cond}\t${rc}\t${wall}" >> "$SUMMARY"
    echo "[$(date '+%F %T')] DONE  $s exit=$rc wall=${wall}s" | tee -a "$PROGRESS"
}

# ---- concurrency-limited job pool ----
CONCURRENCY=3
declare -a CONDS=(control control control AD AD AD)
SAMPLES=(GSM6801751 GSM6801752 GSM6801753 GSM6801754 GSM6801755 GSM6801756)
pids=()
for i in "${!SAMPLES[@]}"; do
    s="${SAMPLES[$i]}"; c="${CONDS[$i]}"
    run_sample "$s" "$c" &
    pids+=($!)
    # throttle to CONCURRENCY
    while [ "$(jobs -r | wc -l)" -ge "$CONCURRENCY" ]; do
        wait -n 2>/dev/null || sleep 5
    done
done
wait

echo "[$(date '+%F %T')] ALL SAMPLES COMPLETE" | tee -a "$PROGRESS"
echo "===== BATCH SUMMARY =====" | tee -a "$PROGRESS"
column -t -s $'\t' "$SUMMARY" | tee -a "$PROGRESS"

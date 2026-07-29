#!/usr/bin/env bash
# Run scAPAtrap on 2 GSE206391 (human skin, psoriasis lesional) Space Ranger BAMs.
# Uses scripts/run_scapatrap_spaceranger.py (HUMAN standard launcher, no modification).
#   GSM6252925 (Patient 4 lesional psoriasis rep1)
#   GSM6252926 (Patient 4 lesional psoriasis rep2)
# Hostname-aware env (CLAUDE.md), per-sample timing, concurrency-limited.
# --tissue "human skin (psoriasis)" --species human. Do NOT git commit.
set -uo pipefail

# ---- hostname-aware environment (CLAUDE.md) ----
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91)
        export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
        RSCRIPT="/home/mengzijun/anaconda3/envs/r442/bin/Rscript"
        PY="/home/mengzijun/anaconda3/envs/spagapa/bin/python"
        export TMPDIR="/s3/mengzijun/tmp"
        export OPENBLAS_NUM_THREADS=8
        ;;
    *)
        echo "ERROR: unknown host $HOSTNAME" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

B=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOGROOT="$B/pipeline_output/gse206391_scapatrap_batch_logs"
mkdir -p "$LOGROOT"
SUMMARY="$LOGROOT/batch_summary.tsv"
PROGRESS="$LOGROOT/batch_progress.log"
echo -e "sample\tcondition\texit\twall_s" > "$SUMMARY"

# sample -> condition mapping (both are psoriasis lesional)
declare -A TITLE
TITLE[GSM6252925]="P4 lesional psoriasis rep1"
TITLE[GSM6252926]="P4 lesional psoriasis rep2"

run_sample () {
    local s="$1"
    local title="${TITLE[$s]}"
    local t0 t1 rc
    local gsm_lc=$(echo "$s" | tr 'A-Z' 'a-z')
    local bam="$B/pipeline_output/gse206391_${s}_sr/outs/possorted_genome_bam.bam"
    local spatial="$B/pipeline_output/gse206391_${s}_sr/outs/spatial"
    local outroot="$B/pipeline_output/gse206391_${s}_scapatrap"
    local processed="$B/data/processed/gse206391_${gsm_lc}_scapatrap"
    local slabel="GSE206391 ${s} Visium ${title} + scAPAtrap"
    local tissue="human skin (psoriasis)"
    local slog="$LOGROOT/${s}.log"

    echo "[$(date '+%F %T')] START $s ($title) -> $slog" | tee -a "$PROGRESS"
    t0=$(date +%s)
    "$PY" "$B/scripts/run_scapatrap_spaceranger.py" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --output-root "$outroot" \
        --processed-dir "$processed" \
        --dataset-name "gse206391_${gsm_lc}_scapatrap" \
        --source-label "$slabel" \
        --tissue "$tissue" \
        --species human \
        --threads 12 \
        > "$slog" 2>&1
    rc=$?
    t1=$(date +%s)
    local wall=$((t1 - t0))
    echo -e "${s}\t${title}\t${rc}\t${wall}" >> "$SUMMARY"
    echo "[$(date '+%F %T')] DONE  $s exit=$rc wall=${wall}s" | tee -a "$PROGRESS"
}

CONCURRENCY=2
SAMPLES=(GSM6252925 GSM6252926)
pids=()
for i in "${!SAMPLES[@]}"; do
    s="${SAMPLES[$i]}"
    run_sample "$s" &
    pids+=($!)
    while [ "$(jobs -r | wc -l)" -ge "$CONCURRENCY" ]; do
        wait -n 2>/dev/null || sleep 5
    done
done
wait

echo "[$(date '+%F %T')] ALL SAMPLES COMPLETE" | tee -a "$PROGRESS"
echo "===== BATCH SUMMARY =====" | tee -a "$PROGRESS"
column -t -s $'\t' "$SUMMARY" | tee -a "$PROGRESS"

#!/usr/bin/env bash
# =============================================================================
# Run scAPAtrap on GSE263303 mouse brain (Nf1+/-) Space Ranger BAMs.
#
# Uses the MOUSE-adapted launcher scripts/run_scapatrap_spaceranger_mouse.py
# (chromosomes chr1-19+XY, mouse GTF) — NOT the human-hardcoded launcher.
#
# Processes 2 selected GSMs sequentially:
#   GSM8189356 (K73-6-FMFC, animal K73)
#   GSM8189359 (K75-2-FMFC, animal K75)
#
# Prereqs (must already exist):
#   pipeline_output/gse263303_<GSM>_sr/outs/possorted_genome_bam.bam  (from spaceranger)
#   pipeline_output/gse263303_<GSM>_sr/outs/spatial/                  (from spaceranger)
#
# Hostname-aware (CLAUDE.md S90/S91/S97/S98). Designed to run in a tmux session.
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md mandatory) ---
HOSTNAME_SHORT=$(hostname -s | tr 'a-z' 'A-Z')
case "$HOSTNAME_SHORT" in
    S90) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         export TMPDIR="/s2/mengzijun/tmp"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    S91) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         export RSCRIPT="$HOME/anaconda3/envs/r442/bin/Rscript"
         export SAMTOOLS="$HOME/anaconda3/envs/samtools/bin/samtools"
         export TMPDIR="/s3/mengzijun/tmp"
         PYTHON_BIN="$HOME/anaconda3/envs/spagapa/bin/python" ;;
    S97) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         export TMPDIR="/s972/mengzijun/tmp"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    S98) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         export RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         export SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         export TMPDIR="/s982/mengzijun/tmp"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python" ;;
    *) echo "错误：未知服务器 $HOSTNAME_SHORT (expected S90/S91/S97/S98)." >&2; exit 1 ;;
esac
export TMPDIR
export OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
LAUNCHER="$ROOT/scripts/run_scapatrap_spaceranger_mouse.py"
# Mouse GRCm39 reference (matches the spaceranger transcriptome used upstream).
GTF="/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"
SR_OUT="$ROOT/pipeline_output"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs/${DATE_TAG}_gse263303_scapatrap.log"
mkdir -p "$(dirname "$SUITE_LOG")"

THREADS="${THREADS:-8}"

# GSM -> "sample_name tissue"
declare -A SAMPLE
SAMPLE[GSM8189356]="K73_6_FMFC"
SAMPLE[GSM8189359]="K75_2_FMFC"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }

run_one() {
    local gsm="$1"
    local sample="${SAMPLE[$gsm]}"
    local sr_id="gse263303_${gsm}_sr"
    local bam="${SR_OUT}/${sr_id}/outs/possorted_genome_bam.bam"
    local spatial="${SR_OUT}/${sr_id}/outs/spatial"
    local out_root="${SR_OUT}/gse263303_${gsm}_scapatrap"
    local proc_dir="$ROOT/data/processed/gse263303_${gsm}_scapatrap"
    local ds_name="gse263303_${gsm}_scapatrap"
    local per_log="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs/${DATE_TAG}_${gsm}_scapatrap.log"

    log "=== ${gsm} (${sample}) ==="
    if [[ ! -f "$bam" ]]; then
        log "ERROR(${gsm}): BAM missing ${bam}; skipping"
        return 2
    fi
    if [[ ! -d "$spatial" ]]; then
        log "ERROR(${gsm}): spatial dir missing ${spatial}; skipping"
        return 2
    fi
    # skip if processed output already present
    if [[ -f "${out_root}/raw_scapatrap/scapatrap_qc.json" && -f "${proc_dir}/qc_summary.json" ]]; then
        log "SKIP(${gsm}): scapatrap already complete (${proc_dir})"
        return 0
    fi

    log "bam=${bam}"
    log "spatial=${spatial}"
    log "gtf=${GTF} (mouse GRCm39)"
    log "out=${out_root}"
    log "processed=${proc_dir}"

    "$PYTHON_BIN" "$LAUNCHER" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --gtf "$GTF" \
        --output-root "$out_root" \
        --processed-dir "$proc_dir" \
        --dataset-name "$ds_name" \
        --source-label "GSE263303 ${gsm} ${sample} mouse brain Nf1+/- Space Ranger BAM + scAPAtrap" \
        --species mouse \
        --tissue "brain" \
        --threads "$THREADS" \
        --readlength 90 \
        > "$per_log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log "FAIL(${gsm}) rc=${rc}; tail of ${per_log}:"
        tail -n 40 "$per_log" | sed 's/^/      /' | tee -a "$SUITE_LOG"
        return $rc
    fi
    log "OK(${gsm}) -> ${proc_dir}"
    # report PAS counts
    if [[ -f "${proc_dir}/qc_summary.json" ]]; then
        python3 -c "
import json
q=json.load(open('${proc_dir}/qc_summary.json'))
print('    n_called_sites=', q.get('n_called_sites'), 'n_gene_annotated=', q.get('n_gene_annotated_sites'), 'n_apa_usage_sites=', q.get('n_apa_usage_sites'), 'n_spots=', q.get('n_spots'))
" 2>/dev/null | tee -a "$SUITE_LOG"
    fi
    return 0
}

# --- main ---
log "=== GSE263303 scAPAtrap (MOUSE) ==="
log "host=$HOSTNAME_SHORT launcher=$LAUNCHER threads=$THREADS tmp=$TMPDIR"
log "gtf=$GTF  species=mouse  chrs=chr1-19+chrX/chrY"

rc_all=0
for gsm in GSM8189356 GSM8189359; do
    run_one "$gsm" || rc_all=1
done

log "=== DONE (overall rc=${rc_all}) ==="
log "processed dirs:"
ls -d "$ROOT"/data/processed/gse263303_*_scapatrap 2>/dev/null | tee -a "$SUITE_LOG"
exit $rc_all

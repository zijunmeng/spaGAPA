#!/usr/bin/env bash
# =============================================================================
# GSE263303 (mouse brain Nf1+/- Visium) -- merge FASTQs -> spaceranger -> scAPAtrap
#
# fasterq-dump is ALREADY COMPLETE for all 4 SRRs:
#   data/raw/gse263303/fastq_split/<GSM>/SRR*_{1,2,3,4}.fastq
# Read layout (verified): _1=10bp I1, _2=10bp I2, _3=28bp barcode+UMI (->R1),
#                         _4=90bp cDNA (->R2)
#
# For each GSM (2 SRR technical splits each):
#   1. cat both SRR _3 files -> one <sample>_S1_L001_R1_001.fastq.gz (28bp)
#   2. cat both SRR _4 files -> one <sample>_S1_L001_R2_001.fastq.gz (90bp)
#      (_1/_2 index reads are NOT needed by spaceranger.)
#   3. spaceranger count (mouse refdata-gex-GRCm39-2024-A, --create-bam)
#   4. scAPAtrap via scripts/run_scapatrap_spaceranger.py (mouse)
#
# Idempotent: skips a step if its done-marker/output already exists.
# Hostname-aware (CLAUDE.md S90/S91/S97/S98). Run inside tmux. Do NOT git commit.
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection (CLAUDE.md mandatory) ---
HOSTNAME_SHORT=$(hostname -s | tr 'a-z' 'A-Z')
case "$HOSTNAME_SHORT" in
    S90) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
         TMPDIR="/s2/mengzijun/tmp" ;;
    S91) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         RSCRIPT="$HOME/anaconda3/envs/r442/bin/Rscript"
         SAMTOOLS="$HOME/anaconda3/envs/samtools/bin/samtools"
         PYTHON_BIN="$HOME/anaconda3/envs/spagapa/bin/python"
         TMPDIR="/s3/mengzijun/tmp" ;;
    S97) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
         TMPDIR="/s972/mengzijun/tmp" ;;
    S98) export R_LIBS="/s1/SHARE/01_software/R_442_SeuratV5/library"
         RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"
         SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools"
         PYTHON_BIN="/s1/mengzijun/anaconda3/envs/spagapa/bin/python"
         TMPDIR="/s982/mengzijun/tmp" ;;
    *) echo "错误：未知服务器 $HOSTNAME_SHORT (expected S90/S91/S97/S98)." >&2; exit 1 ;;
esac
export TMPDIR
export OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

# --- 2. Paths ---
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A"
GTF="${TRANSCRIPTOME}/genes/genes.gtf.gz"
LAUNCHER="${ROOT}/scripts/run_scapatrap_spaceranger.py"
DATA_RAW="${ROOT}/data/raw/gse263303"
SPLIT_BASE="${DATA_RAW}/fastq_split"
MERGED_BASE="${DATA_RAW}/merged_fastq"     # canonical gz R1/R2 per sample
IMAGE_DIR="${DATA_RAW}/images"             # per-GSM tissue H&E TIFFs (GEO suppl)
OUT_ROOT="${ROOT}/pipeline_output"
PROC_ROOT="${ROOT}/data/processed"
LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse263303_merge_sr_scapatrap.log"
mkdir -p "$LOG_ROOT" "$OUT_ROOT" "$MERGED_BASE" "$PROC_ROOT"

# spaceranger resources
SR_CORES="${SR_CORES:-16}"
SR_MEM="${SR_MEM:-120}"
SC_THREADS="${SC_THREADS:-8}"

# GSM -> "SAMPLE_NAME SRR1 SRR2"
declare -A GSM_INFO
GSM_INFO[GSM8189356]="K73_6_FMFC SRR28566759 SRR28566760"
GSM_INFO[GSM8189359]="K75_2_FMFC SRR28566753 SRR28566754"

# GSM -> tissue image filename stem (GEO suppl file: <GSM>_<stem>.tif).
# Full path constructed as ${IMAGE_DIR}/${gsm}_${SAMPLE_IMG[$gsm]}.tif
declare -A SAMPLE_IMG
SAMPLE_IMG[GSM8189356]="K73-6-FMFC"
SAMPLE_IMG[GSM8189359]="K75-2-FMFC"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}"; }
timed() {
    local label="$1"; shift
    local t0 t1 rc
    t0=$(date +%s)
    log "START ${label}"
    "$@"
    rc=$?
    t1=$(date +%s)
    log "END   ${label} (rc=${rc}, elapsed=$(( (t1-t0)/60 ))min)"
    return $rc
}

# Check read length of first sequence in a (possibly gzipped) FASTQ
seq_len() {
    python3 - "$1" <<'PY'
import gzip, sys
fn = sys.argv[1]
opener = gzip.open if fn.endswith(".gz") else open
with opener(fn, "rt", encoding="utf-8", errors="replace") as h:
    h.readline()
    print(len(h.readline().rstrip("\n")))
PY
}

# -----------------------------------------------------------------------------
# Step 1+2: merge the two SRR _3 (R1) and _4 (R2) reads into one FASTQ.gz each.
# -----------------------------------------------------------------------------
merge_fastqs() {
    local gsm="$1" sample="$2" srr1="$3" srr2="$4"
    local split_dir="${SPLIT_BASE}/${gsm}"
    local out_dir="${MERGED_BASE}/${gsm}"
    mkdir -p "$out_dir"

    local r1="${out_dir}/${sample}_S1_L001_R1_001.fastq.gz"
    local r2="${out_dir}/${sample}_S1_L001_R2_001.fastq.gz"
    local done_marker="${out_dir}/merge.done"

    if [[ -f "$done_marker" && -s "$r1" && -s "$r2" ]]; then
        log "SKIP merge ${gsm} (done marker present)"
        return 0
    fi

    # validate source files
    local f
    for f in "${split_dir}/${srr1}_3.fastq" "${split_dir}/${srr1}_4.fastq" \
             "${split_dir}/${srr2}_3.fastq" "${split_dir}/${srr2}_4.fastq"; do
        if [[ ! -s "$f" ]]; then
            log "ERROR(${gsm}): missing source FASTQ $f"
            return 2
        fi
    done

    # sanity-check lengths (28bp R1 / 90bp R2)
    local l3 l4
    l3=$(seq_len "${split_dir}/${srr1}_3.fastq")
    l4=$(seq_len "${split_dir}/${srr1}_4.fastq")
    log "${gsm} verify: ${srr1}_3=${l3}bp (expect 28)  ${srr1}_4=${l4}bp (expect 90)"
    if [[ "$l3" != "28" || "$l4" != "90" ]]; then
        log "ERROR(${gsm}): unexpected read lengths (_3=${l3}, _4=${l4}); aborting"
        return 2
    fi

    # NOTE: do NOT use `... | pigz > out || gzip > out`. If the pipe breaks
    # mid-stream the fallback reads EOF from the closed stdin and silently
    # writes a truncated/empty gzip (this exact bug truncated R2 to 20 bytes on
    # the first run). Capture the pipe status directly under a pipefail scope.
    log "merging R1: cat ${srr1}_3.fastq ${srr2}_3.fastq | pigz -> ${r1##*/}"
    rm -f "$r1"
    set -o pipefail
    cat "${split_dir}/${srr1}_3.fastq" "${split_dir}/${srr2}_3.fastq" \
        | pigz -p 8 -c > "$r1"
    local rc1=$?
    set +o pipefail
    if [[ $rc1 -ne 0 ]]; then log "FAIL(${gsm}) pigz R1 rc=$rc1"; return $rc1; fi

    log "merging R2: cat ${srr1}_4.fastq ${srr2}_4.fastq | pigz -> ${r2##*/}"
    rm -f "$r2"
    set -o pipefail
    cat "${split_dir}/${srr1}_4.fastq" "${split_dir}/${srr2}_4.fastq" \
        | pigz -p 8 -c > "$r2"
    local rc2=$?
    set +o pipefail
    if [[ $rc2 -ne 0 ]]; then log "FAIL(${gsm}) pigz R2 rc=$rc2"; return $rc2; fi

    # Validate merged gzipped FASTQs: gzip integrity + sane size (>100MB after
    # gzip; Visium runs are always GB-scale). Guards against silent truncation.
    local sz1 sz2
    sz1=$(stat -c%s "$r1"); sz2=$(stat -c%s "$r2")
    log "merged sizes: R1=$(numfmt --to=iec ${sz1}) R2=$(numfmt --to=iec ${sz2})"
    if [[ "$sz1" -lt 104857600 || "$sz2" -lt 104857600 ]]; then
        log "FAIL(${gsm}): merged FASTQ too small (R1=${sz1}, R2=${sz2}); aborting before spaceranger"
        return 2
    fi
    log "verifying gzip integrity of merged FASTQs..."
    gunzip -t "$r1" 2>/dev/null || { log "FAIL(${gsm}): R1 gzip corrupt"; return 2; }
    gunzip -t "$r2" 2>/dev/null || { log "FAIL(${gsm}): R2 gzip corrupt"; return 2; }

    log "merged files:"; ls -lh "$r1" "$r2" | tee -a "$SUITE_LOG" >/dev/null
    touch "$done_marker"
    log "OK merge ${gsm}"
    return 0
}

# -----------------------------------------------------------------------------
# Step 3: spaceranger count
# -----------------------------------------------------------------------------
run_spaceranger() {
    local gsm="$1" sample="$2"
    local out_id="gse263303_${gsm}_sr"
    local fastq_dir="${MERGED_BASE}/${gsm}"
    local sr_dir="${OUT_ROOT}/${out_id}"
    local bam="${sr_dir}/outs/possorted_genome_bam.bam"
    local summary="${sr_dir}/outs/web_summary.html"
    local per_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_spaceranger_count.log"

    # Tissue image: spaceranger 4.1.0 REQUIRES one of --image/--darkimage/
    # --colorizedimage/--cytaimage. GSE263303 publishes a per-GSM H&E brightfield
    # TIFF in the GEO supplementary files (downloaded to $IMAGE_DIR).
    local image_tif="${IMAGE_DIR}/${gsm}_${SAMPLE_IMG[$gsm]}.tif"
    if [[ ! -s "$image_tif" ]]; then
        log "ERROR(${gsm}): tissue image missing: ${image_tif}"
        return 2
    fi

    if [[ -s "$bam" && -f "$summary" ]]; then
        log "SKIP spaceranger ${gsm} (BAM + summary present)"
        return 0
    fi
    if [[ -d "$sr_dir" ]]; then
        log "WARN removing prior partial ${sr_dir}"
        rm -rf "$sr_dir"
    fi
    mkdir -p "$sr_dir"

    log "spaceranger count ${gsm}: fastqs=${fastq_dir} sample=${sample}"
    log "  image=${image_tif}"
    log "  --unknown-slide visium-1 --localcores=${SR_CORES} --localmem=${SR_MEM} --create-bam"
    (
        cd "$OUT_ROOT" && \
        "$SPACERANGER" count \
            --id "$out_id" \
            --description "GSE263303 ${gsm} ${sample} (mouse brain Nf1+/- Visium)" \
            --transcriptome "$TRANSCRIPTOME" \
            --fastqs "$fastq_dir" \
            --sample "$sample" \
            --image "$image_tif" \
            --unknown-slide visium-1 \
            --create-bam true \
            --localcores "$SR_CORES" \
            --localmem "$SR_MEM" \
            --disable-ui
    ) > "$per_log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log "FAIL spaceranger ${gsm} rc=${rc}; tail:"
        tail -n 50 "$per_log" | sed 's/^/      /' | tee -a "$SUITE_LOG" >/dev/null
        return $rc
    fi
    if [[ ! -s "$bam" ]]; then
        log "FAIL spaceranger ${gsm}: BAM missing at ${bam}"
        return 2
    fi
    log "OK spaceranger ${gsm} -> ${bam}"
    ls -lh "$bam" | tee -a "$SUITE_LOG" >/dev/null
    return 0
}

# -----------------------------------------------------------------------------
# Step 4: scAPAtrap via launcher (mouse)
# -----------------------------------------------------------------------------
run_scapatrap() {
    local gsm="$1" sample="$2"
    local sr_id="gse263303_${gsm}_sr"
    local bam="${OUT_ROOT}/${sr_id}/outs/possorted_genome_bam.bam"
    local spatial="${OUT_ROOT}/${sr_id}/outs/spatial"
    local out_root="${OUT_ROOT}/gse263303_${gsm}_scapatrap"
    local proc_dir="${PROC_ROOT}/gse263303_${gsm}_scapatrap"
    local ds_name="gse263303_${gsm}_scapatrap"
    local per_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_scapatrap.log"

    if [[ -f "${out_root}/raw_scapatrap/scapatrap_qc.json" && -f "${proc_dir}/qc_summary.json" ]]; then
        log "SKIP scapatrap ${gsm} (qc_summary.json present)"
        return 0
    fi
    if [[ ! -s "$bam" ]]; then
        log "ERROR(${gsm}): BAM missing ${bam}"; return 2
    fi
    if [[ ! -d "$spatial" ]]; then
        log "ERROR(${gsm}): spatial dir missing ${spatial}"; return 2
    fi

    log "scAPAtrap ${gsm}: bam=${bam} species=mouse"
    log "  out=${out_root}  processed=${proc_dir}"
    "$PYTHON_BIN" "$LAUNCHER" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --gtf "$GTF" \
        --output-root "$out_root" \
        --processed-dir "$proc_dir" \
        --dataset-name "$ds_name" \
        --source-label "GSE263303 ${gsm} ${sample} mouse brain Nf1+/- Space Ranger BAM + scAPAtrap" \
        --species mouse \
        --tissue "mouse brain (Nf1+/-)" \
        --threads "$SC_THREADS" \
        --readlength 90 \
        > "$per_log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log "FAIL scapatrap ${gsm} rc=${rc}; tail:"
        tail -n 50 "$per_log" | sed 's/^/      /' | tee -a "$SUITE_LOG" >/dev/null
        return $rc
    fi
    log "OK scapatrap ${gsm} -> ${proc_dir}"
    if [[ -f "${proc_dir}/qc_summary.json" ]]; then
        "$PYTHON_BIN" - "$proc_dir/qc_summary.json" <<'PY' | tee -a "$SUITE_LOG"
import json, sys
q = json.load(open(sys.argv[1]))
print(f"    n_called_sites={q.get('n_called_sites')} "
      f"n_gene_annotated={q.get('n_gene_annotated_sites')} "
      f"n_apa_usage_sites={q.get('n_apa_usage_sites')} "
      f"n_spots={q.get('n_spots')}")
PY
    fi
    return 0
}

# -----------------------------------------------------------------------------
# process one GSM through all 4 steps
# -----------------------------------------------------------------------------
process_gsm() {
    local gsm="$1"
    local info="${GSM_INFO[$gsm]}"
    # shellcheck disable=SC2206
    local parts=($info)
    local sample="${parts[0]}" srr1="${parts[1]}" srr2="${parts[2]}"
    log "==================== ${gsm} (sample=${sample}, runs=${srr1}+${srr2}) ===================="

    merge_fastqs     "$gsm" "$sample" "$srr1" "$srr2" || return 3
    run_spaceranger  "$gsm" "$sample"                  || return 4
    run_scapatrap    "$gsm" "$sample"                  || return 5
    return 0
}

# --- main ---
log "=== GSE263303 merge -> spaceranger -> scAPAtrap (mouse) ==="
log "host=$HOSTNAME_SHORT sr_cores=$SR_CORES sr_mem=$SR_MEM sc_threads=$SC_THREADS tmp=$TMPDIR"
log "transcriptome=$TRANSCRIPTOME  gtf=$GTF"

# Default: process both; if FIRST_ONLY=1, stop after the first GSM completes.
FIRST_ONLY="${FIRST_ONLY:-0}"

rc_all=0
for gsm in GSM8189356 GSM8189359; do
    if ! process_gsm "$gsm"; then
        rc_all=1
        log "GSM ${gsm} FAILED in pipeline; continuing to next GSM"
    else
        log "GSM ${gsm} COMPLETE"
    fi
    if [[ "$FIRST_ONLY" == "1" && "$gsm" == "GSM8189356" ]]; then
        log "FIRST_ONLY=1 -> stopping after GSM8189356"
        break
    fi
done

log "=== DONE (overall rc=${rc_all}) ==="
log "--- spaceranger outputs ---"
ls -d "${OUT_ROOT}"/gse263303_*_sr 2>/dev/null | tee -a "$SUITE_LOG" >/dev/null
log "--- scapatrap processed dirs ---"
ls -d "$PROC_ROOT"/gse263303_*_scapatrap 2>/dev/null | tee -a "$SUITE_LOG" >/dev/null
for gsm in GSM8189356 GSM8189359; do
    qc="${PROC_ROOT}/gse263303_${gsm}_scapatrap/qc_summary.json"
    if [[ -f "$qc" ]]; then
        "$PYTHON_BIN" - "$qc" <<'PY' | tee -a "$SUITE_LOG"
import json, sys
q = json.load(open(sys.argv[1]))
print(f"  {q.get('dataset')}: PAS n_called={q.get('n_called_sites')} "
      f"n_gene_anno={q.get('n_gene_annotated_sites')} "
      f"n_apa_usage={q.get('n_apa_usage_sites')} n_spots={q.get('n_spots')}")
PY
    fi
done
exit $rc_all

#!/usr/bin/env bash
# =============================================================================
# fasterq-dump GSE263303 SRA -> spaceranger count (MOUSE ref GRCm39).
#
# For each selected GSM (2 technical-split SRR runs each):
#   1. fasterq-dump --split-files --include-technical (per SRR)
#   2. auto-detect Visium read structure (which _N.fastq is R1 barcode/UMI vs R2 cDNA)
#      by sequence length: 28bp == barcode/UMI -> R1; 90bp == cDNA -> R2; I1 = sample index
#   3. create spaceranger-compatible gzipped FASTQs with canonical naming:
#        <sample>_S1_L<lane>_R1_001.fastq.gz  (barcode+UMI, 28bp)
#        <sample>_S1_L<lane>_R2_001.fastq.gz  (cDNA, 90bp)
#   4. spaceranger count --transcriptome refdata-gex-GRCm39-2024-A --create-bam true
#
# Hostname-aware (CLAUDE.md S90/S91/S97/S98).
# =============================================================================
set -uo pipefail

# --- 1. Hostname detection ---
HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S90) export TMPDIR="/s2/mengzijun/tmp"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools" ;;
    S91) export TMPDIR="/s3/mengzijun/tmp"; RSCRIPT="$HOME/anaconda3/envs/r442/bin/Rscript"; SAMTOOLS="$HOME/anaconda3/envs/samtools/bin/samtools" ;;
    S97) export TMPDIR="/s972/mengzijun/tmp"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools" ;;
    S98) export TMPDIR="/s982/mengzijun/tmp"; RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript"; SAMTOOLS="/s1/mengzijun/anaconda3/envs/samtools/bin/samtools" ;;
    *)  echo "ERROR unknown host $HOSTNAME" >&2; exit 1 ;;
esac
mkdir -p "$TMPDIR"

# --- 2. Paths ---
export OPENBLAS_NUM_THREADS=8
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SRATOOLKIT="/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin"
FASTERQ="${SRATOOLKIT}/fasterq-dump"
VDB_CONFIG="${SRATOOLKIT}/vdb-config"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCm39-2024-A"
DATA_ROOT="${ROOT}/data/raw/gse263303"
SRA_DIR="${DATA_ROOT}"
SPLIT_BASE="${DATA_ROOT}/fastq_split"
INPUT_BASE="${DATA_ROOT}/spaceranger_inputs"
OUT_ROOT="${ROOT}/pipeline_output"
LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse263303_fasterq_spaceranger.log"
mkdir -p "$LOG_ROOT"

THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-48}"

# GSM -> "SAMPLE_NAME SRR1 SRR2"  (sample_name used as spaceranger --sample)
declare -A GSM_INFO
GSM_INFO[GSM8189356]="K73_6_FMFC SRR28566759 SRR28566760"
GSM_INFO[GSM8189359]="K75_2_FMFC SRR28566753 SRR28566754"

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

# Read the 2nd line (sequence) length of a fastq (handles .gz)
first_seq_len() {
    python3 - "$1" <<'PY'
import gzip, sys
fn = sys.argv[1]
opener = gzip.open if fn.endswith(".gz") else open
with opener(fn, "rt", encoding="utf-8", errors="replace") as h:
    h.readline(); seq = h.readline().rstrip("\n")
    print(len(seq))
PY
}

# Determine which read index (_1,_2,_3,_4) is the 28bp barcode/UMI (R1) and which is cDNA (R2)
# Returns: echoes "R1_IDX R2_LEN R2_IDX R2_LEN"
detect_reads() {
    local split_dir="$1" run="$2"
    local best_r1_idx="" best_r1_len=99999
    local best_r2_idx="" best_r2_len=0
    for idx in 1 2 3 4; do
        local fq="${split_dir}/${run}_${idx}.fastq"
        local gz="${fq}.gz"
        local f=""
        if [[ -s "$fq" ]]; then f="$fq"; elif [[ -s "$gz" ]]; then f="$gz"; else continue; fi
        local len
        len=$(first_seq_len "$f")
        log "    ${run}_${idx}.fastq  seq_len=${len}"
        # R1 = barcode/UMI ~ 28bp (prefer closest to 28 from below or equal)
        if [[ "$len" -le 30 && "$len" -lt "$best_r1_len" ]]; then
            best_r1_idx="$idx"; best_r1_len="$len"
        fi
        # R2 = cDNA, longest read (typically 90bp)
        if [[ "$len" -gt "$best_r2_len" ]]; then
            best_r2_idx="$idx"; best_r2_len="$len"
        fi
    done
    echo "${best_r1_idx} ${best_r1_len} ${best_r2_idx} ${best_r2_len}"
}

process_gsm() {
    local gsm="$1"
    local info="${GSM_INFO[$gsm]}"
    # shellcheck disable=SC2206
    local parts=($info)
    local sample="${parts[0]}"
    local srr1="${parts[1]}"
    local srr2="${parts[2]}"
    log "=== ${gsm} sample=${sample} runs=${srr1} ${srr2} ==="

    local split_dir="${SPLIT_BASE}/${gsm}"
    local input_dir="${INPUT_BASE}/${gsm}"
    local out_id="gse263303_${gsm}_sr"
    mkdir -p "$split_dir" "$input_dir"

    # --- fasterq-dump each SRR ---
    for srr in "$srr1" "$srr2"; do
        # aria2c S3 produces bare file named "$srr" (no .sra ext); fasterq-dump
        # needs either an accession or a .sra-named file, so create a .sra symlink.
        local sra_file="${SRA_DIR}/${srr}.sra"
        local bare_file="${SRA_DIR}/${srr}"
        if [[ ! -e "$sra_file" ]]; then
            if [[ -s "$bare_file" ]]; then
                ln -sf "$bare_file" "$sra_file"
            elif [[ -s "${bare_file}.sra" ]]; then
                :  # already .sra
            else
                log "ERROR: SRA missing ${srr} (looked for ${bare_file} and ${sra_file}); skipping ${gsm}"
                return 2
            fi
        fi
        sra_file="${SRA_DIR}/${srr}.sra"
        if [[ ! -s "$sra_file" ]]; then
            log "ERROR: SRA file empty/missing ${sra_file}; skipping ${gsm}"
            return 2
        fi
        local marker="${split_dir}/${srr}.fasterq.done"
        if [[ -f "$marker" ]]; then
            log "SKIP fasterq ${srr} (done marker present)"
            continue
        fi
        # fasterq-dump writes ${run}_1.fastq, _2.fastq, _3.fastq, _4.fastq
        timed "fasterq-dump ${srr}" \
            "${FASTERQ}" --split-files --include-technical \
                --threads "${THREADS}" \
                --temp "${split_dir}/tmp_${srr}" \
                --outdir "$split_dir" "$sra_file" \
                && touch "$marker"
        local rc=$?
        rm -rf "${split_dir}/tmp_${srr}"
        if [[ $rc -ne 0 ]]; then log "FAIL fasterq ${srr} rc=${rc}"; return $rc; fi
    done

    # --- detect read structure from first SRR ---
    log "--- detect Visium read structure (${srr1}) ---"
    local detected r1_idx r1_len r2_idx r2_len
    detected=$(detect_reads "$split_dir" "$srr1")
    # shellcheck disable=SC2206
    detected=($detected)
    r1_idx="${detected[0]}"; r1_len="${detected[1]}"
    r2_idx="${detected[2]}"; r2_len="${detected[3]}"
    if [[ -z "$r1_idx" || -z "$r2_idx" || "$r1_idx" == "$r2_idx" ]]; then
        log "ERROR(${gsm}): cannot resolve R1/R2 reads (r1_idx=${r1_idx} r2_idx=${r2_idx}); listing split dir:"
        ls -la "$split_dir" | tee -a "$SUITE_LOG"
        return 2
    fi
    log "resolved: R1=${r1_idx}fastq(${r1_len}bp=barcode/UMI)  R2=${r2_idx}fastq(${r2_len}bp=cDNA)"

    # --- build spaceranger-compatible gzipped FASTQs ---
    local fastq_link_dir="${input_dir}/fastqs"
    mkdir -p "$fastq_link_dir"
    rm -f "$fastq_link_dir"/*.fastq.gz
    local lane=1
    for srr in "$srr1" "$srr2"; do
        # gzip if needed, then symlink with canonical name
        for role in R1:${r1_idx} R2:${r2_idx}; do
            local rname="${role%%:*}"
            local ridx="${role##*:}"
            local src="${split_dir}/${srr}_${ridx}.fastq"
            local srcgz="${src}.gz"
            if [[ -s "$srcgz" ]]; then src="$srcgz"; fi
            if [[ ! -s "$src" ]]; then log "ERROR missing ${src}"; return 2; fi
            local dest="${fastq_link_dir}/${sample}_S1_L$(printf '%03d' $lane)_${rname}_001.fastq.gz"
            if [[ "$src" == *.gz ]]; then
                cp "$src" "$dest"
            else
                pigz -p "${THREADS}" -c "$src" > "$dest" || gzip -c "$src" > "$dest"
            fi
        done
        lane=$((lane+1))
    done
    log "spaceranger fastqs:"; ls -lh "$fastq_link_dir" | tee -a "$SUITE_LOG"

    # --- spaceranger count ---
    if [[ -f "${OUT_ROOT}/${out_id}/outs/count_summary.html" || -f "${OUT_ROOT}/${out_id}/outs/web_summary.html" ]]; then
        log "SKIP spaceranger ${gsm} (output complete)"
        return 0
    fi
    if [[ -d "${OUT_ROOT}/${out_id}" ]]; then
        log "WARN removing prior partial ${OUT_ROOT}/${out_id}"
        rm -rf "${OUT_ROOT}/${out_id}"
    fi
    mkdir -p "$OUT_ROOT"
    local sample_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_spaceranger_count.log"
    (
        cd "$OUT_ROOT" && \
        "$SPACERANGER" count \
            --id "$out_id" \
            --description "GSE263303 ${gsm} ${sample} (mouse brain Nf1+/- Visium)" \
            --transcriptome "$TRANSCRIPTOME" \
            --fastqs "$fastq_link_dir" \
            --sample "$sample" \
            --unknown-slide visium-1 \
            --create-bam true \
            --localcores "$THREADS" \
            --localmem "$MEM_GB" \
            --disable-ui
    ) > "$sample_log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log "FAIL spaceranger ${gsm} rc=${rc}; tail:"
        tail -n 40 "$sample_log" | sed 's/^/      /' | tee -a "$SUITE_LOG"
        return $rc
    fi
    log "OK spaceranger ${gsm} -> ${OUT_ROOT}/${out_id}/outs"
    return 0
}

# --- main ---
log "=== GSE263303 fasterq-dump + spaceranger count (mouse GRCm39) ==="
log "host=$HOSTNAME threads=$THREADS mem=$MEM_GB tmp=$TMPDIR"
log "transcriptome=$TRANSCRIPTOME"

rc_all=0
for gsm in GSM8189356 GSM8189359; do
    if ! process_gsm "$gsm"; then rc_all=1; fi
done

log "=== DONE (overall rc=${rc_all}) ==="
log "spaceranger outputs:"
ls -d "${OUT_ROOT}"/gse263303_*_sr 2>/dev/null | tee -a "$SUITE_LOG"
exit $rc_all

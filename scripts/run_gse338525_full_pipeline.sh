#!/usr/bin/env bash
# =============================================================================
# GSE338525 (human NORMAL LIVER, Fresh Frozen polyA Visium) -> APA analysis.
#
# Chains: fasterq-dump -> spaceranger count (HUMAN GRCh38) -> scAPAtrap (human)
#
# 2 selected normal-liver samples (the 2 smallest SRRs by SRA size):
#   SRR39611188  GSM9876374  Normal liver Patient 2   (9.4 GB)
#   SRR39611189  GSM9876373  Normal liver Patient 1   (9.7 GB)
# Each GSM maps to a SINGLE SRR run (no technical split).
# Read structure: R1=28bp (barcode+UMI), R2=150bp (cDNA), I1/I2=10bp. Standard Visium FF.
#
# SRA may store 4 reads/spot (R1,R2,I1,I2). spaceranger only needs R1+R2; we detect
# which split read is the 28bp barcode/UMI (R1) vs the longest cDNA read (R2) by
# sequence length, then build canonical <sample>_S1_L001_R1/R2_001.fastq.gz.
#
# spaceranger 4.1.0 with --unknown-slide visium-1 does NOT require an image
# (validated in run_gse263303_fasterq_spaceranger.sh).
#
# Hostname-aware (CLAUDE.md S90/S91/S97/S98). Run inside tmux on S91.
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
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A"
GTF="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
DATA_ROOT="${ROOT}/data/raw/gse338525"
SRA_DIR="${DATA_ROOT}/sra"
SPLIT_BASE="${DATA_ROOT}/fastq_split"
INPUT_BASE="${DATA_ROOT}/spaceranger_inputs"
OUT_ROOT="${ROOT}/pipeline_output"
PROC_ROOT="${ROOT}/data/processed"
LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"
SUITE_LOG="${LOG_ROOT}/${DATE_TAG}_gse338525_pipeline.log"
mkdir -p "$LOG_ROOT"

THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-48}"
SCAPA_THREADS="${SCAPA_THREADS:-12}"
# Prefer the spagapa conda env python (numpy/pandas), matching run_scapatrap_gse263303.sh
PYTHON_BIN="${PYTHON_BIN:-$HOME/anaconda3/envs/spagapa/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN="python3"

# GSM -> "SAMPLE_NAME SRR"  (single SRR per GSM; sample_name used as spaceranger --sample)
declare -A GSM_INFO
GSM_INFO[GSM9876374]="NL_P2 SRR39611188"   # Normal liver Patient 2 (smallest, 9.4 GB)
GSM_INFO[GSM9876373]="NL_P1 SRR39611189"   # Normal liver Patient 1 (9.7 GB)

# NOTE: log() writes to STDERR (via tee) so it does NOT pollute stdout when a
# function is captured via command substitution (e.g. `detected=$(detect_reads ...)`).
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "${SUITE_LOG}" >&2; }
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

# Length of the 2nd line (sequence) of a FASTQ (gz or plain)
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

# Determine which split read index (_1.._N) is the 28bp barcode/UMI (R1) and which is cDNA (R2).
# GSE338525 layout: _1=I1(10bp), _2=I2(10bp), _3=R1(28bp barcode+UMI), _4=R2(150bp cDNA).
# Robust selection: R1 = read closest to 28bp (min abs distance, tie->shorter idx);
# R2 = longest read (the cDNA). This disambiguates the two 10bp index reads from the 28bp R1.
detect_reads() {
    local split_dir="$1" run="$2"
    local best_r1_idx="" best_r1_dist=99999 best_r1_len=0
    local best_r2_idx="" best_r2_len=0
    local idx len dist
    for idx in 1 2 3 4 5; do
        local fq="${split_dir}/${run}_${idx}.fastq"
        local gz="${fq}.gz"
        local f=""
        if [[ -s "$fq" ]]; then f="$fq"; elif [[ -s "$gz" ]]; then f="$gz"; else continue; fi
        len=$(first_seq_len "$f")
        log "    ${run}_${idx}.fastq  seq_len=${len}"
        # R1 = barcode/UMI closest to 28bp (exclude the long cDNA read)
        if [[ "$len" -lt 60 ]]; then
            dist=$(( len > 28 ? len - 28 : 28 - len ))
            if [[ -z "$best_r1_idx" || "$dist" -lt "$best_r1_dist" ]]; then
                best_r1_idx="$idx"; best_r1_dist="$dist"; best_r1_len="$len"
            fi
        fi
        # R2 = cDNA, longest read
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
    local srr="${parts[1]}"
    log "=== ${gsm} sample=${sample} run=${srr} ==="

    local split_dir="${SPLIT_BASE}/${gsm}"
    local input_dir="${INPUT_BASE}/${gsm}"
    local out_id="gse338525_${gsm}_sr"
    mkdir -p "$split_dir" "$input_dir"

    # --- fasterq-dump ---
    local sra_file="${SRA_DIR}/${srr}.sra"
    local bare_file="${SRA_DIR}/${srr}"
    if [[ ! -e "$sra_file" ]]; then
        if [[ -s "$bare_file" ]]; then
            ln -sf "$bare_file" "$sra_file"
        elif [[ -s "${bare_file}.sra" ]]; then
            :
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
    else
        # fasterq-dump writes ${run}_1.fastq .. _N.fastq
        timed "fasterq-dump ${srr}" \
            "${FASTERQ}" --split-files --include-technical \
                --threads "${THREADS}" \
                --temp "${split_dir}/tmp_${srr}" \
                --outdir "$split_dir" "$sra_file" \
                && touch "$marker"
        local rc=$?
        rm -rf "${split_dir}/tmp_${srr}"
        if [[ $rc -ne 0 ]]; then log "FAIL fasterq ${srr} rc=${rc}"; return $rc; fi
    fi

    # --- detect read structure ---
    log "--- detect Visium read structure (${srr}) ---"
    local detected r1_idx r1_len r2_idx r2_len
    detected=$(detect_reads "$split_dir" "$srr")
    # shellcheck disable=SC2206
    detected=($detected)
    r1_idx="${detected[0]}"; r1_len="${detected[1]}"
    r2_idx="${detected[2]}"; r2_len="${detected[3]}"
    if [[ -z "$r1_idx" || -z "$r2_idx" || "$r1_idx" == "$r2_idx" ]]; then
        log "ERROR(${gsm}): cannot resolve R1/R2 reads (r1_idx=${r1_idx} r2_idx=${r2_idx}); listing split dir:"
        ls -la "$split_dir" | tee -a "$SUITE_LOG"
        return 2
    fi
    log "resolved: R1=${r1_idx}.fastq(${r1_len}bp=barcode/UMI)  R2=${r2_idx}.fastq(${r2_len}bp=cDNA)"

    # --- build spaceranger-compatible gzipped FASTQs ---
    local fastq_link_dir="${input_dir}/fastqs"
    mkdir -p "$fastq_link_dir"
    rm -f "$fastq_link_dir"/*.fastq.gz
    local src dst
    for role in R1:${r1_idx} R2:${r2_idx}; do
        local rname="${role%%:*}"
        local ridx="${role##*:}"
        src="${split_dir}/${srr}_${ridx}.fastq"
        local srcgz="${src}.gz"
        if [[ -s "$srcgz" ]]; then src="$srcgz"; fi
        if [[ ! -s "$src" ]]; then log "ERROR missing ${src}"; return 2; fi
        dest="${fastq_link_dir}/${sample}_S1_L001_${rname}_001.fastq.gz"
        if [[ "$src" == *.gz ]]; then
            cp "$src" "$dest"
        else
            pigz -p "${THREADS}" -c "$src" > "$dest" || gzip -c "$src" > "$dest"
        fi
    done
    log "spaceranger fastqs:"; ls -lh "$fastq_link_dir" | tee -a "$SUITE_LOG"

    # --- spaceranger count ---
    if [[ -f "${OUT_ROOT}/${out_id}/outs/count_summary.html" || -f "${OUT_ROOT}/${out_id}/outs/web_summary.html" ]]; then
        log "SKIP spaceranger ${gsm} (output complete)"
    else
        if [[ -d "${OUT_ROOT}/${out_id}" ]]; then
            log "WARN removing prior partial ${OUT_ROOT}/${out_id}"
            rm -rf "${OUT_ROOT}/${out_id}"
        fi
        mkdir -p "$OUT_ROOT"
        local sample_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_spaceranger_count.log"
        # NOTE: spaceranger 4.1.0 REQUIRES an image even with --unknown-slide (validated
        # the hard way: it errors with "required arguments not provided: <--image...>").
        # GEO suppl images (fluorescent CD68) are 403-blocked from this host, so we use a
        # synthetic placeholder TIFF. With --unknown-slide, spot detection comes from the
        # barcode whitelist (not the image), so a placeholder does not affect the BAM /
        # tissue_positions outputs that scAPAtrap needs. --reorient-images=false prevents
        # spaceranger from trying to rotate a non-tissue synthetic image.
        local placeholder_img="${DATA_ROOT}/images/placeholder_tissue.tif"
        (
            cd "$OUT_ROOT" && \
            "$SPACERANGER" count \
                --id "$out_id" \
                --description "GSE338525 ${gsm} ${sample} (human normal liver Visium FF)" \
                --transcriptome "$TRANSCRIPTOME" \
                --fastqs "$fastq_link_dir" \
                --sample "$sample" \
                --unknown-slide visium-1 \
                --image "$placeholder_img" \
                --reorient-images false \
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
    fi

    # --- scAPAtrap (human) via launcher ---
    local bam="${OUT_ROOT}/${out_id}/outs/possorted_genome_bam.bam"
    local spatial="${OUT_ROOT}/${out_id}/outs/spatial"
    local scapatrap_out="${OUT_ROOT}/gse338525_${gsm}_scapatrap"
    local proc_dir="${PROC_ROOT}/gse338525_${gsm,,}_scapatrap"
    local dataset="gse338525_${gsm,,}_scapatrap"
    local src_label="GSE338525 ${gsm} ${sample} (human normal liver) Space Ranger BAM + scAPAtrap"

    if [[ -f "${proc_dir}/qc_summary.json" ]]; then
        log "SKIP scAPAtrap ${gsm} (processed qc_summary.json present)"
    else
        if [[ ! -s "$bam" ]]; then
            log "ERROR BAM not found ${bam}; skipping scAPAtrap ${gsm}"
            return 3
        fi
        if [[ ! -d "$spatial" ]]; then
            log "ERROR spatial dir not found ${spatial}; skipping scAPAtrap ${gsm}"
            return 3
        fi
        local sc_log="${LOG_ROOT}/${DATE_TAG}_${gsm}_scapatrap.log"
        log "running scAPAtrap launcher -> ${proc_dir}"
        (
            "$RSCRIPT" -e 'suppressMessages(library(scAPAtrap))' >/dev/null 2>&1 || true
            "$PYTHON_BIN" "${ROOT}/scripts/run_scapatrap_spaceranger.py" \
                --bam "$bam" \
                --spatial-dir "$spatial" \
                --gtf "$GTF" \
                --output-root "$scapatrap_out" \
                --processed-dir "$proc_dir" \
                --dataset-name "$dataset" \
                --source-label "$src_label" \
                --species human \
                --tissue "normal liver" \
                --threads "$SCAPA_THREADS" \
                --readlength 90 \
                --force
        ) > "$sc_log" 2>&1
        local rc=$?
        if [[ $rc -ne 0 ]]; then
            log "FAIL scAPAtrap ${gsm} rc=${rc}; tail:"
            tail -n 40 "$sc_log" | sed 's/^/      /' | tee -a "$SUITE_LOG"
            return $rc
        fi
        log "OK scAPAtrap ${gsm} -> ${proc_dir}"
    fi
    return 0
}

# --- main ---
log "=== GSE338525 pipeline: fasterq-dump -> spaceranger (HUMAN GRCh38) -> scAPAtrap (human) ==="
log "host=$HOSTNAME threads=$THREADS mem=${MEM_GB}GB scapa_threads=${SCAPA_THREADS} tmp=$TMPDIR"
log "transcriptome=$TRANSCRIPTOME"
log "samples: GSM9876374 (SRR39611188, NL P2), GSM9876373 (SRR39611189, NL P1)"

rc_all=0
for gsm in GSM9876374 GSM9876373; do
    if ! process_gsm "$gsm"; then rc_all=1; fi
done

log "=== DONE (overall rc=${rc_all}) ==="
log "spaceranger outputs:"
ls -d "${OUT_ROOT}"/gse338525_*_sr 2>/dev/null | tee -a "$SUITE_LOG"
log "scapatrap processed dirs:"
ls -d "${PROC_ROOT}"/gse338525_*_scapatrap 2>/dev/null | tee -a "$SUITE_LOG"
for gsm in GSM9876374 GSM9876373; do
    qc="${PROC_ROOT}/gse338525_${gsm,,}_scapatrap/qc_summary.json"
    if [[ -f "$qc" ]]; then
        python3 -c "
import json
q=json.load(open('$qc'))
print('  ${gsm}: n_called_sites=', q['n_called_sites'], 'n_gene_annotated=', q['n_gene_annotated_sites'], 'n_apa_usage_sites=', q['n_apa_usage_sites'], 'n_spots=', q['n_spots'])
" 2>/dev/null | tee -a "$SUITE_LOG"
    fi
done
exit $rc_all

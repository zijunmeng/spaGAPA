#!/usr/bin/env bash
# =============================================================================
# GSE243981 (human liver, HEALTHY C73 Visium) — full pipeline:
#   1. Link SRA-split FASTQs to Space Ranger-named R1/R2 per sample
#   2. spaceranger count (human GRCh38) -> BAM + spatial/
#   3. scAPAtrap via scripts/run_scapatrap_spaceranger.py -> PAS sites + matrices
#
# Two HEALTHY samples (NOT PSC/PBC):
#   GSM7697868 (C73_A1_VISIUM) -> SRR25581696
#   GSM7697869 (C73_B1_VISIUM) -> SRR25581695
#
# Hostname rules per CLAUDE.md (S90/S91/S97/S98). Designed for a tmux session.
# Usage:
#   ./run_gse243981_spaceranger_and_scapatrap.sh prep     # link fastqs + dry-check
#   ./run_gse243981_spaceranger_and_scapatrap.sh count    # spaceranger count only
#   ./run_gse243981_spaceranger_and_scapatrap.sh scapatrap # scAPAtrap only
#   ./run_gse243981_spaceranger_and_scapatrap.sh all      # prep + count + scapatrap
# =============================================================================
set -uo pipefail

# -----------------------------------------------------------------------------
# 1. Hostname detection (CLAUDE.md mandatory)
# -----------------------------------------------------------------------------
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
export TMPDIR OPENBLAS_NUM_THREADS=8
mkdir -p "$TMPDIR"

# -----------------------------------------------------------------------------
# 2. Paths
# -----------------------------------------------------------------------------
ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse243981"
SPACERANGER="/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger"
TRANSCRIPTOME="/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A"
LAUNCHER="$ROOT/scripts/run_scapatrap_spaceranger.py"
OUT_ROOT="$ROOT/pipeline_output"
LOG_ROOT="/s1/SHARE/mengzijun/01_project/26_spaGAPA/logs"
mkdir -p "$OUT_ROOT" "$LOG_ROOT"

THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-64}"
DATE_TAG="${DATE_TAG:-$(date +%Y%m%d)}"

# Sample table: GSM -> SRR
declare -A SRR
SRR[GSM7697868]=SRR25581696
SRR[GSM7697869]=SRR25581695
GSMS=(GSM7697868 GSM7697869)

log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

# -----------------------------------------------------------------------------
# 3. Link SRA-split FASTQs to Space Ranger R1/R2 layout
#    GSE243981 SRA is PAIRED Visium: read1 = 28 bp (barcode+UMI), read2 = 90 bp cDNA.
#    spaceranger expects <sample>_S1_L001_R1_001.fastq.gz and R2.
# -----------------------------------------------------------------------------
prepare_links() {
    local gsm="$1" srr="${SRR[$1]}"
    local split_dir="$DATA_ROOT/fastq_split/$gsm"
    local input_dir="$DATA_ROOT/spaceranger_inputs/$gsm"
    local fq_link="$input_dir/fastqs"
    mkdir -p "$fq_link"

    # Detect actual FASTQ naming from fasterq-dump split-files
    local r1 r2
    if [[ -s "$split_dir/${srr}_1.fastq" ]]; then
        r1="$split_dir/${srr}_1.fastq"; r2="$split_dir/${srr}_2.fastq"
    elif [[ -s "$split_dir/${srr}_1.fastq.gz" ]]; then
        r1="$split_dir/${srr}_1.fastq.gz"; r2="$split_dir/${srr}_2.fastq.gz"
    else
        log "ERROR: cannot find split FASTQ for ${srr} in ${split_dir}"
        ls -la "$split_dir" >&2 || true
        return 1
    fi
    rm -f "$fq_link"/*.fastq.gz "$fq_link"/*.fastq
    ln -s "$r1" "$fq_link/${gsm}_S1_L001_R1_001.${r1##*.}"
    ln -s "$r2" "$fq_link/${gsm}_S1_L001_R2_001.${r2##*.}"
    log "linked ${gsm}: R1=${r1##*/} R2=${r2##*/}"
    # quick read-length audit
    python - "$r1" "$r2" <<'PY'
import gzip, sys
def len2(p):
    o = gzip.open if str(p).endswith(".gz") else open
    with o(p, "rt", errors="replace") as h:
        h.readline(); return len(h.readline().rstrip("\n"))
r1, r2 = sys.argv[1], sys.argv[2]
print(f"  R1 len={len2(r1)}  R2 len={len2(r2)}")
PY
}

# -----------------------------------------------------------------------------
# 4. spaceranger count (one sample)
# -----------------------------------------------------------------------------
run_count() {
    local gsm="$1"
    local count_id="gse243981_${gsm}_sr"
    local pipestance="$OUT_ROOT/${count_id}"
    local fq_link="$DATA_ROOT/spaceranger_inputs/$gsm/fastqs"
    local sample_log="$LOG_ROOT/${DATE_TAG}_${gsm}_spaceranger_count.log"

    if [[ -f "$pipestance/outs/possorted_genome_bam.bam" ]]; then
        log "SKIP spaceranger ${gsm}: BAM already exists at ${pipestance}/outs/possorted_genome_bam.bam"
        return 0
    fi

    log "START spaceranger count ${gsm} (id=${count_id})"
    local t0 t1 rc
    t0=$(date +%s)
    "$SPACERANGER" count \
        --id "$count_id" \
        --description "GSE243981 ${gsm} C73 healthy liver Visium count" \
        --transcriptome "$TRANSCRIPTOME" \
        --fastqs "$fq_link" \
        --sample "$gsm" \
        --unknown-slide visium-1 \
        --create-bam true \
        --localcores "$THREADS" \
        --localmem "$MEM_GB" \
        --disable-cell-annotation \
        --disable-ui \
        > "$sample_log" 2>&1
    rc=$?
    t1=$(date +%s)
    log "END spaceranger count ${gsm} (rc=${rc}, elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f min",(b-a)/60}'))"
    if [[ $rc -ne 0 ]]; then
        log "ERROR spaceranger ${gsm} failed; see ${sample_log}"
        tail -20 "$sample_log" >&2 || true
        return $rc
    fi
    if [[ ! -f "$pipestance/outs/possorted_genome_bam.bam" ]]; then
        log "ERROR spaceranger ${gsm}: BAM not produced"
        return 1
    fi
    return 0
}

# -----------------------------------------------------------------------------
# 5. scAPAtrap via launcher
# -----------------------------------------------------------------------------
run_scapatrap() {
    local gsm="$1"
    local count_id="gse243981_${gsm}_sr"
    local pipestance="$OUT_ROOT/${count_id}"
    local bam="$pipestance/outs/possorted_genome_bam.bam"
    local spatial="$pipestance/outs/spatial"
    local proc_name="gse243981_${gsm}_scapatrap"
    local sc_out_root="$OUT_ROOT/${proc_name}"
    local proc_dir="$ROOT/data/processed/${proc_name}"
    local sample_log="$LOG_ROOT/${DATE_TAG}_${gsm}_scapatrap.log"

    if [[ ! -f "$bam" ]]; then
        log "ERROR scAPAtrap ${gsm}: BAM missing ${bam}"
        return 1
    fi
    if [[ -f "$proc_dir/apa_sites.csv.gz" ]]; then
        log "SKIP scAPAtrap ${gsm}: processed output exists at ${proc_dir}"
        return 0
    fi

    log "START scAPAtrap ${gsm}"
    local t0 t1 rc
    t0=$(date +%s)
    "$PYTHON_BIN" "$LAUNCHER" \
        --bam "$bam" \
        --spatial-dir "$spatial" \
        --output-root "$sc_out_root" \
        --processed-dir "$proc_dir" \
        --dataset-name "$proc_name" \
        --tissue "human liver (healthy)" \
        --species human \
        --threads "$THREADS" \
        --force \
        > "$sample_log" 2>&1
    rc=$?
    t1=$(date +%s)
    log "END scAPAtrap ${gsm} (rc=${rc}, elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f min",(b-a)/60}'))"
    if [[ $rc -ne 0 ]]; then
        log "ERROR scAPAtrap ${gsm} failed; see ${sample_log}"
        tail -25 "$sample_log" >&2 || true
        return $rc
    fi
    return 0
}

# -----------------------------------------------------------------------------
# 6. Dispatch
# -----------------------------------------------------------------------------
MODE="${1:-all}"
log "GSE243981 healthy-liver pipeline: host=${HOSTNAME_SHORT} mode=${MODE} threads=${THREADS} mem=${MEM_GB}GB"

case "$MODE" in
    prep)
        for gsm in "${GSMS[@]}"; do prepare_links "$gsm" || exit 1; done
        ;;
    count)
        for gsm in "${GSMS[@]}"; do run_count "$gsm" || exit 1; done
        ;;
    scapatrap)
        for gsm in "${GSMS[@]}"; do run_scapatrap "$gsm" || exit 1; done
        ;;
    all)
        for gsm in "${GSMS[@]}"; do prepare_links "$gsm" || exit 1; done
        for gsm in "${GSMS[@]}"; do run_count "$gsm" || exit 1; done
        for gsm in "${GSMS[@]}"; do run_scapatrap "$gsm" || exit 1; done
        log "ALL DONE"
        ;;
    *)
        echo "Usage: $0 {prep|count|scapatrap|all}" >&2
        exit 1
        ;;
esac
log "mode ${MODE} finished"

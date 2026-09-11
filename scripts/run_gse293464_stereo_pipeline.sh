#!/usr/bin/env bash
# =============================================================================
# GSE293464 human retinal-organoid Stereo-seq — full processing orchestrator.
#
# 4 samples (H7 A81.BRN3b.mCherry, RA± × 16/26 wk), DNBSEQ-T7 PE, per-sample
# 8 SRA technical runs; masks deposited in GEO (see logs/20260909_stereo*.md):
#   GSM8882884  D02266B4  RA-/BMS1uM  26wk   SRR32936139..146
#   GSM8882885  D02266C2  RA1/10uM   26wk   SRR32936131..138
#   GSM8882886  D02266D2  RA-/BMS1uM  16wk   SRR32936123..130
#   GSM8882887  D02266D4  RA1/10uM   16wk   SRR32936115..122
#
# Phases (idempotent, run in order): fasterq | merge | saw | retag | scapatrap
# Usage: bash run_gse293464_stereo_pipeline.sh <phase>
# Run inside tmux on S91. Mirrors the proven GSE263789 chain.
# =============================================================================
set -uo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91) : ;;
    *) echo "ERROR: host $HOSTNAME != S91" >&2; exit 1 ;;
esac
export OPENBLAS_NUM_THREADS=64
export TMPDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse293464/tmp
mkdir -p "$TMPDIR"

DATA=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse293464
OUT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse293464_retina
mkdir -p "$OUT"
SAW=/s1/SHARE/01_software/saw-8.2.2/bin/saw
REF=/s1/SHARE/01_software/SAW_refs/Homo_sapiens_index
FASTERQ=/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump
SAMTOOLS=/home/mengzijun/anaconda3/envs/samtools/bin/samtools
RSCRIPT=/home/mengzijun/anaconda3/envs/r442/bin/Rscript
PY=/home/mengzijun/anaconda3/envs/spagapa/bin/python
RETAG=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/scripts/retag_scapatrap_v2.py
LOG="$OUT/pipeline.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

SAMPLES=(
    "s26_bms:D02266B4:GSM8882884:32936139"
    "s26_ra:D02266C2:GSM8882885:32936131"
    "s16_bms:D02266D2:GSM8882886:32936123"
    "s16_ra:D02266D4:GSM8882887:32936115"
)

PHASE="${1:?phase required: fasterq|merge|saw|retag|scapatrap}"

case "$PHASE" in
# ---- fasterq-dump: SRA -> split FASTQ (compressed) --------------------------
fasterq)
    for entry in "${SAMPLES[@]}"; do
        IFS=':' read -r key sn gsm first <<< "$entry"
        for i in 0 1 2 3 4 5 6 7; do
            srr="SRR$((first + i))"
            fqdir="$DATA/fasterq/$srr"
            [ -s "$fqdir/${srr}_1.fastq.gz" ] && [ -s "$fqdir/${srr}_2.fastq.gz" ] && { log "skip $srr (done)"; continue; }
            [ -s "$DATA/sra/$srr" ] || [ -s "$DATA/sra/$srr.sra" ] || { log "WARN: no sra file $srr"; continue; }
            [ -s "$DATA/sra/$srr.sra" ] || ln -sf "$DATA/sra/$srr" "$DATA/sra/$srr.sra"
            mkdir -p "$fqdir"
            log "fasterq-dump $srr ..."
            $FASTERQ -e 16 --split-files --temp "$TMPDIR" -O "$fqdir" "$DATA/sra/$srr.sra" \
                >> "$OUT/fasterq.log" 2>&1 \
              && pigz -f -p 8 "$fqdir/${srr}_1.fastq" "$fqdir/${srr}_2.fastq" \
              && rm -f "$fqdir/${srr}.fastq" \
              && log "  $srr OK $(du -sh "$fqdir" | cut -f1)" \
              || log "  $srr FAIL (see fasterq.log)"
        done
    done
    # read-structure sanity on first available run
    f=$(ls "$DATA"/fasterq/SRR*/SRR*_1.fastq.gz 2>/dev/null | head -1)
    [ -n "$f" ] && { log "R1 head: $(zcat "$f" | head -2 | tail -1 | wc -c) chars";
                     log "R2 head: $(zcat "${f/_1./_2.}" | head -2 | tail -1 | wc -c) chars"; }
    ;;

# ---- merge: per-sample single R1/R2 pair -----------------------------------
merge)
    for entry in "${SAMPLES[@]}"; do
        IFS=':' read -r key sn gsm first <<< "$entry"
        dst="$DATA/fastq_saw/$sn"; mkdir -p "$dst"
        [ -s "$dst/${sn}_1.fq.gz" ] && [ -s "$dst/${sn}_2.fq.gz" ] && { log "skip merge $sn"; continue; }
        r1s=""; r2s=""
        for i in 0 1 2 3 4 5 6 7; do
            srr="SRR$((first + i))"
            f1="$DATA/fasterq/$srr/${srr}_1.fastq.gz"; f2="$DATA/fasterq/$srr/${srr}_2.fastq.gz"
            [ -s "$f1" ] || { log "WARN missing $f1"; continue; }
            r1s="$r1s $f1"; r2s="$r2s $f2"
        done
        [ -n "$r1s" ] || { log "no fastqs for $sn"; continue; }
        log "merging $sn ($(echo $r1s | wc -w) runs)..."
        cat $r1s > "$dst/${sn}_1.fq.gz" && cat $r2s > "$dst/${sn}_2.fq.gz"
        log "  $sn merged: $(du -sh "$dst" | cut -f1)"
    done
    ;;

# ---- SAW count per sample ---------------------------------------------------
saw)
    for entry in "${SAMPLES[@]}"; do
        IFS=':' read -r key sn gsm first <<< "$entry"
        fqs="$DATA/fastq_saw/$sn"
        mask="$DATA/masks/$sn.barcodeToPos.h5"
        [ -s "$fqs/${sn}_1.fq.gz" ] && [ -s "$fqs/${sn}_2.fq.gz" ] || { log "skip saw $sn (no fastq)"; continue; }
        [ -s "$mask" ] || { log "skip saw $sn (no mask)"; continue; }
        [ -d "$OUT/${gsm}_saw" ] && { log "skip saw $sn (output exists)"; continue; }
        log "SAW count $key $sn ..."
        "$SAW" count \
            --id="${gsm}_saw" --sn="$sn" --omics=transcriptomics \
            --kit-version="Stereo-seq T FF V1.3" --sequencing-type="PE75_50+100" \
            --chip-mask="$mask" --organism=human --tissue=retina \
            --fastqs="$fqs" --reference="$REF" --output="$OUT" \
            --threads-num=32 --memory=200 \
            >> "$OUT/saw_${sn}.log" 2>&1 \
          && log "  SAW $sn OK" || log "  SAW $sn FAIL (saw_${sn}.log)"
    done
    ;;

# ---- retag + chr prefix + sort ---------------------------------------------
retag)
    for entry in "${SAMPLES[@]}"; do
        IFS=':' read -r key sn gsm first <<< "$entry"
        sawbam="$OUT/${gsm}_saw/outs/bam/annotated_bam/${sn}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam"
        [ -s "$sawbam" ] || sawbam="$OUT/${gsm}_saw/STEREO_ANALYSIS_WORKFLOW_PROCESSING/ANNOTATION/${sn}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam"
        rd="$OUT/${gsm}_retag"; mkdir -p "$rd"
        [ -s "$sawbam" ] || { log "skip retag $sn (no SAW bam)"; continue; }
        [ -s "$rd/scapatrap_input.retag.chr.bam" ] && { log "skip retag $sn (done)"; continue; }
        log "retag $sn ..."
        $PY "$RETAG" "$sawbam" "$rd/scapatrap_input.retag.bam"
        $SAMTOOLS view -H "$rd/scapatrap_input.retag.bam" | \
            sed 's/SN:\([0-9XY]\)/SN:chr\1/; s/SN:MT/SN:chrM/' > "$rd/chr_header.sam"
        $SAMTOOLS reheader "$rd/chr_header.sam" "$rd/scapatrap_input.retag.bam" > "$rd/scapatrap_input.retag.chr.bam"
        rm -f "$rd/scapatrap_input.retag.bam"
        $SAMTOOLS sort -@ 16 -o "$rd/tmp.sorted.bam" "$rd/scapatrap_input.retag.chr.bam"
        mv "$rd/tmp.sorted.bam" "$rd/scapatrap_input.retag.chr.bam"
        $SAMTOOLS index "$rd/scapatrap_input.retag.chr.bam"
        log "  retag $sn OK $(du -sh "$rd/scapatrap_input.retag.chr.bam" | cut -f1)"
    done
    ;;

# ---- scAPAtrap (human chrs) -------------------------------------------------
scapatrap)
    for entry in "${SAMPLES[@]}"; do
        IFS=':' read -r key sn gsm first <<< "$entry"
        rd="$OUT/${gsm}_retag"; raw="$OUT/${gsm}_scapatrap_raw"
        [ -s "$rd/scapatrap_input.retag.chr.bam" ] || { log "skip scapatrap $sn (no bam)"; continue; }
        [ -s "$raw/peaks_meta.csv.gz" ] && { log "skip scapatrap $sn (done)"; continue; }
        cat > "$rd/run_scapatrap.R" << REOF
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))
.skipwrap <- function(fname, checkfiles, retfun) {
    orig <- get(fname, envir = asNamespace("scAPAtrap"))
    fast <- function(...) {
        if (all(file.exists(checkfiles(...)))) {
            message(fname, ": prebuilt artifacts found, skipping")
            return(retfun(...))
        }
        res <- orig(...)
        if (fname == "dedupByPos") file.create(paste0(gsub(".bam$", "", list(...)\$input), ".dedup.bam.done"))
        res
    }
    assignInNamespace(fname, fast, ns = "scAPAtrap")
}
.fum_files <- function(samtools.path, input, ...) paste0(gsub(".bam$", "", input), c(".UniqSorted.bam", ".UniqSorted.bam.bai"))
.skipwrap("findUniqueMap", .fum_files, function(samtools.path, input, ...) paste0(gsub(".bam$", "", input), ".UniqSorted.bam"))
.dbp_files <- function(umitools.path, input, ...) c(paste0(gsub(".bam$", "", input), ".dedup.bam"), paste0(gsub(".bam$", "", input), ".dedup.bam.done"))
.skipwrap("dedupByPos", .dbp_files, function(umitools.path, input, ...) paste0(gsub(".bam$", "", input), ".dedup.bam"))
.sbs_files <- function(samtools.path, input, ...) paste0(gsub(".bam", "", input), c(".forward.bam", ".forward.bam.bai", ".reverse.bam", ".reverse.bam.bai"))
.skipwrap("separateBamBystrand", .sbs_files, function(samtools.path, input, ...) paste0(gsub(".bam", "", input), c(".forward.bam", ".reverse.bam")))
.fpb_files <- function(bamFile, ...) paste0(bamFile, ".peaks")
.skipwrap("findPeaksByStrand", .fpb_files, function(bamFile, ...) paste0(bamFile, ".peaks"))
tools <- list(
  samtools = "/home/mengzijun/anaconda3/envs/samtools/bin/samtools",
  umitools = "/home/mengzijun/anaconda3/envs/samtools/bin/umi_tools",
  featureCounts = "/home/mengzijun/anaconda3/envs/samtools/bin/featureCounts",
  star = "/home/mengzijun/anaconda3/envs/samtools/bin/STAR")
trap.params <- setTrapParams(print=FALSE)
trap.params\$TenX <- TRUE
trap.params\$chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                      "chr21","chr22","chrX","chrY")
trap.params\$readlength <- 100
trap.params\$cov.cutoff <- 10
trap.params\$min.cells <- 10
trap.params\$min.count <- 10
trap.params\$tails.search <- "peaks"
trap.params\$thread <- 32
input_bam <- "$rd/scapatrap_input.retag.chr.bam"
output_dir <- "$raw"
log_file <- "$rd/scapatrap_internal.log"
if (dir.exists(file.path(output_dir,"peaks_meta.csv.gz"))) stop("already done")
message("Running scAPAtrap $gsm ...")
rda <- scAPAtrap(tools=tools, trap.params=trap.params, inputBam=input_bam,
                 outputDir=output_dir, logf=log_file, verbose=TRUE)
load(rda)
write.csv(scAPAtrapData\$peaks.meta, gzfile(file.path(output_dir,"peaks_meta.csv.gz")), quote=FALSE)
counts <- as(scAPAtrapData\$peaks.count, "TsparseMatrix")
df <- data.frame(site=rownames(counts)[row(counts)], spot=colnames(counts)[col(counts)], count=as.vector(counts))
df <- df[df\$count > 0, ]
gz <- gzfile(file.path(output_dir,"apa_site_counts.csv.gz"), "w")
write.csv(df, gz, quote=FALSE, row.names=FALSE); close(gz)
qc <- list(n_sites=nrow(scAPAtrapData\$peaks.meta), n_barcodes=ncol(scAPAtrapData\$peaks.count))
writeLines(jsonlite::toJSON(qc, auto_unbox=TRUE, pretty=TRUE), file.path(output_dir,"scapatrap_qc.json"))
message("scAPAtrap $gsm export complete!")
REOF
        rm -rf "$raw" "$rd/scapatrap_internal.log"
        log "scAPAtrap $sn ..."
        $RSCRIPT "$rd/run_scapatrap.R" >> "$OUT/scapatrap_${sn}.log" 2>&1 \
          && log "  scAPAtrap $sn OK" || log "  scAPAtrap $sn FAIL (scapatrap_${sn}.log)"
    done
    ;;
*) echo "unknown phase $PHASE" >&2; exit 1 ;;
esac
log "PHASE $PHASE done"

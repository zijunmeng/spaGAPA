#!/usr/bin/env bash
# =============================================================================
# GSE333693 rat thymus Stereo-seq — processing orchestrator (third species).
#
# 1 DNBSEQ-T7 PE run: SRR38888148 (184 Gb), chip SN Y01052GC,
# mask deposited in GEO (GSM9770943_Y01052GC.barcodeToPos.h5).
#
# Phases: fasterq | saw | retag | scapatrap     Usage: <phase> (tmux, S91)
# Rat reference: /s1/SHARE/01_software/SAW_refs/Rattus_norvegicus_index
# (extracted from reference-data-rat.tar.gz — see acquisition script).
# =============================================================================
set -uo pipefail

HOSTNAME=$(hostname -s)
case "$HOSTNAME" in
    S91) : ;;
    *) echo "ERROR: host $HOSTNAME != S91" >&2; exit 1 ;;
esac
export OPENBLAS_NUM_THREADS=64
export TMPDIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse333693/tmp
mkdir -p "$TMPDIR"

DATA=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse333693
OUT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse333693_thymus
mkdir -p "$OUT"
SN=Y01052GC; GSM=GSM9770943; SRR=SRR38888148
SAW=/s1/SHARE/01_software/saw-8.2.2/bin/saw
REF=/s1/SHARE/01_software/SAW_refs/Rattus_norvegicus_index
FASTERQ=/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump
SAMTOOLS=/home/mengzijun/anaconda3/envs/samtools/bin/samtools
RSCRIPT=/home/mengzijun/anaconda3/envs/r442/bin/Rscript
PY=/home/mengzijun/anaconda3/envs/spagapa/bin/python
RETAG=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/scripts/retag_scapatrap_v2.py
LOG="$OUT/pipeline.log"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

PHASE="${1:?phase required: fasterq|saw|retag|scapatrap}"

case "$PHASE" in
fasterq)
    fqdir="$DATA/fasterq/$SRR"
    [ -s "$fqdir/${SRR}_1.fastq.gz" ] && { log "skip fasterq (done)"; exit 0; }
    [ -s "$DATA/sra/$SRR" ] || [ -s "$DATA/sra/$SRR.sra" ] || { log "no sra file yet"; exit 1; }
    [ -s "$DATA/sra/$SRR.sra" ] || ln -sf "$DATA/sra/$SRR" "$DATA/sra/$SRR.sra"
    mkdir -p "$fqdir"
    log "fasterq-dump $SRR ..."
    $FASTERQ -e 32 --split-files --temp "$TMPDIR" -O "$fqdir" "$DATA/sra/$SRR.sra" >> "$OUT/fasterq.log" 2>&1 \
      && pigz -f -p 16 "$fqdir/${SRR}_1.fastq" "$fqdir/${SRR}_2.fastq" && rm -f "$fqdir/${SRR}.fastq" \
      && log "  OK $(du -sh "$fqdir" | cut -f1)" \
      && log "R1: $(zcat "$fqdir/${SRR}_1.fastq.gz" | head -2 | tail -1 | wc -c) chars; R2: $(zcat "$fqdir/${SRR}_2.fastq.gz" | head -2 | tail -1 | wc -c) chars"
    # SAW layout: fastq_saw/<SN>_1.fq.gz / _2
    dst="$DATA/fastq_saw/$SN"; mkdir -p "$dst"
    [ -s "$dst/${SN}_1.fq.gz" ] || ln -sf "$fqdir/${SRR}_1.fastq.gz" "$dst/${SN}_1.fq.gz"
    [ -s "$dst/${SN}_2.fq.gz" ] || ln -sf "$fqdir/${SRR}_2.fastq.gz" "$dst/${SN}_2.fq.gz"
    ;;
saw)
    fqs="$DATA/fastq_saw/$SN"; mask="$DATA/masks/$SN.barcodeToPos.h5"
    [ -s "$fqs/${SN}_1.fq.gz" ] || { log "no fastq"; exit 1; }
    [ -s "$mask" ] || { log "no mask"; exit 1; }
    [ -d "$OUT/${GSM}_saw" ] && { log "skip saw (output exists)"; exit 0; }
    log "SAW count $SN (rat thymus) ..."
    "$SAW" count \
        --id="${GSM}_saw" --sn="$SN" --omics=transcriptomics \
        --kit-version="Stereo-seq T FF V1.3" --sequencing-type="PE75_50+100" \
        --chip-mask="$mask" --organism=rat --tissue=thymus \
        --fastqs="$fqs" --reference="$REF" --output="$OUT" \
        --threads-num=32 --memory=200 \
        >> "$OUT/saw_${SN}.log" 2>&1 \
      && log "SAW OK" || log "SAW FAIL (saw_${SN}.log)"
    ;;
retag)
    sawbam="$OUT/${GSM}_saw/outs/bam/annotated_bam/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam"
    [ -s "$sawbam" ] || sawbam="$OUT/${GSM}_saw/STEREO_ANALYSIS_WORKFLOW_PROCESSING/ANNOTATION/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam"
    rd="$OUT/${GSM}_retag"; mkdir -p "$rd"
    [ -s "$sawbam" ] || { log "no SAW bam"; exit 1; }
    [ -s "$rd/scapatrap_input.retag.chr.bam" ] && { log "skip retag (done)"; exit 0; }
    log "retag $SN ..."
    $PY "$RETAG" "$sawbam" "$rd/scapatrap_input.retag.bam"
    $SAMTOOLS view -H "$rd/scapatrap_input.retag.bam" | \
        sed 's/SN:\([0-9XY]\)/SN:chr\1/; s/SN:MT/SN:chrM/' > "$rd/chr_header.sam"
    $SAMTOOLS reheader "$rd/chr_header.sam" "$rd/scapatrap_input.retag.bam" > "$rd/scapatrap_input.retag.chr.bam"
    rm -f "$rd/scapatrap_input.retag.bam"
    $SAMTOOLS sort -@ 16 -o "$rd/tmp.sorted.bam" "$rd/scapatrap_input.retag.chr.bam"
    mv "$rd/tmp.sorted.bam" "$rd/scapatrap_input.retag.chr.bam"
    $SAMTOOLS index "$rd/scapatrap_input.retag.chr.bam"
    log "retag OK $(du -sh "$rd/scapatrap_input.retag.chr.bam" | cut -f1)"
    ;;
scapatrap)
    rd="$OUT/${GSM}_retag"; raw="$OUT/${GSM}_scapatrap_raw"
    [ -s "$rd/scapatrap_input.retag.chr.bam" ] || { log "no bam"; exit 1; }
    [ -s "$raw/peaks_meta.csv.gz" ] && { log "skip (done)"; exit 0; }
    cat > "$rd/run_scapatrap.R" << REOF
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))
.fum_orig <- scAPAtrap:::findUniqueMap
.fum_fast <- function(samtools.path, input, thread = 12, sort = TRUE, index = TRUE, ...) {
    out.sorted <- paste0(gsub(".bam$", "", input), ".UniqSorted.bam")
    if (file.exists(out.sorted) && file.exists(paste0(out.sorted, ".bai"))) {
        message("findUniqueMap: prebuilt sorted+indexed bam found, skipping: ", out.sorted)
        return(out.sorted)
    }
    .fum_orig(samtools.path = samtools.path, input = input, thread = thread, sort = sort, index = index, ...)
}
assignInNamespace("findUniqueMap", .fum_fast, ns = "scAPAtrap")
.dbp_orig <- scAPAtrap:::dedupByPos
.dbp_fast <- function(umitools.path, input, TenX = TRUE, ...) {
    out <- paste0(gsub(".bam$", "", input), ".dedup.bam")
    if (file.exists(out) && file.exists(paste0(out, ".done"))) {
        message("dedupByPos: prebuilt dedup bam found, skipping: ", out)
        return(out)
    }
    res <- .dbp_orig(umitools.path = umitools.path, input = input, TenX = TenX, ...)
    file.create(paste0(out, ".done"))
    res
}
assignInNamespace("dedupByPos", .dbp_fast, ns = "scAPAtrap")
.sbs_orig <- scAPAtrap:::separateBamBystrand
.sbs_fast <- function(samtools.path, input, thread = 12, ...) {
    outF <- paste0(gsub(".bam", "", input), ".forward.bam")
    outR <- paste0(gsub(".bam", "", input), ".reverse.bam")
    if (file.exists(outF) && file.exists(paste0(outF, ".bai")) &&
        file.exists(outR) && file.exists(paste0(outR, ".bai"))) {
        message("separateBamBystrand: prebuilt strand bams found, skipping")
        return(c(outF, outR))
    }
    .sbs_orig(samtools.path = samtools.path, input = input, thread = thread, ...)
}
assignInNamespace("separateBamBystrand", .sbs_fast, ns = "scAPAtrap")
.fpb_orig <- scAPAtrap:::findPeaksByStrand
.fpb_fast <- function(bamFile, chrs = NULL, strand, L, maxwidth, cutoff = 10, ofile = NULL, ...) {
    pf <- paste0(bamFile, ".peaks")
    if (file.exists(pf) && file.size(pf) > 0) {
        message("findPeaksByStrand: prebuilt peaks file found, skipping: ", pf)
        return(pf)
    }
    .fpb_orig(bamFile = bamFile, chrs = chrs, strand = strand, L = L, maxwidth = maxwidth, cutoff = cutoff, ofile = ofile, ...)
}
assignInNamespace("findPeaksByStrand", .fpb_fast, ns = "scAPAtrap")
tools <- list(
  samtools = "/home/mengzijun/anaconda3/envs/samtools/bin/samtools",
  umitools = "/home/mengzijun/anaconda3/envs/samtools/bin/umi_tools",
  featureCounts = "/home/mengzijun/anaconda3/envs/samtools/bin/featureCounts",
  star = "/home/mengzijun/anaconda3/envs/samtools/bin/STAR")
trap.params <- setTrapParams(print=FALSE)
trap.params\$TenX <- TRUE
hdr_targets <- names(Rsamtools::scanBamHeader("$rd/scapatrap_input.retag.chr.bam")[[1]]\$targets)
trap.params\$chrs <- hdr_targets[grepl("^NC_", hdr_targets)]
message("rat chrs (from BAM header): ", paste(trap.params\$chrs, collapse=","))
trap.params\$readlength <- 100
trap.params\$cov.cutoff <- 10
trap.params\$min.cells <- 10
trap.params\$min.count <- 10
trap.params\$tails.search <- "peaks"
trap.params\$thread <- 32
rda <- scAPAtrap(tools=tools, trap.params=trap.params,
                 inputBam="$rd/scapatrap_input.retag.chr.bam",
                 outputDir="$raw", logf="$rd/scapatrap_internal.log", verbose=TRUE)
load(rda)
write.csv(scAPAtrapData\$peaks.meta, gzfile(file.path("$raw","peaks_meta.csv.gz")), quote=FALSE)
counts <- as(scAPAtrapData\$peaks.count, "TsparseMatrix")
df <- data.frame(site=rownames(counts)[row(counts)], spot=colnames(counts)[col(counts)], count=as.vector(counts))
df <- df[df\$count > 0, ]
gz <- gzfile(file.path("$raw","apa_site_counts.csv.gz"), "w")
write.csv(df, gz, quote=FALSE, row.names=FALSE); close(gz)
message("scAPAtrap rat export complete!")
REOF
    rm -rf "$raw" "$rd/scapatrap_internal.log"
    log "scAPAtrap $SN ..."
    $RSCRIPT "$rd/run_scapatrap.R" >> "$OUT/scapatrap_${SN}.log" 2>&1 \
      && log "scAPAtrap OK" || log "scAPAtrap FAIL"
    ;;
*) echo "unknown phase $PHASE" >&2; exit 1 ;;
esac
log "PHASE $PHASE done"

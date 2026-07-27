#!/usr/bin/env bash
# WT Stereo-seq: retag SAW BAM → chr prefix → sort → scAPAtrap (mouse chrs)
# Mirrors the AD pilot (GSM8199179) retag pipeline for WT control (GSM8199181).
set -uo pipefail
export OPENBLAS_NUM_THREADS=64 TMPDIR=/s3/mengzijun/tmp
SAMTOOLS=/home/mengzijun/anaconda3/envs/samtools/bin/samtools
RSCRIPT=/home/mengzijun/anaconda3/envs/r442/bin/Rscript
PY=/home/mengzijun/anaconda3/envs/spagapa/bin/python

WT_DIR=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control
SAW_BAM="$WT_DIR/wt_control_F5/STEREO_ANALYSIS_WORKFLOW_PROCESSING/ANNOTATION/SS200000745BL.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam"
OUT_DIR="$WT_DIR/scapatrap_pipeline"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[$(date '+%T')] === WT retag+scAPAtrap pipeline ==="

# Step 1: retag (CB=Cx_Cy, UB=UR padded, GX=GI, GN=GS)
echo "[$(date '+%T')] STEP 1: retag..."
$PY /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_stereo_pilot/gsm8199179_full/retag_scapatrap.py \
  "$SAW_BAM" "$OUT_DIR/scapatrap_input.retag.bam"

# Step 2: chr prefix (samtools reheader — instant, only header)
echo "[$(date '+%T')] STEP 2: chr prefix..."
$SAMTOOLS view -H "$OUT_DIR/scapatrap_input.retag.bam" | \
  sed 's/SN:\([0-9XY]\)/SN:chr\1/; s/SN:MT/SN:chrM/' > "$OUT_DIR/chr_header.sam"
$SAMTOOLS reheader "$OUT_DIR/chr_header.sam" "$OUT_DIR/scapatrap_input.retag.bam" > "$OUT_DIR/scapatrap_input.retag.chr.bam"

# Step 3: sort + index
echo "[$(date '+%T')] STEP 3: sort+index..."
$SAMTOOLS sort -@ 16 -o "$OUT_DIR/scapatrap_input.retag.chr.sorted.bam" "$OUT_DIR/scapatrap_input.retag.chr.bam"
$SAMTOOLS index "$OUT_DIR/scapatrap_input.retag.chr.sorted.bam"
mv "$OUT_DIR/scapatrap_input.retag.chr.sorted.bam" "$OUT_DIR/scapatrap_input.retag.chr.bam"
$SAMTOOLS index "$OUT_DIR/scapatrap_input.retag.chr.bam"

# Step 4: scAPAtrap (mouse chrs chr1-19+XY)
echo "[$(date '+%T')] STEP 4: scAPAtrap (mouse)..."
cat > "$OUT_DIR/run_scapatrap.R" << 'REOF'
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))
tools <- list(
  samtools = "/home/mengzijun/anaconda3/envs/samtools/bin/samtools",
  umitools = "/home/mengzijun/anaconda3/envs/samtools/bin/umi_tools",
  featureCounts = "/home/mengzijun/anaconda3/envs/samtools/bin/featureCounts",
  star = "/home/mengzijun/anaconda3/envs/samtools/bin/STAR"
)
trap.params <- setTrapParams(print=FALSE)
trap.params$TenX <- FALSE
trap.params$chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
trap.params$readlength <- 100
trap.params$cov.cutoff <- 10
trap.params$min.cells <- 10
trap.params$min.count <- 10
trap.params$tails.search <- "peaks"
trap.params$thread <- 16
input_bam <- "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control/scapatrap_pipeline/scapatrap_input.retag.chr.bam"
output_dir <- "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control/scapatrap_raw"
log_file <- "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control/scapatrap_pipeline/scapatrap_internal.log"
if (dir.exists(output_dir)) stop(paste("output_dir exists:", output_dir))
message("Running scAPAtrap WT...")
scapatrap_rda <- scAPAtrap(tools=tools, trap.params=trap.params, inputBam=input_bam,
                           outputDir=output_dir, logf=log_file, verbose=TRUE)
load(scapatrap_rda)
write.csv(scAPAtrapData$peaks.meta, gzfile(file.path(output_dir, "peaks_meta.csv.gz")), quote=FALSE)
# sparse long format export (avoid OOM on dense)
counts <- as(scAPAtrapData$peaks.count, "TsparseMatrix")
df <- data.frame(site=rownames(counts)[row(counts)], spot=colnames(counts)[col(counts)], count=as.vector(counts))
df <- df[df$count > 0, ]
gz <- gzfile(file.path(output_dir, "apa_site_counts.csv.gz"), "w")
write.csv(df, gz, quote=FALSE, row.names=FALSE)
close(gz)
qc <- list(n_sites=nrow(scAPAtrapData$peaks.meta), n_barcodes=ncol(scAPAtrapData$peaks.count),
          tails_search="peaks", readlength=100)
writeLines(jsonlite::toJSON(qc, auto_unbox=TRUE, pretty=TRUE), file.path(output_dir, "scapatrap_qc.json"))
message("scAPAtrap WT export complete!")
REOF

$RSCRIPT "$OUT_DIR/run_scapatrap.R" 2>&1 | tee "$OUT_DIR/scapatrap_run.log"
echo "[$(date '+%T')] === WT PIPELINE COMPLETE ==="

#!/usr/bin/env Rscript
# Sierra (v0.99.27, Winnie09/Sierra) peak calling + per-spot UMI counting on a
# 10x Visium Space Ranger BAM. Second-PAS-caller robustness validation for
# spaGAPA. Datasets: GSE183456 (human kidney) / GSE220442 (human AD brain) /
# GSE169749 (mouse colon); paths come from the driver (run_dataset.sh).
#
# Usage: Rscript run_sierra.R <outdir> <bam> <gtf> <junctions.bed> <whitelist.tsv> [ncores]
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressMessages(library(Sierra))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5)
OUT     <- args[1]
BAM     <- args[2]
GTF     <- args[3]
JUNC    <- args[4]
WL      <- args[5]
NCORES  <- as.integer(ifelse(length(args) >= 6, args[6], 16))

message("Sierra ", as.character(packageVersion("Sierra")))
message("BAM: ", BAM)
message("whitelist: ", length(readLines(WL)), " barcodes")
message("junctions: ", length(readLines(JUNC)), " entries")

# ---------------------------- 1. peak calling -------------------------------
peak.file <- file.path(OUT, "sierra_peaks.txt")
if (!file.exists(peak.file) || file.size(peak.file) < 100) {
  t0 <- Sys.time()
  FindPeaks(output.file    = peak.file,
            gtf.file       = GTF,
            bamfile        = BAM,
            junctions.file = JUNC,
            ncores         = NCORES)
  message("FindPeaks done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

peaks <- read.table(peak.file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
message("raw peak rows: ", nrow(peaks))
good <- peaks[!is.na(peaks$Fit.start) & !is.na(peaks$Fit.end) & peaks$Fit.start < peaks$Fit.end, ]
message("usable peaks (Fit.start < Fit.end): ", nrow(good))
message("genes with peaks: ", length(unique(good$Gene)))
print(table(good$`exon.intron`))

# ---------------------------- 2. UMI counting -------------------------------
count.dir <- file.path(OUT, "sierra_counts")
if (!file.exists(file.path(count.dir, "matrix.mtx.gz"))) {
  t0 <- Sys.time()
  CountPeaks(peak.sites.file = peak.file,
             gtf.file        = GTF,
             bamfile         = BAM,
             whitelist.file  = WL,
             output.dir      = count.dir,
             countUMI        = TRUE,
             ncores          = NCORES)
  message("CountPeaks done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
message("ALL DONE")

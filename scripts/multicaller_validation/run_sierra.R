#!/usr/bin/env Rscript
# Sierra (v0.99.27, Winnie09/Sierra) peak calling + per-spot UMI counting on the
# GSE183456 / GSM6047774 human-kidney Visium Space Ranger BAM.
# Second-PAS-caller robustness validation for spaGAPA.
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressMessages(library(Sierra))

BASE <- "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT  <- file.path(BASE, "pipeline_output/multicaller_validation/sierra")
BAM  <- file.path(BASE, "pipeline_output/gse183456_GSM6047774_sr/outs/possorted_genome_bam.bam")
GTF  <- file.path(OUT, "genes.gtf")  # decompressed from refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
JUNC <- file.path(OUT, "junctions.bed")
WL   <- file.path(OUT, "whitelist_barcodes.tsv")
NCORES <- 16

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

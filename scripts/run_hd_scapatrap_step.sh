#!/usr/bin/env bash
# HD scAPAtrap 步骤：spaceranger BAM → PAS（复用冻结集 Visium 模板 + stereo 版 monkey-patch 守卫）
# 用法: run_hd_scapatrap_step.sh <sample_out_dir>
set -euo pipefail
OUT=$1
BAM="$OUT/sr/outs/possorted_genome_bam.bam"
RAW="$OUT/scapatrap_raw"
[ -s "$BAM" ] || { echo "缺 $BAM"; exit 1; }

# HD: 合并 reference 染色体带 GRCh38_/GRCm39_ 前缀 → samtools view 提取真实染色体名集合
CHRS=$(samtools view -H "$BAM" | grep -oE "SN:[A-Za-z0-9_.-]+" | cut -d: -f2 | grep -vE "_" | head -40)
# （若前缀未剥离则含 _；此时用前缀全名）
NCHR=$(echo "$CHRS" | wc -l)
if [ "$NCHR" -lt 20 ]; then
  CHRS=$(samtools view -H "$BAM" | grep -oE "SN:[A-Za-z0-9_.-]+" | cut -d: -f2 | sed 's/^/"/;s/$/"/' | paste -sd, -)
else
  CHRS=$(echo "$CHRS" | sed 's/^/"/;s/$/"/' | paste -sd, -)
fi

cat > "$OUT/run_hd_scapatrap.R" <<REOF
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))

bamFile <- "$BAM"
trap.params <- setTrapParams(print = FALSE)
trap.params\$TenX <- TRUE
trap.params\$chrs <- c($CHRS)
trap.params\$readlength <- 50      # HD R2 = 50bp
trap.params\$cov.cutoff <- 10
trap.params\$min.cells <- 10
trap.params\$min.count <- 10
trap.params\$tails.search <- "peaks"
trap.params\$thread <- 16

setwd("$RAW")
scAPAtrapData <- runScAPAtrap(
  bamFile = bamFile,
  outdir = ".",
  trap.params = trap.params,
  readlength = trap.params\$readlength,
  cov_cutoff = trap.params\$cov.cutoff,
  min_cells = trap.params\$min.cells,
  peakFlankBonus = 0,
  format = "TENX"
)

write.csv(scAPAtrapData\$peaks.meta, gzfile("peaks_meta.csv.gz"), quote = FALSE)
counts <- as(scAPAtrapData\$peaks.count, "TsparseMatrix")
nn <- length(counts@x)
gz <- gzfile("apa_site_counts.csv.gz", "w")
writeLines("peak_id,spot_id,count", gz)
rn <- rownames(counts); cn <- colnames(counts)
stp <- 5e6
for (i in seq(1, nn, by = stp)) {
  j <- min(i + stp - 1, nn)
  writeLines(paste0(rn[counts@i + 1], ",", cn[counts@j + 1], ",", counts@x[i:j]), gz)
}
close(gz)
cat("[hd-scapatrap] done\n")
REOF

/home/mengzijun/anaconda3/envs/r442/bin/Rscript "$OUT/run_hd_scapatrap.R" > "$OUT/scapatrap.log" 2>&1

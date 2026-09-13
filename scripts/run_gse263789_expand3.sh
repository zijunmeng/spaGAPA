#!/usr/bin/env bash
# =============================================================================
# GSE263789 扩展三切片（AD18_E4 + 3M_E1/E2）Phase 2 管线编排
# 复用全链：SAW(mouse ref + per-slice mask) → retag v2 → scAPAtrap(TenX) →
# bin200(向量化) → conformal+spaGAPA 下游
# 前置：Phase1 下载完成（expand3_download.sh 的 ALL DONE 标记）
# 串行三切片，避免内存叠加。
# =============================================================================
set -uo pipefail
HOSTNAME=$(hostname -s); [ "$HOSTNAME" = "S91" ] || { echo "not S91" >&2; exit 1; }
export OPENBLAS_NUM_THREADS=8
ROOT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
RAW=/s1/SHARE/mengzijun/01_project/26_spaGAPA/data/raw/gse263789/expand_3slices
OUT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_expand
LOG=$ROOT/logs/20260913_expand3_pipeline.log
ts(){ date '+%F %T'; }; log(){ echo "[$(ts)] $*" | tee -a "$LOG"; }

SAW=/s1/SHARE/01_software/saw-8.2.2/bin/saw
REF=/s1/SHARE/01_software/SAW_refs/Mus_musculus_index
FASTERQ=/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump
SAMTOOLS=/home/mengzijun/anaconda3/envs/samtools/bin/samtools
RSCRIPT=/home/mengzijun/anaconda3/envs/r442/bin/Rscript
PY=/home/mengzijun/anaconda3/envs/spagapa/bin/python
RETAG=$ROOT/scripts/retag_scapatrap_v2.py
TMPD=/s3/mengzijun/tmp

# 切片清单: key:GSM:SRR:mask尾缀:样本标签
SLICES=(
  "ad18_e4:GSM8199180:SRR28637903:E4:ad18_E4"
  "3m_e1:GSM8199182:SRR28637889:E1:3m_E1"
  "3m_e2:GSM8199183:SRR28637880:E2:3m_E2"
)

run_slice(){
  local key=$1 gsm=$2 srr=$3 mtag=$4 label=$5
  local d=$OUT/$key; mkdir -p "$d"
  local sawdir=$d/saw retag=$d/retag raw=$d/scapatrap_raw binned=$d/binned_200
  local mask=$RAW/masks/${gsm}_SS200000769BR_${mtag}.barcodeToPos.h5

  # ---- 1. fasterq + 布局 ----
  if [ ! -s "$d/fastq/${srr}_1.fastq.gz" ] || [ ! -s "$d/fastq/${srr}_2.fastq.gz" ]; then
    log "$key [1/6] fasterq"
    mkdir -p "$d/fastq"
    [ -s "$RAW/sra/$srr.sra" ] || ln -sf "$RAW/sra/$srr" "$RAW/sra/$srr.sra"
    $FASTERQ -e 24 --split-files --temp "$TMPD" -O "$d/fastq" "$RAW/sra/$srr.sra" >> "$LOG" 2>&1 \
      && pigz -f -p 8 "$d/fastq/${srr}_1.fastq" "$d/fastq/${srr}_2.fastq" >> "$LOG" 2>&1 \
      || { log "$key fasterq FAIL"; return 1; }
    mkdir -p "$d/fastq_saw"
    ln -sf "$d/fastq/${srr}_1.fastq.gz" "$d/fastq_saw/SS200000769BR_1.fq.gz"
    ln -sf "$d/fastq/${srr}_2.fastq.gz" "$d/fastq_saw/SS200000769BR_2.fq.gz"
  fi

  # ---- 2. SAW count ----
  if [ ! -d "$sawdir" ]; then
    log "$key [2/6] SAW"
    "$SAW" count --id="${gsm}_saw" --sn=SS200000769BR --omics=transcriptomics \
      --kit-version="Stereo-seq T FF V1.3" --sequencing-type="PE75_50+100" \
      --chip-mask="$mask" --organism=mouse --tissue=brain \
      --fastqs="$d/fastq_saw" --reference="$REF" --output="$d" \
      --threads-num=32 --memory=200 >> "$d/saw.log" 2>&1 \
      && log "$key SAW OK" || { log "$key SAW FAIL"; return 1; }
    sawdir=$d/${gsm}_saw
  fi
  [ -d "$sawdir" ] || sawdir=$d/${gsm}_saw
  local sawbam
  sawbam=$(find "$sawdir" -name "*.target.bam" 2>/dev/null | head -1)
  [ -n "$sawbam" ] || { log "$key no SAW bam"; return 1; }

  # ---- 3. retag v2 + chr 前缀 + 排序（本地盘）----
  if [ ! -s "$retag/scapatrap_input.retag.chr.bam" ]; then
    log "$key [3/6] retag"
    mkdir -p "$retag"
    $PY "$RETAG" "$sawbam" "$retag/scapatrap_input.retag.bam" >> "$LOG" 2>&1
    $SAMTOOLS view -H "$retag/scapatrap_input.retag.bam" | \
      sed 's/SN:\([0-9XY]\)/SN:chr\1/; s/SN:MT/SN:chrM/' > "$retag/chr_header.sam"
    $SAMTOOLS reheader "$retag/chr_header.sam" "$retag/scapatrap_input.retag.bam" > "$retag/tmp.bam"
    rm -f "$retag/scapatrap_input.retag.bam"
    $SAMTOOLS sort -@ 16 -T "$TMPD/e3s_$key" -o "$retag/tmp2.bam" "$retag/tmp.bam"
    mv "$retag/tmp2.bam" "$retag/scapatrap_input.retag.chr.bam"; rm -f "$retag/tmp.bam"
    $SAMTOOLS index -@ 16 "$retag/scapatrap_input.retag.chr.bam"
  fi

  # ---- 4. scAPAtrap（本地排序/去重 + TenX + 流式导出，全套跳过守卫）----
  if [ ! -s "$raw/peaks_meta.csv.gz" ]; then
    log "$key [4/6] scAPAtrap"
    mkdir -p "$raw"
    cat > "$retag/run_scapatrap.R" << REOF
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))
$(sed -n '/^\.fum_orig/,/^assignInNamespace("findPeaksByStrand"/p' $ROOT/scripts/run_gse333693_stereo_pipeline.sh | sed 's/\\\$/\$/g')
.rp_orig <- scAPAtrap:::reducePeaks
assignInNamespace("reducePeaks", function(countsfile, peaksfile = NULL, min.cells = 10, min.count = 10,
                                          max.cells = NULL, max.count = NULL, suffix = ".reduced", toSparse = FALSE, ...) {
    if (is.character(countsfile) && length(countsfile) == 1L && !is.na(countsfile) &&
        is.character(peaksfile) && length(peaksfile) == 1L && !is.na(peaksfile)) {
        rp <- paste0(peaksfile, suffix); rc <- paste0(countsfile, suffix)
        if (file.exists(rp) && file.exists(rc)) { message("reducePeaks: skipping"); return(list(countsfile = rc, peaksfile = rp)) }
    }
    .rp_orig(countsfile = countsfile, peaksfile = peaksfile, min.cells = min.cells, min.count = min.count,
             max.cells = max.cells, max.count = max.count, suffix = suffix, toSparse = toSparse, ...)
}, ns = "scAPAtrap")
.gfb_orig <- scAPAtrap:::generateFinalBam
assignInNamespace("generateFinalBam", function(featureCounts.path, samtools.path, input, peakfile, thread = 12, ...) {
    fb <- file.path(dirname(input), "final.bam")
    if (file.exists(fb) && file.exists(paste0(fb, ".bai"))) { message("generateFinalBam: skipping"); return(fb) }
    .gfb_orig(featureCounts.path=featureCounts.path, samtools.path=samtools.path, input=input, peakfile=peakfile, thread=thread, ...)
}, ns = "scAPAtrap")
.ftb_orig <- scAPAtrap:::findTailsByPeaks
assignInNamespace("findTailsByPeaks", function(bamfile, peaksfile, d = 200, tailsfile = NULL, ...) {
    tf <- paste0(bamfile, ".peaks.tails")
    if (is.character(bamfile) && file.exists(tf) && file.size(tf) > 0) { message("findTailsByPeaks: skipping"); return(tf) }
    .ftb_orig(bamfile=bamfile, peaksfile=peaksfile, d=d, tailsfile=tailsfile, ...)
}, ns = "scAPAtrap")
tools <- list(
  samtools = "$SAMTOOLS", umitools = "/home/mengzijun/anaconda3/envs/samtools/bin/umi_tools",
  featureCounts = "/home/mengzijun/anaconda3/envs/samtools/bin/featureCounts",
  star = "/home/mengzijun/anaconda3/envs/samtools/bin/STAR")
trap.params <- setTrapParams(print=FALSE)
trap.params\$TenX <- TRUE
hdr_targets <- names(Rsamtools::scanBamHeader("$retag/scapatrap_input.retag.chr.bam")[[1]]\$targets)
trap.params\$chrs <- grep("^chr[0-9XY]+$", hdr_targets, value=TRUE)
trap.params\$readlength <- 100; trap.params\$cov.cutoff <- 10
trap.params\$min.cells <- 10; trap.params\$min.count <- 10
trap.params\$tails.search <- "peaks"; trap.params\$thread <- 32
rda <- scAPAtrap(tools=tools, trap.params=trap.params,
                 inputBam="$retag/scapatrap_input.retag.chr.bam", outputDir="$raw",
                 logf="$retag/scapatrap_internal.log", verbose=TRUE)
load(rda)
write.csv(scAPAtrapData\$peaks.meta, gzfile(file.path("$raw","peaks_meta.csv.gz")), quote=FALSE)
counts <- as(scAPAtrapData\$peaks.count, "TsparseMatrix")
nn <- length(counts@x)
gz <- gzfile(file.path("$raw","apa_site_counts.csv.gz"), "w")
writeLines("peak_id,spot_id,count", gz)
rn <- rownames(counts); cn <- colnames(counts)
stp <- 5e6
for (i in seq(1, nn, by = stp)) {
  j <- min(i + stp - 1, nn)
  writeLines(paste0(rn[counts@i + 1], ",", cn[counts@j + 1], ",", counts@x[i:j]), gz)
}
close(gz)
message("expand slice export complete: $label")
REOF
    # 本地盘预排序去重链（防 NFS 自锁）
    U="$retag/scapatrap_input.retag.chr.UniqSorted.bam"
    if [ ! -s "$U" ]; then
      $SAMTOOLS view -@ 24 -h -F 256 -bS "$retag/scapatrap_input.retag.chr.bam" > "$retag/u.bam"
      $SAMTOOLS sort -@ 24 -T "$TMPD/u_$key" -o "$U" "$retag/u.bam"; rm -f "$retag/u.bam"
      $SAMTOOLS index -@ 24 "$U"
    fi
    DD="$retag/scapatrap_input.retag.chr.UniqSorted.dedup.bam"
    if [ ! -s "$DD" ]; then
      /home/mengzijun/anaconda3/envs/samtools/bin/umi_tools dedup -I "$U" -S "$DD" \
        --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB >> "$LOG" 2>&1
    fi
    UQ="$retag/scapatrap_input.retag.chr.UniqSorted.dedup.sorted.bam"
    [ -s "$UQ" ] || { $SAMTOOLS sort -@ 24 -T "$TMPD/dd_$key" -o "$UQ.tmp" "$DD" && mv "$UQ.tmp" "$UQ"; }
    $SAMTOOLS index -@ 24 "$UQ"
    for st in forward reverse; do
      S="$retag/scapatrap_input.retag.chr.UniqSorted.dedup.$st.bam"
      [ -s "$S" ] || $SAMTOOLS view -@ 24 -h $([ $st = forward ] && echo "-F 0x10" || echo "-f 0x10") -bS "$UQ" > "$S"
      $SAMTOOLS index -@ 24 "$S"
    done
    # 用预处理好的 bam 覆盖主输入（findUniqueMap 跳过链将命中 UniqSorted+ bai）
    rm -f "$retag/scapatrap_input.retag.chr.UniqSorted.bam.tmp"*
    $RSCRIPT "$retag/run_scapatrap.R" >> "$d/scapatrap.log" 2>&1 \
      && log "$key scAPAtrap OK" || log "$key scAPAtrap FAIL"
  fi

  # ---- 5. bin200 ----
  if [ ! -s "$binned/qc_summary.json" ]; then
    log "$key [5/6] bin200"
    $PY "$ROOT/scripts/stereo_bin200.py" --raw-dir "$raw" --out-dir "$binned" \
      --dataset gse263789_expand --sample "$gsm" --bin 200 >> "$LOG" 2>&1 \
      && log "$key bin200 OK" || log "$key bin200 FAIL"
  fi

  # ---- 6. 下游 ----
  if [ ! -s "$OUT/downstream/$gsm/uncertainty/uncertainty_calibration.json" ]; then
    log "$key [6/6] downstream"
    mkdir -p "$OUT/downstream/$gsm"
    $PY "$ROOT/scripts/calibrate_uncertainty.py" \
      --apa-matrix "$binned/apa_matrix.csv" --coordinates "$binned/coordinates.csv" \
      --output "$OUT/downstream/$gsm/uncertainty" >> "$LOG" 2>&1 \
      && log "$key conformal OK" || log "$key conformal FAIL"
  fi
  log "$key SLICE COMPLETE"
}

log "=== expand3 pipeline start ==="
for s in "${SLICES[@]}"; do IFS=':' read -r k g srr mt lb <<< "$s"; run_slice "$k" "$g" "$srr" "$mt" "$lb"; done
log "=== expand3 pipeline ALL DONE ==="

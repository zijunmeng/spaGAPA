#!/usr/bin/env bash
# Visium HD (fresh-frozen, poly-dT) 端到端管线：
#   ENA FASTQ → spaceranger 4.1.0 (HD 自检) → scAPAtrap PAS → 空间聚合 → conformal
# 输入目录约定: data/raw/visium_hd/<study_dir>/<SRX>_SRR*_1/2.fastq.gz
# 输出:       pipeline_output/visium_hd/<study>/<sample>/
#
# 跳过守卫：每阶段产物存在即跳过（与 expand3 相同模式）
set -euo pipefail
ROOT=/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
SR=$ROOT/pipeline_output/visium_hd
SPACERANGER=/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger
REF_HUMAN=/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38_and_GRCm39-2024-A   # 待解压
PY=~/anaconda3/envs/spagapa/bin/python
OPENBLAS_NUM_THREADS=8
mkdir -p "$SR"

usage() { echo "用法: $0 <raw_dir> <study> <sample_glob> <species:human|mouse>"; exit 1; }
[ $# -eq 4 ] || usage
RAW=$1; STUDY=$2; GLOB=$3; SPECIES=$4

ref_for() { case $1 in human) echo $REF_HUMAN ;; mouse) echo $REF_HUMAN ;; *) echo "?" ;; esac; }
# 注：合并 reference 内含两物种；spaceranger 会按 FASTQ 自动适配？
# 否——spaceranger 需单物种 ref。若合并包结构为顶层双物种目录则需 mkref 拆分，
# 首次运行前先检查 $REF_HUMAN 结构（脚本 phase 0 会验证）。

for fq1 in "$RAW"/$GLOB; do
  [ -e "$fq1" ] || continue
  base=$(basename "$fq1" _1.fastq.gz)
  fq2="$RAW/${base}_2.fastq.gz"
  [ -e "$fq2" ] || { echo "[skip] $base: 无 R2"; continue; }
  out=$SR/$STUDY/$base
  mkdir -p "$out"

  # [1/4] spaceranger (HD 自检; 未知 slide 用 --unknown-slide; 16 线程)
  if [ ! -s "$out/sr/outs/possorted_genome_bam.bam" ]; then
    echo "[$base] spaceranger 启动 $(date '+%T')"
    $SPACERANGER count --id=sr --transcriptome="$REF" \
      --fastqs="$RAW" --sample="$(echo $base | sed 's/^SRX[0-9]*_//')" \
      --unknown-slide --localcores=16 --localmem=64 \
      --create-bam=true > "$out/spaceranger.log" 2>&1 \
      && mv sr "$out/" 2>/dev/null || true
  fi

  # [2/4] scAPAtrap（复用 Visium R 模式: 直接吃 spaceranger BAM）
  if [ ! -s "$out/scapatrap_raw/peaks_meta.csv.gz" ]; then
    echo "[$base] scAPAtrap 启动 $(date '+%T')"
    mkdir -p "$out/scapatrap_raw"
    # R heredoc 与冻结集同构（从 gse237183 runner 模板复制；monkey-patch 守卫同 stereo 版）
    bash scripts/run_hd_scapatrap_step.sh "$out" || echo "[warn] $base scAPAtrap 失败，见 $out/scapatrap.log"
  fi

  # [3/4] 空间聚合（HD 2µm 网格 → 8µm bin = 4x4 聚合，对齐 bin200 语义）
  if [ ! -s "$out/binned_8um/qc_summary.json" ]; then
    $PY scripts/hd_bin8um.py --raw-dir "$out/scapatrap_raw" \
      --out-dir "$out/binned_8um" --um 8 2>&1 | tail -3
  fi

  # [4/4] conformal
  if [ ! -s "$out/uncertainty/uncertainty_calibration.json" ]; then
    $PY scripts/calibrate_uncertainty.py \
      --apa-matrix "$out/binned_8um/apa_matrix.csv" \
      --coordinates "$out/binned_8um/coordinates.csv" \
      --output "$out/uncertainty/" 2>&1 | tail -3
  fi
  echo "[$base] 完成 $(date '+%T')"
done
echo "[all] $STUDY 完成"

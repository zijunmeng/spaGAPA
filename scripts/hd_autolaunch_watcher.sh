#!/usr/bin/env bash
# HD 首样本自动验证：等 SLV11 双端 FASTQ 完整（无 .aria2 残留）→ 自动跑端到端管线
set -u
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOG=logs/hd_autolaunch.log
D=data/raw/visium_hd/gse307215_gvhd_gut_hd
echo "[watch] 启动 $(date '+%F %T')" >> $LOG

for i in $(seq 1 720); do  # 最长 30 天
  for srr in SRR35253946; do
    f1=$D/SRX30356448_${srr}_1.fastq.gz
    f2=$D/SRX30356448_${srr}_2.fastq.gz
    if [ -s "$f1" ] && [ -s "$f2" ] && [ ! -e "$f1.aria2" ] && [ ! -e "$f2.aria2" ]; then
      echo "[watch] $srr 完整 $(date '+%T')，启动管线" >> $LOG
      bash scripts/run_visium_hd_pipeline.sh \
        data/raw/visium_hd/gse307215_gvhd_gut_hd \
        gse307215_gvhd "SRX30356448_${srr}*" human >> $LOG 2>&1
      echo "[watch] $srr 管线退出码 $? $(date '+%T')" >> $LOG
      exit 0
    fi
  done
  sleep 600
done

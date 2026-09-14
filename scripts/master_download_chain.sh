#!/usr/bin/env bash
# Visium HD + GSE233208 总下载链（aria2c 分段 + Huntington prefetch）
set -u
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOG=logs/master_download_chain.log
echo "[chain] 启动 $(date '+%F %T')" >> $LOG

# 阶段 1: Visium HD（342GB, 28 文件）
echo "[chain] HD 阶段启动 $(date '+%T')" >> $LOG
aria2c -i /tmp/aria2_visium_hd.txt \
  -x 16 -s 16 -j 3 -k 20M \
  --file-allocation=none --continue=true \
  --max-tries=8 --retry-wait=30 --timeout=120 \
  --summary-interval=600 --console-log-level=warn --download-result=hide >> $LOG 2>&1
echo "[chain] HD 阶段退出码 $? $(date '+%T')" >> $LOG
du -sh data/raw/visium_hd/* >> $LOG 2>&1
touch data/raw/visium_hd/HD_PHASE_DONE

# 阶段 2: GSE233208 人 AD/DS + 5XFAD（734GB, 156 文件）
echo "[chain] GSE233208 阶段启动 $(date '+%T')" >> $LOG
aria2c -i /tmp/aria2_gse233208.txt \
  -x 16 -s 16 -j 3 -k 20M \
  --file-allocation=none --continue=true \
  --max-tries=8 --retry-wait=30 --timeout=120 \
  --summary-interval=600 --console-log-level=warn --download-result=hide >> $LOG 2>&1
echo "[chain] GSE233208 阶段退出码 $? $(date '+%T')" >> $LOG
for d in human mouse_5xfad mouse_wt; do
  echo "  $d: $(ls data/raw/gse233208_adds/$d 2>/dev/null | wc -l) 文件 $(du -sh data/raw/gse233208_adds/$d 2>/dev/null | cut -f1)" >> $LOG
done
touch data/raw/gse233208_adds/DOWNLOAD_COMPLETE

# 阶段 3: Huntington 尾状核 HD（NCBI SRA，无 ENA 镜像 → prefetch）
echo "[chain] Huntington prefetch 启动 $(date '+%T')" >> $LOG
mkdir -p data/raw/visium_hd/gse337131_huntington_caudate_hd
export PATH=/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin:$PATH
for srr in SRR39394365 SRR39394366; do
  prefetch --type sra "$srr" \
    --output-directory data/raw/visium_hd/gse337131_huntington_caudate_hd >> $LOG 2>&1 \
    && echo "[chain] $srr ✓" >> $LOG
done
du -sh data/raw/visium_hd/gse337131_huntington_caudate_hd >> $LOG 2>&1

{
echo "[chain] ✅ 全部完成 $(date '+%F %T')"
echo "===== 汇总 ====="
du -sh data/raw/visium_hd/* data/raw/gse233208_adds/* 2>/dev/null
df -h /s1 | tail -1
} >> $LOG
touch data/raw/MASTER_DOWNLOAD_COMPLETE

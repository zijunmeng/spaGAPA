#!/usr/bin/env bash
# GSE233208 aria2c 下载 + 完成汇总（人 AD 305GB → 5XFAD 343GB → WT 87GB）
set -u
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOG=logs/gse233208_download.log
echo "[dl] 启动 $(date '+%F %T') — 156 文件 734GB" >> $LOG

aria2c -i /tmp/aria2_gse233208.txt \
  -x 16 -s 16 -j 3 -k 20M \
  --file-allocation=none --continue=true \
  --max-tries=8 --retry-wait=30 --timeout=120 \
  --summary-interval=300 --console-log-level=warn \
  --download-result=hide >> $LOG 2>&1
RC=$?

# 汇总
{
echo "[dl] aria2c 退出码 $RC $(date '+%F %T')"
for d in human mouse_5xfad mouse_wt; do
  p="data/raw/gse233208_adds/$d"
  [ -d "$p" ] && echo "  $d: $(ls $p | wc -l) 文件, $(du -sh $p 2>/dev/null | cut -f1)"
done
if [ $RC -eq 0 ]; then
  echo "[dl] ✅ 全部完成 $(date '+%F %T')"
  touch data/raw/gse233208_adds/DOWNLOAD_COMPLETE
else
  echo "[dl] ⚠️ 有失败项——检查 $LOG 中 error 行后可重跑本脚本续传"
fi
} >> $LOG

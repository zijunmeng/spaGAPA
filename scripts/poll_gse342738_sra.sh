#!/usr/bin/env bash
set -u
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOG=logs/gse342738_sra_poll.log
PREFETCH=/s1/mengzijun/pkgs/sratoolkit.3.1.1-centos_linux64/bin/prefetch
OUT=data/raw/gse342738_amygdala_hd
mkdir -p $OUT
echo "[poll] 启动 $(date '+%F %T')" >> $LOG
for i in $(seq 1 56); do
  R=$(curl -s --max-time 30 "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=PRJNA1508859" 2>/dev/null | head -2)
  if ! echo "$R" | grep -q "not found\|ERROR"; then
    echo "[poll] ✅ SRA 已释放 $(date '+%F %T')" >> $LOG
    $PREFETCH --type sra --output-dir $OUT PRJNA1508859 >> $LOG 2>&1 \
      && echo "[poll] ✅ prefetch 完成 $(du -sh $OUT | cut -f1)" >> $LOG \
      && touch $OUT/DOWNLOAD_COMPLETE
    exit 0
  fi
  sleep 21600
done
echo "[poll] 14 天超时" >> $LOG

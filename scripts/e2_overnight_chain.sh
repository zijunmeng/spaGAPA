#!/usr/bin/env bash
# E2 过夜自动链：peaks_meta 出现 → 峰恢复 → {E1,E2} 终版进程分析 → 19 样本总表
set -u
cd /s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA
LOG=logs/e2_overnight_chain.log
PY=~/anaconda3/envs/spagapa/bin/python
export OPENBLAS_NUM_THREADS=8

echo "[chain] 启动 $(date '+%F %T')" >> $LOG

# 阶段 1: 等 E2 peaks_meta（最长 16h）
for i in $(seq 1 96); do
  [ -s pipeline_output/gse263789_expand/3m_e2/scapatrap_raw/peaks_meta.csv.gz ] && break
  sleep 600
done
if [ ! -s pipeline_output/gse263789_expand/3m_e2/scapatrap_raw/peaks_meta.csv.gz ]; then
  echo "[chain] 16h 超时未等到 peaks_meta $(date '+%T')" >> $LOG; exit 1
fi
echo "[chain] E2 peaks_meta 出现 $(date '+%T')" >> $LOG

# 等 expand3 跑完 E2 的 bin200+conformal（SLICE COMPLETE 标志）
for i in $(seq 1 60); do
  grep -q "3m_e2 SLICE COMPLETE" logs/20260913_expand3_pipeline.log 2>/dev/null && break
  sleep 120
done
echo "[chain] E2 切片完成 $(date '+%T')" >> $LOG

# 阶段 2: E2 峰恢复（向量化版）
sed -e 's|ad18_e4|3m_e2|g' -e 's|E4|E2|g' scripts/e4_peak_recovery.py > scripts/e2_peak_recovery.py
$PY scripts/e2_peak_recovery.py >> $LOG 2>&1
if [ ! -f pipeline_output/gse263789_expand/3m_e2/binned_200_raw_recovered/qc_summary.json ]; then
  echo "[chain] E2 峰恢复失败 $(date '+%T')" >> $LOG; exit 2
fi
echo "[chain] E2 峰恢复完成 $(date '+%T')" >> $LOG

# 阶段 3: {E1,E2} vs {E3,E4} 终版进程分析
sed -i 's|"gse263789_expand/3m_e2/binned_200"|"gse263789_expand/3m_e2/binned_200_raw_recovered"|' scripts/gse263789_phase3_ad2v1.py
$PY scripts/gse263789_phase3_ad2v1.py >> $LOG 2>&1
echo "[chain] Phase 3 终版完成 $(date '+%T')" >> $LOG
echo "[chain] 全链完成 $(date '+%F %T')" >> $LOG

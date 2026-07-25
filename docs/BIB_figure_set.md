# spaGAPA — BIB 稿件 Figure Set

每张图的源路径 + 说明。主图（Main）+ 补图（Supplementary）。所有路径相对 `spaGAPA/`。

## Main Figures

### Fig 1 — Stereo-seq 3'-end enrichment (APA feasibility)
polyA-capture signature：read 3' 端富集于基因 TES，证明 Stereo-seq 保留 3' end → APA calling 可行。
- 鼠脑：`pipeline_output/gse263789_stereo_pilot/figures/fig1_3prime_enrichment_near_TES.png`（53.4% within 500bp TES）
- 人脑：`pipeline_output/gse269906_ad_brain/figures/fig1_3prime_enrichment_near_TES.png`（47.0%）
- 脚本：`pipeline_output/{gse263789_stereo_pilot,gse269906_ad_brain}/figures/fig1_3prime_enrichment.py`

### Fig 2 — spaGAPA pipeline overview（概念图，待绘）
空间验证 → 稀疏 GP 插补（+conformal 不确定性）→ APA 定量 → Leiden domain → 差异 APA → SVAPA → 批次校正。

### Fig 3 — Head-to-head imputation（GP 击败命名竞品）
5 方法（spaGAPA-GP / stAPAminer / spvAPA / spatial-KNN / mean）mask-20% 插补，RMSE/Pearson/Spearman + 空间保真度。GP 全指标胜竞品 + 快 12–17×。
- 肾：`pipeline_output/benchmark_stapaminer_headtohead/gse183456/headtohead.png`
- 脑：`pipeline_output/benchmark_stapaminer_headtohead/gse220442/headtohead.png`
- 脚本：`scripts/benchmark_stapaminer_headtohead.py`

### Fig 4 — Scalability（可扩展性硬证据）
spaGAPA-fast 在 42k (162s) + 100k (511s) 完成；stAPAminer TIMEOUT、spvAPA FAILED。
- 时间：`pipeline_output/benchmark_runtime/scaling.png`
- 内存：`pipeline_output/benchmark_runtime/mem_scaling.png`
- 表：`pipeline_output/benchmark_runtime/runtime_table.csv`
- 脚本：`scripts/benchmark_runtime_scalability.py`

### Fig 5 — Uncertainty calibration（独家 UQ）
conformal 校准：90% coverage 从 1.0（过度保守）→ 0.90 精确命中；drop 高不确定 spot → RMSE −13.6%。
- 数据：`pipeline_output/uncertainty_calibration/uncertainty_calibration.json`
- 脚本：`scripts/calibrate_uncertainty.py`、模块 `spagapa/imputation/calibration.py`

### Fig 6 — Domain recovery（无监督空间域）
MOB 无监督 Leiden 恢复 5 层：ARI=0.60 / NMI=0.68（vs spvAPA 监督，无需 label）。
- domain 图：`pipeline_output/mob_domain_recovery/mob_domain_map.png`
- uncertainty 图：`pipeline_output/mob_domain_recovery/mob_uncertainty_map.png`
- 指标：`pipeline_output/mob_domain_recovery/spagapa_metrics.json`
- 脚本：`scripts/mob_domain_recovery.py`

### Fig 7 — AD brain differential APA + miRNA remodeling
3 control vs 3 AD 差异 APA（统一-peak，去混杂）→ 双向 APA 重塑 → miRNA 靶向网络改变。
- 火山图：`pipeline_output/gse220442_differential_apa_unified/volcano.png`
- miRNA 故事：`pipeline_output/mirna_story/report.md`（突触/tau 缩短去抑制；核糖体/OXPHOS 延长抑制）
- 表：`pipeline_output/gse220442_differential_apa_unified/results.csv`、`pipeline_output/mirna_story/gene_mirna_targets.csv`

## Supplementary Figures

### SVAPA + competitor overlap
- 报告：`pipeline_output/svapa/mob_overlap_report.md`（spaGAPA 恢复竞品"梯度型"基因 Myef2/Htra3/Lrrc6，与 stAPAminer 分层型互补）
- 注释 SVAPA：`pipeline_output/svapa/{gse220442,gse263789}_svapa_genelevel.csv`、`pipeline_output/svapa/mob_svapa.csv`

### Bias correction（novel 方法）
- 合成恢复 0.999、真实跨样本 PCC +0.28：`pipeline_output/apa_bias_correction/summary.json`
- 模块：`spagapa/analysis/bias_correction.py`

### Per-sample spaGAPA（GSE220442 6 样本 + GSE263789）
- AD vs control domain 对比：`pipeline_output/gse220442_spagapa_runs/control_vs_ad_domains.png`
- 每样本 domain/uncertainty 图：`pipeline_output/gse220442_spagapa_runs/gse220442_gsm680175{1..6}_scapatrap/{domain_map,uncertainty_map}.png`
- GSE263789 binned：`pipeline_output/gse263789_spagapa_runs/binned_200/{domain_map,uncertainty_map}.png`

### GSE263789 full Stereo-seq（鼠脑高分辨）
- PAS peak summary：`pipeline_output/gse263789_stereo_pilot/figures/fig2_pas_peak_summary.png`
- spatial PAS：`pipeline_output/gse263789_stereo_pilot/figures/fig3_spatial_pas_peaks.png`
- domain spatial accuracy：`pipeline_output/gse263789_stereo_pilot/figures/fig_domain_spatial_accuracy.png`
- scaling curve：`pipeline_output/gse263789_stereo_pilot/figures/fig_scaling_curve.png`

### Known-gene APA validation
- 报告：`pipeline_output/known_gene_validation/report.md`（7/8 AD hits 有文献锚点，ARPP19 直接确认）

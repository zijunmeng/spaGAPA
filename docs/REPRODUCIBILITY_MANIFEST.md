# spaGAPA — Reproducibility Manifest

复现所有数据集 + 分析的环境、路径、命令。日期 2026-07-25。

## 1. 环境

### 服务器（CLAUDE.md hostname-aware）
| Host | R_LIBS | Python | R | TMPDIR |
|------|--------|--------|---|--------|
| **S91**（主） | /s1/SHARE/01_software/R_442_SeuratV5/library | ~/anaconda3/envs/spagapa/bin/python | ~/anaconda3/envs/r442/bin/Rscript | /s3/mengzijun/tmp |
| S90 | 同上 | /s1/mengzijun/anaconda3/envs/r442/bin/... | 同 | /s2/mengzijun/tmp |
| S97/98 | 同上 | /s1/mengzijun/anaconda3/envs/r442/bin/... | 同 | /s9(8)2/mengzijun/tmp |

### 关键 env（**必须**，否则 OpenBLAS segfault / NFS 慢）
```bash
export OPENBLAS_NUM_THREADS=8    # 重度 linalg（factorizer）用 8，不是 64！64/32 会 segfault
export OMP_NUM_THREADS=8
export TMPDIR=/s3/mengzijun/tmp  # S91
# R 竞品 arm：
source ~/anaconda3/etc/profile.d/conda.sh && conda activate r442
export RSCRIPT=/home/mengzijun/anaconda3/envs/r442/bin/Rscript
```

### 软件版本
- Python 3.10（spagapa env），spaGAPA 包（本 repo，16 commits）
- R 4.4.2（r442 env），共享 lib 含：scAPAtrap 0.2.0、stAPAminer 0.1.0、spvAPA 0.1.0、Seurat 5.1.0、movAPA
- SAW 8.2.2（Stereo-seq）、STAR 2.7.10b、Space Ranger 4.1.0、samtools/umi_tools/featureCounts
- 参考基因组：GRCh38 (`/s1/SHARE/00_ref_genecode/04_genecode_GRCh38_STAR_db`, `refdata-gex-GRCh38-2024-A`)、GRCm38 (SAW Mus_musculus_index)

## 2. 数据集

| GSE | 平台/组织 | 位置 | 状态 |
|-----|----------|------|------|
| GSE237183 ×18 | Visium 人胶质瘤 | `data/processed/gse237183_gsm*_scapatrap/` | ✅ |
| GSE183456 | Visium 人肾 | `data/processed/gse183456_gsm6047774_scapatrap/`（APA+表达+coords）| ✅ |
| GSE179572 | Visium 人脑 | `data/processed/gse179572_gsm5420751_scapatrap/` | ✅ |
| GSE263789 | Stereo-seq 鼠脑 | `pipeline_output/gse263789_stereo_pilot/`（SAW BAM + scAPAtrap 21,455 PAS + bin200）| ✅ |
| GSE220442 ×6 | Visium 人 AD 脑 | `data/processed/gse220442_gsm680175{1..6}_scapatrap/`（3 ctrl + 3 AD）| ✅ |
| GSE269906 | Stereo-seq 人 AD 脑 | `pipeline_output/gse269906_ad_brain/`（STAR BAM 61G + 3'-bias）| ⚠️ mask 不可获，仅可行性 |
| GSE311383 | Visium HD 人肝 | `data/raw/gse311383/`（FASTQ）| ⚠️ probe 化学阻塞 |

## 3. 核心命令

### PAS calling（scAPAtrap，标准 10x BAM）
```bash
~/anaconda3/envs/spagapa/bin/python scripts/run_scapatrap_spaceranger.py \
  --bam <BAM> --spatial-dir <spatial> --output-root <out> \
  --processed-dir <processed> --dataset-name <name>
# 批量：scripts/run_scapatrap_gse220442_all.sh（6 样本 concurrency=3）
```

### Stereo-seq（SAW → retag → scAPAtrap）
见 `pipeline_output/gse263789_stereo_pilot/{merge_fastqs_per_sample.sh,run_saw_full.sh,run_scapatrap_full.R}`；retag 流程见 `logs/20260715_数据处理记录.md` §6。

### 3'-bias 验证
```bash
~/anaconda3/envs/spagapa/bin/python pipeline_output/gse269906_ad_brain/figures/fig1_3prime_enrichment.py
# 输入：STAR BAM（subsample 0.2%）+ human GTF；输出 47% within 500bp TES
```

### Head-to-head 插补 benchmark（5 方法）
```bash
OPENBLAS_NUM_THREADS=8 RSCRIPT=...~/r442/bin/Rscript \
~/anaconda3/envs/spagapa/bin/python scripts/benchmark_stapaminer_headtohead.py --dataset gse183456
# 含 stAPAminer/spvAPA R arm（scripts/run_{stapaminer,spvapa}_impute.R）
```

### Scalability benchmark
```bash
OPENBLAS_NUM_THREADS=8 ~/anaconda3/envs/spagapa/bin/python scripts/benchmark_runtime_scalability.py \
  --scales 1000 5000 15000 42000 100000 --gene-cap 500 --wall-cap 1200
```

### Uncertainty calibration
```bash
~/anaconda3/envs/spagapa/bin/python scripts/calibrate_uncertainty.py
# 模块 spagapa/imputation/calibration.py（conformal）
```

### 差异 APA（统一 peak，3v3）
```bash
~/anaconda3/envs/spagapa/bin/python scripts/analyze_gse220442_diff_apa_unified.py
```

### MOB domain recovery
```bash
~/anaconda3/envs/spagapa/bin/python scripts/mob_domain_recovery.py
# 用 spvAPA 自带 ST11RUD+ST11label（scripts/export_spvapa_st11.R 导出）
```

### SVAPA + competitor overlap
```bash
~/anaconda3/envs/spagapa/bin/python scripts/phase3_svapa_annotate_mob_overlap.py
```

### APA bias correction
```bash
~/anaconda3/envs/spagapa/bin/python scripts/run_bias_correction.py
```

### miRNA story
```bash
# 见 pipeline_output/mirna_story/（gene→NM via NCBI EUtils → miRDB v6 查询）
```

## 4. 测试
```bash
OPENBLAS_NUM_THREADS=8 ~/anaconda3/envs/spagapa/bin/python -m pytest tests/ -q
# 420 passed / 1 pre-existing failure (highres_bioml Leiden label count, 无关) / 2 skipped
```

## 5. Git
- Repo: `/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA`（.git 在此；项目根的 .git 是坏空壳）
- 最新 commits：`ec18a32`（Phase 3）、`7ae7cfe`（Phase 2 竞品补强）
- 处理日志：`logs/2026072{3,4,5}_数据处理记录.md`
- Spec：`docs/superpowers/specs/2026-07-2{1,4}-*.md`

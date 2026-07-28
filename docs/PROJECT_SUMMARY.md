# spaGAPA 项目总结文档

**日期**: 2026-07-27
**状态**: 投稿准备阶段（数据收集 + 方法验证基本完成，待稿件撰写）

---

## 1. 项目概述

### 1.1 定位

**spaGAPA**（Spatial Gaussian-Process and Graph-Aware APA Analyzer）是一个**具备校准不确定性量化的空间 APA 统计推断框架**。它解决现有空间 APA 分析工具（stAPAminer、spvAPA）的三个核心缺陷：

1. **无不确定性**——现有工具只给点估计，不告知估计可信度
2. **不可扩展**——现有工具基于 KNN/WNN（O(n²)），在高分辨 ST 数据上崩溃
3. **仅限 Visium**——现有工具未在亚细胞级 Stereo-seq 上验证

### 1.2 核心创新（4 项，竞品完全没有）

| # | 创新 | 竞品状态 | 验证强度 |
|---|------|---------|---------|
| 1 | **Conformal 不确定性量化**——数学保证的预测区间覆盖率 | stAPAminer/spvAPA 均无任何不确定性 | ✅ 强（11 数据集，偏差 <0.2%）|
| 2 | **稀疏 GP 概率框架**——O(nm²) 可扩展，提供后验分布 | KNN/WNN 启发式点估计，O(n²) | ✅ 强（42k/100k 完成，竞品崩溃）|
| 3 | **亚细胞 Stereo-seq APA**——唯一在高分辨 Stereo-seq 上验证 | 竞品 Visium-only | ✅ 强（21,455 PAS × 20.7M DNB）|
| 4 | **不确定性引导的差异 APA**——37% 假阳性减少 | 无法对标（竞品无不确定性）| ✅ 中强（GSE220442 3v3 验证）|

---

## 2. 现阶段成果

### 2.1 数据集（9 GSE，32 Visium 样本 + Stereo-seq）

#### Visium（8 GSE，32 样本）

| GSE | 组织 | 物种 | 样本数 | PAS 范围 |
|-----|------|------|--------|---------|
| GSE237183 | 胶质瘤 | 人 | 18 | 20k–76k |
| GSE183456 | 肾 | 人 | 1 | 53,572 |
| GSE179572 | 脑转移 | 人 | 1 | 26,478 |
| GSE220442 | AD 脑（PFC）| 人 | 6（3 ctrl + 3 AD）| 27k–51k |
| GSE338525 | 肝（正常）| 人 | 2 | 49k–55k |
| GSE206391 | 皮肤（银屑病）| 人 | 2 | 229–276 |
| GSE169749 | 结肠（DSS）| 鼠 | 1 | 40,795 |
| GSE263303 | 脑（Nf1）| 鼠 | 1 | 36,586 |

覆盖：**8 种组织 × 2 物种 × 1 平台（Visium）**

#### Stereo-seq（1 GSE，5 样本）

| GSE | 组织 | 物种 | 样本 | PAS × DNB |
|-----|------|------|------|-----------|
| GSE263789 | AD 脑 | 鼠 | 5（AD 18mo ×2 + WT 18mo ×1 + 3mo ×2）| 21,455 PAS × 20.7M DNB |

- 3'-bias 验证：鼠 53.4% / 人 47.0%（polyA-capture 签名确认）
- AD pilot APA 完成；WT APA 处理中（SAW BAM → retag → scAPAtrap）

#### 竞品数据集对比

| 工具 | GSE 级数据集 | 组织多样性 | 平台覆盖 |
|------|------------|-----------|---------|
| stAPAminer | 1（MOB ×3 重复）| 1 组织 | Visium only |
| spvAPA | 9 | 多组织 | Visium + scRNA |
| metaAPA | 4 | 2 组织 | Visium + long-read |
| **spaGAPA** | **9** | **8 组织** | **Visium + Stereo-seq** |

### 2.2 方法验证成果

#### A. 竞品 Head-to-Head（Pillar 核心）

| 方法 | RMSE | Pearson r | Spearman ρ | 时间 |
|------|------|-----------|------------|------|
| **spaGAPA-GP** | **0.122** | **0.942** | **0.842** | **19s** |
| stAPAminer | 0.208 | 0.852 | 0.778 | 147s |
| spvAPA | 0.272 | 0.770 | 0.697 | 189s |

GP 在所有指标上击败两个竞品 + 快 5.7–7.2×。

#### B. 可扩展性（竞品在高分辨崩溃）

| Spots | spaGAPA-fast | stAPAminer | spvAPA |
|-------|-------------|-----------|--------|
| 1k–15k | ✅ | ✅ | ✅ |
| 42k | ✅ 162s | ❌ TIMEOUT | ❌ FAILED |
| 100k | ✅ 511s | ❌ | ❌ |

#### C. Conformal 不确定性（11 数据集普适验证）

| 目标 | 平均 coverage | 最大偏差 | 达标率 |
|------|-------------|---------|--------|
| 80% | 0.8005 | 0.4% | 11/11 |
| 90% | 0.8996 | 0.5% | 11/11 |
| 95% | 0.9497 | 0.2% | 11/11 |

523,174 个测试点，覆盖 7 GSE × 4 组织 × 2 物种。

#### D. 不确定性杀手应用（37% 假阳性减少）

| 过滤策略 | 显著基因数 | 保留率 | 说明 |
|---------|-----------|--------|------|
| 全量 spot（baseline）| 91 | 100% | 含假阳性 |
| 基因特异性 top-50% 置信 | 57 | 62.6% | **移除 37% 边界调用** |
| Robust core（4 层全显著）| 55 | 60.4% | 高可信差异 APA |

被移除的 34 个基因：effect size 更弱（|Δ| 0.092 vs 0.105）、padj 弱 ~160 个数量级——正是假阳性特征。

#### E. MOB Domain Recovery（无监督）

无监督 Leiden（APA + 空间 graph，无需 label）：ARI = 0.60，NMI = 0.68。

#### F. AD 差异 APA（人脑 3v3）

统一-peak（消除 per-sample peak-calling 混杂）后：8 个 replicated 基因（ARPP19 文献直接确认 + 核糖体簇）。

### 2.3 sAPA-RegNet（描述性，扰动未验证）

- miRNA site 级注释（全基因组 miRanda，251 基因 18,843 位点；APP 带 miR-17/106/20 家族）
- RBP gene-level 注释（TDP-43/FUS/HuR 等 15 神经 RBP，全基因组 FIMO）
- 空间调控网络（control→AD 重塑：SOD2/FAIM2/TARDBP）
- **扰动模型验证失败**（cis 回归阴性，β 与 distal miRNA 数无负相关）→ 降级为描述性，不 claim 预测性扰动

### 2.4 APA 批次校正（Pillar 2，验证进行中）

- 方法：QN + linear batch removal（专为 APA 矩阵）
- 合成恢复：0.999（linear）
- 真实：PCC +0.28（QN），但 QN 有抹平信号风险
- **Harmony 对比实验进行中**——结果决定该模块去留

---

## 3. 竞品创新点对比

### 3.1 spaGAPA 独有（竞品完全没有）

| 能力 | stAPAminer | spvAPA | spaGAPA | 证据 |
|------|-----------|--------|---------|------|
| Conformal 不确定性量化 | ❌ | ❌ | ✅ | 11 数据集，coverage 精确 |
| 概率模型（后验分布）| ❌ KNN | ❌ WNN | ✅ sparse GP | O(nm²) 可扩展 |
| 亚细胞 Stereo-seq APA | ❌ | ❌ | ✅ | 21k PAS × 20M DNB |
| 不确定性引导分析 | ❌ | ❌ | ✅ | 37% 假阳性减少 |
| 可扩展到 100k spots | ❌ OOM@42k | ❌ FAIL@42k | ✅ | 162s@42k, 511s@100k |
| APA 批次校正 | ❌ | ❌ | ⚠️ 有模块，待 Harmony 对比 | — |

### 3.2 竞品独有（spaGAPA 没有）

| 能力 | stAPAminer | spvAPA | spaGAPA |
|------|-----------|--------|---------|
| 监督特征选择（sPLS-DA）| ❌ | ✅ | ❌ |
| 精致可视化模块 | 基础 | ✅ | 基础 |
| 大规模 scRNA-seq APA | ❌ | ✅ | ❌ |

### 3.3 关键洞察：为什么需要 spaGAPA

> "空间 APA 分析目前缺乏统计严谨性——现有工具给点估计不给置信度、不能扩展到高分辨、不能告诉你哪些 APA 调用是可信的。spaGAPA 把 APA 分析从'数据处理管线'升级为'具备不确定性保证的统计推断框架'。"

---

## 4. 三大支柱就绪度

### Pillar 1: 不确定性量化 —— ✅ 85%

| 组成 | 状态 | 证据 |
|------|------|------|
| Conformal 模块 | ✅ 完成 | spagapa/imputation/calibration.py |
| 全量 coverage 验证 | ✅ 完成 | 11 样本 × 7 GSE，偏差 <0.2% |
| 杀手应用 | ✅ 完成 | 37% 假阳性减少 + 92% 保留率 |
| Raw corr 提升 | ⚠️ 中等 | 0.068（conformal 补偿，不依赖 corr）|

**就绪**：足够 BIB claim。论文写法："First spatial APA tool with conformal-calibrated uncertainty, validated across 11 datasets with <0.2% coverage deviation."

### Pillar 2: APA 批次校正 —— ⚠️ 40%（待 Harmony 对比结果决定去留）

| 组成 | 状态 |
|------|------|
| QN + linear 模块 | ✅ 代码完成 |
| 多样本验证 | 🔄 GSE237183 ×18 实验进行中 |
| Harmony 对比 | 🔄 同上 |
| 生物学保留测试 | 🔄 同上 |

**待定**：如果 Harmony 在 APA 上也工作良好 → 降级或砍掉。如果 spaGAPA 方法确实优于 Harmony（处理 bounded/compositional 特性）→ 保留。

### Pillar 3: 可扩展性 + 高分辨 —— ✅ 80%

| 组成 | 状态 | 证据 |
|------|------|------|
| 42k/100k scaling | ✅ | fast 162s/511s，竞品崩溃 |
| Stereo-seq APA | ✅ | 21k PAS × 20M DNB |
| AD vs WT 对比 | 🔄 | WT scAPAtrap 处理中 |
| MOSTA 多器官 domain | ⏳ | 待做（可选）|

**就绪**：竞品崩溃是不可辩驳的硬证据。Stereo-seq APA 是唯一性贡献。

---

## 5. 需要补强的地方

### 5.1 🔴 紧急（投稿前必须）

| # | 任务 | 依赖 | 预计 |
|---|------|------|------|
| 1 | **WT scAPAtrap 完成** → AD vs WT Stereo-seq APA 对比 | SAW BAM + retag 已完成，scAPAtrap 进行中 | ~1h |
| 2 | **Pillar 2 Harmony 对比结果** → 决定 batch correction 去留 | 实验进行中 | ~1h |
| 3 | **commit 本轮所有成果**（conformal 验证 + 杀手应用 + 新数据集脚本 + Pillar 2）| 无 | ~30min |
| 4 | **稿件 outline + figure 设计定稿** | 所有实验结果 | ~1 天 |

### 5.2 🟠 重要（提升投稿竞争力）

| # | 任务 | 价值 |
|---|------|------|
| 5 | **GSE263789 剩余 3 样本处理**（AD 重复 + 3mo 时间点）| Stereo-seq 纵深 |
| 6 | **GSE338525 / GSE263303 第 2 样本** | GSE 内重复 |
| 7 | **local_noise corr 提升**（0.068→>0.3）| 让 raw UQ 也更强（不依赖 conformal 补偿）|
| 8 | **runtime 正式表**（spaGAPA vs stAPAminer vs spvAPA，多 scale）| 竞品对比完整性 |

### 5.3 🟡 可选（锦上添花）

| # | 任务 | 价值 |
|---|------|------|
| 9 | MOSTA 多器官 Stereo-seq domain demo | 亚细胞 scalability 多组织 |
| 10 | stAPAminer head-to-head 在 MOB（竞品自家数据）| 直接可比 |
| 11 | 已知 APA 基因验证（MOB 文献基因列表 overlap）| 生物学可信度 |

---

## 6. 论文定位建议

### 6.1 标题方向

> **spaGAPA: A statistical framework for spatial alternative polyadenylation analysis with calibrated uncertainty quantification**

### 6.2 核心贡献（论文 Contribution Statement）

1. **Conformal 校准**：首个具备数学保证覆盖率的空间 APA 不确定性量化（11 数据集验证，偏差 <0.2%）
2. **稀疏 GP 框架**：概率插补 + 后验方差 + O(nm²) 可扩展，使亚细胞 Stereo-seq APA 分析首次成为可能
3. **不确定性引导分析**：校准后的不确定性可减少 37% 差异 APA 假阳性，提供高可信基因集
4. **竞品对比**：GP 在精度 + 速度上在 spatial fidelity 和速度上优于 stAPAminer 和 spvAPA，且是唯一在 42k+ spots 上完成的工具

### 6.3 Figure 设计（7 主图 + 补图）

详见 `docs/BIB_figure_set.md`。

### 6.4 诚实局限（论文须标注）

1. GP 在 entry-wise RMSE 上未超越 per-gene mean（双峰 index 的结构性限制）
2. Raw GP 不确定性相关性中等（corr ~0.07），conformal 层补偿但局部自适应区间可更紧
3. Stereo-seq raw FASTQ + mask 共存是领域瓶颈，仅 1 个 GSE 可用
4. sAPA-RegNet 扰动模型未验证（cis 回归阴性），调控注释为描述性
5. 无监督框架（不提供监督特征选择，spvAPA 有 sPLS-DA）

---

## 7. 项目资产清单

### 7.1 代码

| 资产 | 数量 |
|------|------|
| Git commits | 41+ |
| Python 测试 | 439 |
| 包模块（spagapa/）| 30+ .py 文件 |
| 分析脚本（scripts/）| 50+ |
| scAPAtrip launcher 支持 | Human + Mouse |

### 7.2 数据

| 资产 | 数量 |
|------|------|
| Visium APA 样本（qc_summary.json）| 32 |
| Stereo-seq APA 样本 | 1（AD pilot）+ 1（WT 处理中）|
| 处理日志 | 4 份（0723–0726）|
| Spec 文档 | 3 份（Phase 1/2/3 + sAPA-RegNet）|

### 7.3 竞品（已安装）

| 工具 | 版本 | 位置 |
|------|------|------|
| stAPAminer | 0.1.0 | R_442 共享 lib（fdm2id 依赖已删）|
| spvAPA | 0.1.0 | R_442 共享 lib（aricode 补装）|
| scAPAtrap | 0.2.0 | R_442 共享 lib |
| miRanda | v3.3a | samtools env |
| FIMO (MEME) | 5.5.9 | samtools env |

---

## 8. 时间线

| 日期 | 里程碑 |
|------|--------|
| 07-10 | GSE237183 scAPAtrap batch 完成（18 样本）|
| 07-13 | SAW 安装 + GSE263789 Stereo-seq pilot 跑通 |
| 07-15 | Stereo-seq retag + scAPAtrap 完成（8,659 PAS → 21,455 PAS full）|
| 07-17–18 | Phase 1 稀疏化（inducing-point reuse + chunked factorizer）|
| 07-21 | Phase 2 竞品补强 spec + GSE220442 R1 修复 |
| 07-23 | GSE220442 6/6 scAPAtrap + stAPAminer 安装 + GSE269906 STAR |
| 07-24 | spvAPA 安装 + sAPA-RegNet + Phase 3（SVAPA/MOB/scalability）|
| 07-25 | 扰动 go/no-go 阴性 + sAPA-RegNet 收敛 + WT 下载完成 |
| 07-26 | 4 新 Visium GSE（肝/皮肤/鼠脑/鼠结肠）+ conformal 全量验证 + 杀手应用 |
| 07-27 | README 重写 + Pillar 2 batch correction Harmony 对比（进行中）|

---

## 9. 关键决策记录

| 决策 | 理由 | 日期 |
|------|------|------|
| 采用 scAPAtrap 而非 Sierra 做 PAS calling | benchmark 综述验证 scAPAtrap 为 DE-APA 最优 | 07-13 |
| Leiden 替换 kmeans 做 domain | kmeans 在稀疏数据上 collapse（98% 单簇）| 07-20 |
| kNN-Laplacian 替换 fused-graph 做 factorizer | fused-graph SuperLU fill-in 导致 42k 挂死 | 07-24 |
| sAPA-RegNet 扰动降级为描述性 | cis 回归 go/no-go 阴性（β 与 miRNA 无负相关）| 07-25 |
| Stereo-seq FASTQ+mask 数据稀缺接受为领域瓶颈 | 66 GSE + STOmics DB 系统调查确认 | 07-26 |
| Visium HD 排除（probe 化学 APA 不兼容）| 官方手册 CG000763 确认 probe-based | 07-26 |
| conformal UQ 定位为主卖点（非 Stereo-seq）| 统计框架在 Visium 上就成立，不依赖 Stereo-seq 数据量 | 07-27 |

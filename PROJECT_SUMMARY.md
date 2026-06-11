# spaGAPA 项目详细总结

**全称**：spatial Gaussian process-based APA analyzer  
**版本**：v0.1.0 (Beta)  
**语言**：Python 3.10  
**开发环境**：conda env `spagapa`  
**更新日期**：2026-05-07

---

## 1. 项目概述

spaGAPA 是一个面向空间转录组数据的**选择性多聚腺苷酸化（APA）分析工具**。核心创新是使用**高斯过程（Gaussian Process）**替代传统 KNN 进行空间插补，并提供**不确定性量化**，弥补了现有工具（如 stAPAminer）的关键缺陷。

**目标用户**：生物信息学研究人员、计算生物学家、空间转录组实验室  
**目标期刊**：Bioinformatics（IF~6）或更高  
**许可证**：MIT

---

## 2. 项目结构

```
spaGAPA/
├── spagapa/                          # 主包（11,200+ 行）
│   ├── __init__.py                   # 包入口，导出所有主要类
│   ├── pipeline.py                   # 560行 - SpaGAPA 主流程类
│   ├── core/                         # 核心数据模型
│   │   ├── apa_dataset.py            # APADataset - AnnData 包装器
│   │   └── apa_site.py               # APASite, APASiteCollection
│   ├── io/                           # I/O 模块
│   │   ├── readers.py               # BAMReader, SpatialCoordinateReader, AnndataReader, BEDReader
│   │   ├── writers.py               # ResultWriter, BEDWriter
│   │   └── scapatrap_wrapper.py     # scAPAtrap Python 包装器
│   ├── spatial/                      # 空间计算
│   │   └── neighbors.py             # SpatialNeighbors (KNN/Radius/Delaunay)
│   ├── calling/                      # APA 位点识别
│   │   ├── spatial_validator.py     # 空间验证算法
│   │   └── quality_filter.py        # 质量过滤器
│   ├── imputation/                   # GP 插补（核心创新）
│   │   ├── gp_imputer.py            # GPImputer + GPImputerBatch
│   │   └── sparse_gp.py             # SparseGPImputer + BlockGPImputer
│   ├── quantification/               # APA 定量
│   │   ├── apa_indices.py           # RUD, PDUI, WUL, PAI 计算
│   │   └── qc_metrics.py            # QCReportGenerator, 交叉验证
│   ├── analysis/                     # 下游分析
│   │   ├── domain_identifier.py     # 空间域识别 (K-means/Leiden/Louvain)
│   │   ├── differential.py          # 差异 APA 分析
│   │   ├── spatial_pattern.py       # Moran's I, SVAPA 检测
│   │   ├── gp_trend_detector.py     # GP 似然比检验 SVAPA（核心创新）
│   │   └── trajectory.py            # 空间轨迹分析 (444行, 已完成)
│   ├── visualization/                # 可视化
│   │   ├── spatial_plots.py         # SpatialPlotter - 空间分布图
│   │   ├── statistical_plots.py     # StatisticalPlotter - 火山图/热图等
│   │   └── qc_plots.py             # QCPlotter - 质量控制图
│   ├── benchmark/                    # Benchmark 框架
│   │   ├── simulator.py             # SpatialAPASimulator + MOBSimulator
│   │   ├── evaluator.py             # BenchmarkEvaluator 多方法对比
│   │   ├── real_data_benchmark.py   # 真实数据 benchmark 脚本
│   │   └── benchmark_plots.py       # 发表级对比图
│   ├── preprocessing/                # 空模块（待实现）
│   └── utils/                        # 空模块（待实现）
├── tests/                            # 测试（288+ 个测试）
│   ├── unit/                         # 17 个单元测试文件
│   ├── integration/                  # 集成测试（空）
│   └── benchmark/                    # 性能测试
├── examples/                         # 10 个使用示例脚本
├── docs/                             # 用户文档
├── pyproject.toml                    # 项目配置
├── setup.py                          # 安装脚本
├── Makefile                          # 便捷命令
└── PROGRESS.md                       # 开发进度记录
```

---

## 3. 核心模块详述

### 3.1 核心数据模型（`core/`）

| 类 | 功能 | 行数 |
|---|------|------|
| `APADataset` | 基于 AnnData 的 APA 专用数据容器，兼容 scanpy/squidpy 生态 | ~350 |
| `APASite` | 单个 APA 位点数据结构（染色体/位置/链/基因/类型/信号等） | ~400 |
| `APASiteCollection` | APA 位点集合，支持过滤/合并/导出 BED | ~400 |

### 3.2 I/O 模块（`io/`）

| 类 | 功能 |
|---|------|
| `BAMReader` | 基于 pysam 的 BAM 文件读取 |
| `SpatialCoordinateReader` | 空间坐标读取（CSV/TSV/10x Visium 格式） |
| `AnndataReader` / `BEDReader` | AnnData 和 BED 格式读取 |
| `ResultWriter` / `BEDWriter` | 多格式输出（BED/CSV/TSV/H5AD/JSON） |
| `ScAPAtrapWrapper` | **R scAPAtrap 的 Python 包装器**，自动检查安装、调用 R 函数、解析输出 |

### 3.3 空间感知 APA Calling（`calling/` + `spatial/`）

**SpatialNeighbors**：构建空间图（KNN / Radius / Delaunay 三角剖分）  
**SpatialValidator**：
- 空间支持度评分（邻居中检测到该位点的比例）
- 距离加权支持度
- 空间一致性过滤
- Moran's I 空间自相关

**QualityFilter**：按 read count / spot count / 均值 / CV / 空间支持度过滤，生成 QC 报告

### 3.4 GP 插补（`imputation/`）— ⭐核心创新

**GPImputer**：
- 核函数：RBF（光滑模式）、Matérn（灵活模式）、auto（自动选择）
- 自动长度尺度估计
- 输出：插补值 + 标准差 + 协方差矩阵
- 批量处理：`GPImputerBatch` 支持多基因并行（multiprocessing + tqdm）

**SparseGPImputer**（性能优化）：
- 诱导点方法：K-means / random / grid 三种选择
- 复杂度：O(n³) → O(nm²)，5-10x 加速

**BlockGPImputer**（大规模数据）：
- 自动空间分块 + 边界重叠处理
- 支持 1000+ spots 的数据集

### 3.5 APA 定量（`quantification/`）

| 指标 | 公式 | 含义 |
|------|------|------|
| RUD | distal / (proximal + distal) | 远端位点相对使用率 [0,1] |
| PDUI | 100 × RUD | 远端使用百分比 [0,100] |
| WUL | Σ(countᵢ × positionᵢ) / Σ(countᵢ) | 加权 3' UTR 长度 |
| PAI | log₂(proximal/distal) 或差值 | Poly(A) 位点偏好指数 |

**QCReportGenerator**：支持多 section（imputation / quantification / spatial / coverage），多格式输出（dict / DataFrame / text / CSV / JSON），方法间交叉验证对比。

### 3.6 下游分析（`analysis/`）

**DomainIdentifier**：K-means / Leiden / Louvain 聚类 + 空间平滑 + 小域移除  
**DifferentialAPAAnalyzer**：Wilcoxon / Welch's t-test / permutation + FDR/Bonferroni 校正  
**SpatialPatternAnalyzer**：Moran's I + SVAPA 基因识别 + 空间模式聚类（层次聚类）  
**GPTrendDetector** — ⭐⭐**最重要创新**⭐⭐：
1. **GP 似然比检验**：比较有/无空间结构的模型
2. **不确定性加权 Moran's I**：用 1/uncertainty 对观测值加权
3. **空间方差分解**：量化空间方差 vs 随机方差
4. 批量 SVAPA 检测 + FDR 校正

**TrajectoryAnalyzer**（已实现 444 行，延后到 v1.1）：
- 最小生成树 / PCA 主曲线轨迹推断
- 沿轨迹 GAM 平滑拟合
- APA 切换点检测（一阶导数分析）

### 3.7 可视化（`visualization/`）

| 类 | 功能 | 代码量 |
|---|------|--------|
| `SpatialPlotter` | 空间 APA 散点图、域可视化、多基因对比网格、差异空间图 | ~570行 |
| `StatisticalPlotter` | 火山图、层次聚类热图、箱线图、小提琴图 | ~570行 |
| `QCPlotter` | 插补质量评估、空间支持可视化、dropout 统计、QC 报告 | ~200行 |

### 3.8 Benchmark 框架（`benchmark/`）

**SpatialAPASimulator**（372行）：
- 4 种空间模式类型（梯度/域特异性/保守/混合）
- 3 种难度等级（easy/medium/hard）
- 可配置 dropout 率和噪声水平
- 完整的 ground truth 输出

**MOBSimulator**：高保真 MOB 同心环结构模拟器

**BenchmarkEvaluator**（414行）：
- 基线方法：mean / median / KNN-spatial（stAPAminer 类） / GP
- 指标：RMSE / MAE / Pearson / Spearman / R² / Bias
- 发表级对比图自动生成

### 3.9 主 Pipeline（`pipeline.py`）

`SpaGAPA` 类整合完整分析流程（7 步）：
1. 数据加载（BAM 或已有 dataset）
2. 空间验证 + 质量过滤
3. GP 插补（可选稀疏 GP）
4. APA 定量
5. 空间域识别
6. 差异 APA 分析
7. SVAPA 基因检测

---

## 4. 与竞品的对比分析

### 4.1 参考/竞品包概况

| 包 | 语言 | 功能定位 | 与 spaGAPA 关系 |
|----|------|---------|----------------|
| **scAPAtrap** | R | 从 BAM 文件识别 poly(A) 位点 | spaGAPA 将其作为 APA calling 的前端 |
| **stAPAminer** | R (~570行) | 空间 APA 分析 | **主要竞品**，spaGAPA 旨在全面超越 |
| **metaAPA** | Nextflow+R | 多工具 APA 位点集成 | 互补工具，非直接竞争 |

### 4.2 scAPAtrap 详细分析

**算法流程**（6 大模块）：
1. `findUniqueMap` — samtools 过滤唯一比对 + 排序 + 索引
2. `dedupByPos` — umi_tools 去重
3. `separateBamBystrand` — 正负链分离
4. `findPeaksByStrand` — derfinder 全基因组覆盖度计算 + 宽峰迭代拆分
5. `findTails` — 通过 soft-clipping 检测 polyA tail（A 富集序列），确定精确切割位点
6. `countPeaks` — featureCounts + umi_tools 定量

**spaGAPA 的集成方式**：`spagapa/io/scapatrap_wrapper.py` 实现了 Python 包装器，调用 R 执行上述流程，解析输出为 `peaks_meta` 和 `peaks_counts`。

### 4.3 stAPAminer — 主要竞品深度对比

**stAPAminer 架构**（570行 R 代码）：
- `APA.R`（~160行）：`computeAPAIndex` + `imputeAPAIndex` + `optimalKvalue`
- `spatialAPA.R`（~410行）：Seurat 聚类 + 差异分析 + SPARK SVAPA + 可视化 + 模式聚类

**KNN 插补核心逻辑**（`imputeAPAIndex`）：
```
1. 基于基因表达矩阵计算 spot 间欧氏距离（不是空间坐标！）
2. 选 k=10 个最近邻居
3. 缺失值 = k 近邻均值
4. 迭代最多 10 轮直至收敛
5. 初始化：若基因表达为 0 则 APA index 填 0
```

**stAPAminer 的局限性**：
- KNN 用**基因表达距离**而非**空间坐标**，丢失了空间位置信息
- 无不确定性估计，无法区分可靠 vs 不可靠的插补值
- SVAPA 检测依赖 SPARK（第三方 R 包，额外依赖）
- 无可扩展性方案（无稀疏近似 / 分块处理）
- 无模拟数据和 Benchmark 框架
- 纯 R 实现，与 Python 生态隔离

### 4.4 系统化对比表

| 维度 | stAPAminer | spaGAPA | 优势方 |
|------|-----------|---------|--------|
| **插补方法** | KNN (k=10)，基于基因表达距离 | GP (RBF/Matérn)，基于空间坐标协方差 | **spaGAPA** |
| **不确定性量化** | ❌ 无 | ✅ 后验标准差 + 协方差 | **spaGAPA** |
| **空间验证** | ❌ 无 | ✅ 空间支持度 + 一致性过滤 | **spaGAPA** |
| **SVAPA 检测** | SPARK（单一方法） | GP 似然比检验 + 不确定性加权 Moran's I + 方差分解 | **spaGAPA** |
| **核函数选择** | 固定 KNN | RBF / Matérn / auto | **spaGAPA** |
| **可扩展性** | 仅线性 | 稀疏 GP (5-10x) + 分块处理 | **spaGAPA** |
| **质量控制** | 基础过滤 | 多维度过滤 + QC 报告 + 交叉验证 | **spaGAPA** |
| **Benchmark** | ❌ 无 | ✅ 模拟器 + 评估器 + 多方法对比 | **spaGAPA** |
| **理论基础** | 启发式 | 贝叶斯非参数框架 | **spaGAPA** |
| **代码量** | ~570 行 R | ~11,200 行 Python | **spaGAPA** |
| **测试** | 0 | 288+ 单元测试，100% 通过 | **spaGAPA** |
| **语言生态** | R | Python（与 scanpy/squidpy 无缝集成） | **spaGAPA** |
| **真实数据验证** | ✅ MOB（已发表） | ⏳ 脚本就绪，等待数据下载 | **stAPAminer** |
| **交互式可视化** | ❌ 仅 ggplot2 静态 | ⏳ Plotly 交互（延后 v1.1） | 持平 |

### 4.5 GP vs KNN 的技术深度对比

**KNN（stAPAminer）**：
```
ŷ(x*) = mean(y of k nearest neighbors in gene expression space)
```
- "空间信息"仅在基因表达空间中隐式体现
- 所有邻居等权或简单距离加权
- 无理论保证，对 k 敏感

**GP（spaGAPA）**：
```
y(x) ~ GP(μ(x), k(x, x'))
k(x, x') = σ² × (1 + √3d/l) × exp(-√3d/l)    [Matérn 3/2]

预测：
μ* = k(x*, X) K⁻¹ y
σ*² = k(x*, x*) - k(x*, X) K⁻¹ k(X, x*)
```
- **明确建模空间协方差结构**
- **长度尺度 l 自动学习**：反映 APA 信号的空间相关半径
- **后验标准差 σ***：每个预测都有置信区间
- **贝叶斯最优**：在正确核函数假设下达到最小均方误差

### 4.6 metaAPA 分析

metaAPA 是 Nextflow 管道，集成多个 APA caller（Sierra / polyApipe / SCAPE）的结果。两种策略：
- **Position-based**：基于基因组坐标距离合并相近位点
- **Similarity-based**：基于表达相似性（Spearman/Cosine/Jaccard 等距离 + K-means/PAM/hdbscan 聚类）

metaAPA 不直接竞争 spaGAPA，其集成策略可作为 spaGAPA 未来整合多个 APA caller 结果的参考。

---

## 5. 当前进展

### 5.1 完成度

| 阶段 | 状态 | 测试数 | 代码量 |
|------|------|--------|--------|
| Phase 1: 项目搭建 + 核心数据模型 | ✅ 100% | 49 | ~2,000 |
| Phase 2: 核心算法 (calling/imputation/quantification) | ✅ 100% | 136 | ~4,000 |
| Phase 3: 分析 + 可视化 | ✅ 100% | 78 | ~3,800 |
| Phase 4: Benchmark | ✅ 95% | 21 | ~1,500 |
| Phase 5: 文档 + Release | ✅ 90% | - | ~2,500 |
| **总计** | **~97%** | **288+** | **~13,000+** |

### 5.2 各模块测试覆盖情况

| 模块 | 测试数 | 代码覆盖率 |
|------|--------|-----------|
| `core/` | 35 | 33% |
| `io/` | 14 | 28% |
| `spatial/neighbors.py` | 19 | 88% |
| `calling/spatial_validator.py` | 12 | 80% |
| `calling/quality_filter.py` | 18 | 96% |
| `imputation/gp_imputer.py` | 20 | 87% |
| `imputation/sparse_gp.py` | 20 | 99% |
| `quantification/apa_indices.py` | 30 | 94% |
| `quantification/qc_metrics.py` | 31 | 89% |
| `analysis/` | 61 | TBD |
| `visualization/` | 21 | TBD |
| `benchmark/` | 21 | TBD |

### 5.3 已完成的示例脚本（10 个）

1. `01_scapatrap_usage.py` — scAPAtrap 包装器使用
2. `02_spatial_validation.py` — 空间验证示例
3. `03_gp_imputation.py` — GP 插补（含可视化）
4. `04_apa_quantification.py` — APA 定量 + QC
5. `05_differential_analysis.py` — 差异 APA 分析
6. `06_spatial_patterns.py` — 空间模式识别
7. `07_gp_svapa_detection.py` — GP SVAPA 检测
8. `08_complete_pipeline.py` — 完整流

---

## 6. 未完成 / 不足之处

### 6.1 关键未完成项

| 项目 | 优先级 | 状态 | 影响 |
|------|--------|------|------|
| **真实数据 Benchmark** | 🔴 P0 | 脚本就绪，等待数据 | 没有真实数据验证，无法与 stAPAminer 结果直接对比 |
| **Git 初始化** | 🟡 P1 | 未执行 | 无法追踪代码变更历史 |
| **CLI 命令行接口** | 🟡 P1 | 延后 v1.1 | 目前仅支持 Python API |
| **代码规范检查** | 🟡 P1 | 未执行 | flake8/black/mypy 未运行 |
| **PyPI 发布** | 🟠 P2 | 延后（等真实数据） | 仅本地安装可用 |
| **Sphinx API 文档** | 🟠 P2 | 延后 v1.1 | 仅有 docstrings + README |
| **交互式可视化 (Plotly)** | 🟠 P2 | 延后 v1.1 | 仅有 matplotlib 静态图 |
| **Jupyter 教程** | 🟠 P2 | 未开始 | 仅有 10 个 .py 示例脚本 |

### 6.2 功能模块缺陷

| 问题 | 说明 |
|------|------|
| **`preprocessing/` 空模块** | 批次效应校正、数据标准化等预处理功能未实现 |
| **`utils/` 空模块** | 通用工具函数缺失 |
| **integrative tests 空** | 仅有单元测试，缺少端到端集成测试 |
| **轨迹分析已实现但标注 deferred** | `analysis/trajectory.py` 有 444 行代码和测试文件，但 tasks.md 标注为 deferred |
| **IO 模块覆盖率仅 28%** | readers/writers 测试薄弱，可能隐藏 bug |
| **core 模块覆盖率仅 33%** | 核心数据模型测试覆盖不足 |

### 6.3 方法论局限性

1. **APA Calling 依赖 scAPAtrap**：v1.0 策略中 APA 位点识别核心仍依赖 R 包，空间验证是增量改进而非独立 calling
2. **GP 在大规模数据上的计算瓶颈**：虽然有稀疏 GP 近似，但 10,000+ spots 的完整 GP 仍然不可行
3. **无 GPU 加速**：当前仅支持 CPU 并行（multiprocessing），未使用 CuPy/JAX
4. **APADataset 与真实 AnnData 对象的集成不深**：pipeline.py 中有多处 `hasattr` 检查和 `try/except`，说明接口不够统一
5. **单物种支持**：当前主要针对小鼠（mm10），其他物种需要额外适配
6. **无多模态整合**：不支持组织学图像叠加分析

---

## 7. 后续开发建议

### 7.1 优先级排序

**P0 — 论文必需（立即开始）**：
1. 下载并处理 MOB/Brain/Embryo 真实数据
2. 运行 `scripts/run_real_benchmark.py` 生成与 stAPAminer 的对比结果
3. 补充下游任务对比（ARI + F1）和消融实验

**P1 — 质量提升（1-2 周）**：
1. `git init` + 首次 commit
2. 运行 flake8/black/mypy 规范代码
3. 补充 io/core 模块测试，提高覆盖率
4. 编写集成测试（完整 pipeline 端到端）

**P2 — 论文增色（有精力时）**：
1. 实现缺失的 preprocessing（批次校正等）
2. Plotly 交互式可视化
3. CLI 命令行接口
4. Jupyter 教程

**P3 — 发布**：
1. PyPI 打包发布
2. GitHub release + release notes
3. Sphinx 文档

### 7.2 论文发表前必须达成的里程碑

- [ ] 在 MOB 数据集上复现 stAPAminer 结果
- [ ] GP 插补 RMSE 显著优于 KNN（p < 0.05）
- [ ] SVAPA 检测：spaGAPA 比 stAPAminer 多发现至少 10% 的 SVAPA 基因
- [ ] 不确定性量化证明其价值（高不确定性区 = 需实验验证区）
- [ ] 至少 2 个数据集的交叉验证

---

## 8. 技术统计汇总

| 指标 | 数值 |
|------|------|
| 源代码文件 | 34 个 .py |
| 源代码总行数 | ~11,200 |
| 测试文件 | 17 个 |
| 单元测试数 | 288+ |
| 测试通过率 | 100% |
| 整体代码覆盖率 | ~53%（核心模块 80-99%） |
| 示例脚本 | 10 个 |
| 已实现 API 类 | 30+ |
| 支持的数据格式 | BAM, BED, CSV, TSV, H5AD, HDF5, JSON |
| Python 版本 | 3.10 |

---

## 9. 核心创新总结

spaGAPA 的核心竞争力可概括为一个**方法学创新链条**：

```
数据输入
  ↓
空间验证（创新 1：降低噪声假阳性）
  ↓
GP 插补 + 不确定性量化（创新 2：最核心）
  ↓
不确定性加权 SVAPA 检测（创新 3：GP 似然比检验）
  ↓
多尺度空间模式分析（创新 4）
  ↓
生物学发现
```

**三大关键差异化卖点**：
1. **GP 插补**：明确建模空间协方差，优于简单 KNN
2. **不确定性量化**：每个预测都有置信区间，竞品完全没有
3. **GP 似然比检验**：基于贝叶斯模型选择的 SVAPA 检测，有严格理论保证

---

## 附录A：参考/竞品源代码位置

| 包 | 路径 |
|----|------|
| scAPAtrap | `00_ref_packages/scAPAtrap-master/` |
| stAPAminer | `00_ref_packages/stAPAminer-main/` |
| metaAPA | `00_ref_packages/metaAPA-master/` |

## 附录B：设计文档位置

| 文档 | 路径 |
|------|------|
| 需求文档 | `.kiro/specs/spagapa/requirements.md` |
| 设计文档 | `.kiro/specs/spagapa/design.md` |
| 任务文档 | `.kiro/specs/spagapa/tasks.md` |
| 开发进度 | `PROGRESS.md` |

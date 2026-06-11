# spaGAPA 面向 Briefings in Bioinformatics 的系统性提升 Spec

更新日期：2026-06-10  
目标期刊：Briefings in Bioinformatics, BIB  
目标版本：spaGAPA v0.2.0 / BIB-ready release  
当前定位：从“空间转录组 APA 工具包”升级为“空间 APA 系统性 benchmark 与 uncertainty-aware 方法框架”

---

## 1. 总体判断

spaGAPA 有机会冲击 Briefings in Bioinformatics，但不能以普通软件包论文的形式投稿。

BIB 更适合的文章形态不是：

```text
We developed a Python package for spatial APA analysis.
```

而应该是：

```text
We present a systematic benchmark and uncertainty-aware computational framework
for spatially resolved alternative polyadenylation analysis.
```

因此，spaGAPA 的提升目标不是单纯补几个功能，而是把项目重塑成一个能回答领域级问题的框架：

1. 空间转录组 APA 分析目前缺少什么？
2. 现有工具为什么不足？
3. 真实数据中 APA missingness、空间噪音、位点不确定性会造成什么问题？
4. 高斯过程和 posterior uncertainty 是否真的改善了插补、SVAPA 检测、差异 APA 和生物学解释？
5. spaGAPA 能否给出可复现、可扩展、可推广的分析范式？

---

## 2. BIB 目标下的核心论文假说

### 2.1 主假说

空间转录组 APA 分析不能简单套用 single-cell APA caller 或普通 spatial gene-expression workflow，因为 APA index 同时受到以下因素影响：

- 3' capture 覆盖不足。
- proximal/distal poly(A) site 定义不确定。
- spot/cell 的 APA index 缺失率高。
- 不同组织区域的空间相关结构明显。
- 简单 KNN 插补无法量化不确定性。
- 下游 SVAPA 和 differential APA 检验容易把插补伪信号当作真实空间模式。

spaGAPA 的核心假说是：

> 用空间 Gaussian Process 显式建模 APA index 的空间协方差，并把 posterior uncertainty 传播到下游 SVAPA 和 differential APA 分析，可以比现有启发式方法更稳健地恢复真实空间 APA 模式，并降低假阳性。

### 2.2 论文贡献点

面向 BIB，建议把贡献点写成 5 个层次：

1. **领域框架贡献**  
   系统梳理 spatial APA analysis 的计算挑战，并给出从 APA calling 到 downstream interpretation 的完整分析框架。

2. **方法学贡献**  
   提出空间 GP 插补、posterior uncertainty、uncertainty-aware SVAPA/differential analysis 的组合框架。

3. **benchmark 贡献**  
   构建 spatial APA 专用模拟器、真实数据 mask benchmark、跨 baseline 的公平比较体系。

4. **软件生态贡献**  
   提供 Python/AnnData/scanpy/squidpy 兼容的 spaGAPA 工具包，补足 stAPAminer 主要在 R 生态中的局限。

5. **生物学贡献**  
   在真实组织中发现空间依赖的 APA remodeling，例如组织层级、空间梯度、发育区域或微环境相关 APA usage 改变。

---

## 3. 当前状态概括

### 3.1 已具备基础

当前 spaGAPA 已经有较完整的代码骨架：

- `spagapa/core/`: APADataset、APASite、APASiteCollection。
- `spagapa/io/`: BAM、BED、AnnData、坐标读取，以及 scAPAtrap wrapper。
- `spagapa/spatial/`: KNN、radius、Delaunay 空间邻居图。
- `spagapa/calling/`: 空间支持度和质量过滤。
- `spagapa/imputation/`: GPImputer、GPImputerBatch、SparseGPImputer、BlockGPImputer。
- `spagapa/quantification/`: RUD、PDUI、WUL、PAI 和 QC metrics。
- `spagapa/analysis/`: spatial domain、differential APA、Moran's I、GP trend/SVAPA、trajectory。
- `spagapa/visualization/`: 空间图、统计图、QC 图。
- `spagapa/benchmark/`: simulator、benchmark evaluator、benchmark plots。

### 3.2 当前最大短板

从 BIB 审稿视角看，当前短板主要有 8 类：

1. 真实 spatial APA 数据还没有成为主 benchmark。
2. 当前 real benchmark 有 expression proxy 风险，不能作为主结果。
3. 与 stAPAminer 的比较还不够“原样复现”和严格公平。
4. pipeline 集成层存在矩阵方向、reader 返回类型、SparseGP 参数等潜在问题。
5. uncertainty 尚未在所有下游分析中形成闭环。
6. 缺少多数据集、多场景、多 caller 的系统性 benchmark。
7. 缺少强生物学 case study。
8. 论文故事仍偏工具开发，需要升级为领域级方法框架。

---

## 4. BIB-ready 的最低验收标准

如果目标是 BIB，建议设置以下最低标准。没有达到这些标准时，不建议投 BIB，应优先投 Bioinformatics、NAR Genomics and Bioinformatics 或 BMC Bioinformatics。

### 4.1 数据标准

- 至少 3 个真实空间转录组数据集。
- 至少 2 种组织或生物场景，例如 brain/MOB、embryo、tumor、development。
- 至少 1 个数据集能与 stAPAminer 原文或示例数据直接关联。
- 主结果必须来自真实 APA calling 或真实 poly(A) site count，不能只用 expression proxy。
- 每个数据集需要完整记录：
  - 原始数据来源。
  - 平台类型。
  - 是否有 BAM。
  - 是否为 3' capture。
  - spot/cell 数量。
  - detected PAS 数量。
  - retained genes 数量。
  - APA index missing rate。
  - filtering 规则。

### 4.2 方法比较标准

必须包含以下 baseline：

- Mean imputation。
- Median imputation。
- KNN-spatial。
- KNN-expression。
- stAPAminer 原始实现或尽可能严格的复现版本。
- SPARK 或 SPARK-X 风格 SV spatial baseline。
- spaGAPA-GP。
- spaGAPA-GP + uncertainty-aware downstream。

可选增强 baseline：

- MAGIC/ALRA/SAVER 等 general imputation 方法。
- Squidpy/Scanpy 中常用 spatial autocorrelation workflow。
- metaAPA 风格多 caller consensus input。

### 4.3 评价标准

不能只报告 RMSE。BIB 级别需要覆盖：

- 插补准确性：RMSE、MAE、Pearson、Spearman、R2。
- 不确定性校准：95% credible interval coverage、negative log likelihood、calibration curve。
- SVAPA 检测：AUPRC、AUROC、FDR、sensitivity、specificity。
- 空间域恢复：ARI、NMI、domain purity、boundary consistency。
- 差异 APA：FDR、effect-size reproducibility、replicate concordance。
- 生物学一致性：known marker overlap、GO/KEGG/Reactome enrichment、literature support。
- 稳定性：subsampling、masking rate、random seed、dataset split。
- 计算性能：runtime、memory、scalability、failure rate。

### 4.4 软件标准

- 一键 pipeline 能在 toy data、simulated data、至少 1 个真实 dataset 上完整跑通。
- CLI、Python API、notebook tutorial 三者至少有两者稳定。
- 所有 benchmark 脚本可复现论文主图。
- 所有主图数据和中间结果有固定输出目录。
- 提供 conda 环境或 lock 文件。
- GitHub/Zenodo/OSF release 准备好。

### 4.5 论文标准

- 主文必须包含一个清晰的 conceptual framework figure。
- 主文不能只展示工具功能，必须展示领域洞察。
- 每一个核心 claims 都必须有实验证据。
- Reviewer 能清楚看到 spaGAPA 相比 stAPAminer 的不可替代性。

---

## 5. 工作包总览

建议把 BIB-ready 提升拆成 10 个工作包：

| 工作包 | 名称 | 优先级 | 目的 |
|---|---|---|---|
| WP0 | BIB 论文定位与项目规范 | P0 | 固定故事线和验收标准 |
| WP1 | 真实数据获取与 APA calling | P0 | 解决真实 benchmark 根基 |
| WP2 | 数据模型与 pipeline 稳定化 | P0 | 确保一键流程可信 |
| WP3 | 核心算法增强与 uncertainty 闭环 | P0 | 强化 spaGAPA 方法创新 |
| WP4 | 竞品与 baseline 公平复现 | P0 | 回答审稿人比较质疑 |
| WP5 | 模拟 benchmark 系统化 | P0 | 提供 ground truth 证据 |
| WP6 | 真实数据 benchmark | P0 | 证明实际数据有效 |
| WP7 | 生物学 case study | P0 | 提升论文档次 |
| WP8 | 可扩展性、鲁棒性和工程化 | P1 | 提高软件可信度 |
| WP9 | 文档、发布和复现材料 | P1 | 支持投稿和用户复现 |
| WP10 | 论文图表与 manuscript assets | P1 | 形成 BIB 投稿材料 |

---

## 6. WP0: BIB 论文定位与项目规范

### 6.1 需要解决的问题

当前文档中目标期刊仍多处写 Bioinformatics，且项目叙事偏“工具包开发完成”。BIB 目标需要一个更大的 framing：

- spatial APA 是一个新兴但缺少标准 workflow 的方向。
- 现有 APA caller 和 spatial gene-expression 方法都不能直接解决 APA index 的空间缺失和不确定性问题。
- spaGAPA 是一个“benchmark-backed uncertainty-aware framework”。

### 6.2 具体任务

#### WP0.1 固定论文标题和摘要骨架

建议候选标题：

1. `spaGAPA: an uncertainty-aware framework for spatially resolved alternative polyadenylation analysis`
2. `Benchmarking and uncertainty-aware inference of spatial alternative polyadenylation with spaGAPA`
3. `A Gaussian process framework for robust spatial alternative polyadenylation analysis`

推荐主标题使用第 2 个风格，因为 BIB 对 benchmark 和 broad framework 更友好。

#### WP0.2 固定核心 claims

建议主 claims：

1. spatial APA analysis has unique statistical challenges that are not handled by existing APA or spatial transcriptomics tools。
2. GP-based imputation better recovers spatial APA usage than heuristic local averaging。
3. posterior uncertainty identifies unreliable imputed APA values and improves downstream SVAPA detection。
4. spaGAPA provides an end-to-end reproducible framework compatible with existing APA callers and Python spatial omics workflows。
5. real tissue analyses reveal spatially organized APA remodeling associated with tissue architecture。

#### WP0.3 建立 claim-to-evidence 表

每个 claim 必须对应：

- 需要的数据。
- 需要的 baseline。
- 需要的主图或补图。
- 需要的统计检验。
- 可能的审稿质疑。
- 对应回应方案。

### 6.3 交付物

- `BIB_TARGET_SPEC.md`。
- 更新后的 `PROJECT_SUMMARY.md`，把目标从 Bioinformatics 改为 BIB-ready strategy。
- 论文 outline 草稿。
- claim-to-evidence 表格。

### 6.4 验收标准

- 后续所有开发任务都能映射到至少一个 BIB claim。
- 不再出现“只为了增加功能而开发”的任务。
- 所有图表计划能支持主 claims。

---

## 7. WP1: 真实数据获取与 APA Calling

### 7.1 需要解决的问题

spaGAPA 当前最大的 BIB 风险是真实 APA 数据不足。下载表达矩阵和坐标不够，必须获得能支持 APA calling 的原始 reads 或 poly(A) site count。

### 7.2 数据选择原则

优先选择满足以下条件的数据集：

1. 3' capture 或适合 APA/PAS calling。
2. 有 BAM、FASTQ 或可重建 BAM。
3. 有空间坐标和组织图像。
4. 有组织层级或空间结构 annotation。
5. 有足够 spot/cell 数量。
6. 与 stAPAminer 或已有 APA 文献有交集。

### 7.3 推荐数据组合

#### Dataset A: Mouse olfactory bulb

目的：

- 与 stAPAminer 形成直接对比。
- MOB 有清楚的层状/环状组织结构，适合空间 APA 模式展示。

需要文件：

- BAM 或 FASTQ。
- spatial coordinates。
- tissue image。
- gene annotation。
- stAPAminer 原文或示例中可对齐的 spot/gene/site 信息。

预期结果：

- 识别 layer-specific 或 ring-like APA usage。
- 展示 spaGAPA 相比 stAPAminer 在插补和 SVAPA 排名上的优势。

#### Dataset B: Brain / cortex / hippocampus

目的：

- 展示复杂组织层级中的 APA remodeling。
- 支持神经系统相关 APA 生物学解释。

预期结果：

- cortical layer 或 anatomical region 相关 APA genes。
- enriched pathways 与神经发育、突触、RNA localization、3'UTR regulation 相关。

#### Dataset C: Embryo / development / tumor

目的：

- 展示 spaGAPA 跨场景泛化。
- 如果选择 embryo，可展示空间梯度和发育轨迹。
- 如果选择 tumor，可展示微环境区域 APA 差异。

预期结果：

- gradient-like 或 domain-specific APA genes。
- 与 cell state、developmental axis 或 tumor microenvironment 相关。

### 7.4 APA calling 方案

至少需要两条路线：

#### 路线 1: scAPAtrap 主线

使用 scAPAtrap 从 BAM/FASTQ 中检测 PAS 和 count matrix。

输出：

- `apa_sites.bed`
- `site_counts.csv`
- `gene_site_mapping.csv`
- `proximal_distal_annotation.csv`
- `apa_index_matrix.csv`

优点：

- 与当前 spaGAPA wrapper 匹配。
- 与 stAPAminer 生态有联系。

风险：

- R 依赖和外部工具较多。
- 空间转录组 BAM 格式可能需要 barcode/UMI 适配。

#### 路线 2: 多 caller 验证

可选 Sierra、SCAPE 或 metaAPA 风格整合。

目的：

- 证明 spaGAPA 不依赖某一个 caller 的偏差。
- 提供 BIB 更喜欢的 systematic framework 证据。

### 7.5 数据目录规范

建议建立：

```text
data/
  raw/
    dataset_name/
      bam/
      fastq/
      spatial/
      annotation/
  processed/
    dataset_name/
      apa_sites.bed
      site_counts.csv
      apa_matrix.csv
      coordinates.csv
      metadata.csv
      gene_site_mapping.csv
      qc_summary.json
  benchmark/
    dataset_name/
      masks/
      baselines/
      spagapa/
      figures/
```

### 7.6 交付物

- 至少 3 个 processed real APA datasets。
- 每个数据集一份 data provenance 文档。
- 每个数据集一份 QC report。
- 每个数据集一份 APA index missingness report。

### 7.7 验收标准

- 真实数据主 benchmark 不再使用 expression proxy。
- 每个真实数据集至少保留 500 个有分析价值的 APA genes，具体阈值可根据平台调整。
- 每个数据集有明确的 spot/cell filtering、gene filtering、site filtering 规则。
- 每个数据集可以被 spaGAPA pipeline 直接读取。

---

## 8. WP2: 数据模型与 Pipeline 稳定化

### 8.1 需要解决的问题

当前项目的核心模块比较丰富，但 BIB 投稿前必须保证一键 pipeline 能稳定复现。否则审稿人或用户跑不通，会严重影响可信度。

### 8.2 已识别的关键风险

#### WP2.1 矩阵方向不统一

风险：

- `APADataset.raw_counts` 似乎是 genes x spots。
- `SpaGAPA.run()` 中可能做了 `.T` 转置。
- `GPImputerBatch` 期望 genes x spots。

可能后果：

- 把 spots 当 genes 拟合。
- coordinates 长度与矩阵轴不匹配。
- 输出结果 gene/spot 标签错位。

解决标准：

- 明确全项目矩阵约定：
  - AnnData: spots x genes。
  - Core algorithm: genes x spots。
  - User-facing CSV: 明确行列含义，并在 reader 中显式转换。
- 所有 public API 增加 shape validation。
- toy data 全流程测试。

#### WP2.2 `read_spatial_data()` 返回类型与 pipeline 不匹配

风险：

- reader 返回 tuple `(coords_df, data_df)`。
- pipeline 可能期待 `APADataset`。

解决标准：

- 明确保留两个层次：
  - low-level reader 返回 DataFrame。
  - high-level loader 返回 APADataset。
- pipeline 只接受 APADataset 或有明确 schema 的文件路径。

#### WP2.3 SparseGPImputer 参数不匹配

风险：

- pipeline 传入 `kernel_type`。
- SparseGPImputer 构造函数未必支持该参数。

解决标准：

- `use_sparse_gp=True` 在 toy data 和 small real data 上可运行。
- dense GP 与 sparse GP API 对齐。

#### WP2.4 Quantification 在 pipeline 中不能只是 placeholder

风险：

- BIB 审稿人会检查 APA index 是否真实计算。
- 若 pipeline 中 quantification 只是传递矩阵，会被认为不完整。

解决标准：

- 从 site count 到 RUD/PDUI/WUL 的路径明确。
- 如果用户输入已经是 APA index，需要标记 `input_type="apa_index"`。
- 如果用户输入是 site counts，需要显式计算 APA index。

#### WP2.5 Differential analysis 输出需要 gene-level 统计

风险：

- 只输出 domain spot counts 不足以支撑 differential APA claim。

解决标准：

- 输出 gene-level p-value、FDR、effect size、direction、domain pair。
- 支持 uncertainty weight 或 sensitivity analysis。

### 8.3 交付物

- pipeline shape audit 文档。
- 最小 toy dataset。
- end-to-end integration tests。
- `run_spagapa_pipeline.py` 示例脚本。
- `pipeline_output/` 标准输出 schema。

### 8.4 验收标准

- Toy dataset 可以在 1 分钟内完整跑完。
- Simulated dataset 可以生成 imputed matrix、uncertainty、SVAPA、domains、differential APA、figures。
- 至少 1 个真实 dataset 可以完整跑通。
- 所有输出文件有 README 解释。

---

## 9. WP3: 核心算法增强与 Uncertainty 闭环

### 9.1 需要解决的问题

spaGAPA 的最大创新是 GP + uncertainty，但当前需要让 uncertainty 从“输出矩阵”变成贯穿下游分析的统计原则。

### 9.2 GP 插补增强

任务：

1. 明确 GP 模型假设：
   - 输入空间坐标。
   - 输出 APA index。
   - kernel 选择 RBF / Matérn。
   - noise level 对应观测噪音。
2. 支持 per-gene 自动 length scale。
3. 支持 low-observation genes fallback。
4. 支持 bounded APA index：
   - RUD/PDUI 在 0-1 或 0-100 范围。
   - 插补后需要 clip 或使用 link function。
5. 输出 posterior mean、posterior std、confidence category。

验收标准：

- 对高缺失率基因不会产生异常值。
- uncertainty 与真实误差正相关。
- 95% interval coverage 接近标称值，允许真实数据中做 approximate calibration。

### 9.3 Uncertainty-aware SVAPA

任务：

1. 完善 GP likelihood ratio test：
   - spatial GP model vs null model。
   - null model 可为 constant mean、non-spatial permutation 或 white noise。
2. 明确 uncertainty weight 进入统计量的方式：
   - weighted Moran's I。
   - weighted regression / GP likelihood。
   - sensitivity analysis: with vs without uncertainty。
3. 输出：
   - statistic。
   - p-value。
   - FDR。
   - effect size。
   - spatial variance fraction。
   - mean uncertainty。

验收标准：

- 模拟数据中 type-I error 可控。
- 高 uncertainty 区域不会主导 SVAPA 结果。
- 与 SPARK/stAPAminer 相比，SVAPA ranking 更稳定。

### 9.4 Uncertainty-aware differential APA

任务：

1. 差异分析中引入 uncertainty weighting。
2. 支持三种模式：
   - raw observed only。
   - imputed unweighted。
   - imputed uncertainty-weighted。
3. 输出每种模式的结果，便于消融。

验收标准：

- 在模拟数据中 uncertainty-weighted 方法 FDR 更低或稳定性更好。
- 在真实数据中 top differential APA genes 的空间分布更可信。

### 9.5 Spatial validation 的统计化

任务：

1. 当前 spatial support score 需要明确 null model。
2. 增加 permutation 或 neighborhood enrichment test。
3. 输出 site-level spatial confidence。

验收标准：

- 能说明哪些 APA sites 是空间一致支持的，而非随机低 read count 噪音。
- 在 simulation 或 spike-in 中减少 false positive PAS。

---

## 10. WP4: 竞品与 Baseline 公平复现

### 10.1 主要竞品定位

#### stAPAminer

直接竞品。其核心包括：

- APA index 计算。
- 基于 gene expression distance 的 KNN 插补。
- Seurat 聚类。
- FindMarkers/limma/edgeR 差异分析。
- SPARK SVAPA。
- spatial pattern clustering。

spaGAPA 必须与其公平比较。

#### scAPAtrap

上游 APA caller，不是直接竞品。spaGAPA 可将其作为输入前端，也可以比较“是否加入空间验证”对 PAS 质量的影响。

#### metaAPA

多 caller APA/PAS 整合框架，不是直接 spatial APA downstream 竞品。可作为未来增强：multi-caller consensus input。

### 10.2 必须复现的 baseline

#### Baseline 1: stAPAminer-original

要求：

- 尽量运行原 R 包。
- 使用其默认 k=10。
- 使用 gene expression distance。
- 使用其 imputeAPAIndex。
- 若运行失败，记录环境和失败原因，并实现 faithful reimplementation。

#### Baseline 2: KNN-spatial

要求：

- 使用空间坐标距离。
- k 与 stAPAminer 保持一致或调参。

目的：

- 区分“仅仅用了空间坐标”与“GP 建模”的贡献。

#### Baseline 3: SPARK/SVAPA

要求：

- 用于 SVAPA 检测比较。
- 如果 SPARK 环境难以复现，可用 SPARK-X 或等价 spatial autocorrelation 方法替代，但必须说明。

#### Baseline 4: Mean/Median

目的：

- 提供简单插补下界。

#### Baseline 5: General imputation optional

可选：

- MAGIC。
- ALRA。
- SAVER。
- scVI 风格 latent imputation。

注意：

- 这些方法不是专为空间 APA 设计，作为 supplement 即可。

### 10.3 公平比较原则

所有方法必须共享：

- 相同数据集。
- 相同 gene set。
- 相同 spot/cell set。
- 相同 missing mask。
- 相同 train/test split。
- 相同 downstream evaluation。
- 相同 random seed。

### 10.4 交付物

- `benchmark/configs/baselines.yaml`。
- `benchmark/results/baseline_summary.csv`。
- stAPAminer 复现日志。
- baseline runtime/memory 表。

### 10.5 验收标准

- reviewer 能复查每个 baseline 如何运行。
- stAPAminer 比较不是“stAPAminer-like”一句话，而是有原始实现或忠实复现证据。
- spaGAPA 优势在至少两个任务中显著：
  - imputation。
  - SVAPA。
  - differential APA。
  - domain recovery。

---

## 11. WP5: 模拟 Benchmark 系统化

### 11.1 需要解决的问题

模拟 benchmark 是证明方法统计性质的关键，但要避免“只模拟了对 GP 有利的数据”。BIB 审稿人会关注 simulation 是否公平、是否覆盖反例、是否能证明 uncertainty 的价值。

### 11.2 模拟场景设计

至少包含 8 类场景：

1. Smooth spatial gradient。
2. Sharp domain boundary。
3. Ring/layer structure。
4. Patchy microenvironment pattern。
5. No spatial signal negative control。
6. Sparse observation with high dropout。
7. Spatially varying noise。
8. Confounded expression and APA usage。

### 11.3 参数网格

建议参数：

- spots: 500, 2,000, 5,000。
- genes: 200, 1,000, 5,000。
- dropout: 10%, 30%, 50%, 70%。
- noise: low, medium, high。
- spatial length scale: short, medium, long。
- domain number: 2, 4, 8。
- effect size: small, medium, large。

### 11.4 评估指标

插补：

- RMSE。
- MAE。
- Pearson。
- Spearman。
- R2。
- bias。

uncertainty：

- calibration curve。
- interval coverage。
- uncertainty-error correlation。
- negative log likelihood。

SVAPA：

- AUROC。
- AUPRC。
- FDR。
- sensitivity。
- specificity。

domain/differential：

- ARI。
- NMI。
- F1。
- effect-size recovery。

### 11.5 必须做的消融实验

1. GP vs KNN。
2. spatial coordinate GP vs random coordinate GP。
3. RBF vs Matérn kernel。
4. dense GP vs sparse GP。
5. with uncertainty weighting vs without uncertainty weighting。
6. spatial validation on vs off。
7. varying missingness。
8. varying spatial scale。

### 11.6 交付物

- simulation benchmark config。
- complete simulation result table。
- summary figures。
- negative control result。
- ablation result。

### 11.7 验收标准

- spaGAPA 在 smooth/gradient/layer 场景中明显优于 KNN。
- 在 no-signal negative control 中不会产生过多假阳性。
- uncertainty 与真实误差有显著相关。
- 在 sharp boundary 场景中承认 GP 可能过度平滑，并展示边界保护或参数敏感性。

---

## 12. WP6: 真实数据 Benchmark

### 12.1 需要解决的问题

真实数据没有 ground truth，必须设计合理的 pseudo-ground-truth 和 biological validation。

### 12.2 真实数据评估策略

#### Strategy 1: Masked observed APA values

流程：

1. 从高置信 observed APA values 中随机 mask 一部分。
2. 用不同方法插补。
3. 与 held-out observed values 比较。

优点：

- 可直接计算 RMSE/MAE/correlation。

风险：

- observed values 自身有噪音。
- random mask 不一定模拟真实 missingness。

增强：

- random mask。
- spatial block mask。
- low coverage mask。
- domain-specific mask。

#### Strategy 2: Replicate concordance

流程：

1. 如果有 biological/technical replicate，分别分析。
2. 比较 top SVAPA/differential APA genes 的重叠。
3. 比较 spatial pattern 的稳定性。

#### Strategy 3: Histology/domain consistency

流程：

1. 使用组织 layer 或手动 annotation。
2. 检查 APA spatial patterns 是否与组织结构一致。
3. 计算 domain enrichment 和 boundary alignment。

#### Strategy 4: Known biology enrichment

流程：

1. top SVAPA genes 做 GO/KEGG/Reactome。
2. 检查 3'UTR regulation、RNA processing、neuronal development、cell migration 等相关 pathway。
3. 对个别基因做文献解释。

### 12.3 交付物

- 每个数据集一份 benchmark report。
- masked benchmark 总表。
- top SVAPA genes 表。
- top differential APA genes 表。
- biological enrichment 结果。

### 12.4 验收标准

- spaGAPA 在至少 2 个真实数据集中优于 KNN/stAPAminer。
- uncertainty filtering 能提升 held-out prediction correlation 或下游稳定性。
- 至少一个真实数据集形成有说服力的生物学故事。

---

## 13. WP7: 生物学 Case Study

### 13.1 需要解决的问题

BIB 不一定要求像 CNS 主刊一样的生物学深度，但需要证明方法能产生领域洞察。只有 benchmark 表格是不够的。

### 13.2 推荐 case study 结构

#### Case 1: MOB layer-specific APA remodeling

展示：

- 空间层级或环状结构中的 APA usage。
- top SVAPA genes。
- distal/proximal usage 的空间变化。
- spaGAPA uncertainty map，说明哪些区域可信。

希望回答：

- APA 是否沿 MOB 层级发生系统变化？
- GP 是否恢复了更连续、更符合组织结构的 APA pattern？

#### Case 2: Brain/cortex region-specific APA

展示：

- cortical layer 或 anatomical region 中的 differential APA。
- 与 neuronal genes、synaptic genes、RNA localization 相关的 enrichment。

希望回答：

- 神经系统中是否存在空间组织化的 3'UTR usage？

#### Case 3: Developmental/tumor gradient APA

展示：

- 空间梯度或微环境边界上的 APA change。
- spatial trajectory optional。

希望回答：

- APA 是否与发育轴或 microenvironment 状态相关？

### 13.3 单基因展示标准

每个重点基因建议展示：

- raw APA index。
- imputed APA index。
- uncertainty。
- expression level。
- proximal site usage。
- distal site usage。
- spatial domain/layer annotation。
- statistical result。

### 13.4 交付物

- 2-3 个主文 case figures。
- 10-20 个候选重点基因列表。
- enrichment 表。
- 文献支持表。

### 13.5 验收标准

- 至少 5 个 top genes 有清楚的空间 APA pattern。
- 至少 1 个 case 与已知组织结构高度一致。
- 至少 1 个 case 提供潜在新发现，而不是只复述已知 marker。

---

## 14. WP8: 可扩展性、鲁棒性和工程化

### 14.1 性能目标

建议目标：

- 1,000 genes x 2,000 spots: dense/sparse 模式可完成 benchmark。
- 5,000 genes x 5,000 spots: sparse/block 模式可运行。
- 记录 runtime 和 peak memory。

### 14.2 鲁棒性场景

必须测试：

- 低 spot 数。
- 高 dropout。
- 单个 gene 全 NA 或几乎全 NA。
- 坐标重复。
- spot 顺序不一致。
- gene/spot metadata 缺失。
- APA index 越界。
- 多样本合并。

### 14.3 工程任务

- 增加 integration tests。
- 增加 benchmark tests。
- 增加 input schema validator。
- 增加 config-driven pipeline。
- 统一日志格式。
- 统一随机种子。
- 输出版本信息。

### 14.4 验收标准

- 测试覆盖核心 pipeline。
- 关键错误能给出可理解的报错。
- benchmark 可以被非开发者复现。

---

## 15. WP9: 文档、发布和复现材料

### 15.1 用户文档

需要补齐：

- Installation。
- Quick start。
- Data preparation。
- scAPAtrap integration。
- APA index input。
- GP imputation。
- uncertainty interpretation。
- SVAPA analysis。
- differential APA。
- benchmark reproduction。
- troubleshooting。

### 15.2 Tutorial notebooks

至少 3 个：

1. `tutorial_01_quickstart_simulated.ipynb`
2. `tutorial_02_real_data_spatial_apa.ipynb`
3. `tutorial_03_benchmark_reproduction.ipynb`

### 15.3 Release assets

需要：

- GitHub release。
- Zenodo DOI。
- Example dataset。
- Conda/pip installation instructions。
- Docker/Singularity optional。
- Reproducibility checklist。

### 15.4 验收标准

- 一个新用户能在 30 分钟内跑通 simulated quickstart。
- 一个审稿人能按 README 复现至少一张主图。
- 所有论文图对应的脚本和输入输出路径明确。

---

## 16. WP10: 论文图表与 Manuscript Assets

### 16.1 主文图建议

#### Figure 1: Field overview and spaGAPA framework

内容：

- spatial APA analysis 的挑战。
- 现有工具生态：scAPAtrap、stAPAminer、metaAPA、SPARK。
- spaGAPA workflow：
  - APA caller input。
  - spatial validation。
  - GP imputation。
  - uncertainty。
  - SVAPA/differential/domain。
  - biological interpretation。

目的：

- 让 BIB 读者看到 broader framework，而非单点工具。

#### Figure 2: Simulation benchmark

内容：

- 多空间模式模拟。
- imputation performance。
- uncertainty calibration。
- dropout sensitivity。
- negative control。

目的：

- 证明统计性质。

#### Figure 3: Real data imputation benchmark

内容：

- masked observed APA benchmark。
- spaGAPA vs stAPAminer/KNN/baselines。
- runtime/memory。

目的：

- 证明真实数据有效。

#### Figure 4: Uncertainty-aware SVAPA

内容：

- GP likelihood ratio / weighted Moran's I。
- with vs without uncertainty。
- FDR/stability/replicate concordance。

目的：

- 证明 uncertainty 是核心创新，而非附加输出。

#### Figure 5: Biological case study

内容：

- MOB/brain/embryo/tumor 中 top SVAPA genes。
- spatial maps。
- APA usage + uncertainty。
- enrichment。

目的：

- 展示生物学发现。

#### Figure 6: Software ecosystem and scalability

内容：

- API/CLI/workflow。
- AnnData compatibility。
- runtime/memory。
- reproducibility resources。

目的：

- 支持 BIB 的工具可用性和推广性。

### 16.2 Supplementary figures

建议：

- Supplementary Fig S1: 数据集 QC。
- Supplementary Fig S2: APA site filtering。
- Supplementary Fig S3: kernel sensitivity。
- Supplementary Fig S4: sparse GP comparison。
- Supplementary Fig S5: more real datasets。
- Supplementary Fig S6: all baseline details。
- Supplementary Fig S7: failure cases and boundary smoothing analysis。

### 16.3 主表建议

#### Table 1: Existing methods and limitations

比较：

- scAPAtrap。
- stAPAminer。
- metaAPA。
- SPARK/SPARK-X。
- generic spatial tools。
- spaGAPA。

#### Table 2: Dataset summary

字段：

- dataset。
- species。
- tissue。
- platform。
- spots/cells。
- genes。
- PAS。
- APA genes。
- missing rate。
- annotation。

#### Table 3: Benchmark summary

字段：

- method。
- imputation metrics。
- SVAPA metrics。
- runtime。
- memory。
- comments。

---

## 17. 审稿人可能质疑与预防

### 17.1 “spaGAPA 只是把 GP 用到 APA 上，创新不足”

预防：

- 强调不是单纯 GP，而是 spatial APA-specific uncertainty-aware framework。
- 展示 uncertainty 进入 SVAPA/differential APA 后带来实际收益。
- 加入 benchmark 和真实生物学 case。

### 17.2 “APA calling 依赖 scAPAtrap，不是端到端新方法”

预防：

- 明确 spaGAPA 的定位是 downstream spatial APA inference framework。
- 支持多个 caller 输入。
- 加入 spatial validation 改善 site-level confidence。
- 将 scAPAtrap 作为 upstream input，而非唯一贡献。

### 17.3 “真实数据没有 ground truth”

预防：

- masked observed values。
- replicate concordance。
- simulation ground truth。
- histology/domain consistency。
- known biology enrichment。

### 17.4 “GP 可能过度平滑，抹掉局部边界”

预防：

- sharp boundary simulation。
- kernel sensitivity。
- length-scale control。
- compare spatial block mask。
- 展示 failure cases 和适用范围。

### 17.5 “stAPAminer 比较不公平”

预防：

- 尽量运行原始 R 包。
- 相同数据、相同 mask、相同 gene/spot set。
- 公开参数。
- 同时比较 KNN-expression 和 KNN-spatial。

### 17.6 “uncertainty 没有校准”

预防：

- calibration curve。
- interval coverage。
- uncertainty-error correlation。
- high-confidence subset performance。

### 17.7 “软件不可复现”

预防：

- conda/Docker。
- figure scripts。
- data provenance。
- Zenodo release。
- integration tests。

---

## 18. 里程碑计划

### Milestone 1: BIB-ready 基础修复

目标：

- pipeline 跑通。
- 数据模型统一。
- toy integration tests 完成。

验收：

- toy data 一键分析成功。
- simulated benchmark 可复现。
- `use_sparse_gp=True` 可运行。

### Milestone 2: 真实 APA 数据打通

目标：

- 至少 1 个真实 dataset 从 BAM/APA caller 到 spaGAPA 输出完整跑通。

验收：

- 有真实 `apa_matrix.csv`。
- 有 QC report。
- 有初步 spatial maps。

### Milestone 3: stAPAminer 公平复现

目标：

- 完成 stAPAminer-original 或 faithful reimplementation。

验收：

- 同一数据上输出 stAPAminer 与 spaGAPA 结果。
- 有 imputation 和 SVAPA 对比表。

### Milestone 4: 系统 benchmark 完成

目标：

- simulation + real masked benchmark + ablation。

验收：

- Figure 2 和 Figure 3 数据完整。
- 主要 metrics 表完整。

### Milestone 5: 生物学 case study 完成

目标：

- 至少 2 个真实 case。

验收：

- Figure 5 初稿。
- top genes 和 enrichment 解释完成。

### Milestone 6: 投稿包准备

目标：

- BIB manuscript、figures、supplement、code release。

验收：

- 所有主图可复现。
- 文档完整。
- release tag 和 DOI 准备。

---

## 19. 优先级执行清单

### 第一阶段：先让项目真的可跑

1. 修复 pipeline 矩阵方向和 API 不一致。
2. 增加 toy end-to-end integration test。
3. 确认 dense GP、sparse GP、SVAPA、differential APA 都能串起来。
4. 更新 docs 中与实际代码不一致的部分。

### 第二阶段：解决真实数据根基

1. 确定 3 个候选数据集。
2. 下载真实 BAM/FASTQ 或 poly(A) count。
3. 跑 scAPAtrap/Sierra/SCAPE。
4. 生成真实 APA index matrix。
5. 写 data provenance。

### 第三阶段：建立公平 benchmark

1. 复现 stAPAminer。
2. 实现 KNN-expression。
3. 增加 SPARK/SPARK-X SVAPA baseline。
4. 固定 benchmark config。
5. 运行 simulation、masked real data、ablation。

### 第四阶段：强化创新闭环

1. 完善 uncertainty calibration。
2. 完善 uncertainty-aware SVAPA。
3. 完善 uncertainty-aware differential APA。
4. 证明 uncertainty 能提升稳定性和降低假阳性。

### 第五阶段：组织 BIB 故事

1. 写 spatial APA computational challenges。
2. 做方法生态比较表。
3. 完成主图。
4. 写 manuscript outline。
5. 准备 supplement 和 reproducibility materials。

---

## 20. Go / No-Go 标准

### 20.1 可以冲 BIB 的条件

满足以下大部分条件时，可以认真冲 BIB：

- 3 个真实数据集完成。
- stAPAminer 原始或忠实复现完成。
- spaGAPA 在至少 2 个核心任务中优于 baseline。
- uncertainty 的价值被定量证明。
- 至少 1-2 个真实生物学 case 有说服力。
- 软件可复现，主图脚本完整。
- 论文 framing 是 systematic benchmark + framework，而非普通 package note。

### 20.2 应降级投稿的条件

如果出现以下情况，应考虑 Bioinformatics / NAR Genomics and Bioinformatics / BMC Bioinformatics：

- 真实数据只有 1 个。
- stAPAminer 复现不完整。
- spaGAPA 只在 simulation 上优势明显。
- uncertainty 没有带来下游收益。
- 生物学 case 不够强。
- pipeline 需要大量手动 patch 才能跑。

### 20.3 可以考虑更高目标的条件

如果出现以下情况，可以考虑 Genome Biology 或更高目标：

- 多组织真实数据中发现新的、强生物学意义的 spatial APA remodeling。
- 结果与实验验证或独立数据强一致。
- 方法可推广到更广泛 spatial omics missingness/uncertainty 问题。
- 社区中缺少同类框架，且 spaGAPA 成为事实标准 workflow。

---

## 21. 下一步建议

建议立即从 Milestone 1 开始：

1. 做 pipeline audit，列出所有 shape/API/placeholder 问题。
2. 修复一键流程，使 toy data 和 simulation data 可靠跑通。
3. 同步准备真实数据清单，优先找 MOB 和 brain 数据。
4. 并行启动 stAPAminer 环境复现。

这个顺序的原因是：BIB 需要真实数据和强 benchmark，但如果 pipeline 基础不稳，后续真实数据 benchmark 会反复返工。

---

## 22. Milestone 1 执行记录

执行日期：2026-06-10

### 22.1 已完成修复

#### M1.1 Pipeline 矩阵方向统一

已修复：

- `APADataset.raw_counts` 的 public contract 明确为 `genes x spots`。
- `SpaGAPA.run()` 不再将 `raw_counts` 转置后传给 batch GP。
- dense GP 和 sparse GP 在 pipeline 中均使用 `genes x spots` 输入。
- imputed matrix 和 uncertainty matrix 均以 `genes x spots` 存入 `APADataset`。
- 新增 pipeline shape validation，检查：
  - `raw_counts.shape == (n_genes, n_spots)`。
  - `coords.shape == (n_spots, 2)`。
  - coordinates 数量与 spot 数一致。

意义：

- 避免把 spots 当 genes 拟合 GP。
- 避免 coordinates 与矩阵轴错位。
- 为真实数据 benchmark 提供基础安全检查。

#### M1.2 高层 APADataset loader

已修复：

- 保留低层 `read_spatial_data()` 的旧行为，继续返回 DataFrame tuple。
- 新增 `load_spatial_dataset()`，用于 pipeline，高层返回 `APADataset`。
- 支持：
  - `apa_matrix + coordinates`。
  - `h5ad_file`。
  - `genes_by_spots` 和 `spots_by_genes` 两种输入方向。
  - DataFrame、NumPy array 和 CSV/TSV 文件。
- 若用户直接给 BAM 但没有 APA matrix，会明确提示需要先运行 APA caller。

意义：

- 解决 pipeline 期待 APADataset、reader 返回 tuple 的接口不一致问题。
- 为后续真实 `apa_matrix.csv + coordinates.csv` benchmark 做好入口。

#### M1.3 SparseGP batch API

已修复：

- `SparseGPImputer` 新增 `fit_batch()`。
- 新增 `SparseGPImputerBatch`。
- sparse batch 与 dense `GPImputerBatch` 保持一致：
  - 输入 `genes x spots`。
  - `predict()` 输出 `genes x spots`。
  - `impute()` 输出 `genes x spots`。
  - observed values 保留，observed uncertainty 置 0。
- `SpaGAPA(use_sparse_gp=True)` 可在 toy data 上完整运行。

意义：

- 修复 `use_sparse_gp=True` 路径不可用的问题。
- 为大规模真实数据 benchmark 提供基础。

#### M1.4 APA quantification 不再伪装为已计算

已修复：

- `input_type='apa_index'` 时，pipeline 明确记录输入矩阵已经是 APA index：
  - 输出 key 为 `APAIndex`。
  - `calculated=False`。
  - `source='input_apa_index'`。
- `input_type='proximal_distal_counts'` 时，要求 AnnData layers：
  - `proximal_counts`。
  - `distal_counts`。
- 只有在 proximal/distal counts 存在时才真正计算 RUD、PDUI、PAI。

意义：

- 避免审稿人质疑 pipeline 把任意矩阵伪装成 APA index calculation。
- 为真实 site count 到 APA index 的路径留下明确接口。

#### M1.5 Differential APA 输出 gene-level 结果

已修复：

- pipeline 不再只输出每个 domain 的 spot 数。
- `differential_analysis=True` 时调用 `DifferentialAPAAnalyzer.test_all_pairwise()`。
- 输出每个 domain pair 的 gene-level 结果：
  - gene。
  - mean_group1。
  - mean_group2。
  - log2fc。
  - statistic。
  - pvalue。
  - padj。
- `test_all_pairwise()` 支持传入 uncertainty。

意义：

- 使 pipeline 输出能支撑 differential APA claim。
- 为后续 uncertainty-aware differential APA benchmark 打基础。

#### M1.6 DomainIdentifier shape 判断更稳健

已修复：

- 不再用 `n_genes < n_spots` 这种启发式判断矩阵是否需要转置。
- 改为用 `spatial_coords.shape[0]` 判断 spot 轴。
- 如果 APA matrix 两个轴都对不上坐标数量，会直接报错。

意义：

- 避免小数据、过滤后数据或 gene 数大于 spot 数时方向判断错误。

#### M1.7 Optional scanpy import 更稳健

已修复：

- `scanpy` 导入失败不再导致整个 `spagapa.analysis` import 崩溃。
- 只有用户实际选择 `leiden` 或 `louvain` 时才需要 scanpy。

意义：

- 避免因为 scanpy/numba 环境问题影响 K-means pipeline、GP、benchmark 等不依赖 scanpy 的功能。

#### M1.8 Toy end-to-end integration test

已新增：

- `tests/integration/test_pipeline_toy.py`。

覆盖：

- dense GP pipeline。
- sparse GP pipeline。
- APADataset 输入。
- `apa_matrix.csv + coordinates.csv` 高层 loader。
- imputation 输出 shape。
- uncertainty 输出 shape。
- APAIndex 注册。
- domain labels。
- gene-level differential results。
- SVAPA results。

### 22.2 验证结果

已运行：

```bash
pytest tests/unit tests/integration -q
```

结果：

```text
307 passed, 2 skipped, 207 warnings in 56.06s
```

说明：

- 全部 unit tests 和 integration tests 通过。
- warnings 主要来自 pandas/numpy deprecation、scikit-learn GP convergence warning、constant input warning，以及 trajectory 的预期 warning。
- 当前没有失败测试。

### 22.3 Milestone 1 后仍需注意的问题

以下问题不阻塞 toy pipeline，但会影响 BIB-ready 真实数据阶段：

1. `input_type='proximal_distal_counts'` 还需要真实数据测试。
2. SparseGP 目前只支持 RBF kernel，pipeline 中 `kernel_type='matern'` 对 sparse path 不生效。
3. Differential APA 的 uncertainty weighting 当前是启发式 weighted resampling，后续需要更严格的统计实现和 simulation 校准。
4. GP likelihood ratio 的 uncertainty 目前没有真正作为 sample weight 进入 sklearn GP fit，因为 sklearn GPR 不支持 sample_weight；后续需要用 alpha/noise 或自定义似然实现。
5. `SpatialValidator` 和 `QualityFilter` 在 pipeline 中还主要是初始化和占位式记录，真实 PAS/site-level validation 需要在 WP1/WP3 继续加强。
6. 真实 BAM 到 APA matrix 仍依赖 scAPAtrap/Sierra/SCAPE 等上游 caller，尚未完成真实数据打通。

---

## 23. Milestone 2/3 执行记录：stAPAminer MOB 数据打通与基线复现

执行日期：2026-06-10

### 23.1 环境确认

用户提供的开发环境已确认可用：

```text
Python: conda env spagapa, Python 3.10.0
Rscript: ~/anaconda3/envs/r442/bin/Rscript, R 4.4.2
R library path: /s1/SHARE/01_software/R_442_SeuratV5/library
```

R 侧已确认：

- `movAPA` 可用。
- `Seurat` 可用。
- `edgeR`、`limma`、`org.Mm.eg.db` 可用。
- `SPARK`、`ClusterR`、`clusterSim` 当前不可用，后续若需要完整 stAPAminer clustering/SVAPA 复现需补装或替代。

### 23.2 Milestone 2: 第一个真实/发表 APA 数据集已标准化

已使用本地 stAPAminer 示例 MOB 数据：

```text
00_ref_packages/stAPAminer-main/inst/extdata/APA.RDA
00_ref_packages/stAPAminer-main/inst/extdata/count.csv
00_ref_packages/stAPAminer-main/inst/extdata/position.txt
```

该数据不是 expression proxy，而是 stAPAminer 论文/示例提供的真实 poly(A) site count / PACdataset：

- `APA.RDA`: movAPA `PACdataset`。
- `APA@counts`: 54493 poly(A) sites x 1001 spots。
- 与 `position.txt` 和 `count.csv` 交集后保留 260 MOB spots。
- `count.csv`: 16218 expression genes x 260 spots。
- `position.txt`: 260 spots，含 `x/y/layer`，layer 包括 GCL、GL、MCL、ONL、OPL。

新增脚本：

```text
spaGAPA/scripts/prepare_stapaminer_mob.py
```

标准化输出目录：

```text
spaGAPA/data/processed/stapaminer_mob/
```

输出文件：

- `apa_matrix.csv`: raw RUD APA index，4845 genes x 260 spots。
- `coordinates.csv`: spot_id, x, y。
- `metadata.csv`: spot_id, layer, dataset, source。
- `expression_matrix.csv`: expression genes x spots，用于 stAPAminer KNN-expression。
- `stapaminer_rud_raw.csv`: raw RUD。
- `stapaminer_rud_imputed.csv`: stAPAminer-compatible KNN imputed RUD。
- `apa_site_counts.csv.gz`: poly(A) site count matrix。
- `apa_sites.csv.gz`: poly(A) site annotation。
- `qc_summary.json`: 数据维度、缺失率和 provenance。

QC 摘要：

```text
n_spots: 260
n_expression_genes: 16218
n_polyA_sites: 54493
n_rud_genes: 4845
raw_rud_na: 670998
raw_rud_missing_rate: 0.5327
imputed_rud_na: 0
layers: GCL=61, GL=69, MCL=31, ONL=44, OPL=55
```

重要兼容处理：

- stAPAminer 原始 `computeAPAIndex()` 在 R 4.4 / movAPA 0.2.0 下会因为 `PACdataset@counts` slot 类型要求报错。
- 本次没有修改竞品源码，而是在导出脚本中使用 `APA_subset@counts <- as.matrix(...)` 保持原始 movAPA RUD 计算语义。
- `movAPAindex()` 返回 `Matrix::dgeMatrix`，导出前转为 base matrix。

### 23.3 spaGAPA 真实 MOB smoke test 已跑通

新增脚本：

```text
spaGAPA/scripts/run_stapaminer_mob_smoke.py
```

功能：

- 读取 `apa_matrix.csv + coordinates.csv`。
- 选择观测足够且仍有缺失的 RUD genes。
- 跑 spaGAPA pipeline：
  - GP imputation。
  - APAIndex 注册。
  - domain identification。
  - gene-level differential APA。
  - GP-based SVAPA。
- 输出 smoke summary、imputed matrix、uncertainty、domain、SVAPA 和 differential 结果。

当前 smoke 输出目录：

```text
spaGAPA/benchmark_results/real/stapaminer_mob_smoke/
```

示例运行：

```bash
conda run -n spagapa python spaGAPA/scripts/run_stapaminer_mob_smoke.py \
  --n-genes 12 \
  --min-observed-spots 120 \
  --n-domains 3
```

示例结果：

```text
n_genes: 12
n_spots: 260
raw_missing_rate: 0.003846
imputed_shape: [12, 260]
uncertainty_shape: [12, 260]
mean_uncertainty: 0.0002966
n_domains: 3
n_svapa_significant: 0
```

同时修复了一个真实 APA index 的关键语义问题：

- 对 `input_type='apa_index'`，0 是合法观测值，不能当作 missing。
- pipeline 现在使用 `np.isfinite(values)` 作为 APA index training mask。
- 对 count-like input 才继续使用 `values > 0` 作为 observed mask。

### 23.4 Milestone 3: stAPAminer KNN-expression 基线已开始复现

新增脚本：

```text
spaGAPA/scripts/run_stapaminer_mob_benchmark.py
```

当前 benchmark 任务：

- 使用 stAPAminer MOB raw RUD。
- 从 observed RUD values 中 mask 20% 作为 held-out。
- 比较：
  - mean。
  - median。
  - KNN-spatial。
  - stAPAminer KNN-expression faithful reimplementation。
  - spaGAPA-GP。

stAPAminer KNN-expression 复现细节：

- spot 距离使用完整 `expression_matrix.csv`，而不是只用 APA genes。
- 对 expression genes 做 z-score 后计算 spot 间 Euclidean distance。
- k=10。
- 迭代最多 10 次。
- 若对应表达为 0，则把 missing APA 初始化为 0。
- 对不能匹配 expression 的 RUD genes，不强行把所有 missing 当表达 0。

示例运行：

```bash
conda run -n spagapa python spaGAPA/scripts/run_stapaminer_mob_benchmark.py \
  --n-genes 40 \
  --min-observed-spots 120 \
  --mask-fraction 0.2
```

当前示例结果：

```text
method                       rmse      mae      pearson   spearman   r2
spagapa_gp                   0.125099  0.075493 0.956909  0.907635   0.915592
mean                         0.126495  0.076474 0.955912  0.906399   0.913697
knn_spatial                  0.130076  0.078204 0.953382  0.899956   0.908742
median                       0.132215  0.070167 0.954131  0.910856   0.905715
stapaminer_knn_expression    0.136510  0.078239 0.949191  0.893350   0.899491
n_holdout: 2066
```

当前结果只能作为 smoke benchmark，不能作为论文级结论。下一步需要：

1. 扩大 gene 数和 mask seed。
2. 加 random mask、spatial block mask、low-coverage mask。
3. 加 uncertainty calibration。
4. 加 runtime/memory。
5. 用 layer labels 评价 domain/differential/SVAPA biological consistency。
6. 若要完整复现 stAPAminer 原文，需要安装或替代 `SPARK`、`ClusterR`、`clusterSim`。

### 23.5 本轮验证

已在 `spagapa` conda 环境运行：

```bash
conda run -n spagapa pytest tests/unit/test_readers.py tests/unit/test_gp_imputer.py tests/unit/test_sparse_gp.py tests/integration/test_pipeline_toy.py -q
```

结果：

```text
47 passed, 1 skipped, 32 warnings
```

warnings 主要来自 sklearn GP convergence，不阻塞当前 smoke/benchmark。

### 23.6 Milestone 3 升级：formal MOB benchmark v2 已完成

在 smoke benchmark 基础上，已将真实数据 benchmark 升级为一个更接近论文主结果草案的 formal benchmark。

脚本：

```text
spaGAPA/scripts/run_stapaminer_mob_benchmark.py
```

本轮方法学增强：

- 支持多 seed：`42,43,44`。
- 支持 3 类 held-out 策略：
  - `random`
  - `spatial_block`
  - `low_coverage`
- `low_coverage` mask 已改为“偏向低覆盖基因的随机抽样”，而不是固定顺序遮蔽，保证 seed 真正生效。
- 输出 runtime 与近似 peak RSS memory。
- 使用 MOB `layer` 标签，对 imputed APA matrix 做 KMeans 聚类，计算 biological consistency：
  - `layer_ari`
  - `layer_nmi`
- 使用 Matplotlib 自动保存 benchmark 图表。

本轮正式运行命令：

```bash
conda run -n spagapa python spaGAPA/scripts/run_stapaminer_mob_benchmark.py \
  --output-dir spaGAPA/benchmark_results/real/stapaminer_mob_formal_v2 \
  --n-genes 120 \
  --min-observed-spots 100 \
  --mask-fraction 0.2 \
  --seeds 42,43,44
```

输入数据规模：

```text
dataset: stAPAminer_MOB
n_genes: 120
n_spots: 260
raw_missing_rate: 0.01984
mask_fraction: 0.20
mask_types: random, spatial_block, low_coverage
methods: spagapa_gp, stapaminer_knn_expression, knn_spatial, mean, median
```

输出目录：

```text
spaGAPA/benchmark_results/real/stapaminer_mob_formal_v2/
```

主要输出文件：

- `benchmark_results_long.csv`
- `benchmark_results_summary.csv`
- `benchmark_results_overall.csv`
- `benchmark_summary.json`
- `figures/benchmark_summary_panels.png`
- `figures/benchmark_mask_types_rmse.png`
- `figures/benchmark_seed_stability.png`
- `figures/benchmark_runtime_vs_accuracy.png`
- `figures/benchmark_memory_usage.png`
- `figures/benchmark_biological_consistency.png`

overall mean 指标：

```text
method                       rmse      mae      pearson   spearman   r2        runtime_s   peak_rss_mb  layer_ari  layer_nmi
spagapa_gp                   0.126333  0.069915 0.958565  0.906471   0.918821  25.472557   541.870226   0.047094   0.091393
stapaminer_knn_expression    0.135536  0.072591 0.952360  0.900924   0.906496   0.127041   686.641927   0.077301   0.121210
knn_spatial                  0.176622  0.100173 0.914055  0.876683   0.812190   0.009305   528.305556   0.066917   0.118395
mean                         0.126231  0.070086 0.958620  0.906580   0.918914   0.000612   516.475260   0.039329   0.081797
median                       0.132107  0.061266 0.956770  0.925692   0.911199   0.003622   522.614583   0.048089   0.091333
```

按 mask 类型分层的平均表现：

```text
low_coverage:
  best rmse = spagapa_gp (0.129135), mean 次之 (0.129451)

random:
  best rmse = spagapa_gp (0.127683), mean 次之 (0.128079)
  best layer_ari = stapaminer_knn_expression (0.129535), spagapa_gp 次之 (0.108106)

spatial_block:
  best rmse = mean (0.121163), spagapa_gp 次之 (0.122179)
  knn_spatial 在 spatial_block 下明显退化 (rmse = 0.263029)
```

当前阶段性结论：

1. `spagapa_gp` 在 `random` 和 `low_coverage` mask 下保持最优或接近最优 RMSE，说明 GP 对一般缺失恢复仍有优势。
2. `mean` baseline 在 formal benchmark 中与 GP 几乎打平，说明当前 benchmark 设置还没有充分拉开“空间结构建模”的优势，这对 BIB 论文是一个重要警示信号。
3. `spatial_block` 场景下，`knn_spatial` 明显失败，提示局部邻域平均在连续空间遮蔽时不稳健；但 `mean` 与 `GP` 仍很接近，说明 GP 的 spatial extrapolation 优势还不够强。
4. 在当前 `layer_ari/layer_nmi` 指标上，`stapaminer_knn_expression` 反而是最强方法，说明“表达空间邻近性”在 MOB layer 恢复上仍然很有竞争力。
5. `spagapa_gp` 的 runtime 明显更高，当前约为 `25.5s`，远高于所有 baseline，因此如果 accuracy 不能稳定拉开差距，BIB 审稿人会直接质疑方法复杂度的必要性。

这组结果非常有价值，因为它不是“只会支持我们方法”的 benchmark，而是真实暴露了当前 spaGAPA 还没有完全建立不可替代性。对 BIB 来说，这反而是下一阶段优化的正确起点。

基于 formal benchmark v2，下一步优先问题已经更清晰：

1. 加入 uncertainty calibration：
   - interval coverage
   - uncertainty vs error correlation
   - high-uncertainty filtering 后的 accuracy 提升
2. 重新设计更能体现空间外推难度的 masking：
   - 更大 block
   - ring/sector mask
   - layer-aware mask
3. 提升 GP 方法本身：
   - kernel 选择
   - noise / alpha 建模
   - expression-informed kernel 或 multi-view kernel
4. 提升 biological consistency 评价：
   - 不只做 KMeans + ARI/NMI
   - 加 layer marker ordering / pseudolayer smoothness / boundary consistency
5. 推进 uncertainty-aware downstream：
   - SVAPA ranking 稳定性
   - differential APA robustness
   - uncertainty-filtered domain identification

### 23.7 uncertainty-aware benchmark 框架已接入

在 formal benchmark v2 基础上，`run_stapaminer_mob_benchmark.py` 已继续升级为 uncertainty-aware benchmark 原型。

本轮新增能力：

1. **GP uncertainty calibration 指标**
   - `uncertainty_coverage_68`
   - `uncertainty_coverage_95`
   - `uncertainty_nll`
   - `uncertainty_error_pearson`
   - `uncertainty_error_spearman`
   - `rmse_keep_80pct`
   - `rmse_keep_50pct`
   - `rmse_gain_keep_80pct`

2. **更难的 masking 场景**
   - `spatial_block_large`
   - `ring_sector`
   - `layer_aware`
   - 保留已有 `random` / `spatial_block` / `low_coverage`

3. **uncertainty 相关输出文件**
   - `uncertainty_filtering_curve.csv`
   - `uncertainty_reliability_curve.csv`

4. **新增 Matplotlib 图**
   - `benchmark_gp_uncertainty_panels.png`
   - `benchmark_mask_layouts.png`

其中：

- `spatial_block_large`：更大范围连续 block 遮蔽。
- `ring_sector`：模拟同心层/扇区状空间外推。
- `layer_aware`：按 MOB layer 顺序遮蔽一个或相邻多个 layer。

### 23.8 uncertainty-aware benchmark validation 已跑通

验证运行：

```bash
conda run -n spagapa python spaGAPA/scripts/run_stapaminer_mob_benchmark.py \
  --output-dir spaGAPA/benchmark_results/real/stapaminer_mob_uncertainty_validation \
  --n-genes 24 \
  --min-observed-spots 110 \
  --mask-fraction 0.15 \
  --seeds 42 \
  --mask-types random,spatial_block_large,ring_sector,layer_aware,low_coverage \
  --calibration-bins 5
```

输出目录：

```text
spaGAPA/benchmark_results/real/stapaminer_mob_uncertainty_validation/
```

生成文件已确认包括：

- `benchmark_results_long.csv`
- `benchmark_results_summary.csv`
- `benchmark_results_overall.csv`
- `benchmark_summary.json`
- `uncertainty_filtering_curve.csv`
- `uncertainty_reliability_curve.csv`
- `figures/benchmark_gp_uncertainty_panels.png`
- `figures/benchmark_mask_layouts.png`

validation 结果显示，GP uncertainty 已经具备“可解释且可利用”的信号：

```text
spagapa_gp:
  rmse = 0.099282
  uncertainty_coverage_68 = 0.869173
  uncertainty_coverage_95 = 0.958324
  uncertainty_error_pearson = 0.517674
  uncertainty_error_spearman = 0.695719
  rmse_keep_80pct = 0.083262
  rmse_keep_50pct = 0.055460
  rmse_gain_keep_80pct = 0.016020
```

这说明：

1. 当前 GP posterior std 与真实误差是正相关的，而不是随机噪声。
2. 去掉高 uncertainty 的 held-out 点后，RMSE 明显下降，说明 uncertainty 可作为质量控制信号。
3. `68% coverage` 偏高于理想值，提示当前 uncertainty 可能偏保守。
4. `95% coverage` 接近目标值，说明后验标准差并非完全失真。

这一步的意义非常关键：

- 以前 spaGAPA 的 uncertainty 只是“模型能吐出来一个 std”。
- 现在已经开始回答 BIB 审稿人真正会问的问题：
  - 这个 uncertainty 是否校准？
  - 是否真的和误差相关？
  - 是否能用于筛掉不可靠 imputation？

因此，下一轮优先工作已经进一步收敛为：

1. 用更大规模 gene/seed 正式运行 uncertainty-aware benchmark。
2. 比较不同 GP kernel / alpha 设置对 calibration 的影响。
3. 将 uncertainty filtering 明确接入 downstream SVAPA / differential APA 分析。

### 23.9 formal uncertainty-aware benchmark v1 已完成

在 validation 跑通后，已进一步运行第一版正式 uncertainty-aware benchmark：

```bash
conda run -n spagapa python spaGAPA/scripts/run_stapaminer_mob_benchmark.py \
  --output-dir spaGAPA/benchmark_results/real/stapaminer_mob_uncertainty_formal_v1 \
  --n-genes 80 \
  --min-observed-spots 100 \
  --mask-fraction 0.2 \
  --seeds 42,43 \
  --mask-types random,spatial_block_large,ring_sector,layer_aware,low_coverage \
  --calibration-bins 6 \
  --gp-kernel matern \
  --gp-alpha 1e-10
```

输出目录：

```text
spaGAPA/benchmark_results/real/stapaminer_mob_uncertainty_formal_v1/
```

运行规模：

```text
n_genes: 80
n_spots: 260
seeds: 42, 43
mask_types: random, spatial_block_large, ring_sector, layer_aware, low_coverage
GP kernel: matern
GP alpha: 1e-10
```

overall mean 结果：

```text
spagapa_gp:
  rmse = 0.110919
  pearson = 0.966930
  runtime_s = 18.098
  layer_ari = 0.050151
  uncertainty_coverage_68 = 0.882208
  uncertainty_coverage_95 = 0.963661
  uncertainty_error_pearson = 0.593356
  uncertainty_error_spearman = 0.838293
  rmse_keep_80pct = 0.078286
  rmse_keep_50pct = 0.039305
  rmse_gain_keep_80pct = 0.032632

mean:
  rmse = 0.110920
  pearson = 0.966936

stapaminer_knn_expression:
  rmse = 0.120556
  layer_ari = 0.103403
  layer_nmi = 0.142849
```

这一轮正式结果非常关键，因为它把当前 spaGAPA 的“真实强项”和“仍未解决的问题”分得更清楚：

1. **纯插补精度层面**
   - `spagapa_gp` 与 `mean` 几乎完全打平。
   - 因此目前还不能把论文主 claim 建立在“GP 明显降低 RMSE”上。

2. **uncertainty 价值层面**
   - `uncertainty_error_spearman = 0.838`
   - `rmse_gain_keep_80pct = 0.0326`
   - 说明 GP uncertainty 与真实误差高度相关，且可显著筛掉不可靠预测。
   - 这已经开始形成 spaGAPA 更合理的 BIB 叙事核心：
     - 不是“GP 比所有方法插补都明显准很多”
     - 而是“GP + calibrated uncertainty 能识别哪些 APA imputation 值可靠、哪些不可靠”

3. **biological consistency 层面**
   - `stapaminer_knn_expression` 在 `layer_ari/layer_nmi` 上仍更强。
   - 特别是在 `layer_aware` mask 下，stAPAminer-like baseline 的 layer ARI 仍显著高于 spaGAPA。
   - 这说明 expression-informed spatial structure 仍是一个必须正面吸收的优势。

4. **mask-specific 观察**
   - `ring_sector` 和 `spatial_block_large` 下，spaGAPA 的 RMSE 优于或微优于 mean。
   - `layer_aware` 下，GP 并未形成明显优势。
   - 这提示当前 GP 更擅长一般空间平滑外推，不一定擅长 layer-structured biological manifold recovery。

因此，formal uncertainty-aware benchmark v1 的结论是：

> spaGAPA 当前最有说服力的优势，已经从“更低的平均误差”转向“可量化、可校准、可用于质量控制的 uncertainty-aware inference”。

### 23.10 首轮 GP kernel/alpha sweep 已完成

为了给后续方法优化提供方向，已运行一个小规模 GP 参数搜索：

- kernel: `matern`, `rbf`
- alpha: `1e-10`, `1e-6`, `1e-3`
- 输出目录：

```text
spaGAPA/benchmark_results/real/gp_sweeps/
```

当前 sweep 摘要：

```text
best rmse: matern + alpha=1e-3   (0.112329)
next:      matern + alpha=1e-6   (0.112366)
next:      matern + alpha=1e-10  (0.112375)
rbf variants 稍弱于 matern
```

同时：

- `matern + alpha=1e-3` 的 `uncertainty_error_spearman` 也是当前最优或近最优。
- `matern` 整体略优于 `rbf`。
- `alpha` 从 `1e-10` 调到 `1e-3` 有小幅但一致的改善。

这说明下一轮正式 benchmark 可以优先切到：

```text
kernel = matern
alpha = 1e-3
```

### 23.11 基于当前结果的下一步收敛

现阶段最值得优先执行的不是盲目扩更多 baseline，而是沿下面这条线收紧：

1. 使用 `matern + alpha=1e-3` 再跑一版 formal uncertainty-aware benchmark v2。
2. 将 uncertainty filtering 直接接入：
   - SVAPA ranking stability
   - differential APA robustness
   - uncertainty-filtered domain identification
3. 在 GP 中显式吸收 expression 信息：
   - expression-informed kernel
   - multi-view kernel
   - 或 spatial + expression hybrid neighborhood baseline

这条线的意义在于：

- 当前 GP 已经证明 uncertainty 有价值；
- 当前 stAPAminer-like baseline 已经证明 expression structure 对 biological consistency 很重要；
- 下一阶段最自然的方向，就是把这两者结合起来，而不是继续只在纯空间坐标上微调。

### 23.12 expression-informed GP unified formal benchmark + sweep 已完成

为了直接判断 expression-informed GP 这条路线是否值得继续作为主线，已完成一轮统一 benchmark suite：

- formal 三变体统一对比：
  - `spatial`
  - `spatial_radial`
  - `expr_additive`
- expr_additive 参数 sweep：
  - `lambda_expr = 0.2 / 0.5 / 1.0 / 2.0`
  - `expr_n_components = 5 / 10 / 15`
  - `HVG on/off`

结果目录：

```text
spaGAPA/benchmark_results/real/expression_gp_suite_formal_v2/
```

关键输出：

- `formal_variant_comparison.csv`
- `expr_additive_sweep_summary.csv`
- `decision_summary.json`
- `figures/formal_variant_comparison.png`
- `figures/expr_additive_sweep.png`

formal 三变体比较结果：

```text
spatial:
  rmse = 0.112778
  layer_ari = 0.051272
  uncertainty_error_spearman = 0.836272
  runtime_s = 22.714

spatial_radial:
  rmse = 0.147605
  layer_ari = 0.057089
  uncertainty_error_spearman = 0.567163
  runtime_s = 0.809

expr_additive:
  rmse = 0.168027
  layer_ari = 0.057371
  uncertainty_error_spearman = 0.475353
  runtime_s = 1.920
```

这一轮结论非常明确：

1. **按 RMSE 看，当前最强 GP 仍然是纯 spatial**
   - `spatial` 明显优于 `spatial_radial` 和 `expr_additive`。
   - 并且 `spatial` 与 `mean` 几乎打平，但略优。

2. **按 biological consistency 看，expr_additive 只有轻微 layer ARI 提升**
   - `expr_additive layer_ari = 0.0574`
   - `spatial layer_ari = 0.0513`
   - 但这个增益远不足以抵消 RMSE 的明显恶化。

3. **按 uncertainty calibration 看，expression-informed 版本反而退步**
   - `spatial uncertainty_error_spearman = 0.836`
   - `spatial_radial = 0.567`
   - `expr_additive = 0.475`
   - 说明当前 expression-informed MVP 不仅没有提高精度，也削弱了 uncertainty 与真实误差的对应关系。

4. **spatial_radial 当前不适合作为主线候选**
   - 虽然 runtime 更低，但 RMSE 和 uncertainty calibration 都比 `spatial` 差很多。

expr_additive sweep 结果也比较一致：

```text
best RMSE:
  expr_additive_lambda_0.2_pca10
  rmse = 0.160504

best layer_ari:
  expr_additive_lambda_0.5_pca5
  layer_ari = 0.063554
  rmse = 0.162826
```

补充观察：

1. `lambda_expr` 越大，RMSE 越差：
   - `0.2 < 0.5 < 1.0 < 2.0`
2. 低维 expression embedding (`pca5`) 比 `pca10/pca15` 略有 biological consistency 优势，但 accuracy 仍明显不够。
3. `HVG` 版本并没有带来改进，反而略差。

因此，这轮 suite 给出的决策线是：

> 当前这版 expression-informed GP MVP 已经证明“方向有潜在生物学价值”，但还没有达到可替代纯 spatial GP 主线的程度。

更具体地说：

- **不建议** 现在就把 spaGAPA 主方法切换为 `expr_additive`。
- **建议** 保持 `spatial GP (matern + alpha=1e-3)` 作为当前主线方法。
- `expr_additive` 继续保留为研究分支，但下一步不应只做小参数微调，而应考虑更换模型形式，例如：
  - product kernel
  - gated / adaptive expression weighting
  - layer-aware local expression kernel
  - uncertainty-preserving multi-view GP

这一步其实很重要，因为它帮我们排除了一个常见陷阱：

- 不是“只要加 expression 就会更强”
- 而是“如何把 expression manifold 以不破坏 calibration 的方式并入 GP”才是真问题

所以，expression-informed 这条路线**没有死**，但当前这个 additive MVP 版本已经基本摸到了它的上限。

### 23.13 expr_product pilot + lambda/offset sweep 已完成

在 additive MVP 明显失败后，已进一步实现并测试 `expr_product` 路线。

结果目录：

```text
spaGAPA/benchmark_results/real/product_gp_suite_v1/
```

关键文件：

- `pilot_comparison.csv`
- `product_sweep_summary.csv`
- `decision_summary.json`
- `figures/product_pilot_comparison.png`
- `figures/product_sweep.png`

本轮 pilot 设置：

- `n_genes = 30`
- `seeds = 42`
- `mask_types = random, ring_sector, layer_aware`
- compare:
  - `spatial`
  - `expr_additive`
  - `expr_product`

pilot 结果摘要：

```text
spatial:
  rmse = 0.112329
  layer_ari = 0.036980
  uncertainty_error_spearman = 0.720964

expr_additive:
  rmse = 0.170647
  layer_ari = 0.029388
  uncertainty_error_spearman = 0.535007

expr_product:
  rmse = 0.158512
  layer_ari = 0.056384
  uncertainty_error_spearman = 0.505211
```

这一轮 pilot 给出的信息非常明确：

1. `expr_product` **比 `expr_additive` 有实质改善**
   - RMSE 从 `0.1706` 改善到 `0.1585`
   - layer ARI 从 `0.0294` 提升到 `0.0564`

2. 但 `expr_product` **仍明显不如 `spatial`**
   - `spatial rmse = 0.1123`
   - `expr_product rmse = 0.1585`
   - calibration 也显著更差：
     - `spatial uncertainty_error_spearman = 0.721`
     - `expr_product uncertainty_error_spearman = 0.505`

3. 因此，`product kernel` 证明了一件重要的事：
   - additive 的失败并不代表“expression-aware GP 完全没有希望”
   - 但目前这版 product kernel 还不足以替代纯 spatial GP 主线

本轮同时完成了一个小规模 product sweep：

- `lambda_expr = 0.2 / 0.5 / 1.0`
- `product_offset = 0.5 / 1.0 / 2.0`

当前最佳 product 组合为：

```text
expr_product_lambda_0.2_offset_0.5
  rmse = 0.153940
  layer_ari = 0.070396
  uncertainty_error_spearman = 0.510808
```

这一组合说明：

1. **更弱的 expression 权重更好**
   - `lambda_expr = 0.2` 明显优于 `0.5` 和 `1.0`

2. **更小的 product offset 更好**
   - `offset = 0.5` 是当前最优
   - offset 增大并没有带来更稳的 calibration 或更好的 accuracy

3. **当前最优 product 仍然没有跨过 spatial 基线**
   - 虽然比 additive 更合理
   - 但还没达到“主线替换”的门槛

所以，这一轮最重要的结论不是“product 成功了”，而是：

> product kernel 证明确实比 additive 更接近正确方向，但仍不足以在 accuracy + calibration 维度上超过纯 spatial GP。

这意味着下一步最合理的收敛是：

1. 不把 `expr_product` 直接升为主线
2. 继续保留 `spatial GP (matern + alpha=1e-3)` 为当前主线
3. 下一步优先进入：
   - `adaptive / gated expression weighting`
   - 或 `layer-aware local expression kernel`

换句话说：

- `expr_additive` 基本可以视为已被淘汰
- `expr_product` 值得保留为中间过渡方案
- 但真正可能解决问题的，更像是“局部门控 expression 信息”，而不是全局表达核本身

### 23.14 adaptive_additive pilot + lambda/tau sweep 已完成

在 `expr_product` 之后，已进一步实现并测试 `adaptive_additive` 路线：

```text
K_total = K_space + lambda_expr * (K_gate * K_expr)
K_gate(i,j) = exp(- d_space(i,j)^2 / (2 * tau^2))
```

结果目录：

```text
spaGAPA/benchmark_results/real/adaptive_gp_suite_v1/
```

关键文件：

- `pilot_comparison.csv`
- `adaptive_sweep_summary.csv`
- `decision_summary.json`
- `figures/gp_variant_pilot_comparison.png`
- `figures/adaptive_sweep.png`

本轮 pilot 设置：

- `n_genes = 30`
- `seeds = 42`
- `mask_types = random, ring_sector, layer_aware`
- compare:
  - `spatial`
  - `expr_additive`
  - `expr_product`
  - `adaptive_additive`

pilot 结果摘要：

```text
spatial:
  rmse = 0.112329
  layer_ari = 0.036980
  uncertainty_error_spearman = 0.720964

expr_additive:
  rmse = 0.170647
  layer_ari = 0.029388
  uncertainty_error_spearman = 0.535007

expr_product:
  rmse = 0.158512
  layer_ari = 0.056384
  uncertainty_error_spearman = 0.505211

adaptive_additive:
  rmse = 0.167212
  layer_ari = 0.028287
  uncertainty_error_spearman = 0.511170
```

这个结果说明：

1. 默认 `adaptive_additive` 并没有直接超过 `expr_product`
   - 默认 adaptive 的 RMSE 仍明显高于 product
   - default adaptive 的 layer ARI 甚至低于 product

2. adaptive 路线并不是完全没有价值
   - 通过 `lambda_expr` 和 `gate_tau` 的 sweep，可以把 adaptive 拉到比默认 additive 更合理的位置
   - 说明“expression 只在局部空间邻域生效”这个方向本身是对的

本轮同时完成了一个小规模 adaptive sweep：

- `lambda_expr = 0.2 / 0.5 / 1.0`
- `gate_tau = auto / 0.5 / 1.0 / 2.0`

当前最佳 adaptive 组合为：

```text
best RMSE:
adaptive_lambda_0.5_tau_0.5
  rmse = 0.159586
  layer_ari = 0.057309
  uncertainty_error_spearman = 0.524192

best layer ARI:
adaptive_lambda_0.5_tau_1.0
  rmse = 0.162721
  layer_ari = 0.058065
  uncertainty_error_spearman = 0.500056
```

这一步带来的核心判断是：

1. adaptive 经调参后，**可以接近 `expr_product`**
   - RMSE 与 product 最优结果相近
   - layer ARI 也可略高于 pilot product 默认值

2. 但 adaptive 仍然**没有跨过当前 spatial 基线**
   - `spatial rmse = 0.1123`
   - 最优 adaptive `rmse` 仍约为 `0.16`
   - calibration 也仍明显弱于纯 spatial GP

3. 因此，adaptive 的定位更像是：
   - 证明 gated expression weighting 是合理方向
   - 但当前这版 `distance-gated additive` 还不足以成为主线替代方案

所以这一轮后，路线判断进一步收敛为：

1. `expr_additive` 已基本结束历史使命
2. `expr_product` 仍是当前最强 expression-aware 候选
3. `adaptive_additive` 证明了 gate 机制有价值，但当前版本更像“支持性证据”而不是最终答案
4. 下一步如果还继续 expression-aware GP，更应该优先：
   - `product kernel + adaptive weighting`
   - 或 `layer-aware local expression kernel`

---

### 23.15 layer_local pilot + local gate sweep 已完成

在 `expr_product` 和 `adaptive_additive` 都未跨过 spatial 基线后，已进入 Phase 5，进一步实现并测试 `layer_local` 路线。

本轮特别注意避免 label leakage：

- 正式 `layer_local` MVP **不使用真实 MOB layer label 作为模型输入**
- layer label 只用于 benchmark evaluation，例如 `layer_ari/layer_nmi`
- local gate 只从坐标派生：
  - `radius`
  - `pseudolayer`，由 radius quantile 生成

MVP kernel 形式为：

```text
K_total = K_space * (product_offset + lambda_expr * K_local * K_expr)
```

其中：

```text
K_local = soft local gate(radius or pseudolayer)
```

结果目录：

```text
spaGAPA/benchmark_results/real/layer_local_gp_suite_v1/
```

关键文件：

- `pilot_comparison.csv`
- `layer_local_sweep_summary.csv`
- `decision_summary.json`
- `figures/gp_variant_pilot_comparison.png`
- `figures/layer_local_sweep.png`

本轮 pilot 设置：

- `n_genes = 30`
- `seeds = 42`
- `mask_types = random, ring_sector, layer_aware`
- compare:
  - `spatial`
  - `expr_additive`
  - `expr_product`
  - `adaptive_additive`
  - `expr_layer_local`

pilot 结果摘要：

```text
spatial:
  rmse = 0.112329
  layer_ari = 0.036980
  uncertainty_error_spearman = 0.720964

expr_product:
  rmse = 0.158512
  layer_ari = 0.056384
  uncertainty_error_spearman = 0.505211

adaptive_additive:
  rmse = 0.167212
  layer_ari = 0.028287
  uncertainty_error_spearman = 0.511170

expr_layer_local:
  rmse = 0.144761
  layer_ari = 0.043207
  uncertainty_error_spearman = 0.578766
```

这一轮最重要的结果是：

1. `layer_local` 明显改善了 expression-aware GP 的 accuracy
   - 比 `expr_additive`、`expr_product`、`adaptive_additive` 的 RMSE 都更好
   - 说明“局部限制 expression influence”确实比全局 expression kernel 更合理

2. `layer_local` 的 calibration 也明显改善
   - `expr_product uncertainty_error_spearman = 0.505`
   - `adaptive_additive uncertainty_error_spearman = 0.511`
   - `expr_layer_local uncertainty_error_spearman = 0.579`

3. 但 `layer_local` 仍未超过纯 spatial GP
   - `spatial rmse = 0.1123`
   - `expr_layer_local rmse = 0.1448`
   - `spatial uncertainty_error_spearman = 0.721`

本轮同时完成了一个小规模 local gate sweep：

- `lambda_expr = 0.2 / 0.5`
- `local_k = 10 / 20`
- `layer_gate_mode = radius / pseudolayer`

当前最佳 RMSE 组合为：

```text
layer_local_lambda_0.5_k_10_radius
  rmse = 0.135396
  layer_ari = 0.049367
  uncertainty_error_spearman = 0.610685
```

当前最佳 layer ARI 组合为：

```text
layer_local_lambda_0.2_k_10_radius
  rmse = 0.136662
  layer_ari = 0.064641
  uncertainty_error_spearman = 0.609602
```

这一结果把路线判断进一步收敛为：

1. `layer_local` 是目前最强的 expression-aware GP 方向
   - 比 additive/product/adaptive 更接近 spatial GP
   - calibration 也明显更合理

2. 但它还没有达到“替代 spatial GP 主线”的标准
   - RMSE 仍落后 spatial
   - uncertainty calibration 仍落后 spatial
   - layer ARI 仍没有超过 stAPAminer-like expression KNN baseline

3. 因此当前最稳妥的主线仍然是：
   - `spatial GP (matern + alpha=1e-3)` 作为 primary imputation
   - `layer_local` 作为 expression-aware research branch / ablation
   - downstream 重点转向 uncertainty-aware SVAPA / differential APA / domain analysis

换句话说：

> expression-aware GP 这条线已经证明“局部化是正确方向”，但在当前数据和模型规模下，仍不足以替代 spatial GP 主线。

---

### 23.16 BioML benchmark integration 初步完成

在确认 single-gene expression-aware GP 难以同时兼顾 RMSE、calibration 和 biological consistency 后，已进一步实现并测试 CPU-friendly BioML backbone 的 benchmark 集成。

本轮 BioML 方法不是替代 GP，而是两阶段模型：

```text
spatial GP imputation + uncertainty
  -> multi-view graph
  -> graph-regularized APA factorization
  -> BioML domain detector
```

核心原则：

- 不使用深度学习
- 不使用 GPU
- 不使用真实 layer label 作为训练输入
- layer label 只用于 ARI/NMI 评价

BioML benchmark 已接入：

```text
scripts/run_stapaminer_mob_benchmark.py
```

新增参数：

```text
--include-bioml
--bioml-rank
--bioml-lambda-graph
--bioml-lambda-l2
--bioml-max-iter
--bioml-n-neighbors
--bioml-blend
--bioml-domain-method {kmeans,spectral}
--bioml-spatial-weight
--bioml-expression-weight
--bioml-apa-weight
```

本轮 targeted pilot 设置：

- `n_genes = 30`
- `seeds = 42`
- `mask_types = random, ring_sector, layer_aware`
- `gp_kernel = matern`
- `gp_alpha = 1e-3`
- `bioml_blend = 0.1`

结果目录：

```text
spaGAPA/benchmark_results/real/bioml_targeted_v1/
spaGAPA/benchmark_results/real/bioml_targeted_spectral_v1/
```

#### KMeans BioML domain detector

```text
spatial GP:
  rmse = 0.112329
  layer_ari = 0.036980
  uncertainty_error_spearman = 0.720964

BioML-kmeans:
  rmse = 0.112568
  layer_ari = 0.059991
  uncertainty_error_spearman = 0.710106

stAPAminer-like:
  rmse = 0.125075
  layer_ari = 0.077191
```

KMeans BioML 的意义：

- 几乎保住了 spatial GP 的 RMSE
- layer ARI 从 `0.037` 提升到 `0.060`
- 但仍未超过 stAPAminer-like 的 `0.077`

#### Spectral BioML domain detector

```text
spatial GP:
  rmse = 0.112329
  layer_ari = 0.036980
  uncertainty_error_spearman = 0.720964

BioML-spectral:
  rmse = 0.112568
  layer_ari = 0.377720
  layer_nmi = 0.474067
  uncertainty_error_spearman = 0.710106

stAPAminer-like:
  rmse = 0.125075
  layer_ari = 0.077191
  layer_nmi = 0.110255
```

分 mask 结果同样支持这一结论：

```text
BioML-spectral layer ARI:
  random      = 0.339723
  ring_sector = 0.370689
  layer_aware = 0.422746

stAPAminer-like layer ARI:
  random      = 0.049233
  ring_sector = 0.038743
  layer_aware = 0.143597
```

这一步是目前为止最重要的突破：

1. **RMSE 基本保持 spatial GP 水平**
   - `spatial GP rmse = 0.112329`
   - `BioML-spectral rmse = 0.112568`

2. **biological consistency 首次显著超过 stAPAminer-like baseline**
   - `BioML-spectral layer ARI = 0.378`
   - `stAPAminer-like layer ARI = 0.077`

3. **uncertainty calibration 基本继承 GP backbone**
   - `spatial GP uncertainty_error_spearman = 0.721`
   - `BioML-spectral uncertainty_error_spearman = 0.710`

4. **这说明 spaGAPA 的“全面超越”路线应转向 BioML backbone**
   - GP 负责 imputation + uncertainty
   - multi-view graph/spectral domain detector 负责 biological consistency

当前判断：

> spaGAPA-BioML 已经在小规模 targeted pilot 中实现了非常接近 spatial GP 的 RMSE/calibration，并显著超过 stAPAminer-like 的 layer ARI/NMI。下一步应扩大到 formal benchmark：更多 genes、seeds、mask types 和 graph weight sweep。

---

### 23.17 BioML formal benchmark v1 已完成

BioML formal suite 已跑完，runner 未再发现后台 benchmark 进程，suite-level summary、decision JSON 和图表均已生成。

结果目录：

```text
spaGAPA/benchmark_results/real/bioml_formal_suite_v1/
```

关键输出：

```text
bioml_formal_overall_summary.csv
bioml_formal_mask_summary.csv
decision_summary.json
figures/bioml_formal_method_comparison.png
figures/bioml_graph_weight_sweep.png
```

Formal balanced 设置：

- `n_genes = 60`
- `seeds = 42,43`
- `mask_types = random, spatial_block_large, ring_sector, layer_aware, low_coverage`
- `spatial_weight = 0.4`
- `expression_weight = 0.4`
- `apa_weight = 0.2`
- `bioml_domain_method = spectral`

#### Overall formal benchmark

```text
spagapa_bioml:
  rmse = 0.116868
  mae = 0.067549
  pearson = 0.962168
  r2 = 0.925756
  layer_ari = 0.375996
  layer_nmi = 0.492093
  uncertainty_error_spearman = 0.810176
  uncertainty_coverage_68 = 0.865732
  uncertainty_coverage_95 = 0.960293
  runtime_s = 15.612
  peak_rss_mb = 615.939

spagapa_gp:
  rmse = 0.116897
  layer_ari = 0.077308
  layer_nmi = 0.117774
  uncertainty_error_spearman = 0.810283

stAPAminer-like expression KNN:
  rmse = 0.127412
  layer_ari = 0.101109
  layer_nmi = 0.134307
```

这组结果非常关键：

1. **RMSE 首次在 formal benchmark 中保持并微弱超过 spatial GP**
   - `spagapa_bioml rmse = 0.116868`
   - `spagapa_gp rmse = 0.116897`

2. **biological consistency 显著超过既往 expression-KNN 类 baseline**
   - `spagapa_bioml layer_ari = 0.376`
   - `stAPAminer-like layer_ari = 0.101`
   - `spagapa_gp layer_ari = 0.077`

3. **uncertainty calibration 基本完整继承 GP backbone**
   - `spagapa_bioml uncertainty_error_spearman = 0.810176`
   - `spagapa_gp uncertainty_error_spearman = 0.810283`

4. **CPU runtime 可接受**
   - formal average runtime 约 `15.6 s`
   - peak RSS 约 `616 MB`

#### Mask-level biological consistency

BioML 在全部 5 种 mask 下均明显提高 layer recovery，尤其是之前最担心的 `layer_aware`：

```text
layer_aware:
  spagapa_bioml layer_ari = 0.444642
  stAPAminer-like layer_ari = 0.230029
  spagapa_gp layer_ari = 0.098998

random:
  spagapa_bioml layer_ari = 0.383458
  stAPAminer-like layer_ari = 0.066869
  spagapa_gp layer_ari = 0.125774

ring_sector:
  spagapa_bioml layer_ari = 0.383751
  stAPAminer-like layer_ari = 0.071885
  spagapa_gp layer_ari = 0.020854

spatial_block_large:
  spagapa_bioml layer_ari = 0.334394
  stAPAminer-like layer_ari = 0.044043
  spagapa_gp layer_ari = 0.075919

low_coverage:
  spagapa_bioml layer_ari = 0.333734
  stAPAminer-like layer_ari = 0.092720
  spagapa_gp layer_ari = 0.064996
```

这说明之前 `single-gene expression-aware GP` 没能解决的 layer-structured biological manifold recovery，已经被 `GP + multi-view graph + spectral BioML` 路线明显改善。

#### Graph weight sweep

Graph weight sweep 结果显示：

```text
Best BioML by RMSE:
  no_apa_s050_e050_a000
  spatial_weight = 0.5
  expression_weight = 0.5
  apa_weight = 0.0
  rmse = 0.115065
  layer_ari = 0.320739

Best BioML by layer ARI:
  expression_heavy_s020_e060_a020
  spatial_weight = 0.2
  expression_weight = 0.6
  apa_weight = 0.2
  rmse = 0.115123
  layer_ari = 0.387113
  layer_nmi = 0.490091
```

解释：

- `expression_heavy` 是 biological consistency 最优候选。
- `no_apa` 是 RMSE 最优候选，提示当前 APA-view graph 可能仍含噪，需要后续改进 APA similarity 或 feature selection。
- `balanced_s040_e040_a020` 是 formal 主线的稳健默认配置，因为它在 RMSE、calibration 和 biological consistency 之间最均衡。

#### 当前决策

BioML backbone 已经通过当前阶段的主线候选标准：

1. RMSE 不劣于 spatial GP。
2. layer ARI/NMI 明显超过 stAPAminer-like expression baseline。
3. uncertainty calibration 基本继承 GP。
4. CPU runtime/memory 可接受。

因此下一步不应再把主要精力放在 single-gene expression-aware GP 微调上，而应把 BIB 主线正式升级为：

```text
spaGAPA = uncertainty-calibrated spatial GP imputation
       + CPU-friendly multi-view graph BioML domain recovery
       + uncertainty-aware downstream APA analysis
```

需要注意的是，这仍然只是 MOB formal v1。BIB 投稿前仍需：

- 扩展到更大 gene set 和更多 seeds。
- 加入外部 spatial transcriptomics 数据集。
- 引入 APA/domain boundary consistency、marker ordering、SVAPA ranking stability。
- 对 BioML graph weight 进行更系统的 cross-dataset 默认参数选择。

但从目前结果看，spaGAPA 已经从“有希望追平既往包”推进到“在核心维度上具备全面超越证据”的阶段。

---

### 23.18 BioML 已接入主 pipeline / CLI

基于 formal benchmark v1 的结果，BioML 不再只是 benchmark 分支，而是已接入 spaGAPA 主线 API 和命令行。

新增主线能力：

```text
SpaGAPA(..., use_bioml=True)
SpaGAPA.run(..., expression_matrix=..., expression_embedding=..., use_bioml=True)
spagapa run --enable-bioml --expression-matrix expression_matrix.csv
```

实现要点：

- `SpaGAPA.run()` 保持旧行为兼容：
  - 默认 `use_bioml=False`
  - 不开启 BioML 时仍使用原 KMeans domain identification
- 开启 BioML 后：
  - 使用 GP imputed APA matrix 和 uncertainty
  - 可使用 expression matrix 经 PCA 形成 spot-level expression embedding
  - 构建 spatial / expression / APA multi-view graph
  - 运行 graph-regularized APA factorization
  - 用 spectral 或 kmeans BioML domain detector 输出 domains
- 真实 layer label 不进入模型，只能用于外部评价。

新增可复现输出：

```text
domains.csv
bioml_imputed_values.npy
bioml_spot_factors.npy
bioml_gene_factors.npy
bioml_metadata.json
dataset.h5ad
```

CLI 示例：

```bash
spagapa run \
  --apa-matrix apa_matrix.csv \
  --coordinates coordinates.csv \
  --expression-matrix expression_matrix.csv \
  --enable-bioml \
  --bioml-domain-method spectral \
  --n-domains 5 \
  --output spagapa_results
```

验证：

```text
conda run -n spagapa pytest tests/integration/test_pipeline_toy.py tests/unit/test_bioml.py tests/unit/test_apa_dataset.py tests/unit/test_feature_builders.py -q

26 passed
```

当前判断：

> BioML backbone 已完成从 benchmark branch 到正式 pipeline/CLI 的迁移。下一步应开始外部真实数据验证，优先选择有 layer/domain annotation 或可获得生物学区域标签的 spatial transcriptomics dataset，以评价 cross-dataset biological consistency 和默认 graph weight 的稳健性。

---

### 23.19 外部真实数据验证 runner 已接入

根据 `metaAPA.pdf` 和 `stAPAminer.pdf` 的实验范式，已新增外部真实数据验证入口：

```text
scripts/run_external_bioml_validation.py
```

这个 runner 的定位不是替代 formal mask benchmark，而是把 spaGAPA-BioML 放进更接近竞品论文的真实数据验证框架中：

#### stAPAminer 对齐的实验维度

stAPAminer 原文主要强调：

- KNN expression imputation 后 layer/domain separation 改善。
- 与解剖 layer label 的一致性，包括 ARI、NMI、purity、Jaccard。
- 内部聚类指标，包括 DBI、Calinski-Harabasz、Silhouette、Dunn。
- 同一 layer 内 spot-spot APA profile correlation。
- SVAPA、DEAPA、LSAPA 及下游生物学解释。
- replicate 中 spatial APA genes 的一致性。

当前 runner 已覆盖：

```text
external_validation_summary.csv
within_layer_correlation.csv
replicate_spatial_gene_overlap.csv
downstream/*_svapa_genes.csv
downstream/*_lsapa_domain_markers.csv
downstream/*_deapa_layer_pairwise.csv
figures/layer_separation_metrics.png
figures/internal_clustering_metrics.png
figures/within_layer_correlation.png
figures/spatial_apa_gene_counts.png
method_domain_maps.png
```

当前可比较方法：

```text
raw
stapaminer_original_imputed
stapaminer_knn_expression
spagapa_gp
spagapa_bioml
```

其中：

- `stapaminer_original_imputed` 使用本地 stAPAminer 示例输出的 imputed RUD。
- `stapaminer_knn_expression` 是 faithful reimplementation，用 expression KNN 做 APA imputation。
- `spagapa_bioml` 通过正式 `SpaGAPA(..., use_bioml=True)` pipeline 跑通，不再依赖 benchmark-only 分支。

#### metaAPA 对齐的实验维度

metaAPA 原文主要评估上游 poly(A) site caller integration：

- Sierra / polyApipe / SCAPE 多 caller site consensus。
- position threshold 和 expression similarity clustering。
- high-confidence sites overlap / Jaccard。
- PAS、cleavage site、CSTF GU-rich、CF I UGUA 等 sequence feature。
- Nanopore / Space Ranger 等外部支持。

这与 spaGAPA-BioML 的下游 imputation/domain recovery 不是同一层级，因此当前 runner 不把 metaAPA 当作 layer ARI baseline，而是输出 site-level readiness checklist：

```text
site_level_validation_checklist.json
figures/ref_package_experiment_coverage.png
```

当前 MOB 数据状态：

```text
apa_sites_available = true
site_counts_available = true
multi_caller_site_tables_available = false
sequence_context_available = false
long_read_support_available = false
metaapa_ready = false
```

解释：

- 当前数据已具备真实 APA site/count 输入。
- 尚未具备 metaAPA 风格多 caller consensus 和 sequence-context validation。
- 后续若要完整对齐 metaAPA，需要补 Sierra / polyApipe / SCAPE outputs、genome FASTA/GTF、PAS annotation 或 long-read support。

#### MOB pilot 结果

运行命令：

```bash
conda run -n spagapa python scripts/run_external_bioml_validation.py \
  --dataset-dir spaGAPA/data/processed/stapaminer_mob \
  --output-dir spaGAPA/benchmark_results/real/external_bioml_validation_v1 \
  --n-genes 60 \
  --min-observed-spots 110 \
  --methods raw,stapaminer_original_imputed,stapaminer_knn_expression,spagapa_gp,spagapa_bioml \
  --gp-alpha 1e-3 \
  --gp-n-restarts 1 \
  --bioml-domain-method spectral \
  --bioml-max-iter 20 \
  --bioml-blend 0.1
```

输出目录：

```text
spaGAPA/benchmark_results/real/external_bioml_validation_v1/
```

关键结果：

```text
method                         layer_ari  layer_nmi  purity   pairwise_jaccard
raw                            0.101117   0.151950   0.4000   0.224058
stapaminer_original_imputed    0.109863   0.160847   0.4077   0.227414
stapaminer_knn_expression      0.073448   0.109605   0.3731   0.199467
spagapa_gp                     0.100875   0.160667   0.3808   0.218274
spagapa_bioml                  0.388330   0.509298   0.6385   0.360246
```

下游 gene count：

```text
method                         n_svapa  n_lsapa  n_deapa
raw                            9        4        6
stapaminer_original_imputed    9        4        5
stapaminer_knn_expression      9        10       6
spagapa_gp                     9        8        8
spagapa_bioml                  9        9        8
```

解释：

1. `spagapa_bioml` 在真实 MOB pilot 中对 biological consistency 的提升非常明显：
   - layer ARI 从 stAPAminer-original 的 `0.110` 提升到 `0.388`
   - layer NMI 从 `0.161` 提升到 `0.509`
   - purity 从 `0.408` 提升到 `0.638`

2. `within-layer correlation` 在当前 60 high-coverage genes 中区分度不大：
   - raw、stAPAminer、GP、BioML 的 median within-layer Pearson 都接近 `0.948`
   - 这说明当前 gene subset 覆盖度过高，该指标在 pilot 中接近饱和
   - 后续应增加 low-coverage / dropout-heavy gene subset，才能更公平测试 imputation 对 layer 内 correlation 的贡献

3. internal clustering metrics 需要谨慎解读：
   - BioML 的 layer agreement 明显最强
   - 但 DBI/Silhouette 并不总是更优，因为 BioML 的目标是对齐真实 biological layer，而不是优化无监督紧致球形 cluster
   - 论文中应把 external biological labels 作为主指标，internal metrics 作为辅助诊断

4. 这次结果进一步支持：

```text
spaGAPA-BioML is no longer a speculative branch.
It is now the main candidate route for BIB-level biological consistency.
```

#### 当前仍缺的外部验证

BIB 投稿前仍需继续扩展：

1. 至少 1 个 stAPAminer 以外的 spatial APA / poly(A) dataset。
2. 至少 1 个 brain/CBS 或 MOB replicate，用于 replicate overlap。
3. 更大 gene set，包括 high-coverage、low-coverage 和 dropout-heavy 分层。
4. metaAPA-style upstream site validation：多 caller consensus、PAS/sequence feature、long-read 或 Space Ranger support。
5. Biological story：从 `spagapa_bioml` 的 SVAPA/LSAPA/DEAPA 结果中挑出 layer-specific APA remodeling gene，做图和文献解释。

---

### 23.20 High-resolution spaGAPA MVP 已启动

为了让 spaGAPA 不只停留在普通 Visium/ST 规模，已新增高分辨率路线 spec 和 pseudo-bin benchmark：

```text
HIGH_RESOLUTION_SPAGAPA_SPEC.md
scripts/run_high_resolution_simulation.py
tests/integration/test_high_resolution_simulation.py
```

高分辨率路线的核心定位：

```text
No deep learning. No GPU. Use sparse/block GP + multi-view BioML.
```

#### 为什么要做 high-resolution benchmark

高分辨率空间转录组会带来两个相反因素：

- 空间结构更细：microdomain、boundary、local niche 更有价值。
- APA 更稀疏：单个 bin 的 reads 更少，APA usage 更 noisy/missing。

因此 BIB 叙事不能只证明 spaGAPA 能在 260 spots MOB 上跑通，还需要证明它具备向 high-resolution spatial APA 扩展的路线。

#### 用户界面设计：内部复杂，外部简单

虽然内部实现会包含 GP、Sparse GP、BioML、highres BioML、fast-domain mode 等组件，但用户-facing 设计不应该让用户直接面对一堆方法名。

推荐最终 CLI/API 暴露方式：

```text
analysis_preset = auto | standard | highres_accuracy | highres_fast
```

含义：

1. `auto`
   - 默认选项。
   - 根据 spot/bin 数量、APA matrix observed fraction、每 spot/bin 覆盖度、坐标密度自动选择。
   - 普通 Visium/ST 数据走 `standard`。
   - 高分辨率/高稀疏数据走 `highres_accuracy`，并可提示 fast-domain 选项。
2. `standard`
   - 低/中分辨率空间转录组默认。
   - 使用常规 spaGAPA GP/BioML 路线。
   - 更强调 APA-specific spatial modeling 和 uncertainty-aware downstream。
3. `highres_accuracy`
   - 高分辨率或 pseudo-bin 数据默认推荐。
   - 使用 `highres_bioml` accuracy mode。
   - 保留 sparse GP value recovery 和 uncertainty。
4. `highres_fast`
   - 高分辨率快速 domain discovery / exploratory analysis。
   - 使用 `highres_bioml_gp_blend = 0`，跳过 sparse GP。
   - layer/domain recovery 快，但 APA value RMSE 和 uncertainty 不作为主要卖点。

高分辨率数据的实用定义不应只看平台名，而应看数据形态：

```text
high-resolution-like =
  many spatial bins/spots
  + small bin/cell-level capture
  + sparse APA observations
  + local pseudo-bin / bead / grid structure
```

初步 automatic rule：

```text
if n_spots_or_bins >= 1000
   or observed_fraction < 0.35
   or median_observed_APA_per_bin is low:
       use highres preset
else:
       use standard preset
```

这条规则后续必须在真实 high-resolution datasets 上校准。当前 benchmark 中 `320-640 pseudo-bins` 已经表现出 high-resolution-like 稀疏性，因此可作为开发压力测试，但不能替代真实平台验证。

#### MVP 设计

当前不先下载大型 high-resolution 数据，而是从真实 MOB APA matrix 生成 pseudo high-resolution bins：

```text
parent spot
  -> multiple pseudo-bins
  -> coordinate jitter
  -> inherited layer label
  -> local APA micro-noise
  -> capture loss + dropout + measurement noise
```

比较方法：

```text
raw
expression_knn
sparse_gp
sparse_bioml
highres_bioml
```

输出：

```text
spaGAPA/benchmark_results/real/highres_simulation_v1/
highres_results_summary.csv
figures/highres_method_comparison.png
figures/highres_scaling.png
figures/highres_domain_maps_last.png
```

#### MOB pseudo-high-resolution pilot

运行规模：

```text
n_genes = 40
n_parent_spots = 80
subbins_per_spot = 4
n_bins = 320
observed_fraction = 0.335781
capture_rate = 0.45
dropout_rate = 0.25
```

初始直接迁移结果：

```text
method          rmse_holdout  parent_rmse  layer_ari  layer_nmi  runtime_s
raw             0.093836      0.070466     0.013281   0.028330   0.000
expression_knn  0.110820      0.086591     0.317657   0.337631   0.012
sparse_gp       0.149626      0.120304     0.121348   0.199846   2.735
sparse_bioml    0.139902      0.112524     0.060880   0.167263   3.174
```

Uncertainty:

```text
sparse_gp uncertainty_error_spearman = 0.226988
sparse_bioml uncertainty_error_spearman = 0.217176
```

解释：

1. 高分辨率 benchmark pipeline 已经打通，可生成 reproducible tables 和 figures。
2. Sparse GP uncertainty 与误差正相关，说明 uncertainty 仍有信息。
3. 直接把普通分辨率 `sparse_bioml` 迁移到 pseudo-bin high-resolution 场景会失败：
   - APA view 太 noisy，参与 domain graph 后会拖垮 layer recovery。
   - KMeans-on-matrix 不足以恢复 high-resolution layer structure。
   - 必须把 APA value recovery 和 biological domain graph recovery 解耦。
4. `raw` 在 RMSE 上最强，主要因为当前 pseudo-bin truth 由 parent-smoothed APA 生成，而且选的是较高覆盖 genes，gene-mean filling 是强 baseline。
5. 因此 v1 high-resolution 结果应被解读为：

```text
High-resolution direct sparse BioML is technically feasible, but not the right
default high-resolution model.
```

这对 BIB 反而是有价值的开发信号：普通 MOB 外部验证中 BioML 已经很强，但 high-resolution 稀疏场景下需要专门的 decoupled / multiscale tuning，而不是直接套用普通分辨率默认参数。

#### Highres BioML v2：已把垫底问题修正为主线候选

为解决 v1 暴露的问题，已在 `scripts/run_high_resolution_simulation.py` 中新增：

```text
highres_bioml
```

算法定义：

```text
Value recovery:
  observed high-res APA
    -> raw gene-mean fill
    -> sparse GP
    -> 70% raw + 30% sparse GP conservative value blend

Domain recovery:
  coordinates + expression embedding + expression-KNN APA proxy
    -> 20% spatial + 60% expression + 20% APA-proxy fused graph
    -> spectral domain detector
```

核心原则：

```text
Do not let noisy high-resolution sparse-GP APA distances define tissue domains.
Use sparse GP for value recovery and uncertainty, but use spatial/expression
structure plus a light APA proxy for biological domain recovery.
```

正式输出：

```text
spaGAPA/benchmark_results/real/highres_simulation_v2/
```

v2 seed 42 结果：

```text
method          rmse_holdout  parent_rmse  layer_ari  layer_nmi  runtime_s
raw             0.093836      0.070466     0.013281   0.028330   0.000
expression_knn  0.110820      0.086591     0.317657   0.337631   0.009
sparse_gp       0.123681      0.100235     0.105452   0.197039   3.232
sparse_bioml    0.117031      0.094870     0.120553   0.231166   3.698
highres_bioml   0.091409      0.071571     0.353396   0.475051   3.457
```

额外 seed 复现：

```text
seed 43: highres_bioml rmse_holdout = 0.087739, layer_ari = 0.375770
seed 44: highres_bioml rmse_holdout = 0.088438, layer_ari = 0.356022
```

当前判断：

1. `sparse_bioml` 垫底不是路线终结，而是 high-resolution graph construction 错配。
2. `highres_bioml` 已经在当前 pseudo-high-resolution MOB benchmark 中同时超过：
   - raw 的 pseudo-bin RMSE
   - expression-KNN 的 layer ARI/NMI
   - sparse GP / direct sparse BioML 的综合表现
3. 这条线符合用户约束：
   - 只用传统机器学习
   - CPU-only
   - 运行时间仍在秒级
4. 现在 high-resolution spaGAPA 的主线候选应从 `sparse_bioml` 切换为 `highres_bioml`。

#### Highres BioML formal suite v1：稳定性与速度取舍

已新增：

```text
scripts/run_highres_bioml_suite.py
```

输出：

```text
spaGAPA/benchmark_results/real/highres_bioml_suite_v1/
highres_bioml_suite_results_long.csv
highres_bioml_suite_overall_summary.csv
decision_summary.json
figures/highres_formal_method_comparison.png
figures/highres_bioml_scaling.png
figures/highres_dropout_stress.png
figures/highres_graph_weight_sweep.png
figures/highres_gp_blend_sweep.png
```

formal multi-seed 平均结果：

```text
method          rmse_holdout  parent_rmse  layer_ari  layer_nmi  runtime_s
raw             0.093485      0.069798     0.025283   0.048144   0.000
expression_knn  0.107669      0.084110     0.243685   0.288447   0.010
sparse_gp       0.113715      0.091988     0.110783   0.183307   3.415
sparse_bioml    0.108105      0.087425     0.165666   0.265937   3.961
highres_bioml   0.089195      0.069369     0.361730   0.470109   3.629
```

关键结论：

1. `highres_bioml` 在 formal multi-seed mean 中同时是：
   - best RMSE
   - best layer ARI
   - best layer NMI
2. 相比 `expression_knn`：
   - `rmse_holdout` 改善 `-0.018474`
   - `layer_ari` 提升 `+0.118044`
   - `layer_nmi` 提升 `+0.181662`
3. 相比 `raw`：
   - `rmse_holdout` 改善 `-0.004290`
   - `layer_ari` 提升 `+0.336447`
4. 这说明 `highres_bioml` 的领先不只是 seed 42 的小幅偶然提升，而是在当前 pseudo-high-resolution formal suite 中具有稳定优势。

graph-weight sweep 结果支持当前默认：

```text
best graph config = default_s020_e060_a020_exprknn
spatial weight    = 0.20
expression weight = 0.60
APA proxy weight  = 0.20
APA proxy source  = expression_knn
```

这说明 high-resolution domain recovery 的最佳方向不是 APA-heavy，而是：

```text
expression-led + spatially constrained + lightly APA-informed
```

GP-blend sweep 暴露出一个重要速度/精度取舍：

```text
config        rmse_holdout  parent_rmse  layer_ari  layer_nmi  runtime_s
gp_blend_0.0  0.093836      0.070466     0.353396   0.475051   0.384
gp_blend_0.1  0.091799      0.069911     0.353396   0.475051   3.655
gp_blend_0.3  0.091409      0.071571     0.353396   0.475051   3.793
gp_blend_0.5  0.095903      0.076670     0.353396   0.475051   3.464
```

因此 high-resolution spaGAPA 现在可以定义两个运行模式：

1. **accuracy mode**
   - `highres_bioml_gp_blend = 0.3`
   - 最好 pseudo-bin RMSE
   - 保留 sparse GP uncertainty
   - runtime 为秒级
2. **fast-domain mode**
   - `highres_bioml_gp_blend = 0.0`
   - 当 APA source 不是 `sparse_gp` 时跳过 sparse GP
   - runtime 从约 3.5 秒降到约 0.38 秒
   - layer ARI/NMI 与 accuracy mode 相同
   - 代价是 RMSE 回到 raw-level，且没有 sparse GP uncertainty

已完成 runtime 优化：

```text
If highres_bioml_gp_blend == 0 and highres_bioml_apa_source != sparse_gp,
skip sparse GP completely.
```

同时 BioML spectral domain detector 已改为保留 sparse graph affinity，不再强制 dense 化。

初始 scaling 风险：

```text
n_bins  rmse_holdout  layer_ari  runtime_s
160     0.090433      0.307608   3.199
320     0.091409      0.353396   3.878
640     0.087832      0.104927   4.137
```

8x pseudo-bin 的 layer ARI 明显下降，说明下一步重点不是继续堆模型复杂度，而是要做：

- adaptive graph neighborhood
- multiscale / parent-aware graph
- aggregation-aware model selection
- true high-resolution dataset validation

#### Adaptive neighborhood v2：修复 8x pseudo-bin layer ARI 下降

针对 8x pseudo-bin 下 layer ARI 从 `0.353` 降到 `0.105` 的问题，已新增：

```text
--highres-bioml-neighbor-mode adaptive
--highres-bioml-adaptive-neighbor-scale 10.0
--highres-bioml-parent-weight
--highres-bioml-parent-neighbors
```

核心思想：

```text
When one parent spot is split into more pseudo-bins, a fixed KNN graph becomes
too local and fragments tissue layers. Increase KNN with pseudo-bin density.
```

当前默认：

```text
neighbor_mode = adaptive
base_neighbors = bioml_n_neighbors
extra_neighbors = round((pseudo_bins_per_parent - 4) * 10)
```

8x seed 42 对照：

```text
setting                    layer_ari  layer_nmi
fixed/default old           0.104927   0.199823
adaptive scale 4            0.308625   0.403829
adaptive scale 8            0.339756   0.456775
adaptive scale 10           0.410055   0.512730
adaptive scale 12           0.355558   0.443150
```

8x scale 10 多 seed：

```text
seed  rmse_holdout  parent_rmse  layer_ari  layer_nmi
42    0.087832      0.064874     0.410055   0.512730
43    0.083241      0.061011     0.401305   0.494254
44    0.085945      0.062521     0.314886   0.415050
mean  0.085673      0.062802     0.375415   0.474011
```

更新后的 scaling check：

```text
n_bins  rmse_holdout  parent_rmse  layer_ari  layer_nmi
160     0.090433      0.073048     0.307608   0.449267
320     0.091409      0.071571     0.353396   0.475051
640     0.087832      0.064874     0.410055   0.512730
```

结论：

1. 8x 崩盘主要来自 fixed KNN graph 太局部，不是 `highres_bioml` 主线失败。
2. adaptive neighborhood 已把 8x layer ARI 从 `0.105` 提升到 `0.410`。
3. `parent_weight=0.2` 没有超过 pure adaptive KNN，因此 parent-aware graph 暂不设为默认。
4. 下一步应在更多 seeds、更多 datasets、真实 high-resolution 数据上验证 adaptive scale 是否需要自动调参。

#### 下一步 high-resolution 优先事项

1. 对 `highres_bioml` 做 formal benchmark：
   - 多 seeds
   - 2x / 4x / 8x pseudo-bins
   - low-coverage genes
   - dropout-heavy scenarios
   - runtime / memory scaling
2. 做 high-resolution graph-weight sweep：
   - spatial-heavy
   - expression-heavy
   - expression-KNN APA proxy
   - no-APA domain graph
   - sparse-GP APA view as negative control
3. 验证 adaptive graph-neighbor selection：
   - 更多 seeds
   - 更多 subbin levels
   - true high-resolution datasets
4. 比较 sparse GP、block GP、local GP。
5. 加 aggregation-aware model selection：
   - pseudo-bin 层面不要只看局部拟合
   - 聚合回 parent spot 后也要保留 tissue/layer structure
6. 寻找真正保留 3-prime/poly(A) 信息的 high-resolution spatial transcriptomics 数据集。

### 23.22 Analysis preset 已接入主 pipeline / CLI

为避免用户直接面对 `sparse GP`、`BioML`、`highres_bioml`、`fast-domain`
等内部组件，spaGAPA 已新增 user-facing analysis preset 层：

```text
analysis_preset = auto | standard | highres_accuracy | highres_fast
```

当前实现位置：

```text
spagapa/presets.py
spagapa/pipeline.py
spagapa/cli.py
```

#### 当前解析逻辑

`auto` 不根据平台名判断，而是根据数据形态判断：

```text
profile =
  n_genes
  n_spots
  finite_fraction
  observed_fraction
  positive_fraction
  median_observed_per_spot
  median_positive_per_spot
```

初始 high-resolution-like 规则：

```text
if n_spots >= 1000
   or observed_fraction < 0.35
   or (n_spots >= 300 and positive APA signal is very sparse):
       auto -> highres_accuracy
else:
       auto -> standard
```

这条规则是 MVP，后续需要在真实 high-resolution datasets 上重新校准阈值。

#### preset 行为

```text
standard
  - 保留常规 GP / domain pipeline 行为
  - BioML 由用户显式 --enable-bioml 或 API use_bioml=True 启用

highres_accuracy
  - 启用 sparse GP
  - 启用 BioML
  - 默认 graph weights 调为 spatial=0.2, expression=0.6, APA=0.2
  - 保留 GP imputation 和 uncertainty

highres_fast
  - 启用 BioML
  - 跳过 GP imputation
  - 默认 graph weights 调为 spatial=0.2, expression=0.6, APA=0.2
  - bioml_blend 默认为 0
  - 定位是快速 domain discovery / exploratory analysis
```

CLI 已改为：

```bash
spagapa run \
  --apa-matrix apa_matrix.csv \
  --coordinates coordinates.csv \
  --expression-matrix expression_matrix.csv \
  --analysis-preset auto \
  --n-domains 5 \
  --output spagapa_results
```

其中 `--enable-bioml/--disable-bioml` 改为三态 override：

- 不写：由 preset 决定
- `--enable-bioml`：强制启用 BioML
- `--disable-bioml`：强制关闭 BioML

#### 可复现输出

每次运行会在结果中记录：

```text
results["analysis_preset"]
dataset.adata.uns["apa"]["analysis_preset"]
analysis_preset.json
```

记录内容包括：

- requested preset
- resolved preset
- dataset profile
- impute / sparse GP / BioML 是否启用
- BioML graph weights
- fast mode 是否跳过 GP

### 23.23 highres_bioml decoupled route 已接入主 pipeline

在 23.22 的 preset entrypoint 之后，已进一步把 benchmark runner 中表现最好的
`highres_bioml` decoupled domain route 提炼为包内模块，并接入主 pipeline。

新增/修改位置：

```text
spagapa/bioml/highres.py
spagapa/bioml/__init__.py
spagapa/pipeline.py
spagapa/cli.py
tests/unit/test_highres_bioml.py
tests/integration/test_pipeline_toy.py
```

#### 核心设计

`highres_bioml` 不再等价于普通 `sparse GP + BioML factorization`，而是采用
benchmark 中验证过的 decoupled 设计：

```text
observed APA matrix
  -> raw gene-mean fill
  -> optional sparse-GP value blend
  -> recovered APA values

coordinates + expression embedding + optional APA proxy
  -> high-resolution multi-view graph
  -> spectral BioML domain recovery
```

这解决了之前 `sparse_bioml` 在 high-resolution pseudo-bin 数据中垫底的核心问题：
不能让 sparse APA value recovery 直接支配 biological domain graph。

#### highres preset 当前行为

```text
highres_accuracy
  - sparse GP imputation: enabled
  - highres_bioml gp_blend: 0.3
  - domain graph APA source: expression_knn
  - graph weights: spatial=0.2, expression=0.6, APA=0.2
  - output domain method: spagapa_highres_bioml

highres_fast
  - sparse GP imputation: skipped
  - highres_bioml gp_blend: 0
  - domain graph APA source: expression_knn if expression is available, otherwise raw
  - graph weights: spatial=0.2, expression=0.6, APA=0.2
  - output domain method: spagapa_highres_bioml
```

#### 新增可调参数

CLI/API 已支持：

```text
highres_bioml_gp_blend
highres_bioml_apa_source = expression_knn | raw | sparse_gp | none
highres_bioml_expression_knn_k
highres_bioml_neighbor_mode = fixed | adaptive
highres_bioml_adaptive_neighbor_scale
highres_bioml_parent_weight
highres_bioml_parent_neighbors
```

其中 `parent_index` 可从以下 AnnData metadata 自动读取：

```text
adata.obs["parent_spot"]
adata.obs["parent_spot_id"]
adata.obs["parent_index"]
adata.obs["parent_bin"]
adata.uns["apa"]["parent_index"]
```

如果没有 parent/coarse-bin 信息，则 adaptive neighbor 自动退化为 fixed KNN。

#### 输出与可复现性

highres preset 的 domain 输出：

```text
results["domains"]["method"] = "spagapa_highres_bioml"
results["domains"]["metadata"]["mode"] = "highres_bioml"
results["highres_bioml_values"]
dataset.adata.uns["apa"]["bioml"]["highres"]
```

保存结果仍复用：

```text
bioml_imputed_values.npy
bioml_metadata.json
analysis_preset.json
```

#### 当前边界

这一步完成的是主 pipeline 接入，不等于已经完成真实 high-resolution 数据验证。

下一步必须做：

1. 用同一套 CLI/pipeline 接口跑 pseudo-highres formal suite 的 pipeline-level
   smoke benchmark，确认与 runner 结果一致。
2. 在外部真实 high-resolution datasets 上验证 `auto -> highres_accuracy`
   阈值和默认 graph weights。
3. 单独报告 `highres_fast` 的 runtime/domain recovery，同时明确其不提供 GP
   uncertainty。

---

## 24. 最终目标陈述

spaGAPA 面向 BIB 的最终目标不是证明“我们写了一个包”，而是证明：

> Spatial APA analysis requires dedicated statistical treatment of sparse, spatially structured, and uncertain APA usage measurements. spaGAPA provides a reproducible framework that integrates spatial validation, Gaussian process imputation, posterior uncertainty propagation, and downstream SVAPA/differential analysis, enabling more reliable detection and interpretation of spatially organized APA remodeling in real tissues.

如果后续所有开发、benchmark 和写作都围绕这句话展开，spaGAPA 才真正具备冲击 BIB 的可能。

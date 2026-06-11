# spaGAPA 当前问题清单与投稿策略评估

更新日期：2026-06-10

本文档用于记录 spaGAPA 在真实数据 benchmark 和论文投稿前需要解决的问题，并从生物信息学专家/顶级期刊审稿人的视角评估其创新性、竞品关系和潜在投稿层级。

本轮只写文档，不修改源码。

---

## 1. 一句话结论

spaGAPA 已经具备一个有潜力的空间转录组 APA 方法学框架：以 scAPAtrap 等 caller 产生的 APA 位点/计数为入口，用空间邻域验证降低噪音，用 Gaussian Process 进行空间插补和不确定性量化，再把不确定性带入 SVAPA、差异 APA、空间模式和 benchmark。

从审稿人视角看，真正决定论文档次的不是代码量，而是以下四件事：

1. 真实 APA 数据而非表达 proxy。
2. 与 stAPAminer 和相关 baseline 的严格、可复现、同数据对比。
3. 证明 GP 与不确定性确实带来下游生物学收益，而不仅是 RMSE 变好。
4. 打通 pipeline 的工程细节，使方法可以被外部用户稳定复现。

若补齐真实数据 benchmark、消融实验、真实组织层级生物学验证，并修复 pipeline 集成问题，spaGAPA 有机会达到 Bioinformatics / PLOS Computational Biology / Genome Biology 边缘到中上水平。Nature Methods 级别需要更强的跨技术、跨组织、跨任务泛化证明，以及明确的生物学发现或广泛方法学影响。

---

## 2. 当前项目定位

### 2.1 软件目标

spaGAPA 的目标是成为空间转录组 APA 分析工具，覆盖：

- APA 位点识别或候选位点接入。
- 空间感知 APA 位点验证。
- APA index 插补。
- 不确定性量化。
- APA 定量指标计算。
- 空间域识别。
- 差异 APA 分析。
- SVAPA 识别。
- 空间模式聚类。
- benchmark 与可视化。

### 2.2 核心方法链条

当前最有价值的主线可以写成：

```text
scAPAtrap / APA caller 输出候选 poly(A) 位点和计数
  -> 空间邻域验证，过滤随机噪音候选位点
  -> APA index 矩阵构建
  -> Gaussian Process 空间插补
  -> posterior uncertainty 量化
  -> 不确定性加权的 SVAPA / 差异 APA / 空间域分析
  -> 真实组织层级和生物学模式解释
```

这个链条比单独说“我们用了 GP”更有说服力，因为它把空间验证、插补、不确定性、SVAPA 检测和生物学解释串成一个完整方法学闭环。

---

## 3. 当前必须解决的问题

下面按照 P0、P1、P2 分优先级。P0 是论文/benchmark 前必须解决，P1 是提高可信度和可维护性的关键工程项，P2 是增强项或 v1.1 方向。

---

## 4. P0 问题：真实 benchmark 与主流程可靠性

### P0.1 真实 APA 数据尚未真正进入 benchmark

**现状**

- 当前 scripts 已能下载 10x Visium 表达矩阵和坐标。
- `prepare_benchmark_data.py` 目前将表达矩阵归一化后模拟为 APA index proxy。
- 这对调试流程可以接受，但不能作为论文中的“真实 APA benchmark”。

**审稿风险**

审稿人会直接问：

- 你比较的是 APA 方法，为什么真实 benchmark 使用 gene expression proxy？
- APA index 的 missingness 是否和表达 dropout 相同？
- distal/proximal poly(A) site 使用变化是否能由基因表达归一化 proxy 代表？
- 如果没有真实 APA calling，如何证明工具真的适用于空间转录组 APA？

**需要解决**

1. 下载或准备真实空间转录组 BAM 文件，而不仅是 filtered feature-barcode matrix。
2. 用 scAPAtrap、Sierra、SCAPE 或 stAPAminer 文档中的 PA extraction 流程获得真实 poly(A) site count matrix。
3. 明确每个 gene 的 proximal/distal site 定义。
4. 构建真实 APA index 矩阵，如 RUD、PDUI 或 WUL。
5. 输出标准化数据：
   - `apa_matrix.csv`: spots x genes 或 genes x spots 必须明确。
   - `coordinates.csv`: spot_id, x, y。
   - `apa_sites.bed` 或 `apa_sites.csv`: poly(A) site 注释。
   - `gene_site_mapping.csv`: gene 到 proximal/distal sites 的映射。
   - `metadata.csv`: 组织层、样本、replicate、技术平台。

**验收标准**

- 至少 2 个真实数据集完成真实 APA calling。
- 至少 1 个数据集与 stAPAminer 原论文或示例数据同源，优先 MOB。
- benchmark 中不再使用 expression proxy 作为主结果。
- expression proxy 仅可放在 supplement 作为 pipeline smoke test。

---

### P0.2 需要复现或严格对齐 stAPAminer

**竞品实际实现**

stAPAminer 主要流程：

1. 用 scAPAtrap/Sierra 等先得到 poly(A) site matrix。
2. 用 movAPA 计算 RUD/WUL/SLR/GPI 等 APA index。
3. 用 gene expression matrix 计算 spot 间欧氏距离。
4. KNN 插补 APA index，默认 k=10。
5. 用 Seurat 做空间聚类/层级标注。
6. 用 Seurat `FindMarkers` / limma / edgeR 做 DEAPA 或 LSAPA。
7. 用 SPARK 做 SVAPA。
8. 用 kmeans 和相关系数做空间 pattern 聚类。

**spaGAPA 的比较优势**

- stAPAminer 的 KNN 距离来自 gene expression space，不是显式空间坐标。
- stAPAminer 没有 posterior uncertainty。
- stAPAminer 的 SVAPA 依赖 SPARK，未把插补不确定性纳入检验。
- stAPAminer 缺少系统模拟框架和消融实验。

**审稿风险**

如果只比较“spaGAPA-GP vs 我们自己写的 KNN-spatial baseline”，审稿人会认为没有真正比较 stAPAminer。

**需要解决**

1. 运行 stAPAminer 原 R 代码，至少在 MOB 数据上复现其 RUD imputation 和 SVAPA。
2. 基线至少包含：
   - stAPAminer 原始 KNN，gene expression distance，k=10。
   - KNN-spatial，使用空间坐标。
   - mean / median imputation。
   - GP-no-spatial，使用表达 PCA 或随机坐标作为消融。
   - spaGAPA-GP，空间坐标 GP。
   - spaGAPA-GP + uncertainty-weighted SVAPA。
   - SPARK-based SVAPA，与 stAPAminer 同类。
3. 统一输入矩阵、gene set、spot set、filtering 规则。
4. 不只汇报 RMSE，还汇报下游结果：
   - domain identification ARI/NMI。
   - DEAPA detection F1/precision/recall。
   - SVAPA reproducibility。
   - spatial pattern stability。

**验收标准**

- 可以用一条脚本复现 `spaGAPA vs stAPAminer` 主要 benchmark 表格。
- 每个 baseline 记录运行时间、内存、失败率。
- 所有比较使用相同 mask、相同 spots、相同 genes。

---

### P0.3 pipeline 集成层有潜在错误

源码中核心算法模块较完整，但 `SpaGAPA.run()` 这类一键 pipeline 有几个需要先修的集成问题。

#### P0.3.1 数据矩阵轴向可能错位

**现象**

- `APADataset.raw_counts` 返回 shape 为 `(n_genes, n_spots)`。
- pipeline 中又执行 `self.dataset_.raw_counts.T`，得到 `(n_spots, n_genes)`。
- 但 `GPImputerBatch` 期望输入是 `(n_genes, n_spots)`。

**风险**

真实数据一跑，可能出现：

- boolean mask 与 coordinates 长度不一致。
- 把 spots 当 genes 拟合。
- imputed matrix 维度错位。
- downstream 结果不可解释。

**需要解决**

- 统一项目内矩阵约定。
- 建议约定：
  - AnnData 内部：spots x genes。
  - 算法 API：genes x spots。
  - benchmark API：spots x genes，若保留则必须显式转换。
- 给每个 public function 加 shape check。

**验收标准**

- 一个 3 genes x 5 spots 的 toy dataset 能完整跑过 pipeline。
- 每一步输出 shape 被测试：
  - raw counts: genes x spots。
  - imputed: genes x spots。
  - uncertainty: genes x spots。
  - domain labels: spots。
  - SVAPA results: genes。

#### P0.3.2 SparseGPImputer 参数不匹配

**现象**

- pipeline 里调用 `SparseGPImputer(kernel_type=..., n_inducing=...)`。
- `SparseGPImputer.__init__` 目前没有 `kernel_type` 参数，只支持 `n_inducing`, `inducing_method`, `length_scale`, `noise_level`。

**风险**

- 一旦 `use_sparse_gp=True`，pipeline 会直接报错。

**需要解决**

- 二选一：
  - 给 SparseGPImputer 增加 `kernel_type` 参数并实现 Matérn/RBF。
  - 或 pipeline 中只传 SparseGPImputer 已支持的参数。

**验收标准**

- `SpaGAPA(use_sparse_gp=True).run(...)` 能在 toy data 上通过。

#### P0.3.3 pipeline 的 APA quantification 步骤目前是占位

**现象**

- pipeline 创建了 `APAIndexCalculator()`。
- 但没有调用 `calculate_all()`。
- `results_['apa_indices']` 只是保存 `{'values': work, 'calculated': True}`。

**风险**

审稿人或用户以为 pipeline 已计算 RUD/PDUI/PAI，但实际上只是把工作矩阵放进去。

**需要解决**

- 明确 pipeline 的输入：
  - 如果输入是 poly(A) site counts，则需要先构建 proximal/distal counts。
  - 如果输入已经是 APA index，则 quantification 不应声称重新计算。
- 将 API 分成：
  - `quantify_from_site_counts()`
  - `use_existing_apa_index()`

**验收标准**

- pipeline 输出真实的 RUD/PDUI/PAI/WUL 或明确标记为 existing APA index。
- 文档中说明每个指标需要的输入。

#### P0.3.4 pipeline 的 differential analysis 目前不是完整统计检验

**现象**

- pipeline 中若开启 differential analysis，只输出每个 domain 的 spot 数。
- 没有调用 `DifferentialAPAAnalyzer.test_differential_apa()` 或 `test_all_pairwise()`。

**风险**

- README 或论文若写“一键差异 APA 分析”，会被复现实验发现不一致。

**需要解决**

- 在 pipeline 中真正调用 pairwise domain differential analysis。
- 输出每个 pair 的 pvalue、padj、log2fc、effect size。
- 支持 uncertainty weighting。

**验收标准**

- `results_['differential']` 是包含 gene-level test results 的 DataFrame 或 dict。
- 至少一个 toy dataset 能检测出人工设置的 domain-specific APA gene。

#### P0.3.5 GP likelihood ratio 中 uncertainty weight 未真正进入 GP fit

**现象**

- `GPTrendDetector.likelihood_ratio_test()` 计算了 inverse uncertainty weights。
- 但 sklearn `GaussianProcessRegressor.fit()` 没有传入这些 weights。
- 因此 LR test 当前不是严格 uncertainty-weighted。

**风险**

如果论文声称“GP likelihood ratio test uses uncertainty weighting”，审稿人可能认为实现不一致。

**需要解决**

可选路径：

1. 改成 heteroscedastic GP，用 per-spot uncertainty 作为 alpha。
2. 用 weighted log likelihood 或 bootstrap。
3. 在论文中把 LR test 与 uncertainty-weighted Moran's I 分开表述：
   - LR test: GP-based spatial trend detection。
   - weighted Moran's I: uncertainty-aware SVAPA screening。

**验收标准**

- 代码和论文描述一致。
- 消融实验展示 uncertainty weighting 的贡献。

---

### P0.4 read_spatial_data 与 pipeline 的输入类型不一致

**现象**

- `read_spatial_data()` 返回 `(coords_df, data_df)` tuple。
- `SpaGAPA.run()` 预期 `self.dataset_` 是 `APADataset`，随后访问 `.n_genes`, `.n_spots`, `.coords`。

**风险**

用 `bam_file` 或 `h5ad_file` 走 pipeline 时可能失败。

**需要解决**

- `read_spatial_data()` 应返回 `APADataset`。
- 或 pipeline 中对 tuple 做转换。
- 明确 BAM 入口是否会调用 scAPAtrap wrapper。

**验收标准**

- `SpaGAPA.run(dataset=APADataset)` 通过。
- `SpaGAPA.run(h5ad_file=...)` 或 CLI dataset 输入通过。
- 若支持 BAM，则必须能调用 scAPAtrap wrapper 并生成 APADataset。

---

### P0.5 需要真实生物学验证，而不仅是算法分数

**审稿人会问**

- 你多发现的 SVAPA genes 是否有生物学意义？
- 是否与 MOB 层级结构、脑区组织结构或发育轴一致？
- 是否能复现 stAPAminer 论文中的空间 APA 模式？
- GP 产生的新模式是否只是平滑过度？

**需要解决**

1. MOB:
   - 复现 ONL、GL、EPL、MCL、GCL 等组织层结构。
   - 比较 stAPAminer 论文报告的层级 APA genes。
   - 展示若干 gene 的空间 APA map。
2. Brain cortex:
   - 验证 cortical layer 或 spatial domain 相关 APA pattern。
3. Embryo:
   - 验证发育轴或组织区特异 APA pattern。
4. GO/KEGG:
   - 对 spaGAPA-only SVAPA genes 做功能富集。
   - 对 consensus genes 与 stAPAminer-only genes 分开分析。
5. 文献支持：
   - 对 top genes 做 APA 或神经/发育相关文献检索。

**验收标准**

- 至少 3 个强 case studies。
- 每个 case 包括：
  - spatial APA map。
  - uncertainty map。
  - stAPAminer 对照。
  - 生物学解释。
  - 可复现脚本。

---

### P0.6 benchmark 设计需要升级为“论文级”

当前 benchmark 框架已经有模拟器和 evaluator，但论文级 benchmark 需要更系统。

**必须有的 benchmark 层次**

1. Simulation benchmark:
   - ground truth APA matrix。
   - domain-specific pattern。
   - gradient pattern。
   - random/no spatial pattern。
   - spatial autocorrelated pattern。
   - dropout rate: 20%, 40%, 60%, 80%。
   - noise level: low/medium/high。
   - spot scale: 500, 2000, 5000。
2. Real-data masking benchmark:
   - 随机 mask observed APA values。
   - 按 gene/spot/dropout strata 分层 mask。
   - 比较恢复能力。
3. Downstream benchmark:
   - domain ARI/NMI。
   - DEAPA F1。
   - SVAPA reproducibility。
   - pattern clustering stability。
4. Ablation study:
   - GP without spatial coordinates。
   - GP without uncertainty。
   - spaGAPA without spatial validation。
   - weighted Moran's I vs unweighted Moran's I。
   - LR test vs Moran's I vs SPARK。
5. Runtime and scalability:
   - full GP vs sparse GP。
   - n_spots scaling。
   - n_genes scaling。
   - memory usage。

**推荐主图**

- Figure 1: spaGAPA workflow。
- Figure 2: Simulation benchmark, imputation accuracy。
- Figure 3: Dropout/noise sensitivity and runtime。
- Figure 4: Real MOB benchmark vs stAPAminer。
- Figure 5: SVAPA detection, uncertainty-weighted improvement。
- Figure 6: Biological case studies。
- Supplementary Figure: Ablation, parameter sensitivity, reproducibility。

---

## 5. P1 问题：代码质量、可复现性和用户信任

### P1.1 测试数量和文档记录不一致

**现状**

- README、PROGRESS、PROJECT_SUMMARY、tasks.md 中测试数量存在 267、268、288、304 等不同说法。

**风险**

- 审稿人或用户会觉得项目状态不清晰。

**需要解决**

- 运行完整 pytest。
- 记录真实测试数量、通过率、覆盖率。
- 更新所有文档。

**验收标准**

- README、PROGRESS、PROJECT_SUMMARY 的测试数一致。
- CI badge 与真实状态一致。

---

### P1.2 缺少端到端集成测试

**现状**

- 单元测试较多。
- integration tests 目录相对空。

**需要解决**

1. Toy dataset end-to-end:
   - APADataset 创建。
   - GP imputation。
   - domain identification。
   - differential APA。
   - SVAPA detection。
   - save results。
2. Realistic simulated dataset end-to-end:
   - 100 genes x 300 spots。
   - dropout 40%。
   - benchmark expected GP > KNN。
3. CLI smoke test:
   - `spagapa impute`。
   - `spagapa diff`。
   - `spagapa run`。

**验收标准**

- 每次修改 pipeline 后，集成测试能捕捉 shape 和 API 错误。

---

### P1.3 CLI 存在但文档仍说 deferred

**现状**

- `spagapa/cli.py` 已存在。
- `pyproject.toml` 中有 `spagapa = "spagapa.cli:main"`。
- 但部分文档仍称 CLI deferred。

**需要解决**

- 确认 CLI 是否实际可运行。
- 更新文档状态。
- 增加 CLI examples。

**验收标准**

- `spagapa --help`、`spagapa impute --help`、`spagapa diff --help` 正常。
- README 和 user guide 有最小可运行 CLI 示例。

---

### P1.4 API 命名和矩阵方向需要统一

**问题**

当前代码中存在两种矩阵方向：

- AnnData/scikit benchmark 常用 spots x genes。
- APA 算法模块多用 genes x spots。

这本身没错，但必须明确转换边界。

**需要解决**

- 在每个函数 docstring 写明 shape。
- 对 public API 输入做 `assert_shape` 或 `validate_matrix_orientation`。
- 对 benchmark 和 pipeline 增加统一 adapter。

**验收标准**

- 用户不会被 `.T` 搞迷糊。
- 错误输入能得到清晰报错，而不是 silent wrong result。

---

### P1.5 代码格式、静态检查和打包

**需要解决**

- 运行 black/isort。
- 运行 flake8。
- 运行 mypy 或至少 pyright/ruff basic checks。
- 确保 `pip install -e ".[dev]"` 成功。
- 确保 examples 不依赖本地相对路径。

**验收标准**

- GitHub Actions 或本地 Makefile 能一键运行：
  - tests。
  - lint。
  - format check。
  - minimal benchmark smoke test。

---

### P1.6 文档需要转成“用户能复现”的结构

**当前文档内容足够多，但需要更像正式软件**

建议文档结构：

1. Installation。
2. Quick start with simulated data。
3. Real data input format。
4. Running scAPAtrap to generate APA matrix。
5. Running spaGAPA imputation。
6. SVAPA detection。
7. Benchmark reproduction。
8. Output explanation。
9. FAQ。
10. Troubleshooting。

**验收标准**

- 一个外部用户只按 user guide 可以跑通 toy example。
- 一个审稿人可以按 benchmark README 复现主图。

---

## 6. P2 问题：未来增强方向

### P2.1 多 caller 整合

metaAPA 的启发：

- position-based integration。
- similarity-based integration。
- 多 caller 输出统一。

spaGAPA 未来可以支持：

- scAPAtrap。
- Sierra。
- SCAPE。
- polyApipe。
- metaAPA-style consensus PAS。

这会降低“依赖单一 caller”的审稿风险。

### P2.2 更先进的 GP 后端

当前 sklearn GP 适合中小规模数据。若目标是 10,000+ spots，建议后续考虑：

- GPyTorch。
- sparse variational GP。
- inducing point learning。
- GPU acceleration。
- per-spot heteroscedastic noise。

### P2.3 组织图像整合

空间转录组常有 H&E 图像。未来可以：

- 用 histology features 改进 kernel。
- 构建 spatial + morphology GP kernel。
- 展示 APA patterns 与组织结构叠加。

### P2.4 trajectory 和多尺度分析正式化

trajectory 模块已实现，但文档中曾标注 deferred。若要作为论文卖点，需要：

- 明确任务。
- 增加真实生物学 case。
- 与现有 trajectory/spatial pseudotime 工具比较。

---

## 7. 竞品分析

### 7.1 scAPAtrap

**定位**

scAPAtrap 是 APA site caller，不是完整空间 APA 分析工具。它从 3' tag-based scRNA-seq/ST BAM 中识别和定量 poly(A) sites。

**核心方法**

- 唯一比对 read 过滤。
- UMI 去重。
- 正负链分离。
- coverage peak calling。
- soft-clipping A/T stretch 识别 poly(A) tail。
- peak count matrix 生成。
- internal priming 过滤和注释。

**优势**

- APA site calling 相对成熟。
- 适合做 spaGAPA upstream。

**限制**

- 不是专为空间下游分析设计。
- 不解决 APA index 插补。
- 不提供空间不确定性建模。

**spaGAPA 与其关系**

不是正面竞争，而是 upstream 依赖/互补。spaGAPA 应该诚实表述为：

> spaGAPA can use scAPAtrap as an upstream PAS detection engine and adds spatial validation, GP imputation, uncertainty-aware downstream analysis, and benchmarking.

### 7.2 stAPAminer

**定位**

stAPAminer 是主要直接竞品，目标是空间转录组 APA usage mining。

**核心方法**

- APA site 输入来自 scAPAtrap 或 Sierra。
- movAPA 计算 RUD/WUL/SLR/GPI 等。
- KNN imputation，默认 k=10。
- KNN 距离来自 gene expression matrix，而不是显式空间坐标。
- Seurat 聚类和 FindMarkers。
- SPARK 检测 SVAPA。
- kmeans 聚类空间 pattern。

**优势**

- 已在 MOB 数据中展示了空间 APA pattern。
- R 生态和 Seurat/SPARK 结合。
- 有真实数据案例。

**限制**

- KNN 插补启发式强，缺少概率模型。
- 无每个插补值 uncertainty。
- imputation error 不进入下游显著性检验。
- SVAPA 依赖 SPARK，和 APA-specific uncertainty 没有结合。
- benchmark 框架较弱。
- Python/scanpy/squidpy 生态不友好。

**spaGAPA 相对创新**

- 用空间坐标直接建模 APA spatial covariance。
- GP posterior std 给出每个 spot/gene 的可信度。
- uncertainty-weighted Moran's I 和 GP-based SVAPA。
- 可做 sparse/block GP scalability。
- Python 生态和 AnnData 数据结构。
- simulation + real benchmark 框架。

### 7.3 metaAPA

**定位**

metaAPA 是多 polyA site caller 输出整合管线，不是空间 APA usage 分析工具。

**核心方法**

- Nextflow pipeline。
- Sierra、polyApipe、SCAPE 等多 caller。
- position-based site clustering。
- similarity-based site clustering。
- 多 distance metrics 和 clustering methods。

**优势**

- 多 caller consensus。
- 对 PAS 集成很有启发。

**限制**

- 不聚焦空间下游 APA。
- 不解决空间插补、不确定性和 SVAPA。

**spaGAPA 可借鉴**

- v1.1 可以引入 multi-caller consensus PAS，减少 scAPAtrap 单一依赖。

---

## 8. 顶级期刊审稿人视角：创新点与质疑点

### 8.1 最强创新点

#### 创新点 1：把空间坐标纳入 APA 插补模型，而不是只做 KNN 平均

stAPAminer 的 KNN imputation 本质上是启发式局部平均。spaGAPA 的 GP 将 APA usage 视为空间连续随机场，用 kernel 描述空间协方差。

审稿人认可度：高。

但前提是必须证明：

- 在空间结构存在时 GP 优于 KNN。
- 在无空间结构时 GP 不产生过强假阳性。
- 参数选择和 kernel choice 稳定。

#### 创新点 2：posterior uncertainty 进入 downstream analysis

这是 spaGAPA 最像“方法学创新”的部分。很多空间工具会做平滑/插补，但不告诉用户哪里不可信。spaGAPA 如果能系统利用 uncertainty，就有较强说服力。

需要强调：

- uncertainty map 可提示低置信区域。
- high uncertainty spots 在 SVAPA 检测中权重降低。
- uncertainty calibration 与真实误差相关。
- uncertainty-aware SVAPA 降低 FDR。

#### 创新点 3：GP-based SVAPA detection

用 GP marginal likelihood 做 spatial trend vs null model 的比较，比单纯 Moran's I 或 SPARK 更贴近 spaGAPA 自身插补模型。

强点：

- 统计框架统一。
- 可以解释空间方差比例。
- 可结合 uncertainty。

需要补强：

- 当前 LR test 中 uncertainty 还没有真正进入 GP fit。
- 需要理论或模拟证明 type I error 控制。
- 需要与 SPARK/Moran's I 比较。

#### 创新点 4：从 APA caller 到 downstream 的完整 Python workflow

空间转录组主流 Python 生态是 AnnData/scanpy/squidpy。spaGAPA 的 APADataset 基于 AnnData，可以成为实际用户选择工具的重要理由。

审稿人认可度：中高。

前提：

- 安装容易。
- examples 可运行。
- CLI/API 文档一致。
- 数据格式清晰。

#### 创新点 5：benchmark framework

如果 benchmark 做得完整，spaGAPA 会比只提供算法函数的工具更有竞争力。

审稿人会喜欢：

- simulation 有 ground truth。
- real data 有复现。
- downstream task benchmark。
- ablation study。

---

### 8.2 审稿人最可能质疑的问题

1. scAPAtrap 已经完成 APA site calling，spaGAPA 的 calling 创新是否只是后处理？
2. GP 是标准方法，应用到 APA 是否足够新？
3. GP 平滑是否会制造不存在的 spatial APA pattern？
4. 不确定性估计是否校准？是否真的反映误差？
5. SVAPA 检测的 type I error 是否受控？
6. 与 stAPAminer 是否公平比较？
7. 是否用了真实 APA 数据，还是表达 proxy？
8. 是否有跨数据集、跨组织、跨平台泛化？
9. 软件是否可复现？
10. 是否有真实生物学发现，还是只有模拟 benchmark？

---

## 9. 投稿期刊层级分析

期刊信息依据官方页面与当前公开信息。投稿前应再次检查最新影响因子、栏目要求和格式。

### 9.1 Nature Methods

**适配度：低到中，属于冲顶选择。**

Nature Methods 明确强调 novel methods、显著改进、强验证、重要 biological application、与现有方法的性能比较，以及 immediate practical relevance。其 scope 也包括 spatial omics 和 computational/statistical methods for biological data analysis。

spaGAPA 想冲 Nature Methods，需要满足：

- 不只是“空间 APA 小众工具”，而是展示一种可推广的 uncertainty-aware spatial transcriptomics modeling framework。
- 至少 3 到 5 个真实数据集，跨技术或跨组织。
- 与 stAPAminer、SPARK、KNN、non-spatial GP 等全面比较。
- 有一个真正重要的 biological question，例如发育/脑区/肿瘤微环境中 APA spatial regulation 的新发现。
- 工具鲁棒、可复现、外部可用。

**当前状态不建议直接投。**

**完善后仍然偏挑战。**

除非真实数据结果非常漂亮，并且从 APA 推广到更广泛 spatial omics uncertainty-aware inference，否则 Nature Methods 概率较低。

### 9.2 Genome Biology

**适配度：中到高，取决于生物学发现。**

Genome Biology 覆盖 genomic/post-genomic biology，包括 new methods and software tools。若 spaGAPA 不仅是软件，还能在 MOB/Brain/Embryo 中发现空间 APA 调控规律，则比较适配。

需要达到：

- 方法学强。
- 真实数据强。
- 有新 biological insight。
- 不是只发工具说明书。

**推荐条件**

- 至少 2 到 3 个真实空间转录组 APA datasets。
- 复现 stAPAminer 的发现，并发现额外高置信 SVAPA。
- 新发现有组织结构、功能富集和文献支持。

### 9.3 PLOS Computational Biology

**适配度：中高。**

PLOS Computational Biology 对方法/软件要求是有广泛采用潜力、能带来新的生物学 insight、严谨可复现。它不要求实验验证必须有，但真实数据应用和 code/data availability 非常重要。

spaGAPA 适合投 PLOS Computational Biology 的条件：

- GP + uncertainty-aware SVAPA 的方法论讲清楚。
- simulation + real benchmark 严谨。
- 生物学案例足够清晰。
- 开源和复现做得好。

**这是一个比较合理的冲刺目标。**

### 9.4 Bioinformatics

**适配度：高。**

Bioinformatics 关注 genome bioinformatics 和 computational biology 的新发展，也有 Application Notes 栏目。官方页面显示其 2024 Impact Factor 为 5.4。

spaGAPA 若定位为软件工具/方法，可以考虑：

- Application Note：篇幅短，适合工具成熟、benchmark 简洁的版本。
- Original Paper：若有完整算法、benchmark 和真实应用，力度更强。

**推荐目标**

如果真实 benchmark 补齐，但生物学发现不够顶刊级，Bioinformatics 是非常现实且合适的目标。

### 9.5 NAR Genomics and Bioinformatics

**适配度：中到高。**

NAR Genomics and Bioinformatics 关注 genomics/bioinformatics large-scale data analysis、novel computational methods、pipelines and workflows。官方页面显示 2024 Impact Factor 为 2.8。

适合条件：

- 工具完整。
- workflow 清晰。
- benchmark 可复现。
- 方法 novelty 中等，但实用性强。

如果 Bioinformatics 被拒，NAR Genomics and Bioinformatics 是合理备选。

### 9.6 BMC Bioinformatics

**适配度：高，保底选择。**

BMC Bioinformatics 接收 novel computational algorithms and software、models and tools，用于 biological data analysis。它更强调 scientifically valid，而不是 perceived impact。

如果目标是稳妥发表软件工具，BMC Bioinformatics 非常合适。但若 spaGAPA 的真实 benchmark 和生物学发现做得好，可以先投更高。

### 9.7 Briefings in Bioinformatics

**适配度：低到中。**

Briefings in Bioinformatics 更偏综述、方法综述和资源性长文。spaGAPA 作为单一新软件工具并不是最自然的目标，除非文章写成空间转录组 APA 方法系统 benchmark + 新工具框架。

### 9.8 GigaScience / Scientific Data / Patterns

**适配度：视数据资源而定。**

如果最终产生了一个高质量空间 APA benchmark dataset 和可复现 workflow，可以考虑这些资源/数据导向期刊。

---

## 10. 推荐投稿策略

### 路线 A：稳健软件论文

**目标**

Bioinformatics 或 NAR Genomics and Bioinformatics。

**最低要求**

- 修复 pipeline。
- 完成真实 MOB benchmark。
- 与 stAPAminer 公平比较。
- 补齐 examples/docs。
- 开源代码和数据处理脚本。

**文章主线**

> spaGAPA is a Python toolkit for uncertainty-aware spatial APA analysis, improving imputation and SVAPA detection over KNN-based workflows.

### 路线 B：方法学加强论文

**目标**

PLOS Computational Biology 或 Genome Biology。

**额外要求**

- 至少 2 到 3 个真实数据集。
- 完整 simulation + real benchmark + ablation。
- 证明 uncertainty improves downstream inference。
- 展示明确 biological findings。

**文章主线**

> Uncertainty-aware spatial Gaussian processes reveal robust spatial APA regulation in tissue architecture.

### 路线 C：冲击 Nature Methods

**目标**

Nature Methods。

**额外要求非常高**

- 方法不仅限于 APA，而是可推广到 spatial omics missingness/inference。
- 多平台、多组织、多任务验证。
- 强生物学问题。
- 外部用户容易使用。
- 结果显著优于现有方法。

**当前不建议作为第一投稿，除非后续真实数据结果远超预期。**

---

## 11. 论文前必须完成的实验清单

### 11.1 数据集

- [ ] MOB real APA matrix。
- [ ] Brain cortex real APA matrix。
- [ ] Embryo or developmental tissue real APA matrix。
- [ ] 至少一个 replicate dataset。
- [ ] metadata: layer/domain labels。
- [ ] APA site annotation。

### 11.2 Baselines

- [ ] Mean。
- [ ] Median。
- [ ] KNN-spatial。
- [ ] stAPAminer KNN expression-distance。
- [ ] stAPAminer full workflow。
- [ ] SPARK SVAPA。
- [ ] Moran's I。
- [ ] non-spatial GP。
- [ ] spaGAPA-GP。
- [ ] spaGAPA weighted SVAPA。

### 11.3 Metrics

- [ ] RMSE。
- [ ] MAE。
- [ ] Pearson。
- [ ] Spearman。
- [ ] R2。
- [ ] Bias。
- [ ] uncertainty calibration。
- [ ] runtime。
- [ ] memory。
- [ ] domain ARI/NMI。
- [ ] DEAPA F1。
- [ ] SVAPA overlap/reproducibility。
- [ ] GO/KEGG enrichment。

### 11.4 Ablation

- [ ] no spatial validation。
- [ ] no uncertainty weighting。
- [ ] RBF vs Matérn。
- [ ] full GP vs sparse GP。
- [ ] different dropout levels。
- [ ] different spot densities。
- [ ] different k in KNN。

### 11.5 Reproducibility

- [ ] `environment.yml`。
- [ ] `requirements.txt` clean。
- [ ] R environment for scAPAtrap/stAPAminer。
- [ ] benchmark Makefile。
- [ ] raw data download script。
- [ ] processed data generation script。
- [ ] figure generation script。
- [ ] random seeds fixed。

---

## 12. 建议的近期执行顺序

### Step 1：先修 pipeline 和数据 shape

目的：让 toy data 能端到端跑通。

输出：

- end-to-end toy test。
- corrected pipeline API。
- README 中最小示例可运行。

### Step 2：完成真实 MOB APA 数据构建

目的：从 expression proxy 切换到真实 APA。

输出：

- MOB `apa_matrix.csv`。
- MOB `coordinates.csv`。
- MOB `metadata.csv`。
- MOB `apa_sites.csv`/BED。

### Step 3：跑 stAPAminer 原始流程

目的：公平比较。

输出：

- stAPAminer RUD imputation。
- stAPAminer SVAPA。
- stAPAminer domain/DEAPA outputs。

### Step 4：统一 benchmark

目的：生成主图需要的表格。

输出：

- `combined_results.csv`。
- `real_cv_results.csv`。
- `ablation_results.csv`。
- `downstream_results.csv`。

### Step 5：生物学 case studies

目的：从工具论文升级到方法/发现论文。

输出：

- 3 到 5 个 top SVAPA genes。
- spatial APA maps。
- uncertainty maps。
- GO enrichment。
- tissue/domain interpretation。

---

## 13. 可用于论文的核心卖点措辞

### 13.1 简短版

spaGAPA introduces an uncertainty-aware Gaussian process framework for spatial alternative polyadenylation analysis, enabling spatially informed imputation, confidence-aware SVAPA detection, and reproducible benchmarking in spatial transcriptomics.

### 13.2 方法创新版

Unlike KNN-based APA recovery methods that impute missing APA indices using fixed local averages, spaGAPA models APA usage as a spatially correlated latent function and estimates posterior uncertainty for each gene-spot pair. This uncertainty is propagated into downstream spatial pattern detection, reducing false positives in SVAPA discovery.

### 13.3 竞品对比版

Compared with stAPAminer, spaGAPA explicitly models spatial covariance using tissue coordinates, quantifies imputation uncertainty, supports uncertainty-weighted spatial statistics, and provides a Python/AnnData-compatible workflow with simulation and real-data benchmarking.

---

## 14. 外部期刊信息来源

- Bioinformatics 官方 About 页面：该刊关注 genome bioinformatics 和 computational biology 新发展，并显示 2024 Impact Factor 为 5.4。  
  https://academic.oup.com/bioinformatics/pages/About
- NAR Genomics and Bioinformatics 官方 About 页面：该刊关注 genomics/bioinformatics large-scale data analysis、novel methods、pipelines/workflows，并显示 2024 Impact Factor 为 2.8。  
  https://academic.oup.com/nargab/pages/About
- Genome Biology 官方 aims and scope：覆盖 genomic/post-genomic biology，包括 research、new methods、software tools 和 reviews。  
  https://link.springer.com/journal/13059/aims-and-scope
- Nature Methods 官方 aims and scope：强调 novel methods、significant improvements、strong validation、biological application、performance comparison，并包括 spatial omics 与 computational/statistical methods。  
  https://www.nature.com/nmeth/aims
- PLOS Computational Biology 官方 journal information：Methods/Software 需要 outstanding importance、biological insight、reproducibility、data/code availability。  
  https://journals.plos.org/ploscompbiol/s/journal-information
- BMC Bioinformatics 官方 aims and scope：接收 biological data analysis 的 novel computational algorithms、software、models and tools。  
  https://link.springer.com/journal/12859/aims-and-scope

---

## 15. 最终判断

spaGAPA 的核心创新不是“写了一个 APA 工具”，而是把空间 GP、不确定性量化和 SVAPA detection 连接到了一个空间 APA 分析任务中。这个方向是有价值的，尤其是在空间组学越来越重视 missingness、dropout 和 uncertainty 的背景下。

但是，顶刊审稿人不会被“代码已经很多”和“模拟数据好看”说服。真正要打动他们，需要真实 APA 数据、严格 baseline、消融实验和生物学发现。当前最关键的短板是：真实数据尚未完成、pipeline 胶水层还需要修、stAPAminer 原流程还没有被严格复现。

如果按本文档的 P0 清单完成，spaGAPA 至少具备投 Bioinformatics / NAR Genomics and Bioinformatics 的条件；若真实数据 case study 足够强，可以冲 PLOS Computational Biology 或 Genome Biology；若想冲 Nature Methods，需要把 spaGAPA 从“空间 APA 工具”提升成“广泛适用于 spatial omics 的 uncertainty-aware inference framework”，并提供非常强的跨场景验证。

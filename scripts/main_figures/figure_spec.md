# spaGAPA 论文图表 Spec v2 — 参照论文校准设计书（NC/NAR 出版级 + Stereo 扩展钩子）

> **版本**：v2（2026-09-12）。本文件是图表集的**设计规格书**，非生成脚本。
> **权威输入**：`scripts/main_figures/figure_index.md`（v1 冻结面板级真值，2026-08-10 重生成）、
> `scripts/supplementary_figures/supp_figure_legends.md`（S1–S16）、
> `scripts/multicaller_validation/report.md`（caller×物种口径）、
> `pipeline_output/stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/uncertainty/uncertainty_calibration.json`（C2 已验证数字）。
> **策略**：演进 + 扩展，不推倒重来。v1 冻结数字**不重算**；本 spec 定义 v2 升级版设计与
> Stereo 扩展集成钩子。图表重生成按本 spec 的 Style Guide 执行，但不在撰写本 spec 的计划内。
> **标注约定**：`[PENDING: pipeline_output/... ]` = 产物未落地，落地后按旁注填数；其余数字均已核实。

---

## Style Guide（出版级硬性规范，重生成时强制执行）

NC（Nature Communications）/ NAR 标准设计。本节 6 条为硬性要求，先于一切设计原则。

1. **字体：全部 Arial。** matplotlib 侧 `rcParams['font.family']='sans-serif'` +
   `rcParams['font.sans-serif']=['Arial']`；R/ggplot2 侧 `theme(text=element_text(family="Arial"))`。
   现有 `_style.py` 的 DejaVu Sans 在图表重生成时替换。数学文本用
   `mathtext.fontset='custom'` + `mathtext.rm='Arial'`（避免 mathtext 回退非 Arial 字形）。
   轴标签/刻度 ≥7 pt（按版面缩放后），面板标题粗体 8–9 pt，图注字母 8 pt 粗体。
2. **版面：Nature 双栏宽 183 mm（7.2 in）为主图标准**；单栏小图 89 mm（3.5 in）。
   NAR 全宽 ~178–190 mm 兼容（取 183 mm 即同时满足两者）。矢量 PDF 主交付 +
   600 dpi PNG/TIFF 预览（v1 为 300 dpi PNG，v2 升级到 600 dpi）。
3. **配色：保留 Okabe-Ito 色盲安全主色**（`#0072B2` / `#D55E00` / `#009E73`，NC/NAR 均要求
   无障碍）。线宽 ≥0.5 pt，散点描边 0.2 pt。沿用 SpliceImpactR 式"事件类型为中心"的
   分布面板配色纪律：同一语义（方法/事件类型/域）在全图集内颜色恒定。
4. **统计图注规范（SCSES/NC 式）**：每个分布面板的图注必须写明：
   N（样本/基因/spot 数）、检验类型（two-sided Wilcoxon rank-sum / t-test，是否多重校正）、
   误差定义（SD / SEM / 95% CI）、箱线定义（median / Q25 / Q75 / 1.5×IQR whisker，单点为止）。
   条形图误差棒统一 SEM-of-replicates 并注明 replicate 数。
5. **Source Data（NC 投稿硬要求）**：每个主图的每个数据面板导出
   `SourceData_FigN[x].csv`（x = 面板字母）。本 spec 逐面板标注导出的数据帧与列
   （见各 Figure 节的 Source Data 表）。纯示意图面板（Fig1 全部、Fig3A、Fig4A、Fig7A）无 Source Data。
6. **作图代码参照仓库**：作图方法不确定时去参照论文的 GitHub 仓库读代码，仓库清单
   见附录 A（`tilgnerlab/Spl-IsoFind_reproducibility`、`fiszbein-lab/SpliceImpactR`、
   SCSES、`algbio/spl-IsoQuant`）。每个 spec 面板行的"参照模板来源"列注明仓库+文件。

---

## 设计原则

从 16 篇参照论文（BIB 5 + NAR/GB/GR/NC 11，清单见附录 B）提炼，逐条注明来源：

1. **叙事主链**：总览 → 可靠性审计 → 基准 → 全景统计 → 空间分解 → 基因案例 → 延伸。
   G4 四篇（Spl-ISO-Seq / Longcell / TUSCO / AF2-isoform）共享此骨架；spaGAPA 图集
   Fig1–Fig7 依次映射到该链的每一环。
2. **诚实基准文化**（TUSCO 式基准叙事强化，沿用仓库现有立场）：
   - per-gene mean 是强逐点 RMSE 基线，获胜时**公开显示**（Fig2B），不藏在补充材料；
   - 部署实现的经验斜率**不裸称 O(N²)****（Fig6A），近似邻域结构可压平观测缩放；
   - conformal 保证是**边际而非条件**覆盖，显式标注（Fig3 / Fig4E / S6D）。
3. **视觉语言**：Okabe-Ito 色盲安全主色 + 全 Arial + 183 mm 版面（v2 从 DejaVu Sans /
   ~178 mm 升级，见 Style Guide）；SpliceImpactR 式"事件类型为中心"的分布面板配色纪律。
4. **每图回答一个审稿人问题**（STIFT / SpaTM 面板模板）：每个主图节首行写明该问题，
   面板不回答该问题的内容移补充图。

---

## 图表总览表

| 编号 | 角色 | 回答的审稿人问题 | 核心指标 | 新数据升级点 | 状态 |
|---|---|---|---|---|---|
| Fig1 | 框架总览 | 方法解决什么问题、输入输出是什么 | 无（纯示意） | B 面板 caller 双通道 + 多物种/多平台条带 | 设计升级 |
| Fig2 | 插补基准 | GP 插补比 mean 好在哪、差在哪 | RMSE、PCC、spatial fidelity、ΔMoran's I、per-gene PCC | 新增 G（Moran's I 前后）+ H（10 折逐基因 PCC） | 冻结 + 扩 2 面板 |
| Fig3 | conformal 覆盖 | 区间是否真有有限样本保证、跨平台物种是否稳健 | coverage±二项 CI、校准偏差、区间宽度 | 11→16 样本（人视网膜 ×4 + 大鼠）+ D 平台/物种不变性面板 | 冻结 + 扩容 |
| Fig4 | 不确定性 | 不确定度可信吗、能做什么 | unc-error r、subgroup coverage、risk-coverage | 新增 G（Longcell 式基因异质性散点） | 冻结 + 扩 1 面板 |
| Fig5 | 域检测 | APA 空间域是否恢复解剖结构 | ARI/NMI、域 Moran's I、混淆矩阵 | C 改双组指标条形 + 新增 F（人视网膜 2×2 域图） | 冻结 + 改造 |
| Fig6 | 可扩展性 | 规模上去还能跑吗、结果稳吗 | runtime/memory 斜率、completion、Pareto、重采样 RMSE | 新增 E（dumbbell）+ F（下采样稳健性） | 冻结 + 扩 2 面板 |
| Fig7 | Stereo-seq showcase | 方法在亚细胞分辨率平台/多物种上成立吗 | PAS 数、observed fraction、3'富集、binning 一致性 | 单鼠 pilot → 三物种 showcase（mouse/human/rat） | 冻结 + 重构 |
| S1–S16 | 补充 | 见补充图节 | 见补充图节 | 沿用，仅重排归属引用 | 冻结 |
| S17 | Stereo 扩展 QC 总览 | 扩展样本质量是否达标 | PAS 数/观察比例/3'富集/域数 | 新增（5 样本） | 新增 |
| S18 | 跨 caller×物种 conformal 全表 | 保证是否独立于 caller 与物种 | coverage 全表 | 新增（multicaller 3 + 新 5 样本合并） | 新增 |

---

## 指标词典

每个指标一段：**定义 / 为何它 / 诚实边界**。

- **RMSE（held-out）**：遮蔽条目上的均方根误差，所有方法共用同一 mask（20% per-gene，
  seed 42）。为何：逐点精度的主指标。诚实边界：per-gene mean 在逐点 RMSE 上获胜
  （GP 仅在 ~15% 基因上胜出，S4），必须公开显示。
- **Pearson / Spearman（held-out）**：遮蔽条目上插补值 vs 真值的相关。为何：RMSE 对尺度
  敏感，相关刻画轮廓恢复。诚实边界：mean 基线同样在 Pearson/Spearman 上不落下风（Fig2B
  公开显示）。
- **spatial fidelity（GP 0.42 vs mean 0.00）**：每基因将插补与真值的空间轮廓（去基因均值后）
  做相关再聚合。mean 插补逐基因输出常数 → 空间轮廓为零向量 → 保真度按构造为 0.00。
  为何：这是 GP 相对 mean 的**结构性**优势所在（恢复空间梯度）。诚实边界：只度量梯度
  恢复，不度量逐点精度；数值依赖去均值口径，图注须写明定义。
- **Moran's I**：变量在 spot 空间图上的自相关。本图集三个用途：(i) Fig2G marker 基因
  插补前后 ΔMoran's I；(ii) Fig5C 域分配的空间连贯度（mean Moran's I of domains）；
  (iii) Fig7E bin 尺寸稳定性。诚实边界：对空间图构造与 bin 尺寸敏感，跨图集比较须同图构造。
- **ARI / NMI**：域检测对 MOB 5 层解剖标注的恢复度。为何：结构恢复的金标准比较。
  诚实边界：mean-impute 在 MOB 上也恢复相当结构（ARI 0.576 vs apa_dominant 展示配置），
  必须诚实报告（S10）。
- **conformal coverage ± 二项 CI**：测试点落入 ŷ±q̂ 的经验比例；每样本每水平的 90%
  二项 CI（Clopper–Pearson）。为何：有限样本边际保证的直接检验。诚实边界：保证是
  边际的，不是条件/子群的（Fig4E 高表达 bin ~0.82 欠覆盖公开显示）。
- **区间宽度**：每水平平均 2q̂（constant 模式）或 Σ2q̂ᵢσᵢ 型（locally adaptive 模式）。
  为何：coverage 之外的有效性维度——宽度无界也能覆盖。诚实边界：A（constant）与
  D（residual-spot）宽度差即模式取舍的量化（Fig3E）。
- **Winkler interval score**（Gneiting & Raftery 2007）：逐观测精确区间得分，越低越好，
  从持久化的 per-observation bounds 计算。诚实边界：早先 q̂+RMSE 近似高估 ~2.3×，已废弃（S6C）。
- **uncertainty–error r（pooled vs within-gene）**：不确定度 vs |误差| 的相关。pooled r ≈ 0.5
  部分由跨基因异质性驱动；within-gene median r ≈ 0.10–0.19 才是同基因内判别力。
  为何：区分"不确定度排序基因"与"排序同一基因内的 spot"两种用途。诚实边界：两者都报告，
  不用 pooled 冒充 within-gene（Fig4C/D）。
- **runtime / memory / 幂律斜率**：实测墙钟与峰值内存（合成网格 1k–100k spots，1200 s 上限）；
  斜率为实测点的幂律拟合。诚实边界：斜率反映**部署实现**（近似邻域结构可压平观测缩放），
  可能偏离教科书复杂度，不裸称 O(N²)。
- **completion**：方法 × 规模的完成矩阵（绿=完成，红=超时/OOM）。为何：可扩展性的
  最直接证据。诚实边界：阈值（1200 s 墙钟 + 内存上限）须在图注写明。
- **ΔPDUI（distal usage shift）**：基因级 distal PAS 使用比例的条件间/域间差，
  Spl-ISO-Seq ΔΠ 的空间版。为何：把 APA 空间变异压缩为可比的基因级量。
  诚实边界：需 ≥10 个共有限 spot 才稳定（multicaller 口径）；n=1 条件下只报告效应量
  不报 p 值（S11/S14）。
- **Jaccard（caller 重叠）**：两 caller PAS 集在 ±50 bp 聚类后的交集/并集。
  为何：量化 caller 依赖性。诚实边界：spaGAPA 实测 Jaccard 仅 0.05–0.18，但统计结论
  （coverage、distal-usage r）跨 caller 不变——重叠低不等于结论不稳（S18/multicaller report）。
- **per-gene PCC（遮蔽 CV）**：10 折遮蔽交叉验证下逐基因 held-out Pearson r 的分布
  （中位数 + IQR）。为何：Fig2B 的数据集级相关会掩盖基因级异质性。诚实边界：折间
  相关非独立，IQR 描述分布不用于推断。

---

## 数据资产映射

逐图列 `pipeline_output/` 数据源（v1 继承 + Stereo 扩展新增）；除表内 `[PENDING: pipeline_output/…]`
标注外，所有路径均已核实存在。

| 图 | 数据源（`pipeline_output/` 下） |
|---|---|
| Fig2 | `benchmark_mean_transparent/transparent_comparison.csv`；`benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json` |
| Fig3 | `conformal_validation/all_samples_coverage.csv`（含 `winkler_*` 列）+ `summary.json`；`per_observation_bounds.npz`；新增：`stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/uncertainty/uncertainty_calibration.json`（已验证） |
| Fig4 | `uncertainty_corr_improvement/per_dataset_corr.csv`；`uncertainty_within_gene_audit/per_gene_corr_distribution.csv`；`risk_coverage_curve/{risk_coverage_data.csv,summary.json}` |
| Fig5 | `mob_domain_recovery/{spagapa_metrics.json,spagapa_domains.csv}` |
| Fig6 | `benchmark_runtime/runtime_table.csv` |
| Fig7 | `gse263789_stereo_pilot/spagapa_downstream_full/`；`stereo_binning_consistency/binning_correlation.csv`；`gse263789_ad_vs_wt_differential/sample_level_effect_sizes.csv` |
| Stereo 扩展 | 已核实：`gse293464_retina/GSM8882885_binned/`（C2：22,762 PAS / 14,905 bins / nnz 12,697,687 / sparsity 0.9626）、`gse293464_retina/GSM8882885_scapatrap_raw_binned/`、`stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/`、`gse293464_retina/GSM8882884_scapatrap_raw/`（B4 scAPAtrap 已完成）、`gse293464_retina/GSM888{86,87}_{retag,saw}/`、`gse333693_thymus/GSM9770943_scapatrap_raw/`（大鼠 23,138 PAS）；待落地：`[PENDING: pipeline_output/gse293464_retina/GSM888{84,86,87}_scapatrap_raw_binned/ 与 pipeline_output/gse333693_thymus/GSM9770943_binned/ —— bin200 矩阵落地后进 Fig3/Fig5/Fig7/S17]` |
| 多 caller | `scripts/multicaller_validation/report.md` 所列 3 数据集 Sierra/scAPAtrap 产物与 coverage 表 |

---

## Figure 1 — 框架总览（`fig1_framework_overview`）

**定位**：一句话讲清 spaGAPA 是什么、解决哪三个缺口、数据从哪来、核心推断与校准两大件。

**审稿人问题**：*“这个方法到底做了什么、输入输出是什么？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | 三缺口：稀疏 APA-usage 矩阵、无校准不确定度、邻域 comparator 在 ≥42k locations 超限 | 概念示意 | 无 | 无（示意） | SCLR-seq 综述 Fig2a 分步总览 |
| B | 输入链：空间 reads → PAS peak counts → APA-usage 矩阵 [gene × spot] → spaGAPA；**caller-agnostic 双通道（scAPAtrap + Sierra）**；**多物种/多平台条带（mouse / human / rat × Visium / Stereo-seq）** | 流程示意 | 无 | 无（示意）；条带内容对应 `multicaller_validation/` 与 Stereo 扩展真实样本 | STIFT Fig1；spl-IsoQuant 流程示意规范（`algbio/spl-IsoQuant`） |
| C | 稀疏 GP：后验均值 + 68/95% CI，M inducing points，O(N·M²) | 概念示意 | 无 | 无（示意） | 沿用 v1 |
| D | split-conformal：train / calibrate / test 划分；sᵢ=\|yᵢ−ŷᵢ\|；分位数 q̂ → ŷ±q̂ 有限样本边际覆盖 | 概念示意 | 无 | 无（示意） | 沿用 v1 |

**指标定义引用**：无（纯示意，无定量面板）。

**诚实边界**：示意图中的矩阵稀疏度/CI 宽度为示意值，图注注明“schematic”；O(N·M²) 写
复杂度而非实测性能（实测见 Fig6）。

**与现有脚本关系**：`fig1_overview.py` 沿用；B 面板重绘（加双通道与物种/平台条带）。

**新数据集成钩子**：无 data-pending——B 面板条带内容（rat / human Stereo）数据已在，
仅示意重绘。

**Source Data**：无（纯示意）。

---

## Figure 2 — 插补基准（`fig2_spatial_vs_mean_benchmark`）

**定位**：5 方法 × 2 数据集（GSE183456 kidney、GSE220442 brain）同 mask（20% per-gene，
seed 42）头对头；诚实框架下展示 GP 的结构性优势与逐点劣势。

**审稿人问题**：*“GP 插补比 per-gene mean 好在哪、差在哪？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | 遮蔽设计：20%-of-observed per-gene holdout，所有方法同条目遮蔽 | 示意+矩阵缩略 | 无 | 无（示意） | 沿用 v1 |
| B | 逐数据集指标（配对点图）：RMSE（左）、Pearson/Spearman（右）；**mean 逐点获胜公开显示** | 配对点图 | RMSE、PCC | `benchmark_mean_transparent/transparent_comparison.csv`；`benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json` | SCSES 分布面板（raincloud 纪律） |
| C | spatial fidelity：GP 0.42 vs Mean 0.00 | 条形+注释 | spatial fidelity（见词典） | `benchmark_mean_transparent/transparent_comparison.csv` | 沿用 v1 |
| D | 精度–空间 2D 散点：每方法每数据集一点 | 散点 | RMSE × spatial fidelity | 同 B | 沿用 v1 |
| E | 代表性空间重构：held-out 真值 / mean / 真实 spaGAPA-GP（`SparseGPImputer`，peak_94938）/ \|error\| | 空间热图四联 | 无（展示） | `_data/fig2E_gse183456_peak_94938_gp_posterior.csv` | Spl-IsoFind_reproducibility `figures/Figure 3.ipynb` 空间散点写法 |
| F | ΔRMSE 分层：真实逐基因 Δ=GP−mean 按 Moran's I 五分位分层；Δ>0 = mean 更好，**差距不随空间信号缩小** | 箱线图（median/Q25/Q75/1.5×IQR） | 逐基因 ΔRMSE | `_cache/s4_stratification.csv` | 沿用 v1；SpaTM 式分层 |
| **G（新）** | **Moran's I 前后对比**：marker 基因插补前（观测）vs GP 后的空间自相关 | 箱线图 + 成对连线 | **ΔMoran's I** | 由 `benchmark_mean_transparent/` 遮蔽后观测矩阵与 GP 后验导出（新增分析：`_cache/fig2G_moran_shift.csv`） | **SpaTM Fig3b** 模板 |
| **H（新）** | **遮蔽 CV 逐基因 PCC 分布**：10 折 per-gene 遮蔽，逐基因 held-out Pearson r | 分布图（raincloud 或小提琴） | **per-gene PCC 中位数 + IQR** | 同 A/B 基准输入重跑 10 折（新增分析：`_cache/fig2H_percgene_pcc.csv`） | **SpaTM Fig2** 式 10 折 |

**指标定义引用**：RMSE、Pearson/Spearman、spatial fidelity、Moran's I、per-gene PCC（见指标词典）。

**诚实边界**：mean 在逐点 RMSE/Pearson/Spearman 获胜（B 面板公开）；spatial fidelity 0.42
vs 0.00 是 GP 的结构性优势但只度量梯度恢复；F 面板最高空间信号五分位反而 Δ 中位数最大
（GP 更差）——保留 v1 的诚实结论。

**与现有脚本关系**：`fig2_benchmark.py` 沿用 A–F；新增 G/H 面板代码（新增分析脚本，
输入为现有基准数据）。

**新数据集成钩子**：G/H 为新增分析（数据已在，导出待做），无待落地产物。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| B | `SourceData_Fig2B.csv` | method, dataset, rmse, pearson, spearman |
| C | `SourceData_Fig2C.csv` | method, dataset, spatial_fidelity |
| D | `SourceData_Fig2D.csv` | method, dataset, rmse, spatial_fidelity |
| E | `SourceData_Fig2E.csv` | x, y, truth, mean, gp_posterior, abs_error |
| F | `SourceData_Fig2F.csv` | gene, moran_quintile, delta_rmse |
| G | `SourceData_Fig2G.csv` | gene, moran_observed, moran_gp, delta_moran |
| H | `SourceData_Fig2H.csv` | gene, fold, heldout_pcc |

---

## Figure 3 — conformal 覆盖（`fig3_conformal_marginal_coverage`）

**定位**：split-conformal 区间的有限样本边际覆盖检验，v2 从 11 Visium 样本扩容到
**16 样本**（11 Visium + 4 人视网膜 Stereo + 1 大鼠 Stereo），并新增平台/物种不变性面板。

**审稿人问题**：*“覆盖保证是真的吗？跨平台、跨物种还成立吗？”*

**已验证新数字（C2 = GSM8882885，s26_ra，人视网膜类器官，Stereo-seq bin200）**：
locally-adaptive coverage **0.7999 / 0.8993 / 0.9498**（80/90/95%），global 0.7996 / 0.8995 /
0.9498；q̂ = 23.31 / 46.47 / 87.64（locally adaptive），14.41 / 28.74 / 54.36（global）；
n_test = 1,265,290（源：`pipeline_output/stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/uncertainty/uncertainty_calibration.json`，与 `pipeline_output/gse293464_retina/GSM8882885_binned/qc_summary.json`
的 22,762 PAS / 14,905 bins / nnz 12,697,687 一致）。

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | split-conformal 流程：calibrate → \|yᵢ−ŷᵢ\| → q̂ → ŷ±q̂ → coverage ≥ 1−α | 流程示意 | 无 | 无（示意） | 沿用 v1 |
| B | 每样本 80/90/95 覆盖点图 + 名义线 + 90% 二项 CI 带（**16 样本**：11 Visium 冻结 + C2 已验证 + 其余 4 落地后并入） | 点图 + CI 带 | coverage ± 二项 CI | `conformal_validation/all_samples_coverage.csv` + `stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/uncertainty/uncertainty_calibration.json`；其余 3 视网膜 + 大鼠：`[PENDING: pipeline_output/stereo_expansion_downstream/gse293464_GSM888{84,86,87}_*/uncertainty/uncertainty_calibration.json 及 pipeline_output/stereo_expansion_downstream/gse333693_GSM9770943_*/uncertainty/uncertainty_calibration.json —— 落地后按 C2 同口径并入，并更新 523,174+1,265,290+… 的总测试点数 ]` | Spl-IsoFind_reproducibility `figures/Figure 4.ipynb` 模拟 PR 面板 |
| C | 校准曲线：nominal vs empirical，y=x 参考线 + 平均偏差标注；11 样本冻结值：mean abs deviation 0.21%（80%）/ 0.16%（90%）/ 0.10%（95%），max deviation 同步报告；**扩容后重算全集合统计，C2 单样本偏差 = \|0.7999−0.80\|=0.01 pp 级，直接写入** | 对角线图 | coverage 偏差 | 同 B | 沿用 v1 + 扩容 |
| **D（新）** | **平台/物种不变性**：coverage 按 Visium-mouse / Stereo-mouse / Stereo-human / Stereo-rat 分组箱线（分组口径参照 multicaller report 的 species×platform 表）。Stereo-mouse = GSE263789（当前不在 conformal_validation 集内，见钩子） | 分组箱线 | coverage 按组 | `conformal_validation/all_samples_coverage.csv` 样本元数据 + Stereo 扩展 5 样本；GSE263789 补跑：`[PENDING: pipeline_output/conformal_validation/ 新增 GSE263789 bin200 行（或 pipeline_output/stereo_expansion_downstream/gse263789_*/uncertainty/）—— 补跑后并入 D 组，全集变 17；若不补跑，D 呈现 Visium 按物种 + Stereo-human + Stereo-rat，并在图注声明缺失组 ]` | multicaller report 分组方式 |
| E | 区间宽度：Constant（A）vs Residual-spot（D） | 配对条形/点图 | 平均 2q̂ | `conformal_validation/all_samples_coverage.csv` + C2 `q_hat` 值 | 沿用 v1（v1-D 移位为本图 E） |
| F | 空间实例：GSE183456 peak_20919 真实后验 μ / σ / \|error\| / covered-vs-uncovered，conformal q̂₉₀=0.284 | 空间热图四联 | 无（展示） | `per_observation_bounds.npz` + 后验导出 | 沿用 v1（v1-E 移位为本图 F） |

**指标定义引用**：conformal coverage ± 二项 CI、区间宽度、Winkler（见指标词典）。

**诚实边界**：80% 水平偏差 0.21% > 0.2%，**不得**笼统写“within 0.2%”（v1 立场保留）；
保证是边际的（Fig4E / S6D 展示子群欠覆盖）；C2 为 16 样本中唯一已验证 Stereo 样本，
扩容数字落地前图注写明 per-sample n。

**与现有脚本关系**：`fig3_conformal.py` 沿用并扩容 B；插入新 D 面板；v1 的 D（宽度）→ E、
E（空间实例）→ F 字母移位；S15 偏差 forest **保持补充图**（主图不膨胀，v1 决策沿用）。

**新数据集成钩子**：见面板表 B/D 的 `[PENDING: pipeline_output/stereo_expansion_downstream/ 与 pipeline_output/conformal_validation/ 增量产物]`（3 视网膜样本 downstream、大鼠 downstream、GSE263789 补跑）。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| B | `SourceData_Fig3B.csv` | sample, platform, species, level, coverage, ci_low, ci_high, n_test |
| C | `SourceData_Fig3C.csv` | level, nominal, empirical_pooled, mean_abs_dev, max_abs_dev |
| D | `SourceData_Fig3D.csv` | sample, group(platform×species), level, coverage |
| E | `SourceData_Fig3E.csv` | sample, level, mode, mean_interval_width |
| F | `SourceData_Fig3F.csv` | x, y, posterior_mean, posterior_sd, abs_error, covered |

---

## Figure 4 — 不确定性（`fig4_noise_and_risk_coverage`）

**定位**：4 种不确定度模型（A constant / B local-gene / C spatial-spot / D residual-spot）
× 5 数据集；不确定度的可信度（与误差的相关）与用途（风险分层、triage）。

**审稿人问题**：*“不确定度可信吗？它能用来做什么？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | 四种不确定度模型示意 | 概念示意 | 无 | 无（示意） | 沿用 v1 |
| B | 多目标方法选择：x=uncertainty–error r，y=区间宽度，bubble=\|coverage 偏差\| | 气泡散点 | r、宽度、\|偏差\| | `uncertainty_corr_improvement/per_dataset_corr.csv` | 沿用 v1 |
| C | pooled r 逐数据集：A vs D；A 近零、D 过 r≥0.3 | 分组条形 | pooled r | 同 B | SpliceImpactR 事件类型分布面板纪律（`fiszbein-lab/SpliceImpactR`） |
| D | per-gene vs pooled r：within-gene median r ≈ 0.10–0.19；pooled ≈ 0.5 部分由跨基因异质性驱动（透明报告） | 分布对比 | within-gene r | `uncertainty_within_gene_audit/per_gene_corr_distribution.csv` | 沿用 v1 |
| E | 经验子群覆盖：左 tissue × uncertainty-quintile 热图；右表达量 bin——高表达 bin 欠覆盖（~0.82），作为诚实局限展示 | 热图 + 条形 | subgroup coverage | `conformal_conditional_coverage/*.csv` | pheatmap 纪律（SpliceImpactR 侧 R 惯例） |
| F | risk-coverage 曲线（RMSE）：retention vs RMSE；不确定度 triage 在 80% retention 下 −23% RMSE（5 数据集配对，p=0.004） | 折线 + 配对点 | RMSE vs retention | `risk_coverage_curve/{risk_coverage_data.csv,summary.json}` | 沿用 v1 |
| **G（新）** | **Longcell 式异质性散点**：x=基因级 distal-usage 均值，y=空间离散度 SD，颜色=GP 平均不确定度——展示“哪些基因的 APA 空间变异性最值得信任” | 三变量散点 | per-gene usage mean × spatial SD × median uncertainty | 基因级 usage/不确定度汇总（由 GP 后验 + `conformal_validation/per_observation_bounds.npz` 导出，新增分析：`_cache/fig4G_gene_heterogeneity.csv`） | **Longcell** 异质性散点模板 |

**指标定义引用**：uncertainty–error r（pooled vs within-gene）、conformal coverage、ΔPDUI（见指标词典）。

**诚实边界**：within-gene r ≈ 0.10–0.19（v1 冻结）——不确定度在同基因内判别力中等，
不得用 pooled r≈0.5 冒充；高表达子群欠覆盖（~0.82）公开显示；triage 收益为 5 数据集
配对结果（p=0.004，来源 `summary.json`）。

**与现有脚本关系**：`fig4_noise.py` 沿用 A–F；新增 G 面板代码。

**新数据集成钩子**：G 为新增分析（数据已在，导出待做），无待落地产物。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| B | `SourceData_Fig4B.csv` | method, dataset, unc_error_r, interval_width, abs_cov_dev |
| C | `SourceData_Fig4C.csv` | method, dataset, pooled_r |
| D | `SourceData_Fig4D.csv` | dataset, gene, within_gene_r, pooled_r |
| E | `SourceData_Fig4E.csv` | subgroup_type, subgroup, coverage, nominal |
| F | `SourceData_Fig4F.csv` | dataset, retention, rmse_triaged, rmse_all |
| G | `SourceData_Fig4G.csv` | gene, usage_mean, spatial_sd, median_uncertainty |

---

## Figure 5 — 域检测（`fig5_domain_recovery`）

**定位**：MOB（小鼠嗅球，5 解剖层，260 spots）上 APA 空间域恢复；v2 把 C 面板升级为
双组指标条形，并新增人视网膜类器官域图作为跨数据集证据。

**审稿人问题**：*“APA 信号定义的空间域恢复解剖结构吗？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | MOB 解剖：金标准 5 层（GCL/GL/MCL/ONL/OPL） | 空间标注图 | 无 | MOB 注释（v1 冻结） | 沿用 v1 |
| B | 域图（公平比较）：同一 Leiden 流程（空间图/表达组件/权重/resolution），**只有 APA 来源不同**——per-gene **mean**-imputed vs spaGAPA **GP**-imputed，两者都是真实运行 | 空间域图 | 无（展示） | `mob_domain_recovery/spagapa_domains.csv` | 沿用 v1 |
| **C（改造）** | **STIFT 式双组指标条形**：结构恢复组（ARI/NMI）+ 空间连贯组（mean Moran's I of domains），各 weight 配置（mean_apa / apa_dominant / spatial_apa / balanced / expression_apa）并排 | 双组条形 | ARI/NMI + 域 Moran's I | `mob_domain_recovery/spagapa_metrics.json` + `spagapa_domains.csv`（域 Moran's I 为新增派生） | **STIFT Fig2c** 双组条形 |
| D | 域–层混淆矩阵：展示配置的域与 5 解剖层的真实对应（取代早先图注与数据矛盾的 Moran's-I 面板） | 混淆热图 | 混淆比例 | 同 B | SpliceImpactR pheatmap 纪律 |
| E | 代表性真实梯度基因：Trim16 / Gdap2 / Tor3a / Pcdhb21——稀疏输入（上）vs 真实 GP 重构（下），逐基因色标 | 空间热图成对 | 无（展示） | v1 冻结导出 | Spl-ISO-Seq Fig1 空间三件套 |
| **F（新）** | **人视网膜类器官域图 2×2**：RA± × 16/26wk 四样本（GSM8882884/B4、GSM8882885/C2-s26_ra、GSM8882886/D2、GSM8882887/D4，条件映射按 GEO 元数据）域着色 + 不确定度叠加 | 2×2 空间域图 | 域数 + 不确定度叠加 | C2：`[PENDING: pipeline_output/stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/domains/domain_map.csv —— 落地后填域数与不确定度分层]`；B4/D2/D4 上游：`[PENDING: pipeline_output/gse293464_retina/GSM888{84,86,87}_scapatrap_raw_binned/ —— bin200 矩阵落地后再跑域检测（B4 scAPAtrap 已完成但 bin200 未落地，D2/D4 在 retag/saw 阶段）]` | **STIFT Fig2c** + Longcell 组织域图 |

**指标定义引用**：ARI/NMI、Moran's I（见指标词典）。

**诚实边界**：mean-impute 在 MOB 上恢复相当结构（ARI 0.576 vs 展示配置；全局最优是
expression_apa ARI 0.597 → S10）——诚实报告；随机种子稳定性未评估（S10 注记沿用）；
F 面板为类器官体外数据，域的生物学对应不做解剖金标准声明。

**与现有脚本关系**：`fig5_domain.py` 沿用 A/B/D/E；C 面板重绘（加域 Moran's I 派生）；
新增 F 面板脚本（Stereo 域检测流程复用 MOB 管线）。

**新数据集成钩子**：F 面板待落地产物见面板表 `[PENDING: pipeline_output/stereo_expansion_downstream/ 域检测产物 与 pipeline_output/gse293464_retina/GSM888{84,86,87}_scapatrap_raw_binned/]`（C2 域检测下游 + B4/D2/D4 bin200）。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| B | `SourceData_Fig5B.csv` | spot_x, spot_y, apa_source, domain_label |
| C | `SourceData_Fig5C.csv` | config, ari, nmi, mean_domain_moran |
| D | `SourceData_Fig5D.csv` | domain_label, layer, fraction |
| E | `SourceData_Fig5E.csv` | gene, x, y, usage_input, usage_gp |
| F | `SourceData_Fig5F.csv` | sample, x, y, domain_label, mean_uncertainty |

---

## Figure 6 — 可扩展性（`fig6_scalability`）

**定位**：合成网格 1k–100k spots、4 方法、1200 s 墙钟上限的规模行为；v2 补两端点对比
与深度稳健性。

**审稿人问题**：*“规模上到 10 万 spot 还能跑吗？结果对测序深度敏感吗？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | 运行时缩放 + 经验斜率：实测点 + 幂律拟合（spaGAPA-fast 0.84 / spaGAPA-accuracy 1.01 / stAPAminer 1.01 / spvAPA 0.49）；斜率反映部署实现，近似邻域结构可压平观测缩放，**不裸称 O(N²)** | log-log 散点+拟合 | 幂律斜率 | `benchmark_runtime/runtime_table.csv` | 沿用 v1 |
| B | 峰值内存缩放：spaGAPA-fast 8.8 GB @100k；R 工具超限 | log-log 折线 | 峰值内存 | 同 A | 沿用 v1 |
| C | 完成矩阵：方法 × 规模，绿=完成 / 红=超时或 OOM | 矩阵热图 | completion | 同 A | 沿用 v1 |
| D | spatial-fidelity–runtime Pareto：Spatial-KNN 保真度最高（0.91）但非概率且不可扩展；spaGAPA-GP 中等保真 + 不确定度 + 可扩展，速度为 stAPAminer/spvAPA 的 6–7×；mean 最快但保真 ~0 | Pareto 散点 | fidelity × runtime | 同 A + 基准保真值 | 沿用 v1 |
| **E（新）** | **TUSCO 式 dumbbell**：1k vs 100k 两端点，runtime 与 memory 双指标，4 方法各一条哑铃线 | 哑铃图 | Δruntime、Δmemory（两端点比） | `benchmark_runtime/runtime_table.csv`（两端点行） | **TUSCO** dumbbell 模板 |
| **F（新）** | **下采样稳健性箱线**：100 次重采样（子采样观测条目）下 RMSE 分布，按方法分面 | 箱线图（median/Q25/Q75/1.5×IQR） | 重采样 RMSE 分布 | 基准数据下采样重跑（新增分析：`_cache/fig6F_downsample_rmse.csv`） | **Spl-ISO-Seq Fig5e** 式深度稳健性 |

**指标定义引用**：runtime/memory/幂律斜率、completion、RMSE（见指标词典）。

**诚实边界**：斜率是部署实现的实测（v1 立场保留）；1200 s 墙钟 + 内存上限阈值写入图注；
F 面板的重采样针对观测条目子采样（深度代理），非文库级别重测序。

**与现有脚本关系**：`fig6_scalability.py` 沿用 A–D；新增 E/F 面板代码（E 用现有
runtime_table 两端点行；F 为新增分析）。

**新数据集成钩子**：F 为新增分析（数据已在，导出待做），无待落地产物。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| A | `SourceData_Fig6A.csv` | method, n_spots, runtime_s, fitted_slope |
| B | `SourceData_Fig6B.csv` | method, n_spots, peak_mem_gb |
| C | `SourceData_Fig6C.csv` | method, n_spots, completed |
| D | `SourceData_Fig6D.csv` | method, spatial_fidelity, runtime_s |
| E | `SourceData_Fig6E.csv` | method, endpoint(1k/100k), runtime_s, peak_mem_gb |
| F | `SourceData_Fig6F.csv` | method, resample_id, rmse |

---

## Figure 7 — Stereo-seq 三物种 showcase（`fig7_stereo_seq`）

**定位**：从单鼠 AD 脑 pilot 升级为 **mouse / human / rat 三物种 showcase**，证明方法在
亚细胞分辨率平台与跨物种上成立。已核实数字：mouse 21,455 PAS；human（C2 视网膜）
22,762 PAS / 14,905 bins / nnz 12,697,687（observed ≈ 3.7%，sparsity 0.9626）；
rat（GSM9770943 胸腺）23,138 PAS（`peaks_meta.csv.gz` 23,139 行含表头）。

**审稿人问题**：*“方法在 Stereo-seq 与多个物种上成立吗？结果依赖 bin 选择吗？”*

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板来源 |
|---|---|---|---|---|---|
| A | 工作流：FASTQ+mask → SAW 8.2.2 → retag → scAPAtrap → PAS 矩阵 → spaGAPA；**mask 获取双途径标注（GEO / STOmics）** | 流程示意 | 无 | 无（示意） | spl-IsoQuant 流程示意规范；沿用 v1 |
| B | 3'QC：mouse 主面板（TES-distance 分布 + 累积 3'富集：**53.4% ≤500 bp，73.5% ≤2 kb**）+ **大鼠/人并排小面板** | 直方图 + 累积曲线 | 3'富集比例 | mouse：`gse263789_stereo_pilot/spagapa_downstream_full/`（v1 冻结）；人/大鼠：`[PENDING: pipeline_output/stereo_expansion_downstream/qc/three_prime_enrichment_all_species.csv —— 落地后填人/大鼠的 ≤500bp 与 ≤2kb 比例]` | 沿用 v1 + 扩展 |
| C | 规模与稀疏度三联：**mouse 20.7M DNB / 21,455 PAS / 15,235 bin200 spots / ~10.3% observed**；**human 22,762 PAS / 14,905 bins / 12,697,687 nnz / ~3.7% observed**；**rat 23,138 PAS** + observed fraction 对比条（nnz/矩阵元素，单一口径，不与原始读段数混同） | 三联条形 | PAS 数、observed fraction | mouse：v1 冻结；human：`gse293464_retina/GSM8882885_binned/qc_summary.json`（已核实）；rat PAS：`gse333693_thymus/GSM9770943_scapatrap_raw/peaks_meta.csv.gz`（已核实）；rat bins/nnz：`[PENDING: pipeline_output/gse333693_thymus/GSM9770943_binned/qc_summary.json —— 落地后填 bins 数与 observed fraction]` | Spl-IsoFind_reproducibility `figures/Figure 1.ipynb` 全景面板 |
| D | 空间三件套 ×3 物种：每物种 UMI 全景 → 域图 → 基因案例。mouse 基因案例（bin200，mm10 注释）：**Cdk8**（peak_69378）/ **Apoe**（peak_312125，AD）/ **Gnb1l**（peak_415285）；human 待选高空间变异 PAS（C2）；rat 待选（GSM9770943） | 空间热图三联 ×3 | 无（展示；选基因按空间变异度排序） | mouse：`gse263789_stereo_pilot/spagapa_downstream_full/`（v1 冻结）；human 候选：`gse293464_retina/GSM8882885_binned/apa_matrix.csv` + 后验；域图三物种与人/大鼠基因案例：`[PENDING: pipeline_output/stereo_expansion_downstream/gse293464_GSM8882885_s26_ra/domains/domain_map.csv 及 pipeline_output/stereo_expansion_downstream/gse333693_GSM9770943_*/domains/domain_map.csv —— 域图与基因案例落地后填位]` | **Spl-ISO-Seq Fig1B–D** 三件套 + **Longcell Fig7A/E** |
| E | binning 稳健性（mouse，v1 冻结：cross-bin Pearson r 50→100、50→200 + Moran's I vs bin size）+ **新增跨物种 PAS 数/深度对齐小面板**（三物种 PAS 数与中位 bin 深度并排，rat 深度 `[PENDING: pipeline_output/gse333693_thymus/GSM9770943_binned/ 落地后填]`） | 折线 + 小条形 | cross-bin r、Moran's I vs bin、PAS/深度对齐 | mouse：`stereo_binning_consistency/binning_correlation.csv`；human：`gse293464_retina/GSM8882885_binned/`；rat：同上 PENDING | 沿用 v1 + 扩展 |

**指标定义引用**：Moran's I、Pearson（cross-bin）、observed fraction（见指标词典与 C 面板单一口径注记）。

**诚实边界**：v1 已移除 41-domain + raw-σ 不确定度面板（200×200 bin-grid × inducing-point
几何伪影、域数过分割）——该决策沿用，v2 域图须通过 bin 与 inducing 几何审计后才入 D 面板；
human observed ~3.7% 低于 mouse ~10.3%（平台/组织差异），对比条图注注明不可解读为方法差异；
AD vs WT 差异分析保持 S14（n=1，仅效应量）。

**与现有脚本关系**：`fig7_stereo_seq.py` 重构为三物种布局：A 沿用；B/C/D 扩展多物种；
E 沿用 mouse 冻结 + 新增对齐小面板。

**新数据集成钩子**：见面板表 B/C/D/E 的 `[PENDING: pipeline_output/…]` 标注（人/大鼠 3'QC、
大鼠 bin200 与 qc_summary、三物种域图与人/大鼠基因案例）。

**Source Data**：

| 面板 | 文件 | 数据帧 / 列 |
|---|---|---|
| B | `SourceData_Fig7B.csv` | species, sample, tes_distance_bin, density, cum_fraction |
| C | `SourceData_Fig7C.csv` | species, sample, n_pas, n_bins, nnz, observed_fraction |
| D | `SourceData_Fig7D.csv` | species, sample, panel_type(umi/domain/gene), x, y, value, gene/peak |
| E | `SourceData_Fig7E.csv` | species, bin_pair, pearson_r；species, n_pas, median_bin_depth |

---

## Figures S1–S18（补充图）

S1–S16 沿用 v1（冻结，数字不重算；标题与数据源以
`scripts/supplementary_figures/supp_figure_legends.md` 为准），仅重排对 v2 主图的归属引用；
S17/S18 为 v2 新增。逐张一行：

| 编号 | 标题 | 指标 | 数据源 | 状态 |
|---|---|---|---|---|
| S1 | 数据集总览（32 样本 Visium 基准） | PAS 数/spot 数 | `supplementary_figures/` 缓存（v1） | 沿用（归 Fig2/fig3 背景） |
| S2 | 稀疏 GP inducing-point 敏感性（GSE183456） | held-out RMSE vs m | v1 缓存 | 沿用（归 Fig1C） |
| S3 | 遮蔽比例敏感性（10/20/30/50%） | median per-gene RMSE | v1 缓存 | 沿用（归 Fig2A） |
| S4 | mean 基线分层（GP 仅 ~15% 基因胜出） | 逐基因 ΔRMSE | v1 缓存 | 沿用（归 Fig2F） |
| S5 | 基准参数表（建议转补充表） | — | v1 缓存 | 沿用（注记保留） |
| S6 | 全样本 conformal 覆盖（11 样本 + Winkler + 子群） | coverage/宽度/Winkler | `conformal_validation/` | 沿用，扩容后同步 16 样本版 |
| S7 | 泄漏审计 + LOOCV 方法选择稳定性 | clean vs leaked r | v1 缓存 | 沿用（归 Fig4） |
| S8 | 逐基因不确定度–误差相关分布 | per-gene r | `uncertainty_within_gene_audit/` | 沿用（归 Fig4D） |
| S9 | 空间块 vs 随机划分 conformal | coverage 差 | v1 缓存 | 沿用（归 Fig3） |
| S10 | 域恢复权重配置稳健性（expression_apa 全局最优 ARI 0.597） | ARI/NMI vs resolution | `mob_domain_recovery/` | 沿用（归 Fig5C） |
| S11 | 伪重复分析（GSE220442；spot 级 91 vs donor 级 0 显著基因） | 显著基因数/效应量一致性 | v1 缓存 | 沿用（归 Fig2 伦理注记） |
| S12 | Stereo QC 与 binning 敏感性（GSE263789） | 染色体 PAS 分布/TES 距离 | v1 缓存 | 沿用（归 Fig7B/E） |
| S13 | 批次校正 QN vs Harmony（APA 矩阵，初探） | PCC vs 残余批次信号 | v1 缓存 | 沿用（可选模块注记） |
| S14 | AD vs WT 效应量（n=1/条件，仅描述） | Δ（效应量） | `gse263789_ad_vs_wt_differential/` | 沿用（归 Fig7D） |
| S15 | conformal 偏差 forest（v1 从 Fig3 移出，**保持补充**） | per-sample 偏差 + CI | `conformal_validation/` | 沿用（v2 决策：主图不膨胀） |
| S16 | MAE risk-coverage（Fig4F 配对单元） | MAE vs retention | `risk_coverage_curve/` | 沿用（归 Fig4F） |
| **S17（新）** | **Stereo 扩展全样本 QC 总览**：每样本 PAS 数 / 观察比例 / 3'富集 / 域数小多联（4 视网膜 + 1 大鼠 + mouse pilot 对照） | PAS 数/observed/3'富集/域数 | `gse293464_retina/*_binned/qc_summary.json`、`[PENDING: pipeline_output/gse333693_thymus/GSM9770943_binned/qc_summary.json —— 落地后填大鼠行]`、`gse263789_stereo_pilot/` | 新增（SCOTCH 式平台基准面板参照） |
| **S18（新）** | **跨 caller × 跨物种 conformal 全表**：multicaller 3 数据集（kidney/brain/colon，Sierra vs scAPAtrap，max 偏差 0.5 pp）+ 新 5 样本合并 coverage 表 | coverage 全表（80/90/95 × global/local） | `scripts/multicaller_validation/report.md` 各表 + Stereo 扩展 5 样本（同 Fig3B PENDING 口径） | 新增 |

---

## 附录 A — 作图代码参照仓库表

作图方法不确定时的权威代码来源（硬性要求：读代码而非猜参数）：

| 仓库 | 论文/用途 | 提供的写法 | 本 spec 引用面板 |
|---|---|---|---|
| `tilgnerlab/Spl-IsoFind_reproducibility` | NC 空间 isoform 论文逐图复现 | `figures/Figure {1,3,4,5}.ipynb` + `figures/fig2/plot_*.py`：空间散点、UMI 全景、模拟 PR 面板的 notebook 写法 | Fig2E、Fig3B、Fig7C/D |
| `fiszbein-lab/SpliceImpactR` | NAR（R + Shiny） | 事件类型分布面板 + pheatmap 纪律、事件类型为中心配色 | Fig4C、Fig5D |
| SCSES 仓库（从其论文 Data availability 节提取；不可得时按其图注统计规范复刻） | 最近缘插补方法 | raincloud + UMAP + read-coverage 验证面板 | Fig2B/H、Fig4D |
| `algbio/spl-IsoQuant` | 流程规范 | 流程示意图规范 | Fig1B、Fig7A |

## 附录 B — 参照论文对照表（简称 → 角色）

16 篇（BIB 5 + NAR/GB/GR/NC 11）的完整清单含 PMID/DOI 以 `02_ref_papers/` 参照清单 docx
（已核对）为准；本表登记 spec 内使用的简称与提供的模板/立场，防歧义：

| 简称 | 在本 spec 中的角色 | 提供的模板/立场 |
|---|---|---|
| Spl-ISO-Seq | G4 叙事骨架；Stereo 基因案例与深度稳健性模板 | Fig1 空间三件套、Fig5E、Fig6F（Fig5e 式）、ΔPDUI=ΔΠ 空间版 |
| Longcell | G4 叙事骨架；基因异质性视图 | Fig4G 异质性散点、Fig7D（Fig7A/E 式）、Fig5F 组织域图 |
| TUSCO | G4 叙事骨架；诚实基准叙事 | Fig6E dumbbell、mean 基线公开报告立场 |
| AF2-isoform | G4 叙事骨架 | 总览→审计→基准→统计→分解→案例→延伸主链 |
| STIFT | 每图一问模板；域检测双组指标 | Fig5C/F（Fig2c 式双组条形）、各图“审稿人问题”字段 |
| SpaTM | 遮蔽 CV 与空间保真模板 | Fig2G（Fig3b 式 ΔMoran's I）、Fig2H（Fig2 式 10 折 PCC） |
| SCSES | 最近缘插补方法；统计图注规范 | Style Guide 第 4 条、Fig2B/H raincloud、附录 A 仓库 |
| SpliceImpactR | NAR 侧 R 实现惯例 | 事件类型为中心配色纪律、pheatmap（Fig4C/E、Fig5D）、附录 A |
| Spl-IsoFind | NC 空间 isoform 最近邻方法 | 空间散点/UMI 全景/模拟 PR notebook 写法（Fig2E/3B/7C/D） |
| spl-IsoQuant | 流程示意规范 | Fig1B、Fig7A |
| SCLR-seq 综述 | 总览图分步模板 | Fig1（Fig2a 式分步总览） |
| SCOTCH | 平台基准面板 | S17 全样本 QC 小多联 |
| scAPAtrap | 主 PAS caller（方法参照） | Fig1B 双通道之一、Fig7 工作流 |
| Sierra | 副 PAS caller（caller-agnostic 论证） | Fig1B 双通道之二、S18、Jaccard 口径 |
| stAPAminer | 比较 comparator | Fig2/6 方法集、S5 参数表 |
| spvAPA | 比较 comparator | Fig2/6 方法集、S5 参数表 |

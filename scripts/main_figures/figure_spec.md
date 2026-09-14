# spaGAPA Figure Specification v2 — NC/NAR-calibrated design

**Status**: design spec (2026-09-13, updated 2026-09-14). Supersedes panel
notes in `figure_index.md` as the *design* document; `figure_index.md` remains
the as-built record of the frozen 2026-07-28 figure set. Target journals:
**Nature Communications / Nucleic Acids Research** standard (BIB-compatible).
All requirements below are mandatory unless marked optional.

Calibrated against the 16 reference papers in `02_ref_papers/` (Appendix B;
PMID/DOI checked against `spaGAPA投稿参照文献清单.docx`).

**Strategy**: 演进 + 扩展，不推倒重来。冻结版 7 主图 + 16 补充图的数字不重算；
NC/NAR 升级点以新面板/新样本进入；stereo 扩展产物已落地（人视网膜类器官 ×4 +
大鼠胸腺 ×1，conformal 11→16 样本），剩余未落地项以 `[PENDING: pipeline_output/…
]` 钩子标注并写明落地后填什么。

---

## 0. Style Guide (硬性规范)

1. **字体：全部 Arial。**
   - matplotlib：`rcParams['font.family']='sans-serif'` +
     `rcParams['font.sans-serif']=['Arial']`；数学文本
     `rcParams['mathtext.fontset']='custom'` + `rcParams['mathtext.rm']='Arial'`。
   - R/ggplot2 面板：`theme(text=element_text(family="Arial"))`。
   - 轴标签/刻度缩放后 ≥7 pt；面板标题粗体 8–9 pt；图注面板字母粗体 8 pt。
   - 现有 `_style.py`（DejaVu Sans）在图表按本 spec 重生成时替换为 Arial。
2. **版面**：主图 Nature 双栏宽 **183 mm**（7.2 in）为标准（同时满足 NAR 全宽
   ~178–190 mm）；单栏小图 89 mm（3.5 in）。矢量 PDF 主交付 + 600 dpi
   PNG/TIFF 预览。
3. **配色**：保留 Okabe-Ito 色盲安全主色（`#0072B2` 蓝 / `#D55E00` 橙 /
   `#009E73` 绿；NC/NAR 均要求无障碍）；线宽 ≥0.5 pt，散点描边 0.2 pt。
4. **统计图注规范**（SCSES/NC 式）：每个分布面板图注必须写明 N、检验类型
   （two-sided Wilcoxon/t-test、是否多重校正）、误差定义（SD/SEM/95% CI）、
   箱线定义（median、Q25/Q75、1.5×IQR whisker）；条形误差棒统一
   SEM-of-replicates 并注明 replicate 数。
5. **Source Data**：每主图每面板导出 `SourceData_FigN<panel>.csv`（Nature 投稿
   硬要求）；§主图 Spec 逐面板标注导出数据帧与列。
6. **作图代码参照仓库**：写任何新面板前先查 Appendix A 仓库表
   （`tilgnerlab/Spl-IsoFind_reproducibility`、`fiszbein-lab/SpliceImpactR`、
   `algbio/spl-IsoQuant`、SCSES 仓库）；每个面板行的"参照模板"列注明
   论文+面板（及仓库文件，若存在）。

## 1. 设计原则

- **叙事主链**：总览 → 可靠性审计 → 基准 → 全景统计 → 空间分解 → 基因案例
  → biology payoff → 延伸。（G4 共性骨架：Spl-ISO-Seq / Longcell / TUSCO /
  AF2-isoform 四篇共享此结构。）
- **诚实基准文化**（TUSCO 式基准叙事强化，沿用仓库现有立场）：per-gene mean
  基线公开报告且公开显示其获胜面板；部署实现的经验斜率不裸称 O(N²)；
  边际覆盖 ≠ 条件覆盖显式标注。
- **视觉语言**：Okabe-Ito 色盲安全配色 + 183 mm 版面（沿用冻结版语言，字体
  换 Arial）+ SpliceImpactR 式"事件类型为中心"的分布面板配色纪律。
- **每图回答一个审稿人问题**（STIFT/SpaTM 面板模板）：§主图 Spec 每图第一行
  写明该问题。
- **冻结数字不重算**：2026-07-28 冻结版 7 主图数字原样引用；新样本/新面板
  数字进入时注明来源产物路径。

## 2. 图表总览

| 编号 | 角色 | 回答的审稿人问题 | 核心指标 | 新数据升级点 | 状态 |
|---|---|---|---|---|---|
| 1 | 框架总览 | spaGAPA 是什么、解决哪三个缺口 | —（纯示意） | B 面板 caller-agnostic 双通道 + 多物种/多平台条带 | 重绘（换 Arial） |
| 2 | 插补基准 | 为什么不用 per-gene mean | RMSE / r / spatial fidelity / ΔRMSE + ΔMoran's I / masked-CV PCC | +G/H 两个面板 | A–F 冻结；G/H 新增 |
| 3 | conformal 覆盖 | 真的 95% 吗 | coverage±二项 CI、偏差、区间宽度 | 11→16 样本 + D 平台/物种不变性面板 | B/C 数据扩容；D 新增 |
| 4 | 不确定性 | σ 有信息量吗 | r (pooled/within-gene)、subgroup coverage、risk–coverage | +G 异质性散点 | A–F 冻结；G 新增 |
| 5 | 域检测 | 能恢复空间结构吗 | ARI/NMI、confusion、domain Moran's I | C 重排双组条形 + F 类器官域图 2×2 | A–E 冻结；C 改造；F 新增（数据已落地） |
| 6 | 可扩展性 | 规模上得去吗 | 斜率 / 内存 / completion / Pareto | +E dumbbell + F 下采样稳健性 | A–D 冻结；E/F 新增 |
| 7 | Stereo showcase | 亚细胞分辨率？跨物种？ | PAS 产出、3′ 富集、binning r | 单鼠 pilot → 三物种 showcase | B/C/D 扩容（数据已落地） |
| 8 | biology payoff（新） | spaGAPA 解锁了什么表达分析拿不到的生物学 | ΔPDUI、program×domain heatmap | MOB 冻结 + 类器官 4/4 完成 | 新图（数据已落地） |
| GA | Graphical abstract（新） | 10 秒讲清故事 | —（纯示意） | — | 新图 |
| S17–S22 | 补充图 | 见 §4 | 见 §4 | 见 §4 | S17/S18 数据就绪；S19/S21 已执行；S20 部分；S22 可选 |

## 3. 指标词典

每个指标一段：定义 / 为何它 / 诚实边界。

- **RMSE** — 定义：20% per-gene 遮蔽（seed 42，所有方法遮同一批 entries）下
  held-out usage 的均方根误差。为何：逐点精度的最直接度量，审稿人第一问。
  诚实边界：per-gene mean 是强 RMSE 基线且在逐点 RMSE 获胜（Fig2B 公开显示）；
  RMSE 不度量空间结构。
- **Pearson / Spearman r** — 定义：同遮蔽协议下 held-out 值与预测值的线性/
  秩相关。为何：与领域基准（SpaTM 式 masked-CV）可比。诚实边界：mean 在
  逐点 r 同样不输；r 对恒定预测无定义（mean 的 per-gene 常数使 r 退化）。
- **Spatial fidelity（空间保真）** — 定义：per-gene 空间梯度相关性
  （真实 vs 重建的跨空间 profile 相关）的聚合；GP 0.42 vs mean 0.00。为何：
  mean 按构造插补 per-gene 常数、恢复零空间梯度——这是 spaGAPA 的核心差异
  量。诚实边界：只度量梯度恢复，不度量逐点精度；两指标必须并排报告。
- **ΔMoran's I** — 定义：marker 基因空间自相关（Moran's I）插补前 vs GP
  插补后的变化。为何：SpaTM Fig3b 式——插补是否"制造或抹平"空间结构的最
  直接审计。诚实边界：只在 marker 基因上计算；GP 可能真实放大自相关也可能
  平滑过度，方向须逐基因展示。
- **Masked-CV per-gene PCC** — 定义：SpaTM Fig2 式 10 折遮蔽交叉验证下
  per-gene PCC 分布（median + IQR）。为何：单次 20% 遮蔽的聚合数可能掩盖
  per-gene 差异；分布让"GP 在哪些基因上赢"可见。诚实边界：折间相关（同基因
  多折非独立），IQR 而非 CI 报告。
- **Conformal coverage ± 90% 二项 CI** — 定义：split-conformal 区间的经验
  覆盖率，逐样本逐水平（80/90/95），附二项 CI 带与平均 |偏差| pp。为何：
  有限样本边际保证是方法核心卖点，必须逐样本验证而非只报池化值。诚实边界：
  保证是边际（marginal）而非条件（conditional）；subgroup 偏差（高表达 bin
  ~0.82）作为诚实限制展示。
- **区间宽度 / Winkler 分数** — 定义：2·q̂ 平均宽度；Winkler interval score
  （Gneiting & Raftery 2007，从持久化的 per-observation bounds 精确计算）。
  为何：coverage 可以靠无限宽区间买到，宽度+严格打分才是校准质量。诚实边界：
  早先 q̂+RMSE 近似高估 Winkler ~2.3×，已替换为精确值。
- **Uncertainty–error r（pooled vs within-gene）** — 定义：σ 与 |error| 的
  Pearson r，池化与 within-gene 两个口径。为何：区分"σ 排序全局可信"与
  "σ 在基因内部排序可信"两种用途。诚实边界：within-gene median r 仅
  ≈0.10–0.19；pooled r≈0.5 部分由跨基因异质性驱动，两个数并排报告。
- **Runtime / memory / 幂律斜率** — 定义：合成网格 1k–100k spots 的实测
  wall time、峰值内存与拟合幂律斜率。为何：部署可扩展性的实测证据。诚实
  边界：斜率反映部署实现（近似邻居结构可压平观测斜率），不裸称理论
  复杂度；1200 s wall cap 下 R 工具超限记 timeout/OOM。
- **Completion** — 定义：method × scale 完成矩阵（绿=完成，红=超时/OOM）。
  为何：可比性前提——没跑完的方法不进精度比较。诚实边界：超限是资源上限
  所致，非方法理论不可行。
- **ARI / NMI** — 定义：检测域 vs MOB 5 层解剖标注的聚类一致性。为何：空间
  结构恢复的社区标准。诚实边界：mean-impute 在 MOB 上也恢复相当结构（诚实
  报告）；随机种子稳定性未评估（S10 注明）。
- **ΔPDUI** — 定义：域间 distal PAS usage shift（domain-level proximal vs
  distal usage index 差），Spl-ISO-Seq ΔΠ 的空间版。为何：把 APA 空间变异
  翻译成生物学可读的 3′UTR 长度程序。诚实边界：依赖 caller 的 PAS 分组与
  基因表达量；类器官 RA-vs-BMS 为 n=1 描述性（无 p 值，同 S14 规则）。
- **Jaccard（caller/平台重叠）** — 定义：±50 bp 窗内两 caller/两平台 PAS
  集的 Jaccard，附随机位点背景对照。为何：结论是否依赖特定 caller 的直接
  检验（multicaller 报告口径）。诚实边界：Jaccard 0.05–0.18 看似低，但
  大多数 Sierra peak 距 scAPAtrap 位点 ≤500 bp——重叠口径必须连同距离分布
  一起报告。

## 4. 数据资产映射

冻结部分逐图继承 `figure_index.md` §Data sources（全部在 `pipeline_output/`
下）：`benchmark_mean_transparent/`、`benchmark_stapaminer_headtohead/`、
`conformal_validation/`（含 `per_observation_bounds.npz`）、
`conformal_conditional_coverage/`、`uncertainty_corr_improvement/`、
`uncertainty_within_gene_audit/`、`risk_coverage_curve/`、
`mob_domain_recovery/`、`benchmark_runtime/`、`gse263789_stereo_pilot/`、
`stereo_binning_consistency/`、`gse263789_ad_vs_wt_differential/`。

stereo 扩展新增（全部 `pipeline_output/` 前缀）：

- 5 样本 conformal 汇总：`stereo_expansion_downstream/conformal_expansion_final5.csv`
  （人 ×4 + 大鼠 ×1，global 与 locally_adaptive 双模式，逐样本 n_test/n_genes）
- 逐样本下游：`stereo_expansion_downstream/gse293464_GSM888288{4..7}_s{16,26}_{ra,bms}/`
  与 `stereo_expansion_downstream/gse333693_GSM9770943_thymus/`（各含
  `spagapa_run/`（domains.csv 等）+ `uncertainty/uncertainty_calibration.json`）
- bin200 矩阵与 QC：`gse293464_retina/GSM888288{4..7}_binned*/qc_summary.json`、
  `gse333693_thymus/binned_200/qc_summary.json`（PAS 数、nnz、sparsity）
- 多 caller 口径：`multicaller_validation/report.md`（species/tissue/platform
  分组约定，Fig3D 与 S18 引用）
- 实验 1–3 产物：`simulation_benchmark/`、`polya_db_overlap/`、
  `longread_ortholog/`、`conformal_transfer/`（见 §6 实验设计）

---

## Figure 1 — Framework overview

**定位**：纯示意概念图——三个缺口、输入链、稀疏 GP、split-conformal 四段
叙事的入口。

**审稿人问题**："spaGAPA 是什么，为什么空间 APA 需要专门的概率插补？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | 空间 APA 三缺口：稀疏 usage 矩阵、无校准不确定度、邻居法基准在 ≥42k 位点超限 | 概念示意 | —（纯示意） | — | SCLR 综述 Fig2a；STIFT Fig1 |
| B | 输入链：空间 reads → PAS peak counts → APA usage 矩阵 [gene × spot] → spaGAPA；新增 caller-agnostic 双通道（scAPAtrap + Sierra）与多物种/多平台条带（mouse/human/rat × Visium/Stereo-seq） | 流程示意 | —（纯示意） | `multicaller_validation/report.md`（双通道依据） | SCLR 综述 Fig2a 分步总览；STIFT Fig1 |
| C | 稀疏 GP：posterior mean + 68/95% CI，M inducing points，O(N·M²) 标注 | 模型示意 | —（纯示意） | — | spl-IsoQuant 流程示意规范（algbio/spl-IsoQuant） |
| D | Split-conformal：train/calibrate/test 划分，sᵢ=\|yᵢ−ŷᵢ\|，q̂ → ŷ±q̂ 有限样本边际覆盖 | 校准示意 | —（纯示意） | — | spl-IsoQuant 流程示意规范 |

**指标定义引用**：无（纯示意图，无定量面板）。

**诚实边界**：O(N·M²) 是部署实现描述；Fig6 用实测斜率，不在此裸称理论复杂度。

**与现有脚本的关系**：`fig1_overview.py` 重绘（Arial + 183 mm + B 面板扩容）。

**新数据集成钩子**：B 面板条带引用已落地的 5 个 stereo 样本（§4 路径）；
无 data-pending 项。

**Source Data**：无。

## Figure 2 — Spatial vs mean benchmark（+G/H 两面板）

**定位**：与 per-gene mean 及 4 个竞品的头对头遮蔽基准，诚实框架的核心图。

**审稿人问题**："既然 per-gene mean 逐点 RMSE 更强，为什么要用 GP？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | 遮蔽设计：20% per-gene holdout，所有方法遮同一批 entries（seed 42） | 示意 | —（设计说明） | `benchmark_mean_transparent/transparent_comparison.csv` | SpaTM Fig2 遮蔽示意 |
| B | 逐数据集 paired metrics：RMSE（左）、Pearson/Spearman（右）；mean 逐点获胜公开显示 | paired dot | RMSE；r | `benchmark_mean_transparent/`、`benchmark_stapaminer_headtohead/{gse183456,gse220442}/results.json` | SpaTM/TUSCO 基准面板 |
| C | 空间保真：GP 0.42 vs mean 0.00（mean 按构造插常数） | 条形 | spatial fidelity | 同 B | — |
| D | 精度–空间 2D 散点（每方法每数据集一点） | 散点 | RMSE × fidelity | 同 B | — |
| E | 代表性空间重建：held-out 真值 / mean / 真实 GP posterior（peak_94938）/ 误差图 | 空间地图 ×4 | 逐点误差 | `_data/fig2E_gse183456_peak_94938_gp_posterior.csv` | Spl-ISO-Seq Fig1B–D 空间三件套 |
| F | ΔRMSE 按 Moran's I 五分位分层（GP−mean；Δ>0 = mean 更好；差距不随空间信号缩小） | 箱线 | ΔRMSE | `_cache/s4_stratification.csv` | S4 镜像 |
| G（新） | **ΔMoran's I 插补前后**：marker 基因 APA 空间自相关 raw vs GP 后，paired dots + 箱线 | paired dot + box | ΔMoran's I | MOB + GSE183456 marker 基因（重算，新产物入 `pipeline_output/` 后引用） | SpaTM Fig3b |
| H（新） | **遮蔽 CV 逐基因 PCC**：SpaTM Fig2 式 10 折 per-gene PCC 分布，GP/KNN/mean 三方法 | 小提琴/箱线 | per-gene PCC median + IQR | GSE183456 + GSE220442（重算） | SpaTM Fig2 |

**指标定义引用**：RMSE、Pearson/Spearman、spatial fidelity、ΔMoran's I、
masked-CV per-gene PCC（§3 词典）。

**诚实边界**：mean 在逐点 RMSE/Pearson/Spearman 获胜（B 面板公开显示）；
fidelity=0.00 是 mean 的构造性质而非 bug；F 面板 GP 仅在 ~15% 基因上胜 mean
且差距不随空间信号缩小——三处诚实显示全部保留。G 只在 marker 基因上计算；
H 折间非独立，报 IQR 不报 CI。

**与现有脚本的关系**：`fig2_benchmark.py` 扩展（+G/H；A–F 冻结数字沿用）。

**新数据集成钩子**：G/H 需要新计算产物（ΔMoran's I 表、10 折 PCC 表），
落地路径 `pipeline_output/`（`mob_markers_morans_i.csv`、
`masked_cv_per_gene_pcc.csv`，命名遵循现有惯例）；面板设计先冻结。

**Source Data**：`SourceData_Fig2{A..H}.csv` 逐面板导出。

## Figure 3 — Conformal marginal coverage（11→16 样本）

**定位**：有限样本边际覆盖保证的逐样本实证 + 平台/物种不变性。

**审稿人问题**："声称的 80/90/95% 覆盖真的成立吗？跨平台跨物种还成立吗？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | Split-conformal 流程：calibrate → \|yᵢ−ŷᵢ\| → q̂ → ŷ±q̂ → coverage ≥ 1−α | 示意 | —（流程） | 冻结版 | spl-IsoQuant 流程示意规范 |
| B | 逐样本 coverage 80/90/95 点图 + nominal 线 + 90% 二项 CI 带（16 样本） | 点图 + CI 带 | coverage ± 二项 CI | `conformal_validation/all_samples_coverage.csv` + `stereo_expansion_downstream/conformal_expansion_final5.csv` | TUSCO 基准覆盖面板 |
| C | 校准曲线 nominal vs empirical + y=x 参考线；冻结 11 样本平均偏差 0.21/0.16/0.10 pp；新 5 样本单独引用（C2 locally-adaptive：**0.7999/0.8993/0.9498**） | 校准曲线 | 平均 \|偏差\| pp | 同 B | — |
| D（新） | **平台/物种不变性**：coverage 按 Visium-human / Visium-mouse / Stereo-human / Stereo-rat 分组箱线（分组口径同 `multicaller_validation/report.md`） | 分组箱线 | 分组 coverage ± CI | 同 B + multicaller 口径 | multi-caller report 分组方式 |
| E | 区间宽度：Constant (A) vs Residual-spot (D) 两模式对比 | 箱线/条形 | 2·q̂ 平均宽度 | `conformal_validation/`（winkler_* 列） | — |
| F | 空间实例（冻结 GSE183456 peak_20919）：posterior mean μ / σ / \|error\| / covered-vs-uncovered @ q̂₉₀=0.284 | 空间地图 | 逐点覆盖 | 冻结版 | Spl-ISO-Seq Fig1B–D |

**指标定义引用**：conformal coverage ± 二项 CI、区间宽度/Winkler、Jaccard
（D 面板分组口径，§3 词典）。

**诚实边界**：边际 ≠ 条件——subgroup 偏差在 Fig4E/S6D 展示；80% 水平冻结
偏差 0.21 pp（"within 0.2%" 不严格成立，原文措辞保留）；新 5 样本与冻结
11 样本数字分开引用不混合平均。Stereo-mouse（GSE263789 pilot）尚无 conformal
产物，D 面板暂为 4 组中 3 组 + Stereo-mouse 预留位 [PENDING: pipeline_output/gse263789_stereo_pilot/ — pilot 补跑 conformal 后加入第 4 组]。

**与现有脚本的关系**：`fig3_conformal.py` 扩展（B/C 数据换成 16 样本；+D；
E/F 沿用）。

**新数据集成钩子**：16 样本数据已落地（§4 路径）；唯一 pending 为
Stereo-mouse 组（见上）。

**Source Data**：`SourceData_Fig3{A..F}.csv`（A 为流程示意无导出）。

## Figure 4 — Uncertainty quality（+G 一面板）

**定位**：σ 是否有信息量的多口径审计 + 风险覆盖收益。

**审稿人问题**："你的不确定度除了好看，能用来做决策吗？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | 四种不确定度模型示意（A constant / B local-gene / C spatial-spot / D residual-spot） | 示意 | —（纯示意） | 冻结版 | — |
| B | 多目标方法选择：x=unc–error r，y=区间宽度，bubble=\|coverage 偏差\| | 气泡散点 | 三目标 | `uncertainty_corr_improvement/per_dataset_corr.csv` | — |
| C | 池化 r 逐数据集：A vs D（A 近零，D 过 r≥0.3） | 条形/点图 | pooled r | 同 B | — |
| D | per-gene vs pooled r：within-gene median ≈0.10–0.19 vs pooled ≈0.5（跨基因异质性驱动，透明报告） | 分布对比 | r 两口径 | `uncertainty_within_gene_audit/per_gene_corr_distribution.csv` | — |
| E | 经验 subgroup coverage：组织 × 不确定度五分位热图 + 表达 bin（高表达 bin ~0.82 欠覆盖，诚实限制） | 热图 + 条形 | subgroup coverage | `conformal_conditional_coverage/*.csv` | — |
| F | 风险覆盖曲线（RMSE）：保留率 vs RMSE；不确定度 triage −23% RMSE @80% 保留（5 数据集配对，p=0.004） | 折线 | RMSE vs retention | `risk_coverage_curve/{risk_coverage_data.csv,summary.json}` | — |
| G（新） | **空间 APA 异质性散点**（Longcell φ-vs-ψ 模板）：x=基因级 distal-usage 均值，y=空间离散度 SD，颜色=GP 平均不确定度——"哪些基因的 APA 空间变异性最值得信任" | 彩色散点 | per-gene usage mean vs spatial SD vs median uncertainty | GSE183456 + MOB（GP posterior 已有） | Longcell Fig7A/E + φ-vs-ψ |

**指标定义引用**：uncertainty–error r（pooled vs within-gene）、区间宽度、
conformal coverage（§3 词典）。

**诚实边界**：within-gene r≈0.10–0.19 与 pooled ≈0.5 并排；subgroup 欠覆盖
（~0.82）作为边际保证的诚实限制展示；F 的 p=0.004 为 5 数据集配对检验。

**与现有脚本的关系**：`fig4_noise.py` 扩展（+G；A–F 冻结）。

**新数据集成钩子**：G 面板数据从现有 GP posterior 派生（无新实验）；
派生表落地 `pipeline_output/` 后引用。

**Source Data**：`SourceData_Fig4{B..G}.csv`。

## Figure 5 — Domain recovery（改造 C + 新增 F）

**定位**：空间结构恢复的公平比较 + 真实生物学梯度 + 类器官延伸。

**审稿人问题**："APA 信号（而非表达信号）真的恢复了空间结构吗？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | MOB 解剖：5 层 ground truth（GCL/GL/MCL/ONL/OPL），260 spots | 空间地图 | —（标注） | 冻结版 | Spl-ISO-Seq Fig1B–D |
| B | 域图公平比较：同一 Leiden 管线，仅 APA 来源不同（mean-imputed vs GP-imputed，双真实运行） | 空间地图 ×2 | 视觉对照 | `mob_domain_recovery/spagapa_domains.csv` | STIFT Fig2c |
| C（改造） | **STIFT 式双组指标条形**：结构恢复组（ARI/NMI）+ 空间连贯组（domain Moran's I），各 weight config 并排（mean_apa / apa_dominant 展示 / spatial_apa / balanced / expression_apa 全局最优→S10） | 双组条形 | ARI/NMI + mean domain Moran's I | `mob_domain_recovery/spagapa_metrics.json` | STIFT Fig2c 双组条形 |
| D | 域–层 confusion matrix（展示 config 的域 vs 5 解剖层真实对应） | 热图 | 混淆矩阵 | 同 C | — |
| E | 代表性真实梯度基因：Trim16 / Gdap2 / Tor3a / Pcdhb21 稀疏输入（上）vs GP 重建（下） | 空间地图 ×8 | 逐基因颜色尺度 | 冻结版（`_run_fig5E_mob_gp.py` 产物） | Spl-ISO-Seq Fig1B–D |
| F（新） | **人视网膜类器官域图 2×2**（RA± × 16/26wk 四样本域着色 + 不确定度叠加） | 空间地图 2×2 | 域数 + 不确定度叠加 | `stereo_expansion_downstream/gse293464_GSM888288{4..7}_s*/spagapa_run/domains.csv`（4/4 已落地） | STIFT 空间域图 + SCOTCH 平台面板 |

**指标定义引用**：ARI/NMI、ΔMoran's I（domain Moran's I 口径，§3 词典）。

**诚实边界**：mean-impute 在 MOB 恢复相当结构（C 面板诚实显示，expression_apa
全局最优 0.597 vs apa_dominant 展示）；随机种子稳定性未评估（S10 注明）；
F 为 n=1 每条件的描述性域图。

**与现有脚本的关系**：`fig5_domain.py` 扩展（C 面板重排双组条形；+F；
A/B/D/E 冻结）。

**新数据集成钩子**：F 数据已全部落地（4 样本 domains.csv，§4 路径）；
无 data-pending 项。

**Source Data**：`SourceData_Fig5{C..F}.csv`（A/B 空间地图导出坐标+标签）。

## Figure 6 — Scalability（+E/F 两面板）

**定位**：部署可扩展性的实测证据 + 稳健性。

**审稿人问题**："在真实规模（10 万 spot）上跑得动吗？结果对测序深度稳吗？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | runtime 标度 + 实测幂律斜率（spaGAPA-fast 0.84 / accuracy 1.01 / stAPAminer 1.01 / spvAPA 0.49） | log-log 折线 | 斜率 ± 拟合 | `benchmark_runtime/runtime_table.csv` | TUSCO Fig4 |
| B | 峰值内存标度（spaGAPA-fast 8.8 GB @100k；R 工具超限） | log-log 折线 | 峰值内存 | 同 A | TUSCO Fig4 |
| C | completion 矩阵（method × scale，绿=完成 红=超时/OOM） | 矩阵热图 | completion | 同 A | — |
| D | 空间保真–runtime Pareto（Spatial-KNN 0.91 最高但不可扩展；GP 中等保真+不确定度+可扩展，6–7× 快于 stAPAminer/spvAPA；mean 最快 ~0 保真） | Pareto 散点 | fidelity × runtime | 同 A | — |
| E（新） | **TUSCO 式 dumbbell**：1k vs 100k 两端点，runtime + memory 双指标，4 方法 | dumbbell | Δruntime / Δmemory | 同 A（两端点行） | TUSCO Fig4 |
| F（新） | **下采样稳健性**：100 次重采样 RMSE 分布 ×3 深度档（仅模拟器种子变化；回答"结果是否深度依赖"） | 箱线 ×3 | RMSE 分布 | `simulation_benchmark/` 下采样子集（Exp1 产物复用） | Spl-ISO-Seq Fig5e |

**指标定义引用**：runtime/memory/幂律斜率、completion、RMSE（§3 词典）。

**诚实边界**：斜率是部署实现的经验值（近似邻居结构可压平观测标度），不裸称
O(N²)；R 工具超限是 1200 s wall cap 下的实测；F 用模拟数据（假设见 §6
实验 1）。

**与现有脚本的关系**：`fig6_scalability.py` 扩展（+E/F；A–D 冻结）。

**新数据集成钩子**：F 复用实验 1 模拟产物（100 重采样 ×3 深度子网格）；
落地路径 `pipeline_output/simulation_benchmark/`。

**Source Data**：`SourceData_Fig6{A..F}.csv`。

## Figure 7 — Stereo-seq 三物种 showcase

**定位**：从单鼠 pilot 升级为 mouse/human/rat 三物种 Stereo-seq 展示。

**审稿人问题**："方法在亚细胞分辨率平台上、跨物种还工作吗？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | 工作流：FASTQ+mask → SAW 8.2.2 → retag → scAPAtrap → PAS 矩阵 → spaGAPA；mask 获取双途径标注（GEO/STOmics） | 流程示意 | —（纯示意） | 冻结版 | spl-IsoQuant 流程示意规范 |
| B | 3′ QC：mouse 冻结（TES-distance 分布 + 累积 3′ 富集 53.4%@500bp / 73.5%@2kb）+ human/rat 并排小面板（各自 TES-distance，从 scapatrap_raw 产物计算） | 直方图 + 累积曲线 ×3 | 3′ 富集比例 | 冻结版 + `gse293464_retina/*_scapatrap_raw*/`、`gse333693_thymus/GSM9770943_scapatrap_raw/` | S12 口径 |
| C | 规模与稀疏度三联：PAS 数 mouse **21,455** / human **22,762**（C2；四样本范围 18,368–32,305）/ rat **23,138** + observed fraction 对比条 | 条形三联 | PAS 数；nnz fraction | `gse263789_stereo_pilot/`（冻结）+ `qc_summary.json` ×5（§4） | SCOTCH 平台规模面板 |
| D | 空间三件套 ×3 物种：UMI 全景 → 域图 → 基因案例；mouse 冻结 Cdk8 (peak_69378) / Apoe (peak_312125) / Gnb1l (peak_415285)；human/rat 各选 1–2 个高空间变异 PAS（从 5 样本 spagapa_run posterior 选） | 空间地图 3×3 | 逐基因空间变异（Moran's I） | 冻结版 + `stereo_expansion_downstream/*/spagapa_run/`（human/rat 基因案例待选型 [PENDING: pipeline_output/stereo_expansion_downstream/ — 选型脚本产出基因清单后填入]） | Spl-IsoFind Fig1B–D 空间三件套；Longcell Fig7A/E |
| E | binning 稳健性：跨 bin Pearson r（50→100、50→200）+ Moran's I vs bin size（mouse 冻结）+ 新增跨物种 PAS 数/深度对齐小面板 | 折线 + 对齐条形 | 跨 bin r；Moran's I | `stereo_binning_consistency/binning_correlation.csv` + `qc_summary.json` ×5 | — |

**指标定义引用**：Pearson/Spearman、ΔMoran's I、Jaccard（跨物种对齐口径，
§3 词典）。

**诚实边界**：mouse 数字冻结引用；human/rat PAS 数逐样本列出（22,762 是
C2 样本值，四样本范围注明）；早期 41 域 + raw-σ 不确定度面板因 bin-grid ×
inducing-point 伪影已删（figure_index.md 记录）；AD vs WT 为 n=1 描述性
（S14）。

**与现有脚本的关系**：`fig7_stereo_seq.py` 扩展（B/C/D 扩容 + E 加对齐
小面板）。

**新数据集成钩子**：5 样本矩阵与下游已落地（§4）；唯 D 面板 human/rat
基因案例选型待选型脚本（见面板表 PENDING 标注）。

**Source Data**：`SourceData_Fig7{B..E}.csv`。

## Figure 8 — Spatial APA programs（biology payoff，新图）

**定位**：方法学的生物学兑现——空间 APA 程序与层/域关联。

**审稿人问题**："spaGAPA 解锁了什么表达分析拿不到的生物学？"

**面板表**：

| 面板 | 内容 | 图型 | 指标 | 数据源 | 参照模板 |
|---|---|---|---|---|---|
| A | MOB 层程序：5 层 × program heatmap（逐层 distal usage）+ 层 marker PAS 列表 | heatmap + 注释 | per-layer distal usage | 冻结 MOB 数据（`mob_domain_recovery/`） | SpliceImpactR 事件类型 heatmap 纪律 |
| B | ΔPDUI 火山：外层（ONL/GL）vs 内层（MCL/GCL）层间 distal usage shift；标注已知基因；FDR 沿用现有差异机器 | 火山 | ΔPDUI + FDR | 同 A | Spl-ISO-Seq ΔΠ 面板 |
| C | 程序梯度基因：2–3 个 PAS 空间图 + 层析 usage profile 曲线（冻结 Fig5E 基因延伸） | 空间地图 + 折线 | 层析 usage | 同 A | Spl-ISO-Seq Fig1B–D |
| D | 类器官域特异程序 4/4 样本（B4/C2/D2/D4 各 n=5 域）：域 × program heatmap + RA-vs-BMS 效应量叠加 | heatmap 2×2 + 叠加 | ΔPDUI（描述性） | `stereo_expansion_downstream/gse293464_GSM888288{4..7}_s*/spagapa_run/domains.csv`（全部落地） | SpliceImpactR heatmap 纪律 |

**指标定义引用**：ΔPDUI、ARI/NMI（§3 词典）。

**诚实边界**：D 面板 RA-vs-BMS 每条件 n=1——只报效应量，无 p 值无 FDR
（同 S14 规则）；A/B 用冻结 MOB 数据不重算。

**与现有脚本的关系**：新增 `fig8_biology.py`（`_style.py` v2 Arial）。

**新数据集成钩子**：类器官 4/4 域数据已落地（§4）；无 data-pending 项。

**Source Data**：`SourceData_Fig8{A..D}.csv`。

## Graphical abstract（新，NAR tools 惯例）

单面板 183 mm：稀疏 APA 矩阵（dropout 网格）→ GP posterior 地图（含区间）→
空间域图 → 物种/平台图标。纯矢量。脚本：新增 `fig0_graphical_abstract.py`。
参照：SpliceImpactR graphical abstract 先例。

## 4. 补充图（S1–S22）

S1–S16 沿用 `supp_figure_legends.md`（标题/面板/数据源以该文件为准），
仅重排归属引用；新增 S17–S22：

| 编号 | 标题 | 指标 | 数据源 | 状态 |
|---|---|---|---|---|
| S1 | Dataset overview（32 样本 Visium 基准总览） | PAS 数/样本 | `supp_figure_legends.md` S1 | 冻结 |
| S2 | 稀疏 GP inducing-point 敏感性 | RMSE vs m | 同 S2 | 冻结 |
| S3 | 遮蔽比例敏感性（10–50%） | median per-gene RMSE | 同 S3 | 冻结 |
| S4 | mean 基线分层（Fig2F 镜像） | ΔRMSE 五分位 | 同 S4 | 冻结 |
| S5 | 基准参数表（转补充表） | 参数表 | 同 S5 | 冻结 |
| S6 | 全样本 conformal 覆盖 11 样本 | coverage/宽度/Winkler | 同 S6 | 冻结 |
| S7 | 泄漏审计 + LOOCV 稳定性 | r 差/选择稳定 | 同 S7 | 冻结 |
| S8 | per-gene 不确定度–误差相关 | per-gene r 分布 | 同 S8 | 冻结 |
| S9 | 空间块 vs 随机划分 conformal | coverage 差 | 同 S9 | 冻结 |
| S10 | MOB 权重配置稳健性 | ARI/NMI 扫描 | 同 S10 | 冻结 |
| S11 | 伪重复分析（GSE220442） | 显著基因数对比 | 同 S11 | 冻结 |
| S12 | Stereo QC + binning 敏感性 | 染色体分布/TES 距离 | 同 S12 | 冻结 |
| S13 | 批次校正 QN vs Harmony | PCC vs 批次信号 | 同 S13 | 冻结 |
| S14 | AD vs WT 效应量（n=1 描述性） | 效应量 | 同 S14 | 冻结 |
| S15 | conformal 偏差 forest（主图 Fig3 移出，保持补充） | 偏差 + CI | 同 S15 | 冻结 |
| S16 | MAE risk–coverage（Fig4F 配对单元） | MAE vs retention | 同 S16 | 冻结 |
| S17（新） | Stereo 扩展全样本 QC 总览：每样本 PAS 数/观察比例/3′ 富集/域数小多联 | PAS 数/nnz/域数 | `qc_summary.json` ×5 + `stereo_expansion_downstream/*/spagapa_run/domains.csv` | 数据就绪，待绘图 |
| S18（新） | 跨 caller × 跨物种 conformal 全表：multicaller 3 数据集 + 新 5 样本合并 coverage 表 | coverage 80/90/95 | `multicaller_validation/report.md` + `stereo_expansion_downstream/conformal_expansion_final5.csv` | 数据就绪，待绘图 |
| S19（新） | 模拟基准（实验 1，A–D 四面板） | 见 §6 实验 1 | `simulation_benchmark/` | 已执行子网格（54 条件）；nPAS/smoothness 轴待补 |
| S20（新） | 长读长正交验证（实验 2） | 见 §6 实验 2 | `polya_db_overlap/` + `longread_ortholog/` | E2a 完成；E2b 待产出 |
| S21（新） | 跨样本 conformal 迁移（实验 3） | 迁移衰减 Δpp | `conformal_transfer/` | 已执行（11 样本） |
| S22（可选） | 通用插补基线（实验 4） | RMSE 同协议 | — | 未启动（仅审稿压力时） |

## 5. Deposition 与合规（投稿前清单）

- GitHub release tag (v1.0) + `docs/` RTD 站点上线（已就绪）。
- Zenodo DOI：代码归档 + 关键中间产物（binned matrices、per-observation
  bounds、benchmark tables；从 `pipeline_output/` 筛选 ~50 GB → 保留
  csv/json 层，弃 BAM）。
- GEO/processed-matrices：binned apa_matrix/coordinates 每样本一份 + README。
- Nature Reporting Summary；NAR 数据可用性声明模板。
- Statistics：全图 two-sided 检验 + 多重校正方法注明（BH）。

## 6. 实验设计：模拟基准与正交验证

### 实验 1 — 参数化模拟自测（→ S19）

- **扩展对象**：`spagapa/benchmark/simulator.py` +
  `scripts/synthetic_apa_generator.py`。
- **设计网格**：测序深度（5 档）× 观察比例（5 档）× PAS 数（3 档）× 空间
  平滑度（3 档）× 噪声模型（4 种，对齐 Fig4A：constant / local-gene /
  spatial-spot / residual-spot）。
- **指标**：RMSE / Pearson / conformal coverage / spatial fidelity / PAS
  检出 precision–recall（Spl-IsoQuant Fig3C–E 式 P–R 曲线）。
- **面板**：S19A 性能–深度功效曲线；S19B 性能–稀疏度曲面；S19C P–R 曲线族；
  S19D 掩蔽"伪新 PAS"恢复率（TUSCO-novel 式）。
- **执行现状（诚实标注）**：已跑子网格 54 条件（depth {5000, 20000, 80000} ×
  observed {0.05, 0.15, 0.35} × noise {constant, residual_spot} × 3 seeds，
  n_genes=50，0 失败；聚合例：depth 5000 → RMSE 0.0703、Pearson 0.9336、
  conformal-90 global 0.8997）。nPAS 与 smoothness 轴及 4 噪声全档待补，
  产物 `pipeline_output/simulation_benchmark/`。
- **诚实边界（模拟器假设显式列出）**：模拟图案族 {autocorrelated 0.3 /
  domain 0.2 / gradient 0.1 / constant 0.4}；GP 核族与模拟核族匹配风险声明
  （RBF GP 拟合 RBF-ish 真值可能高估真实表现——S19 图注必须写明）。

### 实验 2 — 长读长正交验证（→ S20 / Fig7D 选型支持）

- **2a（必做，已完成）PolyA_DB v4 重叠**：PolyA_DB v4.1（PMID 41316728；
  hg38 1,429,829 PAS / mm10 1,346,135 PAS，sha256 记录于
  `polya_db_overlap/overlap_summary.json`）对 18 样本（10 个 GSE）scAPAtrap
  PAS：±50 bp 命中率 **74.7%**（背景 2.5%，富集 **29.7×**），距离中位数
  **10 bp**，88.8% ≤500 bp；分物种 human 73.0% / mouse 83.0%（mm39→mm10
  liftover 修正后）。面板：重叠距离分布 + 命中率条形（含随机位点背景对照）。
- **2b（必做，进行中）Spl-IsoFind 空间长读正交**：从其 Zenodo
  （DOI **10.5281/zenodo.19499423**，reproducibility 仓库 README 记录）取
  人皮层空间 ONT per-gene PAS 使用估计；对重叠基因计算长读 vs spaGAPA
  （Visium 人脑样本）每基因 distal-usage 相关性（SCSES Fig2d 式）。面板：
  每基因散点 + 中位 r + N。参考数据与 ortholog map 已就绪 [PENDING: pipeline_output/longread_ortholog/ — concordance 结果表落地后填中位 r 与 N]。
- **2c（条件做）Longcell MOB ONT**：仅当 2a/2b 相关性中位数 r<0.3 时启用，
  否则不做。
- **对齐口径**：基因级 usage（不做位点级对齐，避免 caller 差异混淆——与
  multicaller 报告同口径）。

### 实验 3 — 跨样本校准迁移（→ S21，已执行）

Leave-one-sample-out：在样本 A 上校准 conformal、在样本 B 上测试，11 样本
全对。结果：80/90/95 水平 within-mean 0.7955/0.8964/0.9496 vs cross-mean
0.7936/0.8932/0.9467；平均迁移衰减 0.19/0.33 pp，最大 11.8 pp（110 跨对）。
诚实边界：跨样本迁移破坏可交换性假设——衰减分布（非均值）才是主信息；
图注必须讨论。数据：`conformal_transfer/`（复用
`conformal_validation/per_observation_bounds.npz`，无新 GP 运行）。

### 实验 4（可选）— 通用插补基线（→ S22）

MAGIC 式扩散插补作用到 APA 矩阵，同遮蔽协议比较；图注说明 usage/AS 语义
错配风险。仅当预期审稿压力时执行。

## Appendix A — 作图代码参照仓库表

写任何新面板前先读本表对应仓库的实现（Step 0.6 硬性规范）：

| 仓库 | 用途 | 关键文件 |
|---|---|---|
| `tilgnerlab/Spl-IsoFind_reproducibility` | NC 空间 isoform 论文逐图复现：空间散点全景、UMI panorama、模拟 P–R 面板的 notebook 写法 | `figures/Figure {1,3,4,5}.ipynb`；`figures/fig2/plot_*.py` |
| `fiszbein-lab/SpliceImpactR` | NAR R+Shiny：事件类型分布面板语法 + pheatmap 纪律 | R 包源码 + Shiny app |
| `algbio/spl-IsoQuant` | 流程示意图规范（Fig1A/Fig7A） | docs/ 示意图 |
| SCSES 仓库（按论文 Data availability 节提取；fallback 按图注复刻） | raincloud、UMAP、read-coverage 验证面板 | — |

本地镜像：`pipeline_output/longread_ortholog/data/Spl-IsoFind_reproducibility-main/`
（实验 2b 已下载，含 SpatialAnalysis.md / Simulation.md 逐图文档）。

## Appendix B — 参照文献对照表（PMID/DOI 已核对）

清单源：`02_ref_papers/spaGAPA投稿参照文献清单.docx`（16 篇 PDF + 2 篇
仅列目录）。

BIB 5 篇：

| # | 论文 | PMID | DOI |
|---|---|---|---|
| 1 | spvAPA（直接竞品：监督式 sc+空间 APA） | 39799000 | 10.1093/bib/bbae720 |
| 2 | APAdeg（APA-seq 差异基因 GLMM） | 42242679 | 10.1093/bib/bbag295 |
| 3 | STIFT（时空转录组整合；面板模板） | 41370630 | 10.1093/bib/bbaf644 |
| 4 | SpaTM（空间主题模型；benchmark 文化+Fig2/Fig3b 模板） | 41359801 | 10.1093/bib/bbaf657 |
| 5 | 单细胞长读长 AS 生信框架（综述；Fig1 分步总览模板） | 41378880 | 10.1093/bib/bbaf655 |

补充竞品（2024-2026 新检索，NAR 2025-09-14）：

| # | 论文 | PMID | DOI | 关系 |
|---|---|---|---|---|
| 19 | stAI（NAR 2025 53(5):gkaf158——深度学习空间转录组缺失基因插补+注释） | 40057378 | 10.1093/nar/gkaf158 | NAR 已发空间插补方法（证明 venue 接收度）；但为表达量插补、无 APA、无覆盖保证——Introduction 引用作差异化锚点 |
| 20 | scASprofiler（BIB 2026 bbag497——深度卷积生成网络恢复 scRNA 剪接 junction counts） | — | 10.1093/bib/bbag497 | 领域内"稀疏恢复"最近缘工作；深度生成网络 vs 我们的 GP+conformal，须在 Related Work 对比 |

NAR / GB / GR / NC 13 篇：

| # | 论文 | PMID | DOI |
|---|---|---|---|
| 6 | SpliceImpactR（NAR；事件类型面板语法） | 42605803 | 10.1093/nar/gkag827 |
| 7 | PolyA_DB v4（NAR；实验 2a 目录源） | 41316728 | 10.1093/nar/gkaf1212 |
| 8 | 3′tag scRNA-seq PAS 方法基准（NAR） | 42120044 | 10.1093/nar/gkag490 |
| 9 | PolyAseqTrap（GB） | 41620776 | 10.1186/s13059-026-03963-w |
| 10 | IFDlong（GB） | 41851882 | 10.1186/s13059-026-04023-z |
| 11 | ScIsoX（GB） | 40983941 | 10.1186/s13059-025-03758-5 |
| 12 | 预测人类 AS 蛋白结构影响（GB） | 40963109 | 10.1186/s13059-025-03744-x |
| 13 | Biosurfer（GR） | 40086882 | 10.1101/gr.279317.124 |
| 14 | SCOTCH（NC；平台基准面板） | 42098110 | 10.1038/s41467-026-72665-5 |
| 15 | SCSES（NC；最近缘插补方法；图注纪律） | 41145486 | 10.1038/s41467-025-64517-5 |
| 16 | TUSCO（NC；基准叙事+dumbbell） | 42026080 | 10.1038/s41467-026-72089-1 |
| 17 | 单细胞+空间剪接 Nanopore（NC；Longcell 系） | 40683866 | 10.1038/s41467-025-60902-2 |
| 18 | 空间长读长近单细胞分辨率（NC；Spl-IsoFind，实验 2b） | 40883294 | 10.1038/s41467-025-63301-9 |

## Appendix C — 参照面板致谢

SpaTM Fig3b（ΔMoran's I）· SpaTM Fig2（masked-CV PCC）· STIFT Fig2c
（双组条形）· Spl-ISO-Seq Fig1B–D & Fig5e（空间三件套、下采样）·
Longcell Fig7A/E + φ-vs-ψ（异质性散点）· TUSCO Fig4 + TUSCO-novel
（dumbbell；伪新 PAS）· SCSES Fig2d（正交一致性）+ 图注纪律 ·
SpliceImpactR（事件类型语法；graphical abstract 先例）· SCOTCH（平台
QC 面板）· APAdeg Fig3（模拟面板组装）· SCLR 综述 Fig2a（分步总览）。

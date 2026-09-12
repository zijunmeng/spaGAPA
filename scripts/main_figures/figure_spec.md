# spaGAPA Figure Specification v2 — NC/NAR-calibrated design

**Status**: design spec (2026-09-13). Supersedes panel notes in `figure_index.md`
as the *design* document; `figure_index.md` remains the as-built record.
Target journals: **Nature Communications / Nucleic Acids Research** standard
(BIB-compatible). All requirements below are mandatory unless marked optional.

Calibrated against 16 reference papers in `02_ref_papers/` (Appendix C).

---

## 0. Style Guide (硬性规范)

1. **Font**: **Arial everywhere.**
   - matplotlib: `rcParams['font.family']='sans-serif'; rcParams['font.sans-serif']=['Arial']; rcParams['mathtext.fontset']='custom'; rcParams['mathtext.rm']='Arial'`
   - R/ggplot2 panels: `theme(text=element_text(family="Arial"))`
   - Axis labels/ticks ≥7 pt after scaling; panel titles bold 8–9 pt; bold
     panel letters 8 pt.
   - Existing `_style.py` (DejaVu Sans) is replaced on regeneration.
2. **Page geometry**: main figures 183 mm (7.2 in) double-column width
   (satisfies Nature 183 mm and NAR full width); small companion figures
   89 mm. Vector PDF primary; 600 dpi PNG/TIFF preview.
3. **Color**: Okabe-Ito colorblind-safe palette retained
   (`#0072B2/#D55E00/#009E73`), line ≥0.5 pt, scatter edge 0.2 pt.
4. **Statistical caption discipline** (SCSES/NC convention — every
   distribution panel caption states): N, test type (two-sided
   Wilcoxon/t-test, adjustment), error definition (SD/SEM/95% CI), boxplot
   definition (median, Q25/Q75, 1.5×IQR whiskers). Bar error bars: SEM
   across stated replicates.
5. **Source Data**: every main-figure panel exports
   `SourceData_FigN<panel>.csv` (Nature submission requirement).
6. **Plot-code reference repos** (consult before writing any new panel):
   `tilgnerlab/Spl-IsoFind_reproducibility` (per-figure NC notebooks,
   spatial panoramas, simulated P–R), `fiszbein-lab/SpliceImpactR` (NAR
   event-type distribution grammar, pheatmap discipline), `algbio/spl-IsoQuant`
   (workflow schematics), SCSES repo per its Data Availability (raincloud,
   read-coverage validation panels; fallback: caption-style reproduction).

## 1. 设计原则

- Narrative spine (Spl-ISO-Seq/Longcell/TUSCO/AF2-isoform): overview →
  reliability audit → benchmark → panorama → spatial decomposition → gene
  cases → **biology payoff** → extension.
- Honest benchmark culture (TUSCO-style): mean baseline openly; deployed
  slopes, no bare O(N²); marginal ≠ conditional flagged.
- Every figure answers one reviewer question (STIFT/SpaTM templates).
- Frozen 2026-07-28 numbers not recomputed; stereo-expansion data enters
  via `[PENDING: path]` hooks.

## 2. 图表总览

| # | Role | Reviewer question | Core metrics | New-data hook | Status |
|---|------|-------------------|--------------|---------------|--------|
| 1 | Overview | What is spaGAPA? | — schematic | 3-species input strip | regenerate |
| 2 | Imputation benchmark | Why not per-gene mean? | RMSE/r/fidelity/ΔRMSE + **ΔMoran's I, masked-CV PCC** | — | +2 panels |
| 3 | Conformal coverage | Really 95%? | coverage±CI, deviation, width | 11→16 samples; invariance panel | extend |
| 4 | Uncertainty | Is σ informative? | r (pooled/gene), subgroup, risk–coverage, heterogeneity | — | +1 panel |
| 5 | Domain recovery | Biology recovered? | ARI/NMI, confusion, domain Moran's I | + organoid 2×2 maps | extend |
| 6 | Scalability | Runs at scale? | slopes/memory/completion/Pareto + dumbbell, downsampling | — | +2 panels |
| 7 | Stereo showcase | Subcellular? Cross-species? | PAS yield, 3′ enrichment, binning r | 3-species | major extend |
| **8** | **Biology payoff (NEW)** | What biology does it unlock? | ΔPDUI programs, program×domain heatmap | MOB frozen + organoid `[PENDING: §7 organoid runs]` | new |
| GA | Graphical abstract (NEW) | 10-second story | — | — | new |
| S19–S22 | see §5–§6 | | | | new |

## 3. 指标词典

RMSE/Pearson/Spearman (20% per-gene mask, seed 42, same entries all
methods) · spatial fidelity (per-gene gradient corr; mean ⇒ 0.00 by
construction; GP 0.42) · ΔMoran's I (marker autocorrelation before/after
imputation) · masked-CV per-gene PCC (10-fold, median+IQR) · conformal
coverage ± 90% binomial CI + mean |deviation| pp · interval width/Winkler ·
unc–error r pooled vs within-gene · ARI/NMI · runtime/memory power-law
slopes + completion · ΔPDUI (domain distal-usage shift; spatial ΔΠ
analogue) · Jaccard ±50 bp (caller/platform overlap).

## 4. 主图 Spec

### Figure 1 — Framework overview
A three gaps · B input pipeline (+caller-agnostic dual channel
scAPAtrap/Sierra; species/platform strip mouse·human·rat ×
Visium·Stereo-seq) · C sparse GP posterior · D split-conformal.
Schematic only (SCLR review Fig2a + STIFT Fig1). Script: `fig1_overview.py`
restyle. No SourceData.

### Figure 2 — Spatial vs mean benchmark (+2 panels)
A–F as built (masking; paired per-dataset metrics — mean wins shown
openly; fidelity 0.42 vs 0.00; accuracy–spatial 2D; representative
reconstruction; ΔRMSE quintiles).
**G — ΔMoran's I before/after** (SpaTM Fig3b): marker APA autocorrelation
raw vs GP; paired dots+boxes. Data: MOB+GSE183456.
**H — masked-CV per-gene PCC** (SpaTM Fig2): 10-fold per-gene PCC
distributions; GP/KNN/mean. Data: GSE183456+GSE220442.
Script: extend `fig2_benchmark.py` (+G,H). SourceData per panel.

### Figure 3 — Conformal coverage (11→16 samples)
A flow · B per-sample dots+binomial CI (16 samples) · C calibration curve
(frozen 0.21/0.16/0.10 pp; new-sample numbers quoted separately —
**C2: 0.7999/0.8993/0.9498**) · **D platform/species invariance**: coverage
by Visium-mouse/Stereo-mouse/Stereo-human/Stereo-rat boxes · E width A vs D
· F spatial instance (frozen peak_20919).
`[PENDING: stereo_expansion_downstream/*_s2*/uncertainty/*.json + rat]`.
Script: extend `fig3_conformal.py` (+D; B/C data swap).

### Figure 4 — Uncertainty quality (+1 panel)
A–F as built (4 noise models; multi-objective bubble; pooled r; per-gene vs
pooled; subgroup coverage incl. ~0.82 high-expression undercoverage;
risk–coverage −23% @80%, p=0.004).
**G — spatial APA heterogeneity scatter** (Longcell φ-vs-ψ template):
x=per-gene mean distal usage, y=spatial SD, color=median GP uncertainty.
Data: GSE183456+MOB. Script: extend `fig4_noise.py`.

### Figure 5 — Domain recovery (extended)
A–E as built (MOB anatomy; fair mean-vs-GP; config sweep; confusion
matrix; gradient genes). **C restyle** — STIFT two-group bars: structure
retention (ARI/NMI) + spatial coherence (domain Moran's I).
**F — organoid domain maps 2×2** (RA±×16/26 wk) + uncertainty overlay.
`[PENDING: pipeline_output/stereo_expansion_downstream/gse293464_GSM888288{4,6,7}_s2*/]`. Script: extend `fig5_domain.py`.

### Figure 6 — Scalability (+2 panels)
A–D as built (slopes 0.84/1.01/1.01/0.49; 8.8 GB @100k; completion;
Pareto). **E — TUSCO dumbbell** (1k vs 100k; runtime+memory; 4 methods).
**F — downsampling robustness** (Spl-ISO-Seq Fig5e): 100× resampling RMSE
boxes ×3 depths; simulator seeds only. Script: extend `fig6_scalability.py`.

### Figure 7 — Stereo-seq three-species showcase
A workflow (+mask dual-path GEO/STOmics) · B 3′-QC (mouse frozen
53.4%@500bp/73.5%@2kb; human/rat flanks `[PENDING: §7 binned dirs]`) · C scale/sparsity
triplet (21,455 / 22,762 / 23,138 PAS; observed-fraction bars) ·
D spatial triplet ×3 species (UMI panorama → domains → gene cases; mouse
Cdk8/Apoe/Gnb1l frozen; human/rat `[PENDING: §7 binned dirs]`) · E binning robustness +
cross-species alignment mini-panel. Script: extend `fig7_stereo_seq.py`.

### Figure 8 — Spatial APA programs (biology payoff; NEW)
**审稿人问题**: "What biology does spaGAPA unlock that expression
analysis cannot?"
- **A — MOB layer programs**: 5-layer × program heatmap (distal-usage per
  layer), layer-marker PAS list; frozen MOB data.
- **B — ΔPDUI volcano** between outer (ONL/GL) vs inner (MCL/GCL) layers;
  label known genes; FDR from existing differential machinery.
- **C — program gradient genes**: 2–3 PAS spatial maps with layer-wise
  usage profile curves (frozen Fig5E genes extended).
- **D — organoid domain-specific programs** `[PENDING: §7 organoid runs]`: domain × program
  heatmap; RA-vs-BMS effect-size overlay (descriptive n=1, no p-values —
  same honesty rule as S14).
- Script: new `fig8_biology.py` (style: `_style.py` v2 Arial).

### Graphical abstract (NEW; NAR tools standard)
Single 183 mm panel: sparse APA matrix (dropout grid) → GP posterior map
with intervals → spatial domain map → species/platform icons. Vector only.
Script: new `fig0_graphical_abstract.py`.

## 5. 补充图

S1–S16 unchanged (`supp_figure_legends.md`). New:
- **S17** stereo-expansion QC overview (per-sample PAS/observed/3′/domains;
  SCOTCH platform-panel style).
- **S18** caller×species conformal table (multicaller 3 + 5 new samples).
- **S19** simulation benchmark (Experiment 1).
- **S20** long-read orthogonal validation (Experiment 2).
- **S21** cross-sample conformal transfer (Experiment 3).
- **S22 (optional)** general-purpose imputation baseline (Experiment 4).

## 6. 实验设计

### 实验 1 — 参数化模拟自测 (→ S19)
Extend `spagapa/benchmark/simulator.py` +
`scripts/synthetic_apa_generator.py`. Grid: depth(5) × observed(5) ×
nPAS(3) × smoothness(3) × noise(4, aligned to Fig4A). Metrics: RMSE,
Pearson, coverage, fidelity, PAS P–R (Spl-IsoQuant Fig3C–E). Panels:
S19A power-vs-depth; S19B vs-sparsity surface; S19C P–R family; S19D
pseudo-novel PAS recovery (TUSCO-novel). Assumptions declared (kernel-family
match risk).

### 实验 2 — 长读长正交验证 (→ S20 / Fig7)
- **2a (mandatory)** PolyA_DB v4 overlap: public catalog; ±50 bp hit rate,
  distance distribution, per-tissue Jaccard vs random-locus background.
- **2b (mandatory)** Spl-IsoFind spatial ONT concordance: Zenodo
  10.5281/zenodo.19499423; per-gene distal-usage correlation long-read vs
  spaGAPA (Visium human brain); SCSES Fig2d style.
- **2c (conditional)** Longcell MOB ONT: ONLY if 2a/2b median r < 0.3.
- Gene-level usage alignment only (multicaller-report convention).

### 实验 3 — 跨样本校准迁移 (→ S21)
Leave-one-sample-out: calibrate conformal on sample A, test on sample B,
across all 16 samples. Metric: coverage transfer decay (Δpp) per level;
compare within- vs cross-sample calibration; per-platform grouping.
Data: existing `per_observation_bounds.npz` + per-sample predictions
(no new GP runs). Honesty: exchangeability violation discussed.

### 实验 4 (optional) — 通用插补基线 (→ S22)
MAGIC-style diffusion imputation on the APA matrix, same mask protocol;
explain usage/AS semantic mismatch in caption. Only if reviewer pressure
anticipates it.

## 7. 数据资产映射

Frozen per `figure_index.md` §Data sources, plus
`stereo_expansion_downstream/*_s2*/uncertainty/*.json` (C2 done),
`gse293464_retina/GSM8882885_binned/` (done),
`GSM888288{4,6,7}_binned/`, `gse333693_thymus/GSM9770943_binned/`
(paths under `pipeline_output/`, pending), `multicaller_validation/report.md` (grouping conventions).

## 8. Deposition 与合规（投稿前清单）

- GitHub release tag (v1.0) + `docs/` RTD 站点上线（已就绪）。
- Zenodo DOI：代码归档 + 关键中间产物（binned matrices、per-observation
  bounds、benchmark tables；从 `pipeline_output/` 筛选 ~50 GB → 保留
  csv/json 层，弃 BAM）。
- GEO/processed-matrices：binned apa_matrix/coordinates 每样本一份 +
  README。
- Nature Reporting Summary；NAR 数据可用性声明模板。
- Statistics：全图 two-sided 检验 + 多重校正方法注明（BH）。

## Appendix C — 参照面板致谢

SpaTM Fig3b (ΔMoran's I) · SpaTM Fig2 (masked-CV PCC) · STIFT Fig2c
(two-group bars) · Spl-ISO-Seq Fig1B–D & Fig5e (spatial triplet,
downsampling) · Longcell Fig7A/E + φ-vs-ψ (heterogeneity scatter) · TUSCO
Fig4 + TUSCO-novel (dumbbell; pseudo-novel) · SCSES Fig2d (orthogonal
concordance) + caption discipline · SpliceImpactR (event-type grammar;
graphical abstract precedent) · SCOTCH (platform QC panel) · APAdeg Fig3
(simulation panel assembly). 16-paper inventory with PMID/DOI:
`02_ref_papers/spaGAPA投稿参照文献清单.docx`.

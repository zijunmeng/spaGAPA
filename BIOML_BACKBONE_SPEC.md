# Milestone 7: CPU-Friendly BioML Backbone

**目标**: 在不使用深度学习、不依赖 GPU、不使用真实 layer label 作为模型输入的前提下，提升 spaGAPA 的 biological consistency，使其在 imputation accuracy、uncertainty calibration 和 domain/layer recovery 上更有机会全面超过既往方法。

**状态**: Draft v1  
**日期**: 2026-06-11

---

## 1. 背景判断

前几轮 benchmark 已经说明：

1. `spatial GP` 在 RMSE 和 uncertainty calibration 上最强。
2. `expr_additive` / `expr_product` / `adaptive_additive` 能引入 expression manifold，但会明显损伤 RMSE 和 calibration。
3. `layer_local` 证明“局部化 expression influence”是正确方向，但仍未超过 `spatial GP`，也未全面超过 stAPAminer-like expression KNN 的 layer ARI。

因此，继续微调 single-gene GP kernel 的收益已经变小。下一阶段需要新增一个 CPU-friendly 的 multi-view machine learning backbone，专门服务 biological consistency 和 domain recovery。

---

## 2. 底线原则

### 2.1 禁止深度学习

本里程碑不使用：

- neural network
- graph neural network
- variational autoencoder
- Transformer
- GPU-only training

### 2.2 禁止 label leakage

真实 MOB layer label 只能用于评价，不能用于模型训练或特征构建。

允许输入：

- spatial coordinates
- coordinate-derived geometry: radius, theta, pseudolayer
- expression matrix / expression PCA
- observed APA matrix
- GP-imputed APA matrix
- GP uncertainty matrix

禁止输入：

- true layer labels
- manually curated layer/domain labels
- marker-defined layer labels as supervised targets

### 2.3 CPU-friendly

核心算法必须可在普通 CPU 上运行：

- sparse graph
- low-rank factorization
- spectral clustering
- graph diffusion
- coordinate descent / alternating least squares
- NMF / PCA / matrix factorization

---

## 3. 总体模型

Milestone 7 不替换现有 GP backbone，而是在其上新增 BioML backbone：

```text
spaGAPA-GP:
  spatial GP / layer_local GP
  -> APA imputation
  -> posterior uncertainty

spaGAPA-BioML:
  multi-view graph
  -> graph-regularized APA representation
  -> biological domain recovery
  -> uncertainty-aware downstream
```

也就是说：

- `spatial GP` 继续作为 primary imputation model
- `BioML` 负责提升 biological consistency / layer recovery / boundary consistency

---

## 4. Multi-View Graph

构建三个 spot-level graph：

```text
W_spatial
W_expr
W_APA
```

### 4.1 Spatial Graph

来源：

- coordinates
- radius / theta
- spatial KNN

形式：

```text
W_spatial(i,j) = exp(-d_spatial(i,j)^2 / sigma_s^2)
```

只保留 KNN sparse edges。

### 4.2 Expression Graph

来源：

- expression PCA
- optional HVG-PCA

形式：

```text
W_expr(i,j) = exp(-d_expr(i,j)^2 / sigma_e^2)
```

只保留 KNN sparse edges。

### 4.3 APA Graph

来源：

- observed APA matrix
- GP-imputed APA matrix
- GP uncertainty

推荐 MVP：

```text
APA_weighted = GP_imputed_APA
confidence = 1 / (uncertainty + eps)
```

构建 uncertainty-weighted spot similarity：

```text
W_APA(i,j) = similarity(APA_i, APA_j; confidence_i, confidence_j)
```

### 4.4 Graph Fusion

统一图：

```text
W = alpha_s * W_spatial + alpha_e * W_expr + alpha_a * W_APA
```

权重选择：

1. 默认手动网格：
   - `alpha_s in {0.3, 0.5, 0.7}`
   - `alpha_e in {0.1, 0.3, 0.5}`
   - `alpha_a in {0.1, 0.3, 0.5}`
2. 约束：
   - `alpha_s + alpha_e + alpha_a = 1`
3. 后续通过 validation mask 和 biological metrics 做 model selection。

---

## 5. Graph-Regularized Matrix Factorization

目标矩阵：

```text
Y: genes x spots APA usage matrix
```

低秩分解：

```text
Y ≈ G @ Z.T
```

其中：

- `G`: gene factors, shape `(n_genes, rank)`
- `Z`: spot factors, shape `(n_spots, rank)`

### 5.1 Objective

```text
min || M * C * (Y - G @ Z.T) ||^2
    + lambda_graph * Tr(Z.T L Z)
    + lambda_g * ||G||^2
    + lambda_z * ||Z||^2
```

其中：

- `M`: observed / reliable mask
- `C`: uncertainty-derived confidence weights
- `L`: fused graph Laplacian
- `Tr(Z.T L Z)`: graph smoothness penalty

### 5.2 直觉

如果两个 spots 在 fused graph 上相近，那么它们的 spot latent factors 应该相近。

这直接针对 biological consistency，因为 layer/domain recovery 评估的是 spot-level latent structure，而不是单 gene 的点预测。

### 5.3 Solver

MVP 使用 CPU-friendly alternating ridge regression：

1. fixed `Z`, update `G`
2. fixed `G`, update `Z` with graph Laplacian regularization
3. repeat until convergence or max_iter

优先实现 sparse operations。

---

## 6. Domain Recovery

从 spot factor `Z` 和 fused graph `W` 得到 biological domains。

候选方法：

1. spectral clustering
2. Leiden / Louvain graph clustering
3. KMeans on graph-diffused factors
4. Gaussian mixture model on `Z`

MVP 推荐：

```text
SpectralClustering(affinity="precomputed")
```

和：

```text
KMeans(Z)
```

并行比较。

---

## 7. Evaluation Metrics

### 7.1 Imputation

- RMSE
- MAE
- Pearson
- Spearman
- R2

### 7.2 Uncertainty

- coverage 68/95
- uncertainty vs error correlation
- uncertainty-filtered RMSE gain

### 7.3 Biological Consistency

- layer ARI
- layer NMI
- boundary consistency
- pseudolayer smoothness
- marker APA ordering

### 7.4 Runtime

- runtime
- peak memory

---

## 8. MVP Tasks

### Phase 7.1: Multi-view graph builder

- Implement `MultiViewGraphBuilder`
- Inputs:
  - coordinates
  - expression embedding
  - APA matrix
  - uncertainty matrix
- Outputs:
  - `W_spatial`
  - `W_expr`
  - `W_APA`
  - `W_fused`
  - graph metadata

### Phase 7.2: Graph-regularized factorization

- Implement `GraphRegularizedAPAFactorizer`
- Support:
  - rank
  - lambda_graph
  - lambda_l2
  - max_iter
  - tolerance
- Outputs:
  - imputed APA matrix
  - gene factors
  - spot factors

### Phase 7.3: Domain recovery

- Implement `BioMLDomainDetector`
- Support:
  - spectral clustering
  - KMeans on spot factors
- Outputs:
  - domain labels
  - factor embedding

### Phase 7.4: Benchmark integration

- Add benchmark variant:

```text
--method spagapa_bioml
```

or:

```text
--gp-variant spatial
--downstream-variant bioml
```

Compare:

- `spatial GP`
- `layer_local GP`
- `stapaminer_knn_expression`
- `spagapa_bioml`

Status:

- [x] `run_stapaminer_mob_benchmark.py` supports `--include-bioml`
- [x] KMeans and spectral BioML domain detectors are supported
- [x] Targeted pilot completed on MOB subset
- [x] Formal benchmark suite completed with multiple masks and seeds
- [x] Graph weight sweep completed

Pilot result summary:

```text
BioML-kmeans:
  rmse = 0.112568
  layer_ari = 0.059991
  uncertainty_error_spearman = 0.710106

BioML-spectral:
  rmse = 0.112568
  layer_ari = 0.377720
  layer_nmi = 0.474067
  uncertainty_error_spearman = 0.710106

stAPAminer-like:
  rmse = 0.125075
  layer_ari = 0.077191
```

Interpretation:

- KMeans on BioML factors improves layer ARI but does not beat stAPAminer-like.
- Spectral clustering on the fused multi-view graph strongly improves layer ARI/NMI while preserving the GP backbone's RMSE/calibration.

Formal suite:

```text
Output:
  spaGAPA/benchmark_results/real/bioml_formal_suite_v1/

Formal balanced:
  n_genes = 60
  seeds = 42,43
  masks = random, spatial_block_large, ring_sector, layer_aware, low_coverage
  graph weights = spatial 0.4, expression 0.4, APA 0.2
```

Formal result summary:

```text
spagapa_bioml:
  rmse = 0.116868
  layer_ari = 0.375996
  layer_nmi = 0.492093
  uncertainty_error_spearman = 0.810176
  runtime_s = 15.612
  peak_rss_mb = 615.939

spagapa_gp:
  rmse = 0.116897
  layer_ari = 0.077308
  layer_nmi = 0.117774
  uncertainty_error_spearman = 0.810283

stAPAminer-like:
  rmse = 0.127412
  layer_ari = 0.101109
  layer_nmi = 0.134307
```

Mask-level highlight:

```text
layer_aware:
  spagapa_bioml layer_ari = 0.444642
  stAPAminer-like layer_ari = 0.230029
  spagapa_gp layer_ari = 0.098998
```

Weight sweep summary:

```text
Best BioML RMSE:
  no_apa_s050_e050_a000
  rmse = 0.115065
  layer_ari = 0.320739

Best BioML layer ARI:
  expression_heavy_s020_e060_a020
  rmse = 0.115123
  layer_ari = 0.387113
  layer_nmi = 0.490091
```

Interpretation:

- BioML now beats the formal benchmark on both RMSE and biological consistency.
- Spectral multi-view graph recovery solves the previous layer consistency bottleneck better than single-gene expression-aware GP variants.
- Expression-heavy graph weighting improves biological consistency, while no-APA graph weighting gives the best RMSE, suggesting APA-view similarity should be revisited before locking the final default.

### Phase 7.5: Decision criteria

BioML 进入主线候选需要满足：

1. [x] RMSE 不显著劣于 `spatial GP`
2. [x] layer ARI/NMI 超过 `stapaminer_knn_expression`
3. [x] calibration 继承或接近 GP backbone
4. [x] runtime 在 CPU 上可接受

Current decision:

```text
BioML is promoted from research branch to main-line candidate.
```

Remaining validation before manuscript freeze:

- Expand to larger gene sets and additional random seeds.
- Validate on at least one external spatial transcriptomics dataset.
- Add boundary consistency / marker ordering / SVAPA ranking stability metrics.
- Choose final graph weights through cross-dataset evidence rather than MOB-only tuning.

### Phase 8: External real-data validation runner

BioML has now been moved from a benchmark-only implementation into the main
pipeline/CLI and into an external real-data validation runner:

```text
scripts/run_external_bioml_validation.py
```

The runner explicitly follows the experimental patterns in the reference
packages:

- stAPAminer-style metrics:
  - layer/domain ARI, NMI, purity, pairwise Jaccard
  - DBI, Calinski-Harabasz, Silhouette, Dunn
  - within-layer APA Pearson correlation
  - SVAPA / LSAPA / DEAPA output tables
  - replicate spatial gene overlap table when multiple datasets are supplied
- metaAPA-style readiness:
  - site count availability
  - APA site annotation availability
  - multi-caller consensus availability
  - sequence-context / long-read support checklist

Important design decision:

```text
metaAPA is treated as upstream site integration evidence, not as a downstream
layer recovery baseline.
```

This avoids an unfair comparison across different layers of the analysis stack.

MOB external pilot:

```text
Output:
  spaGAPA/benchmark_results/real/external_bioml_validation_v1/

n_genes = 60
min_observed_spots = 110
methods:
  raw
  stapaminer_original_imputed
  stapaminer_knn_expression
  spagapa_gp
  spagapa_bioml
```

Layer recovery:

```text
raw:
  layer_ari = 0.101117
  layer_nmi = 0.151950
  purity = 0.400000

stapaminer_original_imputed:
  layer_ari = 0.109863
  layer_nmi = 0.160847
  purity = 0.407692

stapaminer_knn_expression:
  layer_ari = 0.073448
  layer_nmi = 0.109605
  purity = 0.373077

spagapa_gp:
  layer_ari = 0.100875
  layer_nmi = 0.160667
  purity = 0.380769

spagapa_bioml:
  layer_ari = 0.388330
  layer_nmi = 0.509298
  purity = 0.638462
```

Downstream counts:

```text
raw:
  n_svapa = 9
  n_lsapa = 4
  n_deapa = 6

stapaminer_original_imputed:
  n_svapa = 9
  n_lsapa = 4
  n_deapa = 5

stapaminer_knn_expression:
  n_svapa = 9
  n_lsapa = 10
  n_deapa = 6

spagapa_gp:
  n_svapa = 9
  n_lsapa = 8
  n_deapa = 8

spagapa_bioml:
  n_svapa = 9
  n_lsapa = 9
  n_deapa = 8
```

Interpretation:

- BioML is the strongest method in layer ARI/NMI/purity on this external-validation-style MOB pilot.
- The improvement over stAPAminer-original is large enough to support the BIB direction:
  `0.388` vs `0.110` layer ARI.
- Within-layer correlation is not discriminative in this pilot because the selected genes are already high coverage.
- The next validation must include additional datasets and a low-coverage gene stratum.

Implementation notes:

- The runner uses the formal `SpaGAPA(..., use_bioml=True)` API.
- `gp_alpha` and `gp_n_restarts_optimizer` are now exposed in the main pipeline/CLI to align formal benchmarks with reproducible CLI runs.
- Matrix gene and spot IDs are cast to strings during external validation to prevent silent reindexing failures for numeric gene IDs.

---

## 9. 预期论文叙事

如果 BioML 成功，spaGAPA 的 BIB 叙事可以升级为：

```text
spaGAPA combines uncertainty-calibrated spatial GP imputation with
CPU-efficient multi-view graph machine learning to recover biologically
coherent spatial APA domains without deep learning or label supervision.
```

这比单纯宣传一个 GP imputer 更强，也更符合 BIB 对软件方法学完整性的期待。

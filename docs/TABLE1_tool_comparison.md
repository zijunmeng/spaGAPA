# Table 1 — Spatial APA tool capability comparison

**Caption.** Comparison of spaGAPA with the three existing spatial/alternative-polyadenylation analysis frameworks most relevant to this work: stAPAminer (the direct Visium-scale spatial APA competitor), spvAPA (an imputation-plus-supervised-selection tool covering spatial and single-cell APA), and metaAPA (a cross-tool PAS-integration workflow). "N/A" denotes a capability that is outside the tool's intended scope rather than a deficiency. All claims are documented in the manuscript Methods (M8–M9) and Results (R2–R3, R7–R8); version numbers, run conditions, and the reproducibility manifest are given in M12.

| Feature | stAPAminer | spvAPA | metaAPA | spaGAPA |
|---|---|---|---|---|
| **Core task** | Spatial APA imputation / pattern | Imputation + supervised selection | PAS-caller output integration | Probabilistic spatial inference |
| **Input** | APA + expression matrix | APA + expression + labels | Multiple PAS-caller outputs | PAS × spot matrix + coordinates |
| **Spatial modeling** | KNN (heuristic) | WNN (heuristic) | None | Sparse GP (probabilistic) |
| **Uncertainty** | None | None | Site confidence (caller-level) | Conformal-calibrated intervals |
| **Supervised labels** | Optional layer labels | Required (sPLS-DA) | None | None (unsupervised) |
| **High-resolution** | Not validated (Visium) | Not validated (Visium) | N/A | Stereo-seq (21 k PAS) |
| **Scalability** | O(n²), fails @ 42 k | O(n²), fails @ 42 k | N/A | O(n·m²), completes @ 100 k |
| **Batch correction** | None | None | None | QN + linear (optional module) |
| **Primary output** | RUD, SVAPA | Features, visualization | Integrated sites | Posterior mean + interval + domains |
| **Validation data** | MOB × 3 replicates | 9 datasets (sc + ST) | 4 datasets | 9 GSE × 32 samples + Stereo-seq |

**Notes on the comparison.**

- *Scalability and high-resolution rows* report the empirical outcome under the predefined resource limits of M9 (1200 s wall cap, `OPENBLAS_NUM_THREADS = 8`). The failure of stAPAminer/spvAPA at 42 k spots is a timeout / `future.globals.maxSize` outcome under those limits, reported factually; we do not claim the tools are intrinsically unable to scale beyond Visium.
- *Uncertainty row.* spaGAPA is the only tool whose prediction intervals carry a distribution-free marginal coverage guarantee (split-conformal, M5); metaAPA's "site confidence" is a caller-level score from the upstream tools it integrates and is not a calibrated prediction interval on an imputed value.
- *Supervised labels row.* spaGAPA's spatial-domain identification (M7) is fully unsupervised (Leiden on a fused spatial + APA + expression graph); spvAPA's supervised feature selection (sPLS-DA) is a complementary capability that spaGAPA does not provide.
- *Batch correction row.* The QN + linear module is an optional, still-under-evaluation component (Pillar 2; Harmony comparison pending). It is listed for completeness.

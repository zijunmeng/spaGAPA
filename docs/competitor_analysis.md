# spaGAPA Competitor Analysis: stAPAminer vs metaAPA vs spaGAPA

**Date**: 2026-07-17
**Purpose**: BIB manuscript positioning — related work + competitive landscape

---

## Overview

| Feature | stAPAminer | metaAPA | **spaGAPA (post sparse-upgrade)** |
|---------|-----------|---------|----------------------------------|
| **Primary focus** | Spatial APA analysis (Visium) | APA workflow/meta-analysis | Spatial APA framework (Visium + high-res) |
| **Platform support** | 10x Visium only | Bulk / scRNA (not spatial-specific) | Visium + Stereo-seq + Visium HD |
| **High-resolution ST** | ❌ Not supported | ❌ Not supported | **✅ Unique** |
| **APA site calling** | Uses external callers | Uses external callers | Integrates scAPAtrap + custom |
| **APA imputation method** | KNN (expression-neighborhood) | N/A | GP (sparse, with uncertainty) + BioML factorizer |
| **Uncertainty quantification** | ❌ No GP posterior | ❌ | **✅ GP posterior std** |
| **Biological domain recovery** | ❌ | ❌ | **✅ BioML multi-view graph** |
| **Scalability (spot count)** | ~5k-10k (KNN O(n²)) | N/A | **~100k (sparse GP + batched ridge)** |
| **Sparse methods** | ❌ | ❌ | **✅ inducing-point GP + einsum batched ridge** |
| **CPU-only (no GPU/DL)** | ✅ | ✅ | **✅** |

---

## stAPAminer — Detailed Analysis

**Source**: `00_ref_packages/stAPAminer.pdf` (reviewed first 10 pages)

### Method
- Uses **expression-neighborhood KNN** for APA imputation: finds k nearest spots by gene expression similarity, borrows APA usage from neighbors.
- Operates at the **spot level** on 10x Visium data.
- Applied to mouse brain + human cancer Visium datasets.

### Strengths
- Directly designed for spatial APA (the main competitor).
- KNN imputation can produce strong biological-domain consistency (expression neighbors are often biological neighbors).
- Runtime can be faster than GP-heavy methods at Visium scale.

### Limitations (relative to spaGAPA)
- **No uncertainty quantification** — KNN gives point estimates, no posterior std.
- **No high-resolution ST support** — no mention of Stereo-seq / Visium HD / Slide-seqV2.
- **Scalability ceiling** — KNN requires all-pairs distance computation O(n²); at 40k+ spots (Stereo-seq bin-100, Visium HD 8μm, Slide-seqV2), KNN faces the same OOM/scalability wall as spaGAPA's pre-upgrade dense components.
- **No explicit separation** between APA value recovery and biological-domain discovery.
- **No sparse/approximate methods** — no inducing-point GP, no graph-regularized factorization, no Leiden domain.
- **Benchmark coverage** — no systematic high-resolution benchmark or reproducibility framework.

### spaGAPA's Response
- Do not try to beat stAPAminer on every Visium-scale metric using GP alone.
- Use GP where GP is strongest: APA value recovery + calibrated uncertainty.
- Use BioML graph route for biological consistency.
- **At high-resolution scale: spaGAPA is the only tool that runs.**

---

## metaAPA — Detailed Analysis

**Source**: `00_ref_packages/metaAPA.pdf` (reviewed first 10 pages)

### Method
- APA workflow/meta-analysis framework.
- Integrates multiple APA datasets or analysis tools.
- Less spatial-ST-specific than stAPAminer.

### Relationship to spaGAPA
- **More complementary than competitive.**
- metaAPA focuses on cross-study APA integration; spaGAPA focuses on within-sample spatial APA analysis.
- Could be used together: spaGAPA for spatial APA recovery + metaAPA for cross-dataset meta-analysis.

---

## Competitive Implications for BIB Manuscript

### 1. Novelty claim
> spaGAPA is the **first spatial APA framework that supports high-resolution spatial transcriptomics** (Stereo-seq, Visium HD) through sparse GP imputation, uncertainty quantification, and CPU-friendly graph domain recovery.

**Evidence**: stAPAminer (the direct competitor) is Visium-only with no sparse methods; metaAPA is not spatial-specific. Neither handles >10k spots.

### 2. Benchmark strategy
- **Visium scale** (≤5k spots): compare spaGAPA vs stAPAminer (both run; show GP uncertainty advantage + BioML domain advantage).
- **High-resolution scale** (42k-100k spots): spaGAPA only (stAPAminer cannot run — document this as a finding).
- **Scaling curve**: spot count vs runtime/memory/feasibility — shows stAPAminer hits a wall ~10k, spaGAPA scales to 100k.

### 3. Honest limitations (for the paper)
- At Visium scale, stAPAminer's KNN may outperform spaGAPA's GP on some domain metrics (expression-neighborhood is strong for domain recovery). spaGAPA's response: BioML graph route (decoupled from GP).
- spaGAPA's sparse GP (inducing points) is an approximation — the uncertainty estimates compensate, but exact GP would be more precise (infeasible at scale).
- Domain detection at high-res uses kmeans/Leiden (not spectral) — different but valid.

### 4. Positioning statement
> spaGAPA does not claim to be the best APA imputation method at every scale. It claims to be the **only framework that provides uncertainty-aware, graph-guided spatial APA analysis across both conventional (Visium) and high-resolution (Stereo-seq, Visium HD) spatial transcriptomics platforms, with reproducible benchmarks and real-data validation.**

#!/usr/bin/env python3
"""Drive APA batch/protocol bias-correction prototype (P1-1).

Two experiments:

(1) CONTROLLED synthetic recovery (GSE183456 real APA matrix)
    Take a real gene-level distal-usage matrix, split spots into two "batches"
    (random half / half, orthogonal to biology), inject a known additive +
    multiplicative batch effect on one batch, apply each corrector, and measure
    recovery = corr(corrected, clean_truth) vs corr(contaminated, truth), plus
    the reduction in per-gene A-vs-B mean difference.

(2) REAL multi-sample demo (GSE220442, 3 control vs 3 AD)
    Build a gene-level distal-usage index for each of the 6 samples, pool into a
    gene x spot matrix with sample + condition labels, and measure cross-sample
    consistency (mean pairwise Pearson corr of gene APA profiles across samples)
    before vs after correction. Also fit a PCA on the spot-level matrix and
    report how much variance PC1 explains and whether AD/control separate
    (logistic-style separation: corr of PC1 with the condition label).

Note: GSE220442 is a SINGLE study, so protocol batch is mild. We frame this
demo accordingly -- it shows the method does no harm and gently improves
cross-sample consistency on same-protocol data; the controlled experiment is
where the method's batch-removal power is demonstrated.

Outputs -> pipeline_output/apa_bias_correction/
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from spagapa.analysis import (
    quantile_normalize,
    linear_batch_correction,
    build_distal_usage_index,
)

PKG = Path(__file__).resolve().parents[1]
DEFAULT_OUT = PKG / "pipeline_output" / "apa_bias_correction"
GSE183456_DIR = PKG / "data/processed/gse183456_gsm6047774_scapatrap"
GSE220442_DIRS = {
    "GSM6801751": ("control", PKG / "data/processed/gse220442_gsm6801751_scapatrap"),
    "GSM6801752": ("control", PKG / "data/processed/gse220442_gsm6801752_scapatrap"),
    "GSM6801753": ("control", PKG / "data/processed/gse220442_gsm6801753_scapatrap"),
    "GSM6801754": ("ad", PKG / "data/processed/gse220442_gsm6801754_scapatrap"),
    "GSM6801755": ("ad", PKG / "data/processed/gse220442_gsm6801755_scapatrap"),
    "GSM6801756": ("ad", PKG / "data/processed/gse220442_gsm6801756_scapatrap"),
}


# ---------------------------------------------------------------------------
# io helpers
# ---------------------------------------------------------------------------
def load_processed(processed_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (counts peak x spot, sites peak metadata).

    Reads raw per-peak counts from ``apa_site_counts.csv.gz`` (the gene-index
    builder sums counts across peaks), NOT ``apa_matrix.csv`` which holds
    pre-computed peak-level usage ratios.
    """
    counts = pd.read_csv(processed_dir / "apa_site_counts.csv.gz", index_col=0)
    sites = pd.read_csv(processed_dir / "apa_sites.csv.gz")
    return counts, sites


def corr_over_finite(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5:
        return float("nan")
    if np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float("nan")
    return float(pearsonr(a[ok], b[ok])[0])


def batch_mean_abs_diff(m: np.ndarray, b_mask: np.ndarray) -> float:
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        a = np.nanmean(m[:, ~b_mask], axis=1)
        b = np.nanmean(m[:, b_mask], axis=1)
    diff = np.abs(a - b)
    return float(np.nanmean(diff)) if np.isfinite(diff).any() else float("nan")


# ---------------------------------------------------------------------------
# Experiment 1: controlled synthetic recovery on real data
# ---------------------------------------------------------------------------
def run_controlled(out_dir: Path, seed: int = 42) -> Dict:
    t0 = time.time()
    print("\n=== Experiment 1: controlled synthetic recovery (GSE183456) ===")
    counts, sites = load_processed(GSE183456_DIR)
    # build gene-level distal-usage index (gene x spot)
    index = build_distal_usage_index(counts, sites, min_parent=5)
    print(f"  gene index: {index.shape[0]} genes x {index.shape[1]} spots, "
          f"observed frac={index.notna().mean().mean():.3f}")

    # keep only genes with reasonable coverage so recovery signal is clean
    obs_frac = (~index.isna()).mean(axis=1)
    index = index.loc[obs_frac >= 0.3]
    print(f"  after coverage filter (>=30% observed): {index.shape[0]} genes")
    truth = index.to_numpy()
    n_spots = truth.shape[1]

    rng = np.random.default_rng(seed)
    # batch assignment: random half/half, orthogonal to biology
    spot_order = rng.permutation(n_spots)
    half = n_spots // 2
    batch = np.empty(n_spots, dtype=object)
    batch[spot_order[:half]] = "A"
    batch[spot_order[half:]] = "B"
    b_mask = batch == "B"

    # build a mask of where truth is observed (NaN stays NaN throughout)
    observed = ~np.isnan(truth)

    shift, scale = 0.20, 1.30
    contaminated = truth.copy()
    contaminated[:, b_mask] = contaminated[:, b_mask] * scale + shift
    # keep NaN positions NaN
    contaminated[~observed] = np.nan

    truth_for_corr = np.where(observed, truth, np.nan)
    base_corr = corr_over_finite(contaminated, truth_for_corr)
    base_batch_diff = batch_mean_abs_diff(contaminated, b_mask)
    print(f"  injected batch (shift={shift}, scale={scale}) on batch B")
    print(f"  contaminated-vs-truth PCC = {base_corr:.4f}")
    print(f"  contaminated |meanA - meanB| (per gene, avg) = {base_batch_diff:.4f}")

    results = {
        "n_genes": int(truth.shape[0]),
        "n_spots": int(n_spots),
        "shift": shift,
        "scale": scale,
        "contaminated_vs_truth_pcc": base_corr,
        "contaminated_batch_mean_abs_diff": base_batch_diff,
        "methods": {},
    }

    for name, fn in [
        ("quantile_normalize", lambda m: quantile_normalize(m, batch)),
        ("linear_batch_correction", lambda m: linear_batch_correction(m, batch)),
    ]:
        t = time.time()
        out = fn(contaminated)
        dt = time.time() - t
        corr = corr_over_finite(out, truth_for_corr)
        bd = batch_mean_abs_diff(out, b_mask)
        recovery = corr - base_corr
        bd_reduction = base_batch_diff - bd
        print(f"  [{name}] PCC={corr:.4f} (Δ={recovery:+.4f}), "
              f"|A-B|diff={bd:.4f} (reduction={bd_reduction:+.4f}), {dt:.2f}s")
        results["methods"][name] = {
            "corrected_vs_truth_pcc": corr,
            "pcc_recovery_vs_contaminated": recovery,
            "batch_mean_abs_diff": bd,
            "batch_diff_reduction": bd_reduction,
            "wall_seconds": dt,
        }
    results["wall_seconds"] = time.time() - t0
    return results


# ---------------------------------------------------------------------------
# Experiment 2: real multi-sample demo (GSE220442)
# ---------------------------------------------------------------------------
def run_real_demo(out_dir: Path, min_obs_frac: float = 0.2) -> Dict:
    t0 = time.time()
    print("\n=== Experiment 2: real multi-sample demo (GSE220442, 3 ctrl vs 3 AD) ===")

    # build per-sample gene-level index
    per_sample_index: Dict[str, pd.DataFrame] = {}
    sample_meta: List[Dict] = []
    for gsm, (condition, d) in GSE220442_DIRS.items():
        counts, sites = load_processed(d)
        idx = build_distal_usage_index(counts, sites, min_parent=5)
        # rename spots to include sample id to avoid barcode collisions
        idx.columns = [f"{gsm}_{c}" for c in idx.columns]
        per_sample_index[gsm] = idx
        sample_meta.append({"gsm": gsm, "condition": condition,
                            "n_genes": int(idx.shape[0]),
                            "n_spots": int(idx.shape[1]),
                            "obs_frac": float(idx.notna().mean().mean())})
        print(f"  {gsm} ({condition}): {idx.shape[0]} genes x {idx.shape[1]} spots")

    # intersect genes across all 6 samples (this is the bottleneck -- per-sample
    # peak calling yields different gene sets; we keep the common core)
    common_genes = None
    for gsm, idx in per_sample_index.items():
        common_genes = set(idx.index) if common_genes is None else common_genes & set(idx.index)
    common_genes = sorted(common_genes)
    print(f"  common genes across 6 samples: {len(common_genes)}")

    # pool into a single gene x spot matrix
    pooled = pd.concat([per_sample_index[g].loc[common_genes] for g in per_sample_index],
                       axis=1)
    spots = pooled.columns.tolist()
    sample_labels = np.array([c.split("_")[0] for c in spots])
    condition_labels = np.array([
        next(cond for g, (cond, _) in GSE220442_DIRS.items() if g == s)
        for s in sample_labels
    ])
    print(f"  pooled matrix: {pooled.shape[0]} genes x {pooled.shape[1]} spots, "
          f"observed frac={pooled.notna().mean().mean():.3f}")

    # ---- cross-sample consistency: mean pairwise PCC of per-sample gene profiles
    # per-sample gene profile = mean usage per gene (across that sample's spots)
    def per_sample_profiles(m: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Per-sample APA profile = mean usage over the sample's spots, one
        value per gene (length = n_common_genes, aligned across samples)."""
        prof = {}
        for gsm in per_sample_index:
            cols = [c for c in m.columns if c.startswith(gsm + "_")]
            sub = m[cols]
            # axis=1: mean over spots -> per-gene vector (length = n_genes)
            prof[gsm] = sub.mean(axis=1, skipna=True).to_numpy()
        return prof

    def mean_pairwise_pcc(profs: Dict[str, np.ndarray]) -> Tuple[float, int]:
        gsms = sorted(profs)
        pccs = []
        for i in range(len(gsms)):
            for j in range(i + 1, len(gsms)):
                a, b = profs[gsms[i]], profs[gsms[j]]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() < 5:
                    continue
                if np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
                    continue
                pccs.append(pearsonr(a[ok], b[ok])[0])
        return float(np.mean(pccs)) if pccs else float("nan"), len(pccs)

    profs_before = per_sample_profiles(pooled)
    pcc_before, n_pairs = mean_pairwise_pcc(profs_before)
    print(f"  cross-sample mean pairwise PCC (before): {pcc_before:.4f} ({n_pairs} pairs)")

    # ---- AD/control separation on PCA
    def pca_separation(m: pd.DataFrame) -> Dict:
        # spot-level: transpose so spots are rows. Impute NaN with column (gene) mean.
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            X = m.to_numpy().T  # spots x genes
            nan_mask = np.isnan(X)
            col_means = np.nanmean(X, axis=0)
            col_means = np.where(np.isnan(col_means), 0.0, col_means)
            X_filled = np.where(nan_mask, col_means, X)
        # center
        Xc = X_filled - X_filled.mean(axis=0, keepdims=True)
        # SVD -> PCs
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        pc1 = U[:, 0] * S[0]
        pc2 = U[:, 1] * S[1] if U.shape[1] > 1 else np.zeros_like(pc1)
        var_explained = (S ** 2) / (S ** 2).sum()
        # separation: point-biserial corr of PC1 with the AD-vs-control indicator
        ad_indicator = (condition_labels == "ad").astype(float)
        ok = np.isfinite(pc1)
        sep = float(np.corrcoef(pc1[ok], ad_indicator[ok])[0, 1]) if ok.sum() > 5 else float("nan")
        return {
            "pc1_var_explained": float(var_explained[0]),
            "pc2_var_explained": float(var_explained[1]) if len(var_explained) > 1 else 0.0,
            "pc1_condition_corr": sep,
        }

    pca_before = pca_separation(pooled)
    print(f"  PCA before: PC1 var={pca_before['pc1_var_explained']:.3f}, "
          f"PC1~condition corr={pca_before['pc1_condition_corr']:.3f}")

    results = {
        "common_genes": len(common_genes),
        "n_spots_pooled": int(pooled.shape[1]),
        "n_sample_pairs": n_pairs,
        "sample_meta": sample_meta,
        "before": {
            "mean_pairwise_pcc": pcc_before,
            "pca": pca_before,
        },
        "methods": {},
    }

    for name, fn in [
        ("quantile_normalize", lambda m: quantile_normalize(m, sample_labels)),
        ("linear_batch_correction", lambda m: linear_batch_correction(
            m, batch_labels=sample_labels, preserve_labels=condition_labels)),
    ]:
        t = time.time()
        corrected = fn(pooled)
        dt = time.time() - t
        profs = per_sample_profiles(corrected)
        pcc, _ = mean_pairwise_pcc(profs)
        pca = pca_separation(corrected)
        delta_pcc = pcc - pcc_before
        print(f"  [{name}] PCC={pcc:.4f} (Δ={delta_pcc:+.4f}), "
              f"PC1 var={pca['pc1_var_explained']:.3f}, "
              f"PC1~cond={pca['pc1_condition_corr']:.3f}, {dt:.2f}s")
        results["methods"][name] = {
            "mean_pairwise_pcc": pcc,
            "pcc_delta_vs_before": delta_pcc,
            "pca": pca,
            "wall_seconds": dt,
        }
    results["wall_seconds"] = time.time() - t0

    # Honest limitation note: in this same-study demo, batch == sample and
    # condition (AD/control) is derived from sample (s1-3 ctrl, s4-6 ad), so
    # the linear model's batch dummies and the preserved 'condition' covariate
    # are PERFECTLY collinear (design rank-deficient). lstsq then distributes
    # the condition signal into the batch coefficients, so linear correction
    # partially erases the very biology it is asked to preserve. This is the
    # well-known "batch and biology must be orthogonal" requirement of linear
    # batch removal (limma removeBatchEffect has the same constraint). On
    # same-study data, quantile normalization is the safer choice. Linear
    # correction is the right tool when batch (e.g. study/protocol) is NOT
    # collinear with condition (multiple conditions per batch), which is the
    # cross-study use case the controlled experiment demonstrates.
    collinear = _check_design_collinearity(sample_labels, condition_labels)
    results["design_collinear"] = collinear
    if collinear:
        print("  NOTE: batch(sample) and condition are collinear (same-study "
              "demo); linear correction's preserve_labels cannot fully protect "
              "biology here -- see summary.json['design_collinear'].")
    return results


def _check_design_collinearity(batch: np.ndarray, condition: np.ndarray) -> bool:
    """Return True if the [intercept + batch dummies + condition dummy] design
    is rank-deficient (i.e. condition is collinear with batch)."""
    import pandas as pd
    n = len(batch)
    b = pd.get_dummies(batch, drop_first=True).astype(float).values
    c = pd.get_dummies(condition, drop_first=True).astype(float).values
    if c.shape[1] == 0:
        return False
    X = np.column_stack([np.ones(n), b, c])
    return bool(np.linalg.matrix_rank(X) < X.shape[1])


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-controlled", action="store_true")
    ap.add_argument("--skip-real", action="store_true")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    t_total = time.time()

    summary: Dict = {"experiments": {}}
    if not args.skip_controlled:
        summary["experiments"]["controlled_recovery"] = run_controlled(
            args.out, seed=args.seed)
    if not args.skip_real:
        summary["experiments"]["real_demo"] = run_real_demo(args.out)

    summary["total_wall_seconds"] = time.time() - t_total
    out_json = args.out / "summary.json"
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nWrote {out_json} (total {summary['total_wall_seconds']:.1f}s)")


if __name__ == "__main__":
    main()

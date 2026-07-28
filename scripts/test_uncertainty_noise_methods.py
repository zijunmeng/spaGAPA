#!/usr/bin/env python3
"""
Improve raw GP uncertainty-error correlation by testing better noise
estimation methods across representative Visium APA datasets.

THE PROBLEM
-----------
Expert reviewer attack: "Raw GP uncertainty-error correlation is ~0.07 --
the uncertainty is essentially uninformative for ranking prediction quality."

Root cause: a constant ``noise_level=0.1`` makes the posterior std dominated
by distance to the inducing points -> nearly flat across tissue -> corr~0.07.

METHODS TESTED (per dataset, on held-out entries)
------------------------------------------------
A. Constant noise (0.1)             -- baseline (current default)
B. local_noise (per-gene MAD of kNN neighbours)  -- existing Phase 3 addition
C. Spatial noise  (per-spot variance of kNN=8 neighbours' values, per gene)
D. Residual noise  (2-pass empirical-Bayes: fit constant -> per-spot residual
   variance -> re-fit with per-spot noise)

For each method x dataset we report:
  - unc-error Pearson r
  - unc-error Spearman rho
  - n_test (held-out points evaluated)

KEY QUESTION: can we lift raw corr from ~0.07 to >0.3 across multiple datasets?

Outputs (pipeline_output/uncertainty_corr_improvement/)
-------------------------------------------------------
  per_dataset_corr.csv   dataset, method, unc_error_pearson,
                         unc_error_spearman, n_test
  summary.json           mean corr per method across datasets
  figure.png             bar chart: corr by method x dataset
  run.log                execution log
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr

# --- repo-local imports (must run after sys.path has repo root) ---
REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from spagapa.imputation import SparseGPImputer  # noqa: E402

DATA_ROOT = os.path.join(REPO, "data", "processed")
OUT_DIR = os.path.join(REPO, "pipeline_output", "uncertainty_corr_improvement")

# 5 representative datasets spanning tissues / organisms.
# (dataset_label, sample_label, dir_basename)
SAMPLES = [
    ("GSE183456", "kidney",     "gse183456_gsm6047774_scapatrap"),   # kidney 3010
    ("GSE220442", "brain",      "gse220442_gsm6801751_scapatrap"),   # brain 4179
    ("GSE169749", "mouse_colon","gse169749_gsm5213483_scapatrap"),   # mouse colon 2715
    ("GSE338525", "liver",      "gse338525_gsm9876373_scapatrap"),   # liver 4992
    ("GSE237183", "glioma",     "gse237183_gsm7596587_scapatrap"),   # glioma 2541
]

METHODS = ["A_constant", "B_local_gene", "C_spatial_spot", "D_residual_spot"]

# ---- protocol hyperparameters (mirrors calibrate_uncertainty_all_datasets) ----
MASK_FRACTION = 0.20
SEED = 42
MAX_GENES = 300          # subsample top-K most-observed sites per sample
MIN_SPOTS = 500
MIN_GENES_WITH_OBS = 50
K_NEIGHBORS = 8          # kNN for spatial noise estimation (Method C)
RESID_K = 8              # kNN neighbourhood for residual-variance binning (D)
NOISE_SCALE = 1000.0     # kernel-amplitude scaling (mirrors local_noise default)


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


def _safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) >= 2 and np.std(a) > 0 and np.std(b) > 0:
        try:
            r, _ = pearsonr(a, b)
            return float(r)
        except Exception:
            return float("nan")
    return float("nan")


def _safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) >= 3 and np.std(a) > 0 and np.std(b) > 0:
        try:
            r, _ = spearmanr(a, b)
            return float(r)
        except Exception:
            return float("nan")
    return float("nan")


def load_sample(dir_basename: str):
    """Load apa_matrix + coordinates, return (values[n_genes,n_spots], xy)."""
    apa = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "apa_matrix.csv"),
                      index_col=0)
    coords = pd.read_csv(os.path.join(DATA_ROOT, dir_basename, "coordinates.csv"),
                         index_col=0)
    if {"x", "y"} <= set(coords.columns):
        xy_df = coords[["x", "y"]]
    else:
        xy_df = coords.iloc[:, :2].astype(float)
        xy_df.columns = ["x", "y"]
    common = apa.columns.intersection(xy_df.index)
    if len(common) > 0:
        apa = apa.loc[:, common]
        xy_df = xy_df.loc[common]
    elif apa.shape[1] != xy_df.shape[0]:
        raise ValueError(f"spot mismatch: apa={apa.shape}, coords={xy_df.shape}")
    xy = xy_df.values.astype(float)
    values = apa.values.astype(float)
    return values, xy


def prep_sample(dir_basename: str, rng: np.random.Generator):
    """Load, subsample top-K genes, mask 20% per row, return shared tensors."""
    values, xy = load_sample(dir_basename)
    n_genes_full, n_spots = values.shape

    obs_count = (values > 0).sum(axis=1)
    keep_rows = np.where(obs_count >= 5)[0]
    if len(keep_rows) == 0:
        raise RuntimeError("no genes with >=5 observations")
    order = keep_rows[np.argsort(-obs_count[keep_rows])]
    keep = order[:MAX_GENES]
    values_sub = values[keep]
    n_genes = values_sub.shape[0]

    if n_spots < MIN_SPOTS:
        raise RuntimeError(f"too few spots ({n_spots} < {MIN_SPOTS})")
    if n_genes < MIN_GENES_WITH_OBS:
        raise RuntimeError(f"too few genes ({n_genes} < {MIN_GENES_WITH_OBS})")

    # mask 20% per row
    masked = values_sub.copy()
    for g in range(n_genes):
        obs_idx = np.where(values_sub[g] > 0)[0]
        if len(obs_idx) < 5:
            continue
        n_mask = max(1, int(len(obs_idx) * MASK_FRACTION))
        m = rng.choice(obs_idx, size=n_mask, replace=False)
        masked[g, m] = 0.0

    # shared GP hyperparameters (mirror calibrate_uncertainty_all_datasets)
    kdt = cKDTree(xy)
    nn = kdt.query(xy, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    length_scale = nn_dist * 5
    n_inducing = min(500, max(100, n_spots // 100))

    return {
        "values": values_sub,
        "masked": masked,
        "xy": xy,
        "n_genes": n_genes,
        "n_spots": n_spots,
        "nn_dist": nn_dist,
        "length_scale": length_scale,
        "n_inducing": n_inducing,
        "tree": kdt,
    }


def gather_held(values_sub, masked, gp_pred, gp_unc):
    """Pull (truth, pred, unc, abs_err) for held-out entries."""
    truth_obs = values_sub > 0
    masked_zero = masked <= 0
    held = truth_obs & masked_zero
    g_idx, s_idx = np.where(held)
    truth = values_sub[g_idx, s_idx]
    pred = gp_pred[g_idx, s_idx]
    unc = gp_unc[g_idx, s_idx]
    err = np.abs(truth - pred)
    return truth, pred, unc, err


# --------------------------------------------------------------------------
# Noise estimators
# --------------------------------------------------------------------------
def estimate_spatial_noise(values_sub, masked, tree, k=K_NEIGHBORS,
                           scale=NOISE_SCALE):
    """Method C: per-spot noise = variance of k nearest spatial neighbours'
    OBSERVED values, for each gene.

    Returns spot_noise of shape (n_genes, n_spots) in kernel-amplitude units.
    Spots with no observed neighbours fall back to the gene's overall variance.
    """
    n_genes, n_spots = masked.shape
    # neighbour indices shared across genes (depends only on coordinates)
    k_use = min(k, n_spots - 1)
    _, idx = tree.query(tree.data, k=k_use + 1)
    nbr_idx = idx[:, 1:]                  # (n_spots, k)

    spot_noise = np.empty((n_genes, n_spots), dtype=float)
    obs_mask = masked > 0                  # observed-after-masking
    for g in range(n_genes):
        v = masked[g]
        obs = obs_mask[g]
        # variance over each spot's observed neighbours
        nbr_vals = v[nbr_idx]                       # (n_spots, k)
        nbr_obs = obs[nbr_idx]                      # (n_spots, k) bool
        cnt = nbr_obs.sum(axis=1)
        # masked-select: where neighbour unobserved, treat as NaN -> nanvar
        nbr_vals_m = np.where(nbr_obs, nbr_vals, np.nan)
        with np.errstate(all="ignore"):
            with np.testing.suppress_warnings() as sup:
                sup.filter(RuntimeWarning, "Degrees of freedom <= 0 for slice.")
                var = np.nanvar(nbr_vals_m, axis=1)     # (n_spots,)
        # fallback for spots with <2 observed neighbours: gene variance
        gene_var = float(np.var(v[obs])) if obs.any() and np.var(v[obs]) > 0 else 0.01
        bad = (cnt < 2) | ~np.isfinite(var)
        if np.any(bad):
            var[bad] = gene_var
        var = np.clip(var, 1e-6, None)
        spot_noise[g] = var * scale
    return spot_noise


def estimate_residual_noise(values_sub, masked, xy, residual_pred, tree,
                            k=RESID_K, scale=NOISE_SCALE):
    """Method D: per-spot noise estimated from GP residuals on TRAINING spots.

    Two-pass empirical-Bayes:
      pass 1 already done by caller (fit with constant noise -> residual_pred).
      here: residual[g,s] = observed - predicted, on observed entries.
      Bin each spot's residual variance over its k nearest spatial neighbours
      (using only observed spots).  Scale into kernel-amplitude units.

    Returns spot_noise of shape (n_genes, n_spots).
    """
    n_genes, n_spots = masked.shape
    k_use = min(k, n_spots - 1)
    _, idx = tree.query(tree.data, k=k_use + 1)
    nbr_idx = idx[:, 1:]

    spot_noise = np.empty((n_genes, n_spots), dtype=float)
    obs_mask = masked > 0
    for g in range(n_genes):
        v = masked[g]
        obs = obs_mask[g]
        resid = np.where(obs, v - residual_pred[g], np.nan)
        # variance of residuals over each spot's observed neighbours
        nbr_resid = resid[nbr_idx]
        with np.errstate(all="ignore"):
            with np.testing.suppress_warnings() as sup:
                sup.filter(RuntimeWarning, "Degrees of freedom <= 0 for slice.")
                var = np.nanvar(nbr_resid, axis=1)
        cnt = np.sum(np.isfinite(nbr_resid), axis=1)
        # fallback: global residual variance
        glob = np.nanvar(resid)
        if not np.isfinite(glob) or glob <= 0:
            glob = float(np.var(v[obs])) * 0.5 if obs.any() else 0.01
        bad = (cnt < 2) | ~np.isfinite(var) | (var <= 0)
        if np.any(bad):
            var[bad] = max(glob, 1e-6)
        var = np.clip(var, 1e-6, None)
        spot_noise[g] = var * scale
    return spot_noise


# --------------------------------------------------------------------------
# Method runners
# --------------------------------------------------------------------------
def run_method_A(p, rng):
    """Constant noise = 0.1 (baseline)."""
    base = SparseGPImputer(n_inducing=p["n_inducing"], length_scale=p["length_scale"],
                           noise_level=0.1)
    batch = base.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                           verbose=False)
    pred, unc = batch.impute(return_uncertainty=True)
    return pred, unc


def run_method_B(p, rng):
    """local_noise (per-gene kNN MAD)."""
    base = SparseGPImputer(n_inducing=p["n_inducing"], length_scale=p["length_scale"],
                           local_noise=True, noise_scale=NOISE_SCALE)
    batch = base.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                           verbose=False)
    pred, unc = batch.impute(return_uncertainty=True)
    return pred, unc


def run_method_C(p, rng):
    """Spatial per-spot noise (variance of kNN neighbours' observed values)."""
    spot_noise = estimate_spatial_noise(p["values"], p["masked"], p["tree"])
    base = SparseGPImputer(n_inducing=p["n_inducing"], length_scale=p["length_scale"])
    batch = base.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                           verbose=False, spot_noise=spot_noise)
    pred, unc = batch.impute(return_uncertainty=True)
    return pred, unc


def run_method_D(p, rng):
    """Residual-based per-spot noise (2-pass empirical Bayes)."""
    # pass 1: fit with constant noise to get residuals on training spots
    base1 = SparseGPImputer(n_inducing=p["n_inducing"],
                            length_scale=p["length_scale"], noise_level=0.1)
    batch1 = base1.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                             verbose=False)
    pred1, _ = batch1.impute(return_uncertainty=False)
    # estimate per-spot noise from residuals
    spot_noise = estimate_residual_noise(p["values"], p["masked"], p["xy"],
                                         pred1, p["tree"])
    # pass 2: re-fit with per-spot noise
    base2 = SparseGPImputer(n_inducing=p["n_inducing"],
                            length_scale=p["length_scale"])
    batch2 = base2.fit_batch(p["xy"], p["masked"], mask=(p["masked"] > 0),
                             verbose=False, spot_noise=spot_noise)
    pred, unc = batch2.impute(return_uncertainty=True)
    return pred, unc


METHOD_RUNNERS = {
    "A_constant":     run_method_A,
    "B_local_gene":   run_method_B,
    "C_spatial_spot": run_method_C,
    "D_residual_spot":run_method_D,
}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    rng = np.random.default_rng(SEED)

    log_path = os.path.join(OUT_DIR, "run.log")
    log_fh = open(log_path, "w")

    def both(msg):
        log(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    both(f"Uncertainty-noise-method comparison over {len(SAMPLES)} datasets")
    both(f"Methods: {METHODS}")
    both(f"MASK_FRACTION={MASK_FRACTION} MAX_GENES={MAX_GENES} "
         f"K_NEIGHBORS={K_NEIGHBORS} NOISE_SCALE={NOISE_SCALE} SEED={SEED}")
    both("=" * 70)

    rows = []
    for dataset, label, dir_basename in SAMPLES:
        both("")
        both(f">>> {dataset}/{label} ({dir_basename})")
        try:
            t0 = time.time()
            p = prep_sample(dir_basename, rng)
            both(f"  loaded: {p['n_genes']} genes x {p['n_spots']} spots, "
                 f"obs_frac={(p['values']>0).mean():.3f}, "
                 f"nn_dist={p['nn_dist']:.1f}, "
                 f"n_inducing={p['n_inducing']}, length_scale={p['length_scale']:.1f}")
            for method in METHODS:
                tm0 = time.time()
                try:
                    pred, unc = METHOD_RUNNERS[method](p, rng)
                    truth, _, u, err = gather_held(p["values"], p["masked"],
                                                   pred, unc)
                    pear = _safe_pearson(u, err)
                    spear = _safe_spearman(u, err)
                    rows.append({
                        "dataset": dataset,
                        "label": label,
                        "method": method,
                        "unc_error_pearson": pear,
                        "unc_error_spearman": spear,
                        "n_test": int(len(err)),
                        "wall_s": float(time.time() - tm0),
                    })
                    both(f"  {method:18s} pearson={pear:+.4f}  "
                         f"spearman={spear:+.4f}  n={len(err)}  "
                         f"[{time.time()-tm0:.0f}s]")
                except Exception as e:
                    both(f"  {method:18s} ERROR: {e}")
                    rows.append({
                        "dataset": dataset, "label": label, "method": method,
                        "unc_error_pearson": float("nan"),
                        "unc_error_spearman": float("nan"),
                        "n_test": 0, "wall_s": float(time.time() - tm0),
                    })
            both(f"  dataset wall: {time.time()-t0:.0f}s")
        except Exception as e:
            tb = traceback.format_exc()
            both(f"  DATASET ERROR: {e}\n{tb}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUT_DIR, "per_dataset_corr.csv")
    df.to_csv(csv_path, index=False)
    both(f"\nWrote {csv_path}")

    # ---- summary: mean corr per method across datasets ----
    valid = df.dropna(subset=["unc_error_pearson"])
    summary = {
        "n_datasets": int(df["dataset"].nunique()),
        "methods": METHODS,
        "mean_corr_per_method": {},
        "median_corr_per_method": {},
        "mean_spearman_per_method": {},
        "target_corr": 0.30,
        "baseline_constant_mean": None,
        "baseline_constant_median": None,
        "achieves_target": {},
    }
    for m in METHODS:
        sub = valid[valid["method"] == m]
        if len(sub):
            summary["mean_corr_per_method"][m] = float(sub["unc_error_pearson"].mean())
            summary["median_corr_per_method"][m] = float(sub["unc_error_pearson"].median())
            summary["mean_spearman_per_method"][m] = float(sub["unc_error_spearman"].mean())
            summary["achieves_target"][m] = bool(
                (sub["unc_error_pearson"] >= 0.30).mean() >= 0.5)
    base_sub = valid[valid["method"] == "A_constant"]
    if len(base_sub):
        summary["baseline_constant_mean"] = float(base_sub["unc_error_pearson"].mean())
        summary["baseline_constant_median"] = float(base_sub["unc_error_pearson"].median())

    # best method by mean pearson
    means = summary["mean_corr_per_method"]
    if means:
        best = max(means, key=means.get)
        summary["best_method"] = best
        summary["best_mean_pearson"] = means[best]

    json_path = os.path.join(OUT_DIR, "summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    both(f"Wrote {json_path}")

    # ---- figure: bar chart corr by method x dataset ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ds_order = [d for (d, _, _) in SAMPLES]
        method_colors = {
            "A_constant":     "#999999",
            "B_local_gene":   "#4C72B0",
            "C_spatial_spot": "#55A868",
            "D_residual_spot":"#C44E52",
        }
        n_ds = len(ds_order)
        n_m = len(METHODS)
        x = np.arange(n_ds)
        w = 0.8 / n_m
        fig, ax = plt.subplots(figsize=(max(10, n_ds * 1.6), 6))
        for i, m in enumerate(METHODS):
            vals = []
            for d in ds_order:
                row = df[(df["dataset"] == d) & (df["method"] == m)]
                vals.append(float(row["unc_error_pearson"].iloc[0])
                            if len(row) else float("nan"))
            ax.bar(x + i * w - 0.4 + w / 2, vals, w,
                   label=m, color=method_colors[m])
        ax.hlines(0.30, x[0] - 0.5, x[-1] + 0.5,
                  colors="k", linestyles="--", linewidths=1.0,
                  label="target (0.30)")
        ax.hlines(summary.get("baseline_constant_mean", 0.07) or 0.07,
                  x[0] - 0.5, x[-1] + 0.5,
                  colors="gray", linestyles=":", linewidths=1.0,
                  label="baseline mean (A)")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{d}\n{lab}" for (d, lab, _) in SAMPLES],
                           rotation=0, fontsize=9)
        ax.set_ylabel("Uncertainty-error Pearson r")
        ax.set_title("Raw GP uncertainty-error correlation by noise method\n"
                     "(higher = uncertainty more informative for ranking quality)")
        ax.set_ylim(0, max(0.7, df["unc_error_pearson"].max() + 0.05))
        ax.legend(loc="upper right", ncol=3, fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        png_path = os.path.join(OUT_DIR, "figure.png")
        fig.savefig(png_path, dpi=140)
        both(f"Wrote {png_path}")
    except Exception as e:
        both(f"figure error: {e}")

    # ---- console summary ----
    both("")
    both("=" * 70)
    both("SUMMARY (mean Pearson r per method across datasets)")
    both("=" * 70)
    for m in METHODS:
        both(f"  {m:18s} mean={summary['mean_corr_per_method'].get(m, float('nan')):+.4f}  "
             f"median={summary['median_corr_per_method'].get(m, float('nan')):+.4f}  "
             f"spearman_mean={summary['mean_spearman_per_method'].get(m, float('nan')):+.4f}  "
             f">0.30 on >50% datasets: {summary['achieves_target'].get(m)}")
    both("")
    both(f"baseline (A constant) mean corr: {summary['baseline_constant_mean']}")
    if "best_method" in summary:
        both(f"BEST method: {summary['best_method']} "
             f"(mean r={summary['best_mean_pearson']:+.4f})")
    both("")
    both("Verdict: " + (
        "Raw GP uncertainty CAN be made informative (>0.3) via better noise "
        "estimation -> strengthens conformal story."
        if any(summary["achieves_target"].values())
        else "Raw GP uncertainty remains weakly informative; conformal "
             "provides distribution-free coverage regardless."))

    log_fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

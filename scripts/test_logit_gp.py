"""
Logit-space GP vs raw-space GP benchmark on MOB Visium (st11).

Problem
-------
The sparse GP currently regresses directly on APA usage proportions in
[0, 1] with a Gaussian likelihood.  Gaussian assumptions on bounded data
cause predictions outside [0, 1] (observed on Gdap2: 39/260 spots out of
bounds, range -0.298..1.177).  This benchmark evaluates the fix -- a GP
in latent logit space with a sigmoid inverse (``transform='logit'``) --
against the historical raw-space GP (``transform=None``).

Protocol
--------
* Data: ``data/processed/mob_st11`` (260 spots x 4845 genes).
* Genes: top-100 "multi-PAS" genes, ranked by the number of spots with
  genuinely fractional usage (0 < v < 1) -- genes where >1 PAS competes.
* For each seed: mask 20% of the observed (>0) entries per gene, fit the
  sparse GP batch on the masked matrix, predict at all spots.
* Conformal calibration is done in the ORIGINAL space: split-conformal
  scores |y_true - y_pred| (original space) standardised by the GP
  posterior std (latent space for the logit variant), calibration split
  pooled across genes, coverage/width evaluated on pooled held-out eval
  entries.

Reports (per variant x seed, CSV) and a comparison figure:
  a) original-space RMSE on held-out entries
  b) out-of-bound spot counts (<0 or >1): held-out and full field
  c) conformal coverage at 80 / 90 / 95%
  d) median interval width at each level
  e) spatial gradient recovery: per-gene Pearson r between the predicted
     and true spatial fields (mean over genes)
  f) runtime

Usage:
    python scripts/test_logit_gp.py \
        [--data-dir data/processed/mob_st11] \
        [--output benchmark_results/logit_gp] \
        [--mask-fraction 0.2] [--seeds 42 43 44 45 46] \
        [--n-top-genes 100]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def load_dataset(data_dir: Path):
    apa = pd.read_csv(data_dir / "apa_matrix.csv", index_col=0)
    coords = pd.read_csv(data_dir / "coordinates.csv", index_col=0)
    # Align spot order between matrix and coordinates.
    common = [s for s in apa.columns if s in coords.index]
    apa = apa[common]
    xy = coords.loc[common, ["x", "y"]].to_numpy(dtype=float)
    values = apa.to_numpy(dtype=float)
    genes = apa.index.to_numpy()
    return values, xy, genes


def select_multi_pas_genes(values: np.ndarray, n_top: int = 100) -> np.ndarray:
    """Indices of the top ``n_top`` multi-PAS genes.

    Ranked by the number of spots with fractional usage (0 < v < 1) --
    i.e. spots where more than one poly(A) site is actively competing.
    Ties are broken by observed-count then usage variance.
    """
    finite = np.isfinite(values)
    fractional = finite & (values > 0) & (values < 1)
    observed = finite & (values > 0)
    frac_counts = fractional.sum(axis=1)
    obs_counts = observed.sum(axis=1)
    variances = np.nanvar(np.where(finite, values, np.nan), axis=1)

    order = np.lexsort((-variances, -obs_counts, -frac_counts))
    return order[:n_top]


# ---------------------------------------------------------------------------
# Masking
# ---------------------------------------------------------------------------
def mask_observed(values: np.ndarray, mask_fraction: float, seed: int):
    """Mask ``mask_fraction`` of observed (>0) entries per gene.

    Returns the masked matrix, the post-masking training mask, and a list
    of (gene_index, held_out_spot_idx, held_out_truth) records.
    """
    rng = np.random.default_rng(seed)
    masked = values.copy()
    train_mask = np.zeros_like(values, dtype=bool)
    records = []
    for g in range(values.shape[0]):
        obs_idx = np.where(values[g] > 0)[0]
        if len(obs_idx) < 10:
            train_mask[g] = values[g] > 0
            continue
        n_mask = max(1, int(len(obs_idx) * mask_fraction))
        held = rng.choice(obs_idx, size=n_mask, replace=False)
        masked[g, held] = 0.0
        keep = np.setdiff1d(obs_idx, held)
        train_mask[g, keep] = True
        # Split held-out into calibration / evaluation halves (per gene),
        # for split-conformal calibration pooled across genes.
        perm = rng.permutation(len(held))
        n_cal = len(held) // 2
        cal_idx, eval_idx = held[perm[:n_cal]], held[perm[n_cal:]]
        records.append({
            "gene": g,
            "held": held,
            "cal": cal_idx,
            "eval": eval_idx,
            "truth": values[g, held],
        })
    return masked, train_mask, records


# ---------------------------------------------------------------------------
# GP run + metrics
# ---------------------------------------------------------------------------
def run_variant(masked, train_mask, xy, *, transform, n_inducing,
                length_scale, noise_level):
    """Fit the sparse-GP batch and return pure GP predictions + std."""
    from spagapa.imputation.sparse_gp import SparseGPImputer

    t0 = time.perf_counter()
    base = SparseGPImputer(
        n_inducing=n_inducing,
        inducing_method='kmeans',
        length_scale=length_scale,
        noise_level=noise_level,
        transform=transform,
    )
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    pred, std = batch.predict(return_std=True)
    elapsed = time.perf_counter() - t0
    return pred, std, elapsed


def conformal_metrics(records, values, pred, std, alpha):
    """Pooled split-conformal coverage/width at miscoverage ``alpha``.

    Scores: |y_true - y_pred| in the ORIGINAL space, standardised by the
    GP posterior std (latent space when transform='logit' -- the std is a
    conditioning variable only, so interval semantics stay original-space).
    """
    from spagapa.imputation.calibration import ConformalCalibrator

    cal_err, cal_std = [], []
    for rec in records:
        if len(rec["cal"]) == 0:
            continue
        g = rec["gene"]
        cal_err.append(np.abs(values[g, rec["cal"]] - pred[g, rec["cal"]]))
        cal_std.append(std[g, rec["cal"]])
    cal_err = np.concatenate(cal_err)
    cal_std = np.concatenate(cal_std)

    cal = ConformalCalibrator(alpha=alpha, mode="locally_adaptive")
    cal.fit(cal_err, cal_std)

    ev_true, ev_pred, ev_std = [], [], []
    for rec in records:
        if len(rec["eval"]) == 0:
            continue
        g = rec["gene"]
        ev_true.append(values[g, rec["eval"]])
        ev_pred.append(pred[g, rec["eval"]])
        ev_std.append(std[g, rec["eval"]])
    ev_true = np.concatenate(ev_true)
    ev_pred = np.concatenate(ev_pred)
    ev_std = np.concatenate(ev_std)

    lower, upper = cal.predict(ev_pred, ev_std)
    coverage = float(np.mean((lower <= ev_true) & (ev_true <= upper)))
    width = float(np.median(upper - lower))
    return coverage, width


# ---------------------------------------------------------------------------
def spatial_gradient_r(values, pred, records):
    """Mean per-gene Pearson r between predicted and true spatial fields."""
    rs = []
    for rec in records:
        g = rec["gene"]
        truth = values[g]
        p = pred[g]
        if np.std(truth) < 1e-12 or np.std(p) < 1e-12:
            continue
        rs.append(float(np.corrcoef(truth, p)[0, 1]))
    return float(np.mean(rs)) if rs else float("nan")


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path,
                        default=Path("data/processed/mob_st11"))
    parser.add_argument("--output", type=Path,
                        default=Path("benchmark_results/logit_gp"))
    parser.add_argument("--mask-fraction", type=float, default=0.2)
    parser.add_argument("--seeds", type=int, nargs="+",
                        default=[42, 43, 44, 45, 46])
    parser.add_argument("--n-top-genes", type=int, default=100)
    parser.add_argument("--alphas", type=float, nargs="+",
                        default=[0.2, 0.1, 0.05])
    parser.add_argument("--length-scale-multiplier", type=float, default=2.0)
    parser.add_argument("--noise-level", type=float, default=0.08)
    args = parser.parse_args()

    from scipy.spatial import cKDTree

    args.output.mkdir(parents=True, exist_ok=True)

    values_all, xy, genes = load_dataset(args.data_dir)
    print(f"Loaded APA matrix: {values_all.shape[0]} genes x "
          f"{values_all.shape[1]} spots")

    top = select_multi_pas_genes(values_all, args.n_top_genes)
    values = values_all[top]
    gene_names = genes[top]
    frac_spots = ((values > 0) & (values < 1)).sum(axis=1)
    print(f"Selected top-{len(top)} multi-PAS genes "
          f"(fractional spots: min={frac_spots.min()}, "
          f"median={int(np.median(frac_spots))}, max={frac_spots.max()})")

    # Spatial scale: median nearest-neighbour distance (shared by both arms).
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    median_nn = float(np.median(nn))
    length_scale = median_nn * args.length_scale_multiplier
    n_inducing = min(500, max(100, xy.shape[0] // 100))
    print(f"median_nn={median_nn:.3f}, length_scale={length_scale:.3f}, "
          f"n_inducing={n_inducing}, noise_level={args.noise_level}")

    variants = {"raw": None, "logit": "logit"}
    rows = []
    per_seed = {}
    for seed in args.seeds:
        masked, train_mask, records = mask_observed(
            values, args.mask_fraction, seed)
        n_held = sum(len(r["held"]) for r in records)
        print(f"\n=== seed {seed}: {len(records)} genes, {n_held} held-out "
              f"entries ===")
        per_seed[seed] = {}
        for name, transform in variants.items():
            pred, std, elapsed = run_variant(
                masked, train_mask, xy, transform=transform,
                n_inducing=n_inducing, length_scale=length_scale,
                noise_level=args.noise_level)

            # (a) original-space RMSE on held-out entries
            truths = np.concatenate([r["truth"] for r in records])
            preds = np.concatenate([pred[r["gene"], r["held"]]
                                    for r in records])
            rmse = float(np.sqrt(np.mean((truths - preds) ** 2)))

            # (b) out-of-bound counts
            oob_held = int(((preds < 0) | (preds > 1)).sum())
            oob_field = int(((pred < 0) | (pred > 1)).sum())

            row = {
                "seed": seed,
                "variant": name,
                "rmse_held": rmse,
                "oob_held": oob_held,
                "oob_field": oob_field,
                "n_held": len(truths),
                "gradient_r": spatial_gradient_r(values, pred, records),
                "runtime_s": elapsed,
            }
            # (c, d) conformal coverage / width
            for alpha in args.alphas:
                cov, width = conformal_metrics(records, values, pred, std,
                                               alpha)
                tag = f"{int(round((1 - alpha) * 100))}"
                row[f"coverage_{tag}"] = cov
                row[f"width_{tag}"] = width
            rows.append(row)
            per_seed[seed][name] = (pred, std, records)
            print(f"  [{name:5s}] RMSE={rmse:.4f}  OOB(held)={oob_held:4d}  "
                  f"OOB(field)={oob_field:4d}  "
                  f"grad_r={row['gradient_r']:.3f}  "
                  f"cov80/90/95="
                  f"{row['coverage_80']:.3f}/{row['coverage_90']:.3f}/"
                  f"{row['coverage_95']:.3f}  t={elapsed:.1f}s")

    df = pd.DataFrame(rows)
    csv_path = args.output / "logit_vs_raw_metrics.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path}")

    summary = df.groupby("variant").agg(["mean", "std"])
    summary_path = args.output / "logit_vs_raw_summary.csv"
    summary.to_csv(summary_path)
    print(f"Wrote {summary_path}")

    # ---- comparison figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    colors = {"raw": "#4C72B0", "logit": "#DD8452"}

    ax = axes[0, 0]
    for name in variants:
        sub = df[df.variant == name]
        ax.scatter(sub.seed, sub.rmse_held, color=colors[name],
                   label=name, s=60, zorder=3)
        ax.plot(sub.seed, sub.rmse_held, color=colors[name], alpha=0.6)
    ax.set_xlabel("seed")
    ax.set_ylabel("RMSE (original space, held-out)")
    ax.set_title("(a) Accuracy: RMSE")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    width = 0.35
    for k, (name) in enumerate(variants):
        sub = df[df.variant == name]
        ax.bar(np.arange(len(args.alphas)) + (k - 0.5) * width,
               [sub[f"coverage_{int(round((1-a)*100))}"].mean()
                for a in args.alphas],
               width, color=colors[name], label=name)
    ax.axhline(0.80, color="grey", ls="--", lw=1, alpha=0.7)
    ax.axhline(0.90, color="grey", ls="--", lw=1, alpha=0.7)
    ax.axhline(0.95, color="grey", ls=":", lw=1, alpha=0.7)
    ax.set_xticks(range(len(args.alphas)))
    ax.set_xticklabels([f"{int(round((1-a)*100))}%" for a in args.alphas])
    ax.set_ylim(0.5, 1.02)
    ax.set_ylabel("empirical coverage")
    ax.set_title("(c) Conformal coverage (target dashed)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1, 0]
    for k, name in enumerate(variants):
        sub = df[df.variant == name]
        meds = [sub[f"width_{int(round((1-a)*100))}"].mean()
                for a in args.alphas]
        errs = [sub[f"width_{int(round((1-a)*100))}"].std()
                for a in args.alphas]
        ax.bar(np.arange(len(args.alphas)) + (k - 0.5) * width, meds, width,
               yerr=errs, capsize=3, color=colors[name], label=name)
    ax.set_xticks(range(len(args.alphas)))
    ax.set_xticklabels([f"{int(round((1-a)*100))}%" for a in args.alphas])
    ax.set_ylabel("median interval width")
    ax.set_title("(d) Interval width (mean +/- sd over seeds)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1, 1]
    for name in variants:
        sub = df[df.variant == name]
        ax.scatter(sub.seed, sub.gradient_r, color=colors[name],
                   label=name, s=60, zorder=3)
        ax.plot(sub.seed, sub.gradient_r, color=colors[name], alpha=0.6)
    ax.set_xlabel("seed")
    ax.set_ylabel("mean Pearson r (pred field vs true field)")
    ax.set_title("(e) Spatial gradient recovery")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Raw-space GP vs logit-space GP -- MOB st11, "
        f"top-{len(top)} multi-PAS genes, {args.mask_fraction:.0%} masked, "
        f"{len(args.seeds)} seeds", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig_path = args.output / "logit_vs_raw_comparison.png"
    fig.savefig(fig_path, dpi=150)
    print(f"Wrote {fig_path}")

    # ---- acceptance snapshot ----
    logit = df[df.variant == "logit"]
    raw = df[df.variant == "raw"]
    rmse_ratio = (logit.rmse_held.mean() / raw.rmse_held.mean())
    acceptance = {
        "logit_all_within_unit_interval":
            bool((logit.oob_field == 0).all()),
        "rmse_ratio_logit_over_raw": float(rmse_ratio),
        "rmse_within_10pct": bool(rmse_ratio <= 1.10),
        "coverage_90_mean_logit": float(logit.coverage_90.mean()),
        "coverage_90_mean_raw": float(raw.coverage_90.mean()),
        "gradient_r_mean_logit": float(logit.gradient_r.mean()),
        "gradient_r_mean_raw": float(raw.gradient_r.mean()),
        "oob_field_mean_raw": float(raw.oob_field.mean()),
        "oob_field_mean_logit": float(logit.oob_field.mean()),
    }
    acc_path = args.output / "acceptance.json"
    acc_path.write_text(json.dumps(acceptance, indent=2))
    print(f"Wrote {acc_path}")
    print(json.dumps(acceptance, indent=2))


if __name__ == "__main__":
    main()

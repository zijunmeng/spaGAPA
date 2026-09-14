#!/usr/bin/env python3
"""Parametric simulation benchmark (S19) for spaGAPA sparse-GP imputation.

Sweeps a parameter grid over synthetic spatial APA datasets with known ground
truth and measures how imputation accuracy, conformal interval coverage, and
spatial-pattern fidelity scale with sequencing depth (number of spots),
observed fraction, and the observation-noise model.

Grid (default)
--------------
* depth        : 5,000 / 20,000 / 80,000 spots  ("sequencing depth" proxy)
* observed frac: 0.05 / 0.15 / 0.35
* noise model  : ``constant`` (iid Gaussian, sigma=0.05 on observed entries)
                 ``residual_spot`` (per-spot heteroscedastic sigma drawn from
                 a smooth spatial log-Gaussian field -- noisy tissue regions)
* seeds        : 42 / 123 / 456            -> 3*3*2*3 = 54 conditions

Data generation
---------------
The existing ``benchmark/simulator.py`` interface has no ``observed_fraction``
masking and no per-spot noise model, so truth fields are constructed directly
with numpy: near-square integer grid (Visium-like), per-gene spatial patterns
(``autocorrelated`` random-Fourier-feature GP fields, 4-quadrant ``domain``
blocks, ``gradient``, and ``constant`` genes) whose spatial scales are
proportional to the tissue side, micro-heterogeneity sigma=0.02,
clipped to [0,1] like a distal-usage APA index.  Truth is seeded by
``(seed, n_spots)`` only, so conditions sharing depth+seed see the *same*
truth field (paired comparisons across fraction / noise model).

Imputation & calibration
------------------------
``SparseGPImputer`` with ``length_scale = 5 x median NN distance``,
``n_inducing = min(500, max(100, n_spots // 100))`` and ``noise_level=0.1``
(identical to the real-data calibration protocol).  The GP is fit on ALL
observed entries; prediction / metric evaluation uses at most
``PREDICT_CAP`` = 20,000 spots (seeded uniform subsample) -- ~1M held-out
entries give negligible metric noise while capping the runtime of the
80k-spot conditions.  Missing entries are split 50/50 into a conformal
calibration set and an independent test set; ``ConformalCalibrator``
(global + locally_adaptive) produces 80/90/95% intervals.

Outputs (in --out-dir, default pipeline_output/simulation_benchmark)
--------------------------------------------------------------------
* ``simulation_grid_results.csv``  -- one row per condition
* ``power_vs_depth.png``           -- performance vs depth (lines = observed
                                      fraction, linestyle = noise model)
* ``simulation_summary.json``      -- config, aggregates, key findings

Usage
-----
OPENBLAS_NUM_THREADS=8 ~/anaconda3/envs/spagapa/bin/python \
    scripts/simulation_benchmark.py [--workers 8] [--smoke]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

PREDICT_CAP = 20_000  # GP prediction / metric evaluation subsample (see run_condition)
# Make the repo root importable when run as `python scripts/<this>.py`.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from spagapa.imputation import SparseGPImputer  # noqa: E402
from spagapa.imputation.calibration import (  # noqa: E402
    ConformalCalibrator,
    evaluate_coverage,
)
from spagapa.analysis.svapa import morans_i  # noqa: E402

import numpy as np  # noqa: E402

# Pattern mixture for the truth fields (fractions of n_genes).
PATTERN_FRACTIONS = {
    "autocorrelated": 0.30,  # smooth GP field (random Fourier features)
    "domain": 0.20,          # 4-quadrant tissue blocks
    "gradient": 0.10,        # linear gradient, random direction
    "constant": 0.40,        # spatially flat genes
}
CONSTANT_SIGMA = 0.02   # micro-heterogeneity added to every truth field
NOISE_SIGMA = 0.05      # constant-model observation noise
NOISE_LEVELS = (0.80, 0.90, 0.95)
Z_LEVELS = {0.80: 1.2816, 0.90: 1.6449, 0.95: 1.96}


# --------------------------------------------------------------------------
# Data generation
# --------------------------------------------------------------------------
def _rff_gp_field(coords: np.ndarray, rng: np.random.Generator,
                  length_scale: float, n_features: int = 200) -> np.ndarray:
    """One draw from a zero-mean GP with RBF covariance via random Fourier
    features: f(x) = sqrt(2/M) * sum_k cos(w_k . x + b_k).  Exact in the
    M -> inf limit; unit marginal variance."""
    d = coords.shape[1]
    W = rng.normal(0.0, 1.0 / length_scale, size=(d, n_features))
    b = rng.uniform(0.0, 2.0 * np.pi, size=n_features)
    return np.sqrt(2.0 / n_features) * np.cos(coords @ W + b).sum(axis=1)


def _sigmoid(z: np.ndarray) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh(0.5 * z))


def generate_truth(n_spots: int, n_genes: int, seed: int):
    """Ground-truth APA fields on a near-square integer grid.

    Seeded by (seed, n_spots) ONLY, so every condition at the same depth and
    seed shares the identical truth field (paired grid comparisons).
    Returns (truth (n_genes, n_spots) in [0,1], coords (n_spots, 2),
    patterns (n_genes,) str array).
    """
    rng = np.random.default_rng([int(seed), int(n_spots), 1])
    side = int(np.ceil(np.sqrt(n_spots)))
    rows, cols = np.unravel_index(np.arange(n_spots), (side, side))
    coords = np.column_stack([cols, rows]).astype(float)
    center = np.array([(side - 1) / 2.0, (side - 1) / 2.0])
    quadrant = ((cols < side // 2).astype(int)
                + (rows < side // 2).astype(int))  # 0..3

    # pattern assignment (shuffle for a random layout)
    names = list(PATTERN_FRACTIONS)
    counts = {k: int(round(n_genes * f)) for k, f in PATTERN_FRACTIONS.items()}
    counts["constant"] = n_genes - sum(v for k, v in counts.items()
                                       if k != "constant")
    patterns = np.concatenate([[k] * counts[k] for k in names])
    rng.shuffle(patterns)

    truth = np.empty((n_genes, n_spots))
    for gi, pat in enumerate(patterns):
        if pat == "autocorrelated":
            # spatial scale proportional to tissue side -> the field's
            # effective DOF is depth-invariant, so depth isolates sampling
            # density rather than inducing-point capacity
            ls = side * rng.uniform(0.05, 0.12)
            z = _rff_gp_field(coords, rng, ls)
            field = 0.15 + 0.70 * _sigmoid(1.6 * z)
        elif pat == "domain":
            means = rng.uniform(0.1, 0.9, size=4)
            field = means[quadrant] + 0.05 * _rff_gp_field(
                coords, rng, max(2.0, 0.03 * side))
        elif pat == "gradient":
            theta = rng.uniform(0.0, 2.0 * np.pi)
            u = np.array([np.cos(theta), np.sin(theta)])
            proj = (coords - center) @ u
            t = 0.5 * (1.0 + proj / (np.abs(proj).max() + 1e-9))
            field = 0.15 + 0.70 * t
        else:  # constant
            field = np.full(n_spots, rng.uniform(0.2, 0.8))
        field = field + rng.normal(0.0, CONSTANT_SIGMA, n_spots)
        truth[gi] = np.clip(field, 0.0, 1.0)
    return truth, coords, patterns


def _spot_sigma_field(coords: np.ndarray, n_spots: int, seed: int) -> np.ndarray:
    """Per-spot heteroscedastic noise scale for the ``residual_spot`` model:
    smooth log-Gaussian field -> sigma roughly in [0.01, 0.12].  Seeded by
    (seed, n_spots) so the noisy tissue regions are shared across fractions."""
    rng = np.random.default_rng([int(seed), int(n_spots), 777])
    log_sigma = np.log(0.03) + 0.9 * _rff_gp_field(coords, rng, 8.0)
    return np.exp(log_sigma)


def generate_observation(truth: np.ndarray, coords: np.ndarray,
                         observed_fraction: float, noise_model: str, seed: int):
    """Sparse, noisy observation of the truth fields.

    Returns (values (n_genes, n_spots) with 0.0 at missing entries,
    obs_mask bool, sigma_spot (n_spots,) or None).
    """
    n_genes, n_spots = truth.shape
    noise_id = 0 if noise_model == "constant" else 1
    rng = np.random.default_rng(
        [int(seed), int(n_spots), int(round(observed_fraction * 1000)), noise_id])
    obs_mask = rng.random((n_genes, n_spots)) < observed_fraction

    if noise_model == "constant":
        sigma_spot = None
        noisy = truth + rng.normal(0.0, NOISE_SIGMA, truth.shape)
    elif noise_model == "residual_spot":
        sigma_spot = _spot_sigma_field(coords, n_spots, seed)
        noisy = truth + sigma_spot[None, :] * rng.normal(0.0, 1.0, truth.shape)
    else:
        raise ValueError(f"unknown noise_model {noise_model!r}")

    values = np.where(obs_mask, np.clip(noisy, 0.0, 1.0), 0.0)
    return values, obs_mask, sigma_spot


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------
def _safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) > 0 and np.std(b) > 0:
        return float(np.corrcoef(a, b)[0, 1])
    return float("nan")


def _spatial_fidelity(truth, imputed_full, values_sparse, obs_mask, coords,
                      patterns, seed, max_spots=20000):
    """Moran's-I fidelity for structured genes: correlation between the
    truth and reconstructed spatial-autocorrelation profiles.  Spots are
    subsampled above ``max_spots`` (identically for all three views) to keep
    the kNN weight matrix cheap.  ``morans_i`` drops NaN entries, so the
    sparse view (missing -> NaN) is the no-imputation baseline."""
    n_spots = coords.shape[0]
    if n_spots > max_spots:
        sub = np.random.default_rng([int(seed), 555]).choice(
            n_spots, size=max_spots, replace=False)
        coords_s, truth_s = coords[sub], truth[:, sub]
        imp_s, sparse_s = imputed_full[:, sub], values_sparse[:, sub]
        mask_s = obs_mask[:, sub]
    else:
        coords_s, truth_s = coords, truth
        imp_s, sparse_s, mask_s = imputed_full, values_sparse, obs_mask

    structured = np.where(patterns != "constant")[0]
    i_truth, i_imp, i_obs = [], [], []
    for gi in structured:
        i_truth.append(morans_i(truth_s[gi], coords_s, k=8, n_perm=0)[0])
        i_imp.append(morans_i(imp_s[gi], coords_s, k=8, n_perm=0)[0])
        obs_vals = np.where(mask_s[gi], sparse_s[gi], np.nan)
        i_obs.append(morans_i(obs_vals, coords_s, k=8, n_perm=0)[0])
    i_truth, i_imp, i_obs = map(np.asarray, (i_truth, i_imp, i_obs))
    return {
        "morans_r_imputed": _safe_pearson(i_truth, i_imp),
        "morans_r_observed": _safe_pearson(i_truth, i_obs),
        "morans_mean_abs_delta": float(np.nanmean(np.abs(i_imp - i_truth))),
        "morans_mean_truth": float(np.nanmean(i_truth)),
        "n_structured_genes": int(len(structured)),
    }


# --------------------------------------------------------------------------
# One grid condition
# --------------------------------------------------------------------------
def run_condition(cfg: dict) -> dict:
    t0 = time.time()
    n_spots = cfg["n_spots"]
    frac = cfg["observed_fraction"]
    noise_model = cfg["noise_model"]
    seed = cfg["seed"]
    n_genes = cfg["n_genes"]

    truth, coords, patterns = generate_truth(n_spots, n_genes, seed)
    values, obs_mask, _sigma = generate_observation(
        truth, coords, frac, noise_model, seed)

    # ---- sparse GP imputation (same recipe as scripts/calibrate_uncertainty.py).
    # The GP is FIT on every observed entry; the prediction / metric
    # evaluation is capped at ``PREDICT_CAP`` spots (seeded uniform
    # subsample) because predicting at all 80k coordinates dominates
    # runtime while ~1M held-out entries already give negligible metric
    # noise.
    t_gp = time.time()
    from scipy.spatial import cKDTree
    nn = cKDTree(coords).query(coords, k=2)[0][:, 1]
    nn_dist = float(np.median(nn))
    length_scale = nn_dist * 5
    n_inducing = min(500, max(100, n_spots // 100))
    base = SparseGPImputer(n_inducing=n_inducing, length_scale=length_scale,
                           noise_level=0.1)
    batch = base.fit_batch(coords, values, mask=obs_mask, verbose=False)

    if n_spots > PREDICT_CAP:
        sub = np.sort(np.random.default_rng(
            [int(seed), int(n_spots), 999]).choice(
                n_spots, PREDICT_CAP, replace=False))
        coords_ev, truth_ev = coords[sub], truth[:, sub]
        values_ev, mask_ev = values[:, sub], obs_mask[:, sub]
    else:
        coords_ev, truth_ev, values_ev, mask_ev = (
            coords, truth, values, obs_mask)
    pred_ev, unc_ev = batch.predict(coords_ev, return_std=True)
    # complete-map view: stamp observed values back (as batch.impute does)
    pred_ev[mask_ev] = values_ev[mask_ev]
    gp_s = time.time() - t_gp

    # ---- gather missing entries; split 50/50 calibration / test
    g_idx, s_idx = np.where(~mask_ev)
    truths = truth_ev[g_idx, s_idx]
    preds = pred_ev[g_idx, s_idx]
    stds = unc_ev[g_idx, s_idx]
    errors = np.abs(truths - preds)

    rng = np.random.default_rng([int(seed), int(n_spots),
                                 int(round(frac * 1000)),
                                 0 if noise_model == "constant" else 1, 7])
    perm = rng.permutation(len(truths))
    n_cal = int(len(truths) * 0.5)
    cal_idx, test_idx = perm[:n_cal], perm[n_cal:]
    cal_err, cal_std = errors[cal_idx], stds[cal_idx]
    t_truth, t_pred, t_std, t_err = (truths[test_idx], preds[test_idx],
                                     stds[test_idx], errors[test_idx])
    n_test_spots = len(coords_ev)

    res = {
        "n_spots": n_spots, "observed_fraction": frac,
        "noise_model": noise_model, "seed": seed, "n_genes": n_genes,
        "n_observed_per_gene": float(obs_mask.sum(axis=1).mean()),
        "gp_length_scale": length_scale, "gp_n_inducing": n_inducing,
        "n_eval_spots": int(n_test_spots),
        "n_calibration": int(n_cal), "n_test": int(len(test_idx)),
        "unc_error_pearson_cal": _safe_pearson(cal_std, cal_err),
        "rmse": float(np.sqrt(np.mean(t_err ** 2))),
        "pearson": _safe_pearson(t_pred, t_truth),
        "gp_fit_impute_s": float(gp_s),
    }

    # ---- coverage: raw GP std vs conformal (global + locally_adaptive)
    for level, z in Z_LEVELS.items():  # raw GP intervals
        half = z * t_std
        res[f"coverage_raw_{int(level*100)}"] = evaluate_coverage(
            t_pred - half, t_pred + half, t_truth)
    for mode in ("global", "locally_adaptive"):
        for alpha, level in [(0.20, 80), (0.10, 90), (0.05, 95)]:
            cal = ConformalCalibrator(alpha=alpha, mode=mode)
            cal.fit(cal_err, cal_std if mode == "locally_adaptive" else None)
            lo, hi = cal.predict(t_pred, t_std if mode == "locally_adaptive"
                                 else None)
            res[f"coverage_{mode}_{level}"] = evaluate_coverage(lo, hi, t_truth)
            if level == 90:
                res[f"q_hat_{mode}_90"] = cal.interval_.q_hat
                res[f"width_{mode}_90"] = float(np.mean(hi - lo))

    t_mor = time.time()
    fid = _spatial_fidelity(truth_ev, pred_ev, values_ev, mask_ev,
                            coords_ev, patterns, seed)
    res.update(fid)
    res["morans_s"] = float(time.time() - t_mor)
    res["total_s"] = float(time.time() - t0)
    return res


def make_power_plot(df, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    depths = sorted(df["n_spots"].unique())
    fracs = sorted(df["observed_fraction"].unique())
    noises = sorted(df["noise_model"].unique())
    frac_colors = {f: c for f, c in zip(fracs, ["#1f77b4", "#2ca02c", "#d62728"])}
    noise_ls = {n: s for n, s in zip(noises, ["-", "--"])}

    panels = [
        ("rmse", "RMSE (test missing entries)", False),
        ("pearson", "Pearson r (pred vs truth)", True),
        ("coverage_locally_adaptive_90", "90% conformal coverage (locally adaptive)", True),
        ("morans_r_imputed", "Spatial fidelity (Moran's I corr)", True),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    for ax, (metric, title, _) in zip(axes.ravel(), panels):
        for noise in noises:
            for frac in fracs:
                sub = df[(df["noise_model"] == noise)
                         & (df["observed_fraction"] == frac)]
                g = sub.groupby("n_spots")[metric].mean().reindex(depths)
                ax.plot(depths, g.values, marker="o",
                        ls=noise_ls[noise], color=frac_colors[frac],
                        label=f"frac={frac}, {noise}")
        ax.set_xscale("log")
        ax.set_xticks(depths)
        ax.set_xticklabels([f"{d//1000}k" for d in depths])
        ax.set_xlabel("Sequencing depth (spots)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, ncol=2)
    # target line on the coverage panel
    axes[0, 1].axhline(0.90, color="k", lw=1, ls=":", alpha=0.8)
    axes[0, 1].text(depths[0], 0.902, " target 0.90", fontsize=8)
    fig.suptitle("spaGAPA simulation benchmark: performance vs sequencing depth\n"
                 "(mean over 3 seeds; solid = constant noise, dashed = residual-spot noise)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def build_summary(df, cfg_grid: dict, out_files: dict, wall_s: float) -> dict:
    def agg(keys):
        g = df.groupby(keys)
        cols = ["rmse", "pearson", "coverage_raw_90",
                "coverage_global_90", "coverage_locally_adaptive_90",
                "morans_r_imputed", "morans_r_observed"]
        return json.loads(g[cols].mean().round(4).to_json(orient="index"))

    by_depth = agg(["n_spots"])
    by_frac = agg(["observed_fraction"])
    by_noise = agg(["noise_model"])
    by_depth_frac = agg(["n_spots", "observed_fraction"])

    depths_sorted = sorted(df["n_spots"].unique())
    d_lo, d_hi = depths_sorted[0], depths_sorted[-1]
    rmse_gain = (1.0 - by_depth[str(d_hi)]["rmse"] / by_depth[str(d_lo)]["rmse"])
    pearson_gain = (by_depth[str(d_hi)]["pearson"]
                    - by_depth[str(d_lo)]["pearson"])
    cov_dev = float((df["coverage_locally_adaptive_90"] - 0.90).abs().mean())
    fid_gain = float((df["morans_r_imputed"]
                      - df["morans_r_observed"]).mean())

    return {
        "experiment": "S19 parametric simulation benchmark (sparse GP + conformal)",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "grid": cfg_grid,
        "n_conditions": int(len(df)),
        "pattern_fractions": PATTERN_FRACTIONS,
        "constant_noise_sigma": NOISE_SIGMA,
        "method_config": {
            "imputer": "SparseGPImputer (RBF, kmeans inducing)",
            "length_scale": "5 x median NN distance (per condition)",
            "n_inducing": "min(500, max(100, n_spots // 100))",
            "noise_level": 0.1,
            "fit_on": "all observed entries",
            "prediction_eval_subsample_cap": PREDICT_CAP,
            "conformal": {
                "method": "split conformal (50/50 cal/test of missing entries)",
                "modes": ["global", "locally_adaptive"],
                "levels": [80, 90, 95],
            },
        },
        "aggregates": {
            "by_depth": by_depth, "by_observed_fraction": by_frac,
            "by_noise_model": by_noise, "by_depth_fraction": by_depth_frac,
        },
        "key_findings": {
            "rmse_improvement_low_to_high_depth": round(float(rmse_gain), 4),
            "pearson_gain_low_to_high_depth": round(float(pearson_gain), 4),
            "mean_abs_coverage_deviation_from_90pct_locally_adaptive":
                round(cov_dev, 4),
            "mean_morans_fidelity_gain_over_observed_baseline":
                round(fid_gain, 4),
        },
        "outputs": {k: str(v) for k, v in out_files.items()},
        "total_wall_s": round(float(wall_s), 1),
        "mean_condition_wall_s": round(float(df["total_s"].mean()), 1),
    }


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default="pipeline_output/simulation_benchmark")
    ap.add_argument("--n-genes", type=int, default=50)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--depths", default="5000,20000,80000")
    ap.add_argument("--fracs", default="0.05,0.15,0.35")
    ap.add_argument("--seeds", default="42,123,456")
    ap.add_argument("--noise", default="constant,residual_spot")
    ap.add_argument("--resume", action="store_true",
                    help="skip conditions already recorded (without error) in "
                         "the output dir's conditions_progress.jsonl")
    ap.add_argument("--smoke", action="store_true",
                    help="single small condition for pipeline validation")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.smoke:
        grid = [{"n_spots": 5000, "observed_fraction": 0.15,
                 "noise_model": "constant", "seed": 42,
                 "n_genes": args.n_genes}]
    else:
        depths = [int(x) for x in args.depths.split(",")]
        fracs = [float(x) for x in args.fracs.split(",")]
        seeds = [int(x) for x in args.seeds.split(",")]
        noises = args.noise.split(",")
        grid = [{"n_spots": d, "observed_fraction": f, "noise_model": n,
                 "seed": s, "n_genes": args.n_genes}
                for d in depths for f in fracs for n in noises for s in seeds]
    t0 = time.time()
    rows, done = [], 0
    jsonl = out_dir / "conditions_progress.jsonl"
    todo = grid
    if args.resume and jsonl.exists():
        finished = set()
        prior = []
        for line in jsonl.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # torn final line from a killed run
            key = (r.get("n_spots"), r.get("observed_fraction"),
                   r.get("noise_model"), r.get("seed"))
            if "error" in r:
                continue
            if key not in finished:
                finished.add(key)
                prior.append(r)
        rows = prior
        done = len(prior)
        todo = [c for c in grid
                if (c["n_spots"], c["observed_fraction"],
                    c["noise_model"], c["seed"]) not in finished]
        print(f"[resume] {len(prior)} conditions already done, "
              f"{len(todo)} to run", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_condition, cfg): cfg for cfg in todo}
        for fut in as_completed(futs):
            cfg = futs[fut]
            try:
                row = fut.result()
            except Exception as exc:  # keep going; record the failure
                row = {**cfg, "error": f"{type(exc).__name__}: {exc}"}
                print(f"[FAIL] {cfg}: {exc}", flush=True)
            rows.append(row)
            with open(jsonl, "a") as f:
                f.write(json.dumps(row) + "\n")
            done += 1
            if "error" not in row:
                print(f"[{done}/{len(grid)}] spots={cfg['n_spots']} "
                      f"frac={cfg['observed_fraction']} "
                      f"noise={cfg['noise_model']} seed={cfg['seed']} | "
                      f"RMSE={row['rmse']:.4f} r={row['pearson']:.3f} "
                      f"cov90(local)={row['coverage_locally_adaptive_90']:.3f} "
                      f"MoranR={row['morans_r_imputed']:.3f} "
                      f"({row['total_s']:.0f}s)", flush=True)
    wall_s = time.time() - t0


    import pandas as pd
    df = pd.DataFrame(rows)
    csv_path = out_dir / "simulation_grid_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"[out] CSV -> {csv_path}", flush=True)

    if args.smoke:
        print(json.dumps(rows[0], indent=2, default=str))
        return

    ok = df[df.get("error").isna()] if "error" in df.columns else df
    png_path = out_dir / "power_vs_depth.png"
    make_power_plot(ok, png_path)
    print(f"[out] PNG -> {png_path}", flush=True)

    cfg_grid = {"depths": sorted(ok["n_spots"].unique().tolist()),
                "observed_fractions": sorted(ok["observed_fraction"].unique().tolist()),
                "noise_models": sorted(ok["noise_model"].unique().tolist()),
                "seeds": sorted(ok["seed"].unique().tolist()),
                "n_genes": args.n_genes,
                "n_conditions_ok": int(len(ok)),
                "n_conditions_failed": int(len(df) - len(ok))}
    summary = build_summary(ok, cfg_grid,
                            {"grid_results_csv": csv_path, "power_plot": png_path},
                            wall_s)
    json_path = out_dir / "simulation_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[out] JSON -> {json_path}", flush=True)
    print(f"[done] {len(ok)}/{len(grid)} conditions in {wall_s/60:.1f} min",
          flush=True)


if __name__ == "__main__":
    main()

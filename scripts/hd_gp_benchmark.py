#!/usr/bin/env python3
"""HD-scale sparse-GP + conformal benchmark on Visium HD adrenal scAPAtrap data.

External review asked for a REAL sparse-GP (+ conformal) benchmark at
HD scale, replacing the per-peak-mean-only numbers in
``hd_conformal_quick.json``.

Data
----
``pipeline_output/visium_hd/test_adrenal/scapatrap_raw/apa_site_counts_full.csv.gz``
long table (peak_id, spot_id, count): 887 peaks x 2,006,707 barcodes,
6.1M nonzero rows (median 2 APA counts per barcode).

Protocol (per scale N in 5k / 20k / 50k / 100k barcodes, independent
sampling, seed 42)
---------------------------------------------------------------
1. Sample N barcodes uniformly from the barcodes present in the table.
2. Build the top-300-peaks (by global total count) x N count matrix,
   normalise to per-barcode usage fractions
   (usage = peak_count / barcode_total over all 887 peaks).
3. Keep peaks with >= 10 observed (>0) entries (mask/eval viability).
4. Mask 20% of observed entries per peak (seed 42); split held-out
   entries 50/50 into conformal calibration / evaluation (Phase-1.1
   protocol from ``scripts/test_logit_gp.py``).
5. Predict masked entries with
     * spaGAPA sparse GP (n_inducing=100, transform='logit'), the
       pipeline default length scale (median NN dist x 5);
     * the same GP with the length scale matched to the inducing-point
       spacing (2x median distance-to-nearest-inducing) -- diagnoses
       whether the default degenerates at HD density;
     * per-peak training mean (the reviewer's baseline).
   HD barcodes have no spatial coordinates yet (barcode bridging is
   unresolved), so the GP runs on random 2D pseudo-coordinates.  This
   preserves the RMSE / coverage / runtime / memory question (they do
   not need real geometry) but means NO spatial-skill claim is made:
   under pseudo-coordinates the GP can at best match the mean.
6. Split-conformal intervals at 80/90/95% (pooled across peaks);
   GP uses the locally-adaptive mode (posterior std), mean uses the
   global mode (it has no uncertainty).

Outputs -> ``pipeline_output/hd_gp_benchmark/``
    hd_gp_vs_mean.csv    long format: n_barcodes x method x metric x value
    hd_gp_benchmark.png  comparison figure
    hd_gp_summary.json   config + full results + conclusions
"""
from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
SRC = (PKG / "pipeline_output" / "visium_hd" / "test_adrenal" / "scapatrap_raw"
       / "apa_site_counts_full.csv.gz")
OUT = PKG / "pipeline_output" / "hd_gp_benchmark"

SCALES = [5_000, 20_000, 50_000, 100_000]
N_TOP_PEAKS = 300
MASK_FRACTION = 0.2
SEED = 42
MIN_OBS = 10
ALPHAS = [(0.05, 95), (0.10, 90), (0.20, 80)]  # (alpha, level %)


def rss_mb() -> float:
    """Process peak RSS in MB (Linux ru_maxrss is in KB)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def load_long_table():
    df = pd.read_csv(SRC)
    peak_codes, peaks = pd.factorize(df["peak_id"])
    spot_codes, spots = pd.factorize(df["spot_id"])
    cnt = df["count"].to_numpy(np.int64)
    del df
    peak_tot = np.bincount(peak_codes, weights=cnt, minlength=len(peaks))
    spot_tot = np.bincount(spot_codes, weights=cnt, minlength=len(spots))
    return peak_codes, spot_codes, cnt, peaks, spots, peak_tot, spot_tot


def build_usage_matrix(peak_codes, spot_codes, cnt, peak_tot, spot_tot, top,
                       n_barcodes):
    """Sample n_barcodes barcodes (seed 42) -> usage matrix (n_peaks x N)."""
    rng = np.random.default_rng(SEED)
    sampled = rng.choice(len(spot_tot), size=n_barcodes, replace=False)
    keep = np.isin(spot_codes, sampled)
    sc, pc, cc = spot_codes[keep], peak_codes[keep], cnt[keep]

    col_of = np.full(len(spot_tot), -1, dtype=np.int64)
    col_of[sampled] = np.arange(n_barcodes)
    row_of = np.full(peak_tot.shape[0], -1, dtype=np.int64)
    row_of[top] = np.arange(len(top))

    M = np.zeros((len(top), n_barcodes), dtype=np.float64)
    sel = row_of[pc] >= 0  # keep only rows belonging to the top peaks
    M[row_of[pc[sel]], col_of[sc[sel]]] = cc[sel]
    usage = M / spot_tot[sampled][None, :]
    observed = M > 0
    # pseudo-coordinates (HD barcodes have no mapped coordinates yet)
    xy = np.random.default_rng(SEED).random((n_barcodes, 2))
    return usage, observed, xy, sampled


def mask_matrix(observed, min_obs=MIN_OBS):
    """Phase-1.1 protocol: 20% mask of observed per peak, 50/50 cal/eval."""
    rng = np.random.default_rng(SEED)
    n_peaks = observed.shape[0]
    train_mask = np.zeros_like(observed)
    records = []
    for g in range(n_peaks):
        obs_idx = np.where(observed[g])[0]
        if len(obs_idx) < min_obs:
            continue
        n_mask = max(1, int(len(obs_idx) * MASK_FRACTION))
        held = rng.choice(obs_idx, size=n_mask, replace=False)
        keep = np.setdiff1d(obs_idx, held)
        train_mask[g, keep] = True
        perm = rng.permutation(len(held))
        n_cal = len(held) // 2
        records.append({"gene": g, "cal": held[perm[:n_cal]],
                        "eval": held[perm[n_cal:]]})
    return train_mask, records




def run_gp(usage, train_mask, xy, *, length_scale=None, multiplier=None):
    """Fit sparse-GP batch on train entries; predict the FULL surface.

    Returns (pred, std, fit_s, predict_s).
    """
    from spagapa.imputation.sparse_gp import SparseGPImputer

    masked = np.where(train_mask, usage, np.nan)
    base = SparseGPImputer(
        n_inducing=100, inducing_method="kmeans",
        length_scale=1.0 if length_scale is None else length_scale,
        noise_level=0.1, transform="logit",
        **({"length_scale_multiplier": multiplier} if multiplier is not None else {}),
    )
    t0 = time.perf_counter()
    batch = base.fit_batch(xy, masked, mask=train_mask, verbose=False)
    t1 = time.perf_counter()
    pred, std = batch.predict(return_std=True)
    t2 = time.perf_counter()
    return pred, std, t1 - t0, t2 - t1


def run_mean(usage, train_mask):
    """Per-peak mean of training entries, broadcast to all spots."""
    t0 = time.perf_counter()
    n_peaks = usage.shape[0]
    pred = np.zeros_like(usage)
    for g in range(n_peaks):
        pred[g] = usage[g, train_mask[g]].mean() if train_mask[g].any() else np.nan
    return pred, time.perf_counter() - t0


def conformal(usage, records, pred, std, mode):
    """Pooled split-conformal coverage + median width at 80/90/95%.

    ``mode='locally_adaptive'`` requires the per-entry GP ``std`` (score =
    |y - pred| / max(std, floor)); ``mode='global'`` is plain absolute
    residual conformal (used for the mean baseline, which has no
    uncertainty).
    """
    from spagapa.imputation.calibration import ConformalCalibrator

    if mode == "locally_adaptive":
        cal_t, cal_p, cal_s = gather(usage, records, pred, std, "cal")
        ev_t, ev_p, ev_s = gather(usage, records, pred, std, "eval")
    else:
        cal_t, cal_p = gather(usage, records, pred, std, "cal")
        ev_t, ev_p = gather(usage, records, pred, std, "eval")
    out = {}
    for alpha, level in ALPHAS:
        cal = ConformalCalibrator(alpha=alpha, mode=mode)
        if mode == "locally_adaptive":
            cal.fit(np.abs(cal_t - cal_p), cal_s)
            lo, hi = cal.predict(ev_p, ev_s)
        else:
            cal.fit(np.abs(cal_t - cal_p))
            lo, hi = cal.predict(ev_p)
        out[f"coverage_{level}"] = float(np.mean((lo <= ev_t) & (ev_t <= hi)))
        out[f"width_median_{level}"] = float(np.median(hi - lo))
    return out


def gather(usage, records, pred, std, split):
    """Pooled (truth, pred[, std]) over the cal or eval split."""
    t, p, s = [], [], []
    for rec in records:
        idx = rec[split]
        if len(idx) == 0:
            continue
        g = rec["gene"]
        t.append(usage[g, idx])
        p.append(pred[g, idx])
        if std is not None:
            s.append(std[g, idx])
    out = [np.concatenate(t), np.concatenate(p)]
    if std is not None:
        out.append(np.concatenate(s))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"loading {SRC.name} ...")
    t0 = time.perf_counter()
    peak_codes, spot_codes, cnt, peaks, spots, peak_tot, spot_tot = load_long_table()
    top = np.argsort(peak_tot)[::-1][:N_TOP_PEAKS]
    print(f"loaded in {time.perf_counter() - t0:.0f}s | {len(peaks)} peaks | "
          f"{len(spots)} barcodes | top-peak totals "
          f"{peak_tot[top[0]]:.0f}..{peak_tot[top[-1]]:.0f}")

    rows = []
    summary = {"config": {
        "source": str(SRC), "n_peaks_total": int(len(peaks)),
        "n_barcodes_total": int(len(spots)), "n_nonzero_rows": int(len(cnt)),
        "n_top_peaks": int(N_TOP_PEAKS), "scales": SCALES,
        "mask_fraction": MASK_FRACTION, "seed": SEED, "min_obs_per_peak": MIN_OBS,
        "gp": {"n_inducing": 100, "transform": "logit", "noise_level": 0.1,
               "default_length_scale": "median_NN_dist * 5 (pipeline default)",
               "adaptive_length_scale": "2 * median_dist_to_nearest_inducing"},
        "coordinates": ("random 2D pseudo-coordinates in [0,1]^2 (HD barcode "
                        "bridging unresolved); RMSE/coverage/runtime/memory are "
                        "coordinate-independent, no spatial-skill claim"),
        "conformal": "pooled split-conformal; GP=locally_adaptive(std), mean=global",
    }, "scales": {}}

    rss_base = rss_mb()
    for N in SCALES:
        print(f"\n=== scale {N} barcodes ===")
        t_scale = time.perf_counter()
        usage, observed, xy, sampled = build_usage_matrix(
            peak_codes, spot_codes, cnt, peak_tot, spot_tot, top, N)
        # restrict to evaluable peaks
        obs_per_peak = observed.sum(axis=1)
        keep_peaks = obs_per_peak >= MIN_OBS
        usage, observed = usage[keep_peaks], observed[keep_peaks]
        train_mask, records = mask_matrix(observed)
        n_eval = int(sum(len(r["eval"]) for r in records))
        print(f"peaks evaluable: {keep_peaks.sum()}/{N_TOP_PEAKS} | observed "
              f"{observed.sum()} ({observed.mean() * 100:.2f}%) | eval entries {n_eval}")

        scale_res = {
            "n_peaks_evaluable": int(keep_peaks.sum()),
            "n_observed_entries": int(observed.sum()),
            "density": float(observed.mean()),
            "n_eval_entries": n_eval,
            "n_cal_entries": int(sum(len(r["cal"]) for r in records)),
            "peak_obs_per_peak_median": float(np.median(obs_per_peak[keep_peaks])),
        }

        # ---- GP, pipeline-default length scale ----
        pred, std, t_fit, t_pred = run_gp(usage, train_mask, xy, multiplier=5.0)
        metrics = score(usage, records, pred)
        metrics.update(conformal(usage, records, pred, std, "locally_adaptive"))
        metrics.update({"time_fit_s": t_fit, "time_predict_full_s": t_pred})
        # Nyström degeneracy detector: with a kernel much shorter than the
        # inducing spacing, K_nm -> 0 and every prediction collapses to the
        # prior mean sigmoid(0) = 0.5 (logit transform).
        metrics["predictions_collapse_to_prior"] = bool(
            np.std(pred) < 1e-6
            or (abs(float(np.mean(pred)) - 0.5) < 1e-3 and np.std(pred) < 1e-2))
        add_rows(rows, N, "gp_default", metrics)
        scale_res["gp_default"] = metrics
        print(f"gp_default: rmse={metrics['rmse']:.4f} "
              f"cov95={metrics['coverage_95']:.3f} w95={metrics['width_median_95']:.4f} "
              f"fit={t_fit:.1f}s pred={t_pred:.1f}s "
              f"collapse={metrics['predictions_collapse_to_prior']}")

        # ---- GP, inducing-spacing-adapted length scale ----
        from scipy.spatial import cKDTree
        from spagapa.imputation.sparse_gp import SparseGPImputer
        probe = SparseGPImputer(n_inducing=100)
        inducing = probe._select_inducing_points(xy)
        d_ind = cKDTree(inducing).query(xy, k=1)[0]
        ls_adaptive = 2.0 * float(np.median(d_ind))
        med_nn = float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
        print(f"length scales: default={5.0 * med_nn:.5f} adaptive={ls_adaptive:.5f} "
              f"(inducing spacing ratio {(5.0 * med_nn) / ls_adaptive:.2f})")
        pred, std, t_fit, t_pred = run_gp(usage, train_mask, xy,
                                          length_scale=ls_adaptive)
        metrics = score(usage, records, pred)
        metrics.update(conformal(usage, records, pred, std, "locally_adaptive"))
        metrics.update({"time_fit_s": t_fit, "time_predict_full_s": t_pred})
        add_rows(rows, N, "gp_inducing_scale", metrics)
        scale_res["gp_inducing_scale"] = metrics
        scale_res["length_scale_default"] = 5.0 * med_nn
        scale_res["length_scale_adaptive"] = ls_adaptive
        print(f"gp_inducing_scale: rmse={metrics['rmse']:.4f} "
              f"cov95={metrics['coverage_95']:.3f} w95={metrics['width_median_95']:.4f} "
              f"fit={t_fit:.1f}s pred={t_pred:.1f}s")

        # ---- per-peak mean ----
        pred, t_mean = run_mean(usage, train_mask)
        metrics = score(usage, records, pred)
        metrics.update(conformal(usage, records, pred, None, "global"))
        metrics.update({"time_fit_s": t_mean, "time_predict_full_s": 0.0})
        add_rows(rows, N, "mean", metrics)
        scale_res["mean"] = metrics
        print(f"mean: rmse={metrics['rmse']:.4f} cov95={metrics['coverage_95']:.3f} "
              f"w95={metrics['width_median_95']:.4f} fit={t_mean:.2f}s")

        scale_res["peak_rss_mb"] = rss_mb()
        scale_res["matrix_mb"] = float(usage.nbytes / 1024**2)
        scale_res["total_scale_time_s"] = time.perf_counter() - t_scale
        summary["scales"][str(N)] = scale_res
        # checkpoint
        (OUT / "hd_gp_summary_partial.json").write_text(json.dumps(summary, indent=2))
        pd.DataFrame(rows).to_csv(OUT / "hd_gp_vs_mean.csv", index=False)
        print(f"scale done in {scale_res['total_scale_time_s']:.0f}s | "
              f"peak RSS {scale_res['peak_rss_mb']:.0f} MB")

    summary["baseline_rss_mb_after_load"] = rss_base
    summary["conclusions"] = conclusions(summary)
    (OUT / "hd_gp_summary.json").write_text(json.dumps(summary, indent=2))
    pd.DataFrame(rows).to_csv(OUT / "hd_gp_vs_mean.csv", index=False)
    make_figure(rows, summary)
    print(f"\nwrote {OUT / 'hd_gp_vs_mean.csv'}")
    print(f"wrote {OUT / 'hd_gp_summary.json'}")
    print(f"wrote {OUT / 'hd_gp_benchmark.png'}")


def score(usage, records, pred):
    from scipy.stats import pearsonr
    t, p = gather(usage, records, pred, None, "eval")
    return {
        "rmse": float(np.sqrt(np.mean((t - p) ** 2))),
        "mae": float(np.mean(np.abs(t - p))),
        "pearson": float(pearsonr(t, p)[0]) if np.std(p) > 0 else float("nan"),
    }


def add_rows(rows, n, method, metrics):
    for k, v in metrics.items():
        rows.append({"n_barcodes": n, "method": method, "metric": k, "value": v})


def conclusions(s):
    lines = []
    for N in SCALES:
        r = s["scales"][str(N)]
        gp, gpa, mn = r["gp_default"], r["gp_inducing_scale"], r["mean"]
        lines.append(
            f"N={N}: RMSE gp_default={gp['rmse']:.4f} "
            f"gp_inducing_scale={gpa['rmse']:.4f} mean={mn['rmse']:.4f}; "
            f"95% coverage {gp['coverage_95']:.3f}/{gpa['coverage_95']:.3f}/"
            f"{mn['coverage_95']:.3f}; GP fit+predict "
            f"{gp['time_fit_s'] + gp['time_predict_full_s']:.0f}s.")
    n_big = s["scales"][str(SCALES[-1])]
    lines.append(
        "Default length scale (5x median NN) degenerates at HD density: "
        "inducing spacing grows as sqrt(N) while the kernel shrinks as "
        "1/sqrt(N), so K_nm -> 0 and predictions collapse to the prior mean "
        "(sigmoid(0)=0.5) unless the length scale tracks the inducing spacing."
        if n_big["gp_default"].get("predictions_collapse_to_prior")
        else "Default length scale remains functional at the largest scale.")
    return lines


def make_figure(rows, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.DataFrame(rows)
    df["value"] = pd.to_numeric(df.value, errors="coerce")  # bool metric -> object col
    methods = ["gp_default", "gp_inducing_scale", "mean"]
    labels = {"gp_default": "GP (default ls)",
              "gp_inducing_scale": "GP (inducing-scaled ls)",
              "mean": "per-peak mean"}
    colors = {"gp_default": "#e74c3c", "gp_inducing_scale": "#3498db",
              "mean": "#7f8c8d"}
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    def get(metric, method):
        d = df[(df.metric == metric) & (df.method == method)]
        return d.n_barcodes.values, d.value.values

    ax = axes[0, 0]
    for m in methods:
        x, y = get("rmse", m)
        ax.plot(x, y, "o-", label=labels[m], color=colors[m])
    ax.set(xscale="log", xlabel="barcodes sampled", ylabel="RMSE (eval)",
           title="Accuracy: RMSE")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[0, 1]
    for m in methods:
        d = df[(df.metric == "time_fit_s") & (df.method == m)].sort_values("n_barcodes")
        d2 = df[(df.metric == "time_predict_full_s") & (df.method == m)].sort_values("n_barcodes")
        ax.plot(d.n_barcodes.values, d.value.values + d2.value.values, "o-",
                label=labels[m], color=colors[m])
    ax.set(xscale="log", yscale="log", xlabel="barcodes sampled",
           ylabel="fit + full-surface predict (s)", title="Runtime")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[0, 2]
    for m in methods:
        x, y = get("pearson", m)
        ax.plot(x, y, "o-", label=labels[m], color=colors[m])
    ax.set(xscale="log", xlabel="barcodes sampled", ylabel="Pearson r (eval)",
           title="Accuracy: Pearson")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # coverage vs nominal at each scale (GP adaptive + mean)
    ax = axes[1, 0]
    for m in ("gp_inducing_scale", "gp_default", "mean"):
        for N in SCALES:
            d = df[(df.method == m) & (df.n_barcodes == N)]
            lev = [80, 90, 95]
            cov = [d[d.metric == f"coverage_{l}"].value.iloc[0] for l in lev]
            style = "-" if m != "mean" else "--"
            ax.plot(lev, cov, style + "o", color=colors[m], alpha=0.35 + 0.65 * (N == SCALES[-1]),
                    label=labels[m] if N == SCALES[-1] else None, ms=4)
    ax.plot([80, 95], [0.80, 0.95], "k:", label="nominal")
    ax.set(xlabel="nominal coverage (%)", ylabel="empirical coverage",
           title="Conformal coverage (faint->dark = 5k->100k)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1, 1]
    for m in methods:
        x, y = get("width_median_95", m)
        ax.plot(x, y, "o-", label=labels[m], color=colors[m])
    ax.set(xscale="log", xlabel="barcodes sampled",
           ylabel="median interval width (95%)", title="Conformal interval width")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1, 2]
    rss = [summary["scales"][str(N)]["peak_rss_mb"] for N in SCALES]
    ax.plot(SCALES, rss, "o-", color="#2c3e50")
    ax.set(xscale="log", xlabel="barcodes sampled", ylabel="process peak RSS (MB)",
           title="Memory (incl. data load)")
    ax.grid(alpha=0.3)

    fig.suptitle("Visium HD sparse-GP vs per-peak mean — adrenal scAPAtrap, "
                 "top-300 peaks, 20% masking, seed 42", fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "hd_gp_benchmark.png", dpi=150, bbox_inches="tight")




if __name__ == "__main__":
    main()

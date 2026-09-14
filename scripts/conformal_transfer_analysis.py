#!/usr/bin/env python
"""Cross-sample conformal transfer analysis (leave-one-sample-out).

For each frozen sample (src) and each held-out sample (tgt), ask whether the
conformal interval half-width calibrated on src still covers tgt's errors:

    q_hat_src   = median over src observations of (hi_X - lo_X) / 2
    mid_tgt     = (hi_X + lo_X) / 2          (prediction proxy, see notes)
    err_tgt     = |truth_tgt - mid_tgt|
    cross cov   = P(err_tgt <= q_hat_src)

Within-sample coverage (src == tgt) is computed identically as the reference.

Outputs (pipeline_output/conformal_transfer/):
    transfer_results.csv   one row per (src, tgt, level)
    transfer_summary.json  aggregated within vs cross statistics
    run.log                console capture
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

NPZ = Path("pipeline_output/conformal_validation/per_observation_bounds.npz")
OUT_DIR = Path("pipeline_output/conformal_transfer")
LEVELS = (80, 90, 95)


def load_samples(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    names = sorted({k.rsplit("__", 1)[0] for k in d.files})
    data = {}
    for s in names:
        data[s] = {a: d[f"{s}__{a}"].astype(np.float64)
                   for a in ("truth", "lo_80", "hi_80", "lo_90", "hi_90",
                             "lo_95", "hi_95")}
    return names, data


def per_sample_stats(arr, level):
    lo, hi = arr[f"lo_{level}"], arr[f"hi_{level}"]
    half = (hi - lo) / 2.0
    q_hat = float(np.median(half))
    mid = (hi + lo) / 2.0
    err = np.abs(arr["truth"] - mid)
    within_cov = float(np.mean(err <= q_hat))
    return q_hat, err, within_cov


def main():
    names, data = load_samples(NPZ)
    n = len(names)
    print(f"Loaded {n} frozen samples from {NPZ}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    summary = {
        "analysis": "cross_sample_conformal_transfer_loo",
        "input": str(NPZ),
        "n_samples": n,
        "samples": names,
        "levels": list(LEVELS),
        "per_level": {},
        "notes": (
            "Prediction proxy: the true point predictions were not stored in "
            "per_observation_bounds.npz, so the interval midpoint (hi_X + lo_X)/2 "
            "is used as the prediction proxy. For symmetric (centered) conformal "
            "intervals the midpoint is an unbiased proxy for the point prediction; "
            "any residual bias only affects |truth - midpoint| mildly. "
            "q_hat_src = median half-width of src's intervals at level X; "
            "cross coverage = fraction of tgt observations with "
            "|truth - midpoint| <= q_hat_src. Within-sample rows (src == tgt) use "
            "the identical estimator and serve as the reference; the small gap "
            "between within coverage and the nominal level reflects the median-"
            "half-width summary of a heterogeneous width distribution."
        ),
    }

    for level in LEVELS:
        stats = {s: per_sample_stats(data[s], level) for s in names}
        within = {s: stats[s][2] for s in names}

        level_rows = []
        for src in names:
            q_src = stats[src][0]
            for tgt in names:
                err_t = stats[tgt][1]
                cov = float(np.mean(err_t <= q_src))
                rows.append({
                    "level": level,
                    "pair_type": "within" if src == tgt else "cross",
                    "src_sample": src,
                    "tgt_sample": tgt,
                    "n_obs_tgt": int(err_t.size),
                    "q_hat_src": round(q_src, 6),
                    "coverage": round(cov, 6),
                    "within_coverage_tgt": round(within[tgt], 6),
                    "decay_vs_within_pp": round((within[tgt] - cov) * 100.0, 3),
                })
                if src != tgt:
                    level_rows.append((src, tgt, cov))

        w = np.array([within[s] for s in names])
        c = np.array([cov for _, _, cov in level_rows])
        decays = np.array([within[t] - cov for _, t, cov in level_rows])
        worst_i = int(np.argmax(decays))
        nominal = level / 100.0

        summary["per_level"][str(level)] = {
            "nominal": nominal,
            "within_mean": round(float(w.mean()), 6),
            "within_std": round(float(w.std(ddof=1)), 6),
            "cross_mean": round(float(c.mean()), 6),
            "cross_std": round(float(c.std(ddof=1)), 6),
            "transfer_decay_pp": round(float((w.mean() - c.mean()) * 100.0), 3),
            "decay_std_pp": round(float(decays.std(ddof=1)) * 100.0, 3),
            "max_decay_pp": round(float(decays.max()) * 100.0, 3),
            "max_decay_pair": {
                "src": level_rows[worst_i][0],
                "tgt": level_rows[worst_i][1],
            },
            "cross_max_dev_from_nominal_pp": round(
                float(np.max(np.abs(c - nominal)) * 100.0), 3),
            "within_max_dev_from_nominal_pp": round(
                float(np.max(np.abs(w - nominal)) * 100.0), 3),
            "n_cross_pairs": len(level_rows),
        }

        lv = summary["per_level"][str(level)]
        print(f"level {level}: within {lv['within_mean']:.4f}±{lv['within_std']:.4f}  "
              f"cross {lv['cross_mean']:.4f}±{lv['cross_std']:.4f}  "
              f"decay {lv['transfer_decay_pp']:+.2f} pp  "
              f"max decay {lv['max_decay_pp']:.2f} pp  "
              f"max |cross-nominal| {lv['cross_max_dev_from_nominal_pp']:.2f} pp")

    csv_path = OUT_DIR / "transfer_results.csv"
    with open(csv_path, "w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(rows)

    json_path = OUT_DIR / "transfer_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {csv_path} ({len(rows)} rows) and {json_path}")


if __name__ == "__main__":
    sys.exit(main())

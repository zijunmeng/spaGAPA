#!/usr/bin/env python
"""GP-vs-mean RMSE benchmark run identically on both caller inputs
(scAPAtrap baseline and Sierra), mirroring the supplementary-fig S4 protocol:
mask 20% of observed entries per gene, impute with
SparseGPImputer(n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1),
per-gene held-out RMSE for GP vs per-gene mean. 200 genes per input, seed 42.
Usage: python gp_vs_mean.py --dataset gse220442_gsm6801751
Writes <outdir>/metrics_d.json.
"""
import os, sys, json, time, argparse
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spagapa import SparseGPImputer
from datasets import get as get_dataset

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--dataset", required=True,
                    help="dataset key from datasets.py; runs both caller inputs")
parser.add_argument("--out-json", default=None)
args = parser.parse_args()
ds = get_dataset(args.dataset)
BASE = ds["baseline"]
SIERRA = ds["outdir"]
OUTJSON = args.out_json or os.path.join(SIERRA, "metrics_d.json")
N_GENES = 200
MASK_FRAC = 0.2
SEED = 42


def run_one(name, apa_path, coords_path):
    apa = pd.read_csv(apa_path, index_col=0)
    coords = pd.read_csv(coords_path)
    X = coords[["x", "y"]].values.astype(float)

    obs_frac = apa.notna().mean(axis=1)
    eligible = apa[obs_frac.between(0.05, 0.55)]
    n = min(N_GENES, len(eligible))
    keep = eligible.sample(n=n, random_state=SEED)
    V = keep.values.astype(float)
    rng = np.random.RandomState(SEED)
    mask = np.zeros_like(V, dtype=bool)
    for g in range(V.shape[0]):
        obs = np.where(np.isfinite(V[g]))[0]
        m = max(1, int(round(MASK_FRAC * len(obs))))
        mask[g, rng.choice(obs, size=m, replace=False)] = True

    imp = SparseGPImputer(n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1)
    rows = []
    t0 = time.time()
    for g in range(V.shape[0]):
        vals = V[g]
        tm = np.isfinite(vals) & (~mask[g])
        pred_mean = np.full(len(vals), np.nanmean(vals[tm]))
        try:
            pred_gp, _ = imp.impute(X, vals, mask=tm, return_uncertainty=True)
        except Exception:
            pred_gp = pred_mean.copy()
        hg = mask[g] & np.isfinite(vals)
        rmse_gp = np.sqrt(np.nanmean((vals[hg] - pred_gp[hg]) ** 2))
        rmse_mean = np.sqrt(np.nanmean((vals[hg] - pred_mean[hg]) ** 2))
        rows.append({"rmse_gp": rmse_gp, "rmse_mean": rmse_mean,
                     "delta": rmse_gp - rmse_mean, "n_obs": int(tm.sum())})
        if (g + 1) % 50 == 0:
            print(f"[{name}] {g+1}/{V.shape[0]} ({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    return {
        "n_genes": int(len(df)),
        "rmse_gp_median": float(df.rmse_gp.median()),
        "rmse_gp_mean": float(df.rmse_gp.mean()),
        "rmse_mean_median": float(df.rmse_mean.median()),
        "rmse_mean_mean": float(df.rmse_mean.mean()),
        "frac_gp_better": float((df.delta < 0).mean()),
        "delta_median": float(df.delta.median()),
        "delta_mean": float(df.delta.mean()),
    }


def main():
    res = {"protocol": {
        "n_genes": N_GENES, "mask_fraction": MASK_FRAC, "seed": SEED,
        "imputer": "SparseGPImputer(n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1)",
        "gene_filter": "obs_frac in [0.05, 0.55]",
    }}
    res["scapatrap"] = run_one("scapatrap",
                               os.path.join(BASE, "apa_matrix.csv"),
                               os.path.join(BASE, "coordinates.csv"))
    res["sierra"] = run_one("sierra",
                            os.path.join(SIERRA, "apa_matrix.csv"),
                            os.path.join(BASE, "coordinates.csv"))
    with open(OUTJSON, "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()

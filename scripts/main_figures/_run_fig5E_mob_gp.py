"""Pick layer-differentiated real MOB genes + run real spaGAPA sparse-GP imputation.

This is the data-generation script for Figure 5 Panel E. It is run ONCE to produce
the cached imputed maps consumed by fig5_domain.py, so the figure script itself stays
pure-plotting (fast, deterministic).

The MOB apa_matrix.csv is a dense per-spot APA-usage matrix (values in [0,1], no NaN).
To honestly demonstrate spaGAPA's spatial-field reconstruction we mask 50% of spots
per gene (simulating sparse observation), fit the GP on the remaining 50%, and
reconstruct the full map. The full dense matrix is kept as ground truth.

Outputs (consumed by fig5_domain.py Panel E):
  pipeline_output/main_figures/_data/fig5E_mob_gp_imputed.csv
  pipeline_output/main_figures/_data/fig5E_mob_gp_imputed_meta.json
"""
import os, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")
import numpy as np
import pandas as pd

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
DATA = f"{REPO}/data/processed/mob_st11"
OUT  = f"{REPO}/pipeline_output/main_figures/_data/fig5E_mob_gp_imputed.csv"

apa  = pd.read_csv(f"{DATA}/apa_matrix.csv", index_col=0)        # genes x spots
coord = pd.read_csv(f"{DATA}/coordinates.csv", index_col=0)      # spots ; x,y
lab  = pd.read_csv(f"{DATA}/labels.csv", index_col=0)            # spots ; label
spots = apa.columns
apa = apa[spots]; coord = coord.loc[spots]; lab = lab.loc[spots]
X = coord[["x", "y"]].values.astype(float)
layers = ["GCL", "GL", "MCL", "ONL", "OPL"]
labv = lab["label"].values

# ---- rank genes by across-layer APA-usage spread (real signal) ----
means_all = {L: apa.loc[:, labv == L].mean(axis=1) for L in layers}
mean_df = pd.DataFrame(means_all)
spread = (mean_df.max(axis=1) - mean_df.min(axis=1)).sort_values(ascending=False)
nz = (apa.values.astype(float) > 0).mean(axis=1)
spread = spread[pd.Series(nz, index=apa.index) >= 0.10]

# ---- pick one gene per layer argmax, prefer outer->inner for visual contrast ----
picked, seen = [], set()
for L in ["ONL", "GCL", "GL", "MCL", "OPL"]:
    cands = spread.index[~spread.index.isin(seen)]
    am = mean_df.loc[cands].idxmax(axis=1)
    cL = am[am == L].index
    if len(cL):
        g = cL[0]; picked.append((L, g)); seen.add(g)
picked = picked[:4]
print("[select] PICKED genes:")
for L, g in picked:
    r = mean_df.loc[g]
    print(f"  {L:3s}  {g:12s}  spread={spread[g]:.3f}  "
          + "  ".join(f"{x}={r[x]:.3f}" for x in layers), flush=True)
genes = [g for _, g in picked]

# ---- real spaGAPA sparse-GP imputation (50% mask) ----
from spagapa import SparseGPImputer
imp = SparseGPImputer(n_inducing=150, length_scale_multiplier=5.0, noise_level=0.1)
rng = np.random.RandomState(0)
MASK_FRAC = 0.50

out = pd.DataFrame({"spot_id": spots, "x": X[:, 0], "y": X[:, 1], "label": labv})
for gname in genes:
    v = apa.loc[gname].values.astype(float)
    n = len(v)
    n_obs = max(20, int(round((1.0 - MASK_FRAC) * n)))
    obs_idx = rng.choice(n, size=n_obs, replace=False)
    mask_obs = np.zeros(n, dtype=bool); mask_obs[obs_idx] = True
    try:
        pred, unc = imp.impute(X, v, mask=mask_obs, return_uncertainty=True)
    except Exception as e:
        print(f"[gp] {gname} FAILED ({e}); fallback mean", flush=True)
        pred = np.where(mask_obs, v, np.nanmean(v[mask_obs]))
    out[f"{gname}__raw"] = v
    out[f"{gname}__observed_mask"] = mask_obs.astype(int)
    out[f"{gname}__gp"] = pred
    held = ~mask_obs
    rmse_gp = np.sqrt(np.nanmean((v[held] - pred[held]) ** 2))
    rmse_mean = np.sqrt(np.nanmean((v[held] - np.nanmean(v[mask_obs])) ** 2))
    print(f"[gp] {gname}: obs={mask_obs.sum()}/{n}  RMSE_gp={rmse_gp:.3f}  RMSE_mean={rmse_mean:.3f}", flush=True)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.to_csv(OUT, index=False)
meta = {"genes": genes, "layer_argmax": {g: L for L, g in picked},
        "method": "spaGAPA SparseGPImputer; 50% spots masked, GP reconstructs the rest",
        "mask_fraction": MASK_FRAC, "random_seed": 0,
        "n_inducing": 150, "length_scale_multiplier": 5.0, "noise_level": 0.1,
        "n_spots": int(len(spots)), "source_apa": f"{DATA}/apa_matrix.csv"}
json.dump(meta, open(OUT.replace(".csv", "_meta.json"), "w"), indent=2)
print(f"\n[done] wrote {OUT}")
print(f"[done] wrote {OUT.replace('.csv', '_meta.json')}")

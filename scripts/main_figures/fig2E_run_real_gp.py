#!/usr/bin/env python
"""Figure 2E provenance: run REAL spaGAPA SparseGPImputer on a single
high-spatial-signal gene from GSE183456 and save per-spot posterior.

This is NOT a benchmark (that is Panel A-D / transparent_comparison.csv).
Panel E is a *visual* demonstration on one representative gene: truth /
mean / real-spaGAPA-GP / |GP-truth|. The per-spot GP posterior is saved
for provenance so the figure can be regenerated deterministically.

Gene selection: compute Moran's I on every gene with enough observations,
then pick the gene with the highest Moran's I among those whose observed
fraction is in a visually informative range (10-30%). This yields a gene
with strong spatial structure where the difference between a flat per-gene
mean and the spatial GP posterior is visually obvious.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
from spagapa import SparseGPImputer

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT = os.path.join(ROOT, "pipeline_output/main_figures/_data")
APA = os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap/apa_matrix.csv")
COORD = os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap/coordinates.csv")

# Match the SparseGPImputer call used in supp_fig04 (fidelity to the method actually published).
SEED = 42
MASK_FRAC = 0.20  # 20% held-out, same as Panel A-D masking design.


def morans_i(values, W):
    finite = np.isfinite(values)
    if finite.sum() < 5:
        return np.nan
    x = values[finite]
    z = x - x.mean()
    denom = (z ** 2).sum()
    if denom == 0:
        return 0.0
    Wf = W[finite][:, finite]
    s0 = Wf.sum()
    if s0 == 0:
        return np.nan
    n = finite.sum()
    cross = (Wf @ z) @ z
    return (n / s0) * (cross / denom)


def main():
    os.makedirs(OUT, exist_ok=True)
    apa = pd.read_csv(APA, index_col=0)            # rows=genes, cols=spots
    coord = pd.read_csv(COORD)                     # cols: spot_id, x, y
    # align spots to coord order
    apa = apa[coord["spot_id"].astype(str).values]
    X = coord[["x", "y"]].values.astype(float)

    # spatial weights (8-NN, symmetric, row-normalised) for Moran's I -- same as supp_fig04
    K = 8
    knn = NearestNeighbors(n_neighbors=K + 1).fit(X)
    _, idx = knn.kneighbors(X)
    rows = np.repeat(np.arange(len(X)), K)
    cols = idx[:, 1:].ravel()
    Wraw = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(len(X), len(X)))
    Wsym = (Wraw + Wraw.T) * 0.5
    deg = np.asarray(Wsym.sum(axis=1)).ravel(); deg[deg == 0] = 1.0
    W = Wsym.multiply(1.0 / deg[:, None]).tocsr()

    # candidate genes: observed fraction in an informative visual range.
    # apa values are floats (APA usage); treat nonzero as "observed".
    obs_frac = (apa > 0).mean(axis=1)
    cand = apa[(obs_frac >= 0.10) & (obs_frac <= 0.30)]
    print(f"[2E] candidate genes (obs_frac in [0.10,0.30]): {len(cand)}")

    # compute Moran's I per candidate on the raw observed signal
    print("[2E] computing Moran's I per candidate gene ...", flush=True)
    mis = []
    for gene, row in cand.iterrows():
        v = row.values.astype(float)
        v[~np.isfinite(v)] = 0.0
        mis.append((gene, morans_i(v, W), obs_frac.loc[gene], int((v > 0).sum())))
    mi_df = pd.DataFrame(mis, columns=["gene", "morans_i", "obs_frac", "n_obs"]).dropna(subset=["morans_i"])
    mi_df = mi_df.sort_values("morans_i", ascending=False).reset_index(drop=True)
    print("[2E] top-10 candidate genes by Moran's I:")
    print(mi_df.head(10).to_string(index=False))

    # pick the top gene (highest spatial autocorrelation)
    gene = mi_df.iloc[0]["gene"]
    mi_val = mi_df.iloc[0]["morans_i"]
    print(f"[2E] SELECTED gene = {gene}  (Moran's I = {mi_val:+.4f}, obs_frac = {mi_df.iloc[0]['obs_frac']:.3f})")

    row = apa.loc[gene].values.astype(float)
    rng = np.random.RandomState(SEED)
    obs_idx = np.where(np.isfinite(row))[0]
    nmask = int(round(MASK_FRAC * len(obs_idx)))
    held = rng.choice(obs_idx, size=nmask, replace=False)
    mask_held = np.zeros(len(row), dtype=bool)
    mask_held[held] = True
    tm = np.isfinite(row) & (~mask_held)   # training mask = observed & not held-out

    # --- REAL spaGAPA-GP ---
    imp = SparseGPImputer(n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1)
    pred_gp, unc_gp = imp.impute(X, row, mask=tm, return_uncertainty=True)

    # --- per-gene mean baseline ---
    pred_mean = np.full(len(row), np.nanmean(row[tm]))

    # held-out RMSE on the held-out entries only
    rmse_gp_held = float(np.sqrt(np.nanmean((row[mask_held] - pred_gp[mask_held]) ** 2)))
    rmse_mean_held = float(np.sqrt(np.nanmean((row[mask_held] - pred_mean[mask_held]) ** 2)))

    # save per-spot posterior table
    out = pd.DataFrame({
        "spot_id": coord["spot_id"].values,
        "x": X[:, 0], "y": X[:, 1],
        "truth": row,
        "held_out": mask_held,
        "train": tm,
        "pred_mean": pred_mean,
        "pred_gp": pred_gp,
        "unc_gp": unc_gp,
        "abs_err_gp": np.abs(pred_gp - row),
        "abs_err_mean": np.abs(pred_mean - row),
    })
    out_path = os.path.join(OUT, f"fig2E_gse183456_{gene}_gp_posterior.csv")
    out.to_csv(out_path, index=False)
    print(f"[2E] wrote -> {out_path}")

    meta = {
        "gene": gene,
        "morans_i": float(mi_val),
        "obs_frac": float(mi_df.iloc[0]["obs_frac"]),
        "n_obs": int(mi_df.iloc[0]["n_obs"]),
        "n_heldout": int(nmask),
        "rmse_gp_held": rmse_gp_held,
        "rmse_mean_held": rmse_mean_held,
        "delta_rmse_gp_minus_mean": rmse_gp_held - rmse_mean_held,
        "seed": SEED,
        "sparse_gp_kwargs": "n_inducing=100, length_scale_multiplier=5.0, noise_level=0.1",
        "data": "GSE183456 GSM6047774 (scapatrap)",
    }
    import json
    meta_path = os.path.join(OUT, f"fig2E_gse183456_{gene}_meta.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[2E] meta -> {meta_path}")
    print(f"[2E] held-out RMSE: GP={rmse_gp_held:.4f}  mean={rmse_mean_held:.4f}  delta={rmse_gp_held-rmse_mean_held:+.4f}")
    print("[2E] DONE")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""MOB domain-recovery validation (Phase 3 Task Y, quantitative).

Recovers MOB olfactory-bulb layers (GCL/MCL/ONL/GL/OPL) **unsupervised** with
spaGAPA on spvAPA's ST11RUD (gene x spot RUD/APA matrix, 260 spots) and scores
the recovery against spvAPA's ST11label expert layer labels with ARI/NMI.

Pipeline per spaGAPA:
  1. GP imputation (exact GP -- 260 spots is tiny) of the RUD matrix with
     uncertainty.
  2. Build the fused spatial + APA (+ optional expression) affinity graph via
     spaGAPA's MultiViewGraphBuilder.
  3. Leiden community detection on the fused graph (the high-res route in
     spagapa/bioml/highres.py) plus a spectral k=5 baseline (what the standard
     pipeline uses for small slides).
  4. ARI / NMI of recovered domains vs MOB layers.

Inputs : data/processed/mob_st11/{apa_matrix.csv, coordinates.csv,
         labels.csv, expression.csv}   (produced by export_spvapa_st11.R)
Outputs: pipeline_output/mob_domain_recovery/
           - spagapa_domains.csv        spot_id, x, y, true_label, domain_<run>
           - spagapa_metrics.json       ARI/NMI per run + params + wall time
           - mob_domain_map.png         spots colored by domain + true label
           - mob_uncertainty_map.png    per-spot GP uncertainty
           - fused_graph_weights.json   effective multi-view weights

Competitor comparison (spvAPA makeCluster, stAPAminer makeStCluster+findLabels)
is run from R via the matching wrappers; this script only does the spaGAPA
side and writes its metrics. The R wrappers append competitor ARI/NMI to a
sibling JSON.

Usage (S91):
  ~/anaconda3/envs/spagapa/bin/python scripts/mob_domain_recovery.py
  # optional competitor runs (R, shared lib):
  Rscript scripts/run_mob_competitor_clusters.R
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

# GP imputation parallelizes across genes via joblib (n_jobs). Set BLAS /
# OpenMP threads to 1 so the per-gene sklearn GPR workers don't oversubscribe
# the CPU (16 jobs x 64 BLAS threads = 1024 contended threads). Gene-level
# GPR problems are small enough that single-threaded BLAS per worker is best.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.makedirs("/tmp/matplotlib-spagapa", exist_ok=True)
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib-spagapa"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

# GP hyperparameter optimizers emit many ConvergenceWarnings; they don't
# affect domain recovery. Silence them to keep logs readable.
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.exceptions import ConvergenceWarning  # noqa: E402
warnings.filterwarnings("ignore", category=ConvergenceWarning)
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from spagapa.bioml import MultiViewGraphBuilder
from spagapa.imputation import GPImputer, SparseGPImputer, ExpressionFeatureBuilder

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
DATA = ROOT / "data" / "processed" / "mob_st11"
OUT = ROOT / "pipeline_output" / "mob_domain_recovery"
OUT.mkdir(parents=True, exist_ok=True)

N_MOB_LAYERS = 5  # GCL / GL / MCL / ONL / OPL


# ── helpers ───────────────────────────────────────────────────────────
def _leiden_on_graph(fused_graph, resolution: float, seed: int = 42):
    """Leiden community detection on a sparse affinity graph (highres.py path)."""
    import igraph as ig
    import leidenalg

    coo = fused_graph.tocoo()
    g = ig.Graph(
        n=fused_graph.shape[0],
        edges=list(zip(coo.row.tolist(), coo.col.tolist())),
        directed=False,
    )
    g.es["weight"] = coo.data.tolist()
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=resolution,
        seed=seed,
    )
    return np.array(partition.membership, dtype=int)


def _spectral_on_graph(fused_graph, k: int, seed: int = 42):
    from sklearn.cluster import SpectralClustering

    sc = SpectralClustering(
        n_clusters=k, affinity="precomputed", assign_labels="kmeans",
        random_state=seed,
    )
    return sc.fit_predict(fused_graph.tocsr())


def _mask(values: np.ndarray) -> np.ndarray:
    """APA-index mask: observed where finite and > 0 (qc_metrics convention)."""
    return np.isfinite(values) & (values > 0)


# ── load ──────────────────────────────────────────────────────────────
def load_st11():
    apa = pd.read_csv(DATA / "apa_matrix.csv", index_col=0)        # gene x spot
    coords = pd.read_csv(DATA / "coordinates.csv")                 # spot_id,x,y
    labels = pd.read_csv(DATA / "labels.csv")                      # spot_id,label
    expr_path = DATA / "expression.csv"
    expr = pd.read_csv(expr_path, index_col=0) if expr_path.exists() else None

    spots = list(apa.columns)
    assert list(coords["spot_id"]) == spots, "coords spot order != apa cols"
    assert list(labels["spot_id"]) == spots, "labels spot order != apa cols"
    if expr is not None:
        expr = expr[spots]  # align columns

    return apa, coords, labels, expr


# ── main ──────────────────────────────────────────────────────────────
def main():
    t0 = time.time()
    apa_df, coords_df, labels_df, expr_df = load_st11()
    spots = list(apa_df.columns)
    n_spots = len(spots)

    coords = coords_df[["x", "y"]].to_numpy(dtype=float)
    true_labels = labels_df["label"].to_numpy()
    apa = apa_df.to_numpy(dtype=float)  # gene x spot, RUD in [0,1]

    # Keep informative genes: observed in >= 5% of spots (enough spatial
    # signal to fit a GP and contribute to domain structure). ST11 RUD is
    # sparse; this drops all-zero / single-spot genes that only add GP cost.
    obs_mask = _mask(apa)
    gene_obs_frac = obs_mask.mean(axis=1)
    keep = gene_obs_frac >= 0.05
    n_genes_full = apa.shape[0]
    apa = apa[keep, :]
    apa_df = apa_df.loc[keep]
    n_genes = apa.shape[0]

    print(f"ST11: {n_genes_full} -> {n_genes} informative genes "
          f"(obs>=5%) x {n_spots} spots; "
          f"observed_fraction={_mask(apa).mean():.3f}")
    print(f"MOB layers: {dict(zip(*np.unique(true_labels, return_counts=True)))}")

    metrics = {
        "dataset": "spvAPA_ST11_MOB",
        "n_genes_full": int(n_genes_full),
        "n_genes": int(n_genes),
        "gene_filter": "observed_fraction >= 0.05",
        "n_spots": int(n_spots),
        "n_mob_layers": N_MOB_LAYERS,
        "true_layer_counts": {
            str(k): int(v) for k, v in zip(*np.unique(true_labels, return_counts=True))
        },
        "observed_fraction": float(_mask(apa).mean()),
        "runs": {},
    }

    # ── Step 1: GP imputation ─────────────────────────────────────────
    # Sparse GP is used (exact GP's per-gene hyperparameter optimizer is the
    # bottleneck across thousands of genes). Inducing points capped below
    # n_spots (260); length scale resolved from local coordinate density
    # (the same auto-scale the highres pipeline uses).
    t_gp = time.time()
    from sklearn.neighbors import NearestNeighbors
    k_ls = min(max(2, 6), n_spots - 1)
    nn_d = NearestNeighbors(n_neighbors=k_ls).fit(coords).kneighbors(coords)[0]
    pos_d = nn_d[:, -1]
    pos_d = pos_d[np.isfinite(pos_d) & (pos_d > 0)]
    length_scale = float(np.median(pos_d)) if pos_d.size else 1.0
    n_inducing = min(150, max(20, n_spots - 1))
    train_mask = _mask(apa)
    imputer = SparseGPImputer(
        n_inducing=n_inducing, inducing_method="kmeans",
        length_scale=length_scale, noise_level=0.1,
    )
    # fit_batch returns a *Batch imputer; its impute() takes no coordinates
    # (it predicts over the stored training coordinates).
    batch = imputer.fit_batch(coords, apa, mask=train_mask, n_jobs=32, verbose=False)
    imputed, uncertainty = batch.impute(return_uncertainty=True)
    imputed = np.clip(imputed, 0.0, 1.0)
    gp_wall = time.time() - t_gp
    obs = train_mask
    mae = float(np.mean(np.abs(apa[obs] - imputed[obs]))) if obs.sum() else float("nan")
    print(f"GP impute (sparse, n_ind={n_inducing}, ls={length_scale:.2f}) done in "
          f"{gp_wall:.1f}s (MAE={mae:.4f}); median uncertainty={np.median(uncertainty):.4f}")
    metrics["gp_impute"] = {
        "mode": "sparse_gp", "n_inducing": int(n_inducing),
        "length_scale": length_scale, "wall_s": gp_wall, "mae": mae,
        "median_uncertainty": float(np.median(uncertainty)),
        "mean_uncertainty": float(np.mean(uncertainty)),
    }

    # ── Step 2: optional expression embedding ─────────────────────────
    expr_embedding = None
    if expr_df is not None:
        # expression matrix is gene x spot -> ExpressionFeatureBuilder handles genes_by_spots
        emb = ExpressionFeatureBuilder(n_components=10, orientation="genes_by_spots")
        expr_embedding = emb.fit_transform(expr_df)
        print(f"Expression embedding: {expr_embedding.shape}")

    # ── Step 3: fused graph + Leiden / spectral domain recovery ───────
    # ST11 is tiny (260 spots). Use a small spatial kNN and weight APA heavily
    # since the task is APA-domain recovery. Try a few weight schemes.
    work = imputed
    apa_view = np.clip(work, 0.0, 1.0)
    apa_mean_view = None  # lazily built in the loop (cached after first use)

    # Schemes: (name, spatial_w, expression_w, apa_w, n_neighbors, apa_source).
    #   apa_source = "gp"  -> GP-imputed APA (spaGAPA; imputed/uncertainty above)
    #   apa_source = "mean"-> per-gene-mean-imputed APA (fair baseline; same
    #                        graph build + same Leiden pipeline, NO GP). Missing
    #                        entries filled with the observed per-gene mean, the
    #                        same convention MultiViewGraphBuilder uses
    #                        internally (spagapa/bioml/multiview_graph.py:147).
    # The mean_apa config carries the SAME weights as the main display config
    # (apa_dominant) so the only difference is the APA imputation source.
    # apa_mean_view is built lazily inside the loop (cached after first use).
    if expr_embedding is not None:
        schemes = [
            ("apa_dominant",  0.2, 0.2, 0.6, 6, "gp"),
            ("balanced",      0.4, 0.4, 0.2, 6, "gp"),
            ("spatial_apa",   0.5, 0.0, 0.5, 6, "gp"),
            ("expression_apa", 0.1, 0.5, 0.4, 6, "gp"),
            ("mean_apa",      0.2, 0.2, 0.6, 6, "mean"),
        ]
    else:
        schemes = [
            ("apa_dominant", 0.2, 0.0, 0.6, 6, "gp"),
            ("balanced",     0.4, 0.0, 0.2, 6, "gp"),
            ("spatial_apa",  0.5, 0.0, 0.5, 6, "gp"),
            ("mean_apa",     0.2, 0.0, 0.6, 6, "mean"),
        ]

    domains_out = {"spot_id": spots, "x": coords[:, 0], "y": coords[:, 1],
                   "true_label": true_labels}

    best = {"ari": -1.0, "name": None}
    for name, sw, ew, aw, nn, apa_src in schemes:
        t_g = time.time()
        if apa_src == "mean":
            # Fair baseline: per-gene-mean-imputed APA (no GP). Built once and
            # cached. Same fill convention as MultiViewGraphBuilder
            # (np.where(missing, gene_mean)); uncertainty=None so the graph
            # builder does NOT apply confidence weighting (mean impute has no
            # meaningful uncertainty to weight by).
            if apa_mean_view is None:
                apa_raw_mean = np.where(
                    obs, apa, np.nanmean(np.where(obs, apa, np.nan), axis=1,
                                         keepdims=True),
                )
                apa_raw_mean = np.nan_to_num(apa_raw_mean, nan=0.0)
                apa_mean_view = np.clip(apa_raw_mean, 0.0, 1.0)
            this_apa = apa_mean_view
            this_unc = None
        else:
            this_apa = apa_view
            this_unc = uncertainty
        graph = MultiViewGraphBuilder(
            n_neighbors=min(nn, n_spots - 1),
            spatial_weight=sw,
            expression_weight=ew,
            apa_weight=aw,
        ).build(
            coords,
            expression_embedding=expr_embedding,
            apa_matrix=this_apa,
            uncertainty=this_unc,
        )
        weights_eff = graph.weights
        fused = graph.fused

        # Spectral at k = N_MOB_LAYERS (fixed-k baseline, comparable to competitors)
        labels_sp = _spectral_on_graph(fused, k=N_MOB_LAYERS, seed=42)

        # Leiden at a few resolutions; pick the run closest to 5 communities,
        # and also record best-ARI Leiden run for honesty.
        leiden_runs = {}
        for res in (0.3, 0.5, 0.8, 1.0, 1.5):
            lab = _leiden_on_graph(fused, resolution=res, seed=42)
            leiden_runs[res] = lab

        g_wall = time.time() - t_g

        # Score everything
        ari_sp = adjusted_rand_score(true_labels, labels_sp)
        nmi_sp = normalized_mutual_info_score(true_labels, labels_sp)
        domains_out[f"{name}__spectral_k5"] = labels_sp

        run_rec = {
            "weights_requested": {"spatial": sw, "expression": ew, "apa": aw},
            "weights_effective": weights_eff,
            "n_neighbors": nn,
            "apa_source": apa_src,
            "graph_wall_s": g_wall,
            "spectral_k5": {"n_domains": int(len(np.unique(labels_sp))),
                            "ari": ari_sp, "nmi": nmi_sp},
            "leiden": {},
        }
        best_leiden_ari = -1.0
        for res, lab in leiden_runs.items():
            a = adjusted_rand_score(true_labels, lab)
            n = normalized_mutual_info_score(true_labels, lab)
            run_rec["leiden"][f"res_{res}"] = {
                "n_domains": int(len(np.unique(lab))), "ari": a, "nmi": n,
            }
            domains_out[f"{name}__leiden_res{res}"] = lab
            if a > best_leiden_ari:
                best_leiden_ari = a
                run_rec["leiden_best"] = {
                    "resolution": res, "n_domains": int(len(np.unique(lab))),
                    "ari": a, "nmi": n,
                }
            if a > best["ari"]:
                best = {"ari": a, "nmi": n, "name": f"{name}/leiden_res{res}",
                        "n_domains": int(len(np.unique(lab)))}

        # also track best spectral
        if ari_sp > best["ari"]:
            best = {"ari": ari_sp, "nmi": nmi_sp, "name": f"{name}/spectral_k5",
                    "n_domains": int(len(np.unique(labels_sp)))}

        metrics["runs"][name] = run_rec
        print(f"  [{name}] spectral_k5 ARI={ari_sp:.3f} NMI={nmi_sp:.3f} | "
              f"best Leiden ARI={best_leiden_ari:.3f}")

    metrics["best"] = best
    metrics["wall_total_s"] = time.time() - t0

    # ── Save domains + metrics ────────────────────────────────────────
    pd.DataFrame(domains_out).to_csv(OUT / "spagapa_domains.csv", index=False)
    (OUT / "spagapa_metrics.json").write_text(json.dumps(metrics, indent=2))
    print(f"\nBest spaGAPA recovery: {best['name']} "
          f"ARI={best['ari']:.3f} NMI={best['nmi']:.3f} "
          f"(n_domains={best['n_domains']}, MOB has 5)")

    # ── Figures ───────────────────────────────────────────────────────
    _plot_domain_map(domains_out, metrics, OUT / "mob_domain_map.png")
    _plot_uncertainty(coords, uncertainty, OUT / "mob_uncertainty_map.png")
    print(f"Saved -> {OUT}")


def _plot_domain_map(domains: dict, metrics: dict, path: Path):
    """Domain map: true labels + best spaGAPA run, side by side."""
    x = domains["x"]; y = domains["y"]
    true_lab = domains["true_label"]

    best_name = metrics["best"]["name"]
    # find the column matching best name
    scheme, run = best_name.split("/")
    col = f"{scheme}__{run}"
    pred = np.asarray(domains[col])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    layers = sorted(set(true_lab))
    layer_cmap = {l: i for i, l in enumerate(layers)}
    tc = [layer_cmap[l] for l in true_lab]
    sc1 = axes[0].scatter(x, y, c=tc, cmap="tab10", s=45, edgecolor="k", linewidth=0.2)
    axes[0].set_title(f"MOB expert layers (true)\n{', '.join(layers)}")
    axes[0].set_aspect("equal"); axes[0].set_xlabel("x"); axes[0].set_ylabel("y")

    n_dom = len(np.unique(pred))
    axes[1].scatter(x, y, c=pred, cmap="tab10", s=45, edgecolor="k", linewidth=0.2)
    ari = metrics["best"]["ari"]; nmi = metrics["best"]["nmi"]
    axes[1].set_title(f"spaGAPA unsupervised domains\n{best_name} "
                      f"({n_dom} domains, ARI={ari:.3f}, NMI={nmi:.3f})")
    axes[1].set_aspect("equal"); axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_uncertainty(coords: np.ndarray, uncertainty: np.ndarray, path: Path):
    """Per-spot median GP uncertainty map."""
    if uncertainty.ndim == 2:
        per_spot = np.median(uncertainty, axis=0)
    else:
        per_spot = uncertainty
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=per_spot, cmap="magma",
                    s=50, edgecolor="k", linewidth=0.2)
    plt.colorbar(sc, label="median GP uncertainty")
    ax.set_title("ST11 per-spot GP imputation uncertainty")
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

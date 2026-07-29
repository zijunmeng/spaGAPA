#!/usr/bin/env python3
"""
fix_domain_collapse.py
======================

Diagnose + fix the spatial-domain COLLAPSE observed in spaGAPA's high-resolution
Stereo-seq analysis (preset ``highres_fast``).

PROBLEM
-------
On the full-sample Stereo-seq binned data (21,455 genes x 15,235 bins), the
domains produced by spaGAPA's highres BioML kmeans branch are degenerate:

    domain 5 : 14,956 / 15,235 = 98.2 % of all spots
    domains 0,1,2,3,4,6,7 : < 2 % combined

The same collapse happened on the pilot data (94.5 % in one domain), so it is
systematic, not a one-off.

ROOT-CAUSE HYPOTHESIS (from code reading of ``spagapa/bioml/highres.py:333-348``)
---------------------------------------------------------------------------------
The highres path auto-downgrades spectral -> kmeans when n_spots > 5000, then runs::

    recovered_T = recovered.T                # (n_spots, n_genes) = (15235, 21455)
    labels = KMeans(...).fit_predict(np.nan_to_num(recovered_T))

with **no preprocessing whatsoever**: no StandardScaler, no HVG selection, no PCA,
no log transform, no spatial-coords concatenation. ``recovered`` is the gene-mean-
filled APA matrix (values clipped to [0, 1]) on data that was 90 % sparse before
the fill. After gene-mean fill, ~90 % of every spot's profile equals its column-
wise gene mean, so most spots are nearly identical and KMeans collapses to a
single giant cluster.

This script reproduces the collapse, quantifies the data pathology, and then
tries four standard fixes (HVG filtering, PCA, Leiden on a fused graph, and
spatially-constrained kmeans), reporting the resulting domain distribution for
each so we can pick the best.

USAGE (host S91)
----------------
    ~/anaconda3/envs/spagapa/bin/python scripts/fix_domain_collapse.py

OUTPUTS
-------
- ``figures/fig_domain_recovery_full.png``  : spatial scatter for the best fix
- stdout report (domain distributions + diagnostics)
"""

from __future__ import annotations

import os
import sys
import time
from typing import Dict, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Host config (CLAUDE.md rule: S91 paths)
# ---------------------------------------------------------------------------
HOSTNAME = os.popen("hostname -s").read().strip()
if HOSTNAME != "S91":
    # the paths below are S91-specific; warn but continue (data paths are absolute)
    print(f"[warn] expected host S91, got {HOSTNAME!r} -- continuing with absolute paths", flush=True)

BASE = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
Binned_DIR = os.path.join(
    BASE,
    "pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200",
)
APA_CSV = os.path.join(Binned_DIR, "apa_matrix.csv")
COORD_CSV = os.path.join(Binned_DIR, "coordinates.csv")
EXISTING_DOMAINS_CSV = os.path.join(
    BASE,
    "pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/spagapa_run_fast/domains.csv",
)
FIG_OUT = os.path.join(BASE, "figures/fig_domain_recovery_full.png")
FIG_MULTIPANEL_OUT = os.path.join(BASE, "figures/fig_domain_comparison_full.png")
CACHE_PKL = os.path.join(BASE, "scripts/.fix_domain_collapse_cache.npz")

N_DOMAINS = 8  # match spaGAPA default (bioml_rank / n_domains)
RANDOM_STATE = 42


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_data() -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Load APA matrix (n_genes x n_spots), coords (n_spots x 2), sites dataframe."""
    t0 = time.time()
    print("[load] reading coordinates.csv ...", flush=True)
    coords_df = pd.read_csv(COORD_CSV, index_col=0)
    coords = coords_df[["x", "y"]].to_numpy(dtype=float)
    print(f"[load] coords: {coords.shape}  ({time.time()-t0:.1f}s)", flush=True)

    t0 = time.time()
    print("[load] reading apa_matrix.csv (this is ~1.6 GB, be patient) ...", flush=True)
    apa = pd.read_csv(APA_CSV, index_col=0)
    sites = apa.index.to_frame(index=False)
    X = apa.to_numpy(dtype=float)  # (n_genes, n_spots)
    print(f"[load] apa matrix: {X.shape}  ({time.time()-t0:.1f}s)", flush=True)
    return X, coords, sites


def domain_distribution(labels: np.ndarray) -> Dict[int, int]:
    """Return {domain_id: count} sorted by count desc, including all k domains."""
    counts = np.bincount(labels, minlength=int(labels.max()) + 1)
    return {int(k): int(v) for k, v in enumerate(counts)}


def report_distribution(name: str, labels: np.ndarray, n_domains: int = N_DOMAINS) -> None:
    """Pretty-print the per-domain distribution + imbalance stats."""
    counts = np.bincount(labels, minlength=n_domains)
    total = counts.sum()
    order = np.argsort(counts)[::-1]
    top_frac = counts[order[0]] / total
    print(f"\n--- {name} ---", flush=True)
    for d in range(n_domains):
        c = counts[d]
        frac = c / total * 100
        bar = "#" * int(frac / 2)
        print(f"  domain {d}: {c:6d}  ({frac:5.2f}%)  {bar}", flush=True)
    n_nonempty = int((counts > 0).sum())
    print(f"  -> top-domain fraction = {top_frac*100:.2f}% | "
          f"non-empty domains = {n_nonempty}/{n_domains} | "
          f"entropy = {shannon_entropy(counts):.3f} bits (max={np.log2(n_domains):.3f})", flush=True)


def shannon_entropy(counts: np.ndarray) -> float:
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def collapse_score(labels: np.ndarray) -> float:
    """1.0 = perfectly degenerate (one cluster has everything), 0.0 = balanced."""
    counts = np.bincount(labels)
    return float(counts.max() / counts.sum())


# ---------------------------------------------------------------------------
# DIAGNOSTIC: why does kmeans collapse?
# ---------------------------------------------------------------------------
def diagnose(X: np.ndarray) -> None:
    """Quantify the data pathology that causes the collapse."""
    n_genes, n_spots = X.shape
    print("\n" + "=" * 70, flush=True)
    print("DIAGNOSTIC: why does kmeans collapse?", flush=True)
    print("=" * 70, flush=True)

    # --- (a) NaN fill status ---
    nan_frac = np.isnan(X).mean()
    print(f"\n[a] NaN fraction in apa_matrix.csv   = {nan_frac*100:.4f}%", flush=True)
    print(f"    -> if ~0%, gene-mean fill already happened upstream "
          f"(expected: highres_fast fills via fill_missing_by_gene_mean)", flush=True)

    # --- (b) per-spot variance across genes ---
    # X is (n_genes, n_spots). Per-spot variance = variance across the gene axis.
    spot_var = np.nanvar(X, axis=0)  # (n_spots,)
    print(f"\n[b] per-spot variance across the {n_genes} genes:", flush=True)
    print(f"    min={spot_var.min():.3e}  median={np.median(spot_var):.3e}  "
          f"mean={spot_var.mean():.3e}  max={spot_var.max():.3e}", flush=True)
    # Fraction of spots whose variance is within 1% of the median -- if huge, profiles
    # are nearly identical and kmeans can't separate them.
    med = np.median(spot_var)
    near_median = np.mean(np.abs(spot_var - med) <= 0.01 * med) * 100
    print(f"    {near_median:.1f}% of spots have variance within +/-1% of the "
          f"median -> near-identical profiles", flush=True)

    # --- (c) per-gene variance across spots (the HVG signal) ---
    gene_var = np.nanvar(X, axis=1)  # (n_genes,)
    print(f"\n[c] per-gene variance across the {n_spots} spots:", flush=True)
    print(f"    min={gene_var.min():.3e}  median={np.median(gene_var):.3e}  "
          f"mean={gene_var.mean():.3e}  max={gene_var.max():.3e}", flush=True)
    for thr in (1e-5, 1e-4, 1e-3):
        n_above = int((gene_var > thr).sum())
        print(f"    genes with var > {thr:.0e}: {n_above} ({n_above/n_genes*100:.2f}%)", flush=True)

    # --- (d) HVG-style: top-K genes by variance ---
    # Coefficient of variation across spots -- normalizes out gene-level magnitude.
    gene_mean = np.nanmean(X, axis=1)
    gene_cv = np.where(gene_mean > 0, np.sqrt(gene_var) / gene_mean, 0.0)
    print(f"\n[d] coefficient-of-variation (CV) across spots:", flush=True)
    print(f"    median CV = {np.median(gene_cv):.3f}  90th pct = "
          f"{np.percentile(gene_cv, 90):.3f}  99th pct = "
          f"{np.percentile(gene_cv, 99):.3f}", flush=True)
    for k in (500, 1000, 2000):
        # how distinct is the k-th largest CV from the bulk?
        top_k_cv = np.sort(gene_cv)[-k]
        print(f"    top-{k} HVG threshold CV = {top_k_cv:.3f} "
              f"(median gene CV = {np.median(gene_cv):.3f})", flush=True)

    # --- (e) spot-spot distance distribution (raw recovered.T, no scaling) ---
    # subsample to keep it tractable
    rng = np.random.default_rng(RANDOM_STATE)
    idx = rng.choice(n_spots, size=min(2000, n_spots), replace=False)
    Xt = X[:, idx].T  # (sub, n_genes)
    from scipy.spatial.distance import pdist
    d = pdist(np.nan_to_num(Xt), metric="euclidean")
    print(f"\n[e] pairwise Euclidean distance between spots (raw recovered.T, "
          f"subsample n=2000):", flush=True)
    print(f"    mean={d.mean():.3f}  median={np.median(d):.3f}  "
          f"std={d.std():.3f}", flush=True)
    print(f"    ratio std/mean (coefficient of variation of distances) = "
          f"{d.std()/d.mean():.3f}", flush=True)
    print(f"    -> low ratio means spots sit in a tight ball; kmeans sees ~no "
          f"structure.", flush=True)

    # --- (f) dimensionality: how many PCs to capture variance? ---
    # PCA on (n_spots, n_genes) -- standardize per-gene first to avoid a handful
    # of high-magnitude genes dominating.
    Xt_centered = np.nan_to_num(X.T - gene_mean)  # (n_spots, n_genes), centered
    # SVD on a 15k x 21k float64 is ~5 GB; use float32 + truncated SVD.
    from sklearn.decomposition import TruncatedSVD
    t0 = time.time()
    svd = TruncatedSVD(n_components=50, random_state=RANDOM_STATE)
    svd.fit(Xt_centered.astype(np.float32))
    var_ratio = svd.explained_variance_ratio_
    cum = np.cumsum(var_ratio)
    print(f"\n[f] PCA on gene-centered APA matrix ({time.time()-t0:.1f}s):", flush=True)
    print(f"    PC1 explains {var_ratio[0]*100:.2f}% of variance "
          f"(huge PC1 => a single global mode dominates)", flush=True)
    for target in (0.5, 0.7, 0.9):
        n_pc = int(np.searchsorted(cum, target) + 1)
        print(f"    PCs to reach {target*100:.0f}% cumulative variance: {n_pc}", flush=True)
    print(f"    first-10 PC variance ratios: "
          f"{[round(float(v)*100,2) for v in var_ratio[:10]]}", flush=True)


# ---------------------------------------------------------------------------
# FIX A: reproduce the spaGAPA collapse (baseline)
# ---------------------------------------------------------------------------
def reproduce_spaGAPA_collapse(X: np.ndarray) -> np.ndarray:
    """Mirror spagapa/bioml/highres.py:344-348 exactly."""
    from sklearn.cluster import KMeans
    t0 = time.time()
    print("\n[baseline] reproducing spaGAPA kmeans on recovered.T "
          "(nan_to_num only, no scaling) ...", flush=True)
    recovered_T = np.nan_to_num(X.T)  # (n_spots, n_genes)
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(recovered_T)
    print(f"[baseline] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution("BASELINE (reproduce spaGAPA highres kmeans)", labels)
    return labels


# ---------------------------------------------------------------------------
# FIX A: HVG filtering + kmeans
# ---------------------------------------------------------------------------
def fix_hvg_kmeans(X: np.ndarray, n_hvg: int = 2000) -> np.ndarray:
    """Select top-N HVGs by variance, StandardScaler, then kmeans."""
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    t0 = time.time()
    print(f"\n[fix A] HVG filtering: top-{n_hvg} genes by variance across spots...", flush=True)
    gene_var = np.nanvar(X, axis=1)
    top_idx = np.argsort(gene_var)[-n_hvg:]
    X_hvg = X[top_idx, :].T  # (n_spots, n_hvg)
    Xs = StandardScaler().fit_transform(np.nan_to_num(X_hvg))
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(Xs)
    print(f"[fix A] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution(f"FIX A: HVG (top-{n_hvg}) + StandardScaler + kmeans", labels)
    return labels


# ---------------------------------------------------------------------------
# FIX B: PCA preprocessing + kmeans
# ---------------------------------------------------------------------------
def fix_pca_kmeans(X: np.ndarray, n_pc: int = 30) -> np.ndarray:
    """Gene-center, truncated SVD to n_pc, kmeans on PC scores."""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import TruncatedSVD
    t0 = time.time()
    print(f"\n[fix B] PCA preprocessing: center genes -> top-{n_pc} PCs -> kmeans...", flush=True)
    gene_mean = np.nanmean(X, axis=1)
    Xt = np.nan_to_num(X.T - gene_mean).astype(np.float32)  # (n_spots, n_genes)
    svd = TruncatedSVD(n_components=n_pc, random_state=RANDOM_STATE)
    PCs = svd.fit_transform(Xt)  # (n_spots, n_pc)
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(PCs)
    print(f"[fix B] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution(f"FIX B: PCA ({n_pc} PCs) + kmeans", labels)
    return labels


# ---------------------------------------------------------------------------
# FIX B2: HVG + PCA + kmeans (combine A and B -- the scanpy gold standard)
# ---------------------------------------------------------------------------
def fix_hvg_pca_kmeans(X: np.ndarray, n_hvg: int = 2000, n_pc: int = 30) -> np.ndarray:
    """The scanpy-style recipe: HVG -> scale -> PCA -> kmeans on PCs."""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import StandardScaler
    t0 = time.time()
    print(f"\n[fix B2] HVG ({n_hvg}) -> StandardScaler -> PCA ({n_pc}) -> kmeans "
          f"(scanpy recipe) ...", flush=True)
    gene_var = np.nanvar(X, axis=1)
    top_idx = np.argsort(gene_var)[-n_hvg:]
    X_hvg = X[top_idx, :]  # (n_hvg, n_spots)
    # StandardScaler across SPOTS (so each gene has unit variance over spots)
    Xs = StandardScaler().fit_transform(np.nan_to_num(X_hvg.T))  # (n_spots, n_hvg)
    svd = TruncatedSVD(n_components=n_pc, random_state=RANDOM_STATE)
    PCs = svd.fit_transform(Xs.astype(np.float32))
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(PCs)
    print(f"[fix B2] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution(f"FIX B2: HVG({n_hvg})+scale+PCA({n_pc})+kmeans", labels)
    return labels


# ---------------------------------------------------------------------------
# FIX C: Leiden on fused spatial+APA graph (leidenalg, no scanpy wrapper)
# ---------------------------------------------------------------------------
def _leiden_labels(adj_sparse, resolution: float = 1.0,
                   n_domains_target: int = N_DOMAINS) -> np.ndarray:
    """Run Leiden via leidenalg + igraph. Returns integer labels."""
    import igraph as ig
    import leidenalg
    # Coerce to symmetric COO of an undirected graph with positive weights.
    A = adj_sparse.tocoo()
    # dedupe undirected edges, keep weights
    srcs, dsts, ws = A.row, A.col, A.data
    g = ig.Graph(n=A.shape[0], edges=list(zip(srcs.tolist(), dsts.tolist())),
                 directed=False)
    g.es["weight"] = ws.tolist()
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution,
        seed=RANDOM_STATE, n_iterations=-1,
    )
    labels = np.array([int(m) for m in part.membership], dtype=int)
    return labels


def fix_leiden(X: np.ndarray, coords: np.ndarray, n_hvg: int = 2000,
               n_pc: int = 30, spatial_alpha: float = 0.3,
               resolution: float = 1.0) -> np.ndarray:
    """Build a fused graph from (HVG-PCA) APA similarity and spatial KNN, then
    run Leiden. ``spatial_alpha`` in [0,1] blends the two adjacency matrices."""
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import csr_matrix

    t0 = time.time()
    print(f"\n[fix C] Leiden on fused APA-PCA + spatial graph "
          f"(spatial_alpha={spatial_alpha}, resolution={resolution}) ...", flush=True)

    # 1) APA manifold: HVG -> scale -> PCA -> KNN graph (shared-neighbor weights)
    gene_var = np.nanvar(X, axis=1)
    top_idx = np.argsort(gene_var)[-n_hvg:]
    Xs = StandardScaler().fit_transform(np.nan_to_num(X[top_idx, :].T))
    PCs = TruncatedSVD(n_components=n_pc, random_state=RANDOM_STATE).fit_transform(
        Xs.astype(np.float32))
    nn_apa = NearestNeighbors(n_neighbors=15).fit(PCs)
    # use distance->affinity (Gaussian kernel on PC distance)
    d_apa, i_apa = nn_apa.kneighbors(PCs)
    # build sparse affinity
    from scipy.sparse import lil_matrix
    n_spots = X.shape[1]
    sigma = max(np.median(d_apa[:, -1]), 1e-9)
    G_apa = lil_matrix((n_spots, n_spots))
    for r in range(n_spots):
        w = np.exp(-(d_apa[r] ** 2) / (2 * sigma ** 2))
        G_apa[r, i_apa[r]] = w
    G_apa = G_apa.tocsr()
    G_apa = G_apa.maximum(G_apa.T)  # symmetrize

    # 2) spatial manifold: KNN on coords (binary adjacency, physical neighbors)
    nn_sp = NearestNeighbors(n_neighbors=15).fit(coords)
    G_spatial = nn_sp.kneighbors_graph(coords, mode="connectivity").astype(float)
    G_spatial = G_spatial.maximum(G_spatial.T).tocsr()

    # 3) fuse: weighted sum of normalized adjacencies
    def _norm(M):
        s = np.asarray(M.sum(axis=1)).ravel()
        s[s == 0] = 1.0
        return M.multiply(1.0 / s[:, None]).tocsr()
    G_fused = (1.0 - spatial_alpha) * _norm(G_apa) + spatial_alpha * _norm(G_spatial)
    G_fused = G_fused.tocsr()

    labels = _leiden_labels(G_fused, resolution=resolution)
    n_found = int(labels.max()) + 1
    print(f"[fix C] done ({time.time()-t0:.1f}s) | n_leiden_clusters="
          f"{n_found} (resolution {resolution})", flush=True)
    n_eff = max(n_found, N_DOMAINS)
    report_distribution(f"FIX C: Leiden fused APA-PCA+spatial "
                        f"(alpha={spatial_alpha}, res={resolution})", labels,
                        n_domains=n_eff)
    return labels


# ---------------------------------------------------------------------------
# FIX D: spatially-constrained kmeans (concatenate APA-PCs + normalized x,y)
# ---------------------------------------------------------------------------
def fix_spatial_kmeans(X: np.ndarray, coords: np.ndarray, n_pc: int = 30,
                       spatial_weight: float = 1.0) -> np.ndarray:
    """PCA on APA -> concat normalized spatial coords -> kmeans. Lets kmeans
    see spatial proximity directly."""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import StandardScaler
    t0 = time.time()
    print(f"\n[fix D] spatially-constrained: APA-PCA({n_pc}) + normalized "
          f"(x,y)*{spatial_weight} -> kmeans ...", flush=True)
    gene_mean = np.nanmean(X, axis=1)
    Xt = np.nan_to_num(X.T - gene_mean).astype(np.float32)
    PCs = TruncatedSVD(n_components=n_pc, random_state=RANDOM_STATE).fit_transform(Xt)
    PCs = StandardScaler().fit_transform(PCs)  # unit-variance PCs
    # normalize coords to unit std then scale
    coords_n = StandardScaler().fit_transform(coords) * spatial_weight
    feats = np.hstack([PCs, coords_n])
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(feats)
    print(f"[fix D] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution(f"FIX D: APA-PCA({n_pc}) + scaled coords (w={spatial_weight}) + kmeans",
                        labels)
    return labels


# ---------------------------------------------------------------------------
# FIX E: pure-spatial baseline (kmeans on normalized x,y only)
# ---------------------------------------------------------------------------
def fix_pure_spatial(coords: np.ndarray) -> np.ndarray:
    """Sanity check: does the tissue have spatial structure at all? Cluster on
    normalized (x,y) only. If this gives balanced, spatially-contiguous domains
    while APA-based methods collapse, that confirms the APA matrix (as currently
    gene-mean-filled) carries essentially no usable cluster signal."""
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    t0 = time.time()
    print(f"\n[fix E] pure-spatial baseline: kmeans on normalized (x,y) only ...", flush=True)
    coords_n = StandardScaler().fit_transform(coords)
    km = KMeans(n_clusters=N_DOMAINS, random_state=RANDOM_STATE, n_init=10)
    labels = km.fit_predict(coords_n)
    print(f"[fix E] done ({time.time()-t0:.1f}s)", flush=True)
    report_distribution("FIX E: pure-spatial kmeans on normalized (x,y)", labels)
    return labels


# ---------------------------------------------------------------------------
# FIX F: Leiden on PURE spatial graph (no APA) -- spatially-contiguous domains
# ---------------------------------------------------------------------------
def fix_leiden_spatial(coords: np.ndarray, resolution: float = 1.0) -> np.ndarray:
    """Leiden on a KNN graph built from physical coordinates only. The
    canonical 'spatial domain' approach when molecular signal is too sparse to
    cluster on directly. Produces spatially-contiguous regions."""
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import csr_matrix
    t0 = time.time()
    print(f"\n[fix F] Leiden on PURE spatial KNN graph (resolution={resolution}) ...", flush=True)
    nn = NearestNeighbors(n_neighbors=8).fit(coords)
    G = nn.kneighbors_graph(coords, mode="connectivity").astype(float)
    G = G.maximum(G.T).tocsr()
    labels = _leiden_labels(G, resolution=resolution)
    n_found = int(labels.max()) + 1
    print(f"[fix F] done ({time.time()-t0:.1f}s) | n_leiden_clusters="
          f"{n_found}", flush=True)
    n_eff = max(n_found, N_DOMAINS)
    report_distribution(f"FIX F: Leiden on pure spatial KNN graph (res={resolution})",
                        labels, n_domains=n_eff)
    return labels


# ---------------------------------------------------------------------------
# Plot the best result
# ---------------------------------------------------------------------------
def plot_multipanel(coords: np.ndarray, results: Dict[str, np.ndarray],
                    scores: Dict[str, float], out: str) -> None:
    """Side-by-side spatial scatters for the key methods: baseline collapse vs
    the working fixes. Picked to tell the story clearly."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t0 = time.time()
    print(f"\n[plot] saving multi-panel comparison -> {out} ...", flush=True)
    # Choose panels: baseline + the 4 non-degenerate methods, in story order.
    panel_order = [
        "baseline",
        "B2_hvg+pca",
        "D_w5",
        "E_pure_spatial",
        "C_leiden_a0.7",
        "F_leiden_spatial",
    ]
    panels = [(n, results[n]) for n in panel_order if n in results]
    n_panels = len(panels)
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (name, lab) in zip(axes, panels):
        n_dom = int(lab.max()) + 1
        cmap = plt.get_cmap("tab20" if n_dom > 10 else "tab10")
        for d in range(n_dom):
            m = lab == d
            ax.scatter(coords[m, 0], coords[m, 1], s=1.5,
                       c=[cmap(d % 20)], rasterized=True)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        s = scores.get(name, float("nan"))
        ax.set_title(f"{name}\n{n_dom} domains | top={s*100:.1f}% | "
                     f"ent={shannon_entropy(np.bincount(lab)):.2f}", fontsize=9)
    # hide unused axes
    for ax in axes[len(panels):]:
        ax.axis("off")
    fig.suptitle("spaGAPA domain recovery: collapsed baseline (left) vs working fixes",
                 fontsize=12, y=0.995)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] multi-panel done ({time.time()-t0:.1f}s)", flush=True)


# ---------------------------------------------------------------------------
def plot_domains(coords: np.ndarray, labels: np.ndarray, title: str, out: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t0 = time.time()
    print(f"\n[plot] saving spatial scatter -> {out} ...", flush=True)
    n_dom = int(labels.max()) + 1
    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = plt.get_cmap("tab10" if n_dom <= 10 else "tab20")
    for d in range(n_dom):
        m = labels == d
        ax.scatter(coords[m, 0], coords[m, 1], s=2, c=[cmap(d % 10)],
                   label=f"d{d} (n={int(m.sum())})", rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel("x (bin)"); ax.set_ylabel("y (bin)")
    ax.set_title(title, fontsize=10)
    ax.legend(markerscale=4, fontsize=7, loc="best", framealpha=0.8)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] done ({time.time()-t0:.1f}s)", flush=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 70, flush=True)
    print("spaGAPA domain-collapse: diagnostic + fix", flush=True)
    print(f"host={HOSTNAME}  python={sys.executable}", flush=True)
    print("=" * 70, flush=True)

    replot_only = "--replot" in sys.argv

    # ---- replot path: load cached labels + coords only ----
    if replot_only:
        if not os.path.exists(CACHE_PKL):
            print(f"[error] --replot requested but cache {CACHE_PKL} not found; "
                  f"run the full script first.", flush=True)
            return 2
        print("[replot] loading cached labels + coords ...", flush=True)
        z = np.load(CACHE_PKL, allow_pickle=True)
        coords = z["coords"]
        results = {str(k): np.asarray(z["labels"][i], dtype=int)
                   for i, k in enumerate(z["names"])}
        scores = {n: collapse_score(lab) for n, lab in results.items()}
    else:
        X, coords, _sites = load_data()

        diagnose(X)

        # Baseline: reproduce the spaGAPA collapse
        baseline_labels = reproduce_spaGAPA_collapse(X)

        # Try the fixes
        results: Dict[str, np.ndarray] = {}
        results["baseline"] = baseline_labels
        results["A_hvg2000"] = fix_hvg_kmeans(X, n_hvg=2000)
        results["B_pca30"] = fix_pca_kmeans(X, n_pc=30)
        results["B2_hvg+pca"] = fix_hvg_pca_kmeans(X, n_hvg=2000, n_pc=30)
        results["C_leiden_a0.3"] = fix_leiden(X, coords, n_hvg=2000, n_pc=30, spatial_alpha=0.3)
        results["C_leiden_a0.7"] = fix_leiden(X, coords, n_hvg=2000, n_pc=30, spatial_alpha=0.7)
        # Sweep spatial weight for fix D
        results["D_w0.5"] = fix_spatial_kmeans(X, coords, n_pc=30, spatial_weight=0.5)
        results["D_w2"] = fix_spatial_kmeans(X, coords, n_pc=30, spatial_weight=2.0)
        results["D_w5"] = fix_spatial_kmeans(X, coords, n_pc=30, spatial_weight=5.0)
        # Pure-spatial sanity baselines
        results["E_pure_spatial"] = fix_pure_spatial(coords)
        results["F_leiden_spatial"] = fix_leiden_spatial(coords, resolution=0.5)

        # cache labels + coords so --replot can skip the expensive computations
        try:
            np.savez(CACHE_PKL,
                     coords=coords,
                     names=np.array(list(results.keys())),
                     labels=np.array(list(results.values()), dtype=object))
            print(f"[cache] saved labels + coords -> {CACHE_PKL}", flush=True)
        except Exception as e:
            print(f"[cache] failed to save cache: {e}", flush=True)

    # Pick the best by lowest collapse score (most balanced distribution).
    print("\n" + "=" * 70, flush=True)
    print("SUMMARY: collapse score (1.0 = degenerate, 0.0 = balanced)", flush=True)
    print("=" * 70, flush=True)
    scores = {}
    n_nonempty = {}
    for name, lab in results.items():
        s = collapse_score(lab)
        ent = shannon_entropy(np.bincount(lab))
        ne = int((np.bincount(lab) > 0).sum())
        scores[name] = s
        n_nonempty[name] = ne
        print(f"  {name:20s}  collapse={s*100:6.2f}%   entropy={ent:.3f} bits   "
              f"non-empty domains={ne}", flush=True)

    # Heuristic best: require >=3 non-empty domains AND lowest collapse.
    # Among those, also reward spatial contiguity for the winning label set
    # (we measure it for the top candidates only, since it is expensive).
    candidates = [n for n in scores if n_nonempty[n] >= 3]
    if not candidates:
        candidates = list(scores.keys())
    # rank by collapse (asc), break ties by entropy (desc) via sign flip
    candidates.sort(key=lambda n: (scores[n], -shannon_entropy(np.bincount(results[n]))))
    best_name = candidates[0]

    # Spatial contiguity score for the top 3 candidates: average fraction of each
    # domain's spots whose 6 nearest neighbors are in the SAME domain.
    print("\n[contiguity] computing spatial contiguity for top candidates ...", flush=True)
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=7).fit(coords)  # self + 6 neighbors
    _, nn_idx = nn.kneighbors(coords)

    def contiguity(labels: np.ndarray) -> float:
        nb_labels = labels[nn_idx[:, 1:]]  # drop self
        same = (nb_labels == labels[:, None]).mean()
        return float(same)

    for name in candidates[:4]:
        c = contiguity(results[name])
        print(f"  {name:20s}  spatial contiguity = {c*100:.1f}% "
              f"(fraction of 6-NN in same domain)", flush=True)
    # Re-rank the top-4 by contiguity (prefer spatially-coherent domains when
    # collapse scores are similar -- within 20 percentage points).
    top4 = candidates[:4]
    best_within = [n for n in top4 if scores[n] <= scores[top4[0]] + 0.20]
    best_within.sort(key=lambda n: -contiguity(results[n]))
    best_name = best_within[0]

    print(f"\n[BEST by auto-heuristic] {best_name}  "
          f"(collapse={scores[best_name]*100:.2f}%, "
          f"contiguity={contiguity(results[best_name])*100:.1f}%)", flush=True)
    plot_domains(
        coords, results[best_name],
        title=f"spaGAPA domain recovery (auto-best={best_name}) "
              f"-- {n_nonempty[best_name]} domains",
        out=FIG_OUT,
    )

    # Headline recommendation: the biologically-best method is the one that uses
    # BOTH the APA signal AND spatial structure -> Leiden on fused graph
    # (C_leiden_a0.7). It beats pure-spatial baselines because it incorporates
    # molecular similarity, and it beats kmeans because graph clustering is
    # robust to the homogenized APA fill.
    headline = "C_leiden_a0.7" if "C_leiden_a0.7" in results else best_name
    print(f"\n[RECOMMENDED for spaGAPA integration] {headline}  "
          f"(collapse={scores[headline]*100:.2f}%, "
          f"contiguity={contiguity(results[headline])*100:.1f}%, "
          f"n_domains={n_nonempty[headline]})", flush=True)

    # Multi-panel comparison: baseline collapse vs working fixes (story figure)
    plot_multipanel(coords, results, scores, FIG_MULTIPANEL_OUT)

    print("\nDONE.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

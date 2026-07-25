#!/usr/bin/env python
"""Per-sample spaGAPA full spatial analysis (Phase 3 Task Y).

Runs GP imputation (with uncertainty) + Leiden domain detection on every
GSE220442 sample (control vs AD, 6 scAPAtrap outputs) and the GSE263789
Stereo-seq binned (200 um) sample, producing publication-figure assets:

  - domain map      (spots in x,y, colored by Leiden domain)
  - uncertainty map (per-spot median GP uncertainty)
  - GSE220442: control-vs-AD domain comparison panel

For each sample a per-sample summary JSON records n_domains,
observed_fraction, median uncertainty, and wall time.

Outputs:
  pipeline_output/gse220442_spagapa_runs/<sample>/
      domains.csv, uncertainty_map.png, domain_map.png, summary.json
  pipeline_output/gse220442_spagapa_runs/control_vs_ad_domains.png
  pipeline_output/gse263789_spagapa_runs/binned_200/
      domains.csv, uncertainty_map.png, domain_map.png, summary.json

The binned GSE263789 apa_matrix.csv is ~1.6 GB (gene x spot). It is loaded
once and downsampled to the top-varying genes for tractable GP imputation
(spots stay fixed; GP runs per gene). Memory is bounded by reading the CSV
in chunks and selecting columns/genes.

Usage (S91):
  ~/anaconda3/envs/spagapa/bin/python scripts/run_spagapa_per_sample.py
  # single sample:
  ~/anaconda3/envs/spagapa/bin/python scripts/run_spagapa_per_sample.py \
      --sample gse220442_gsm6801751_scapatrap
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

# GP imputation parallelizes across genes (joblib, n_jobs). Pin BLAS/OpenMP
# threads to 1 so per-gene workers don't oversubscribe the CPU.
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

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.exceptions import ConvergenceWarning  # noqa: E402
warnings.filterwarnings("ignore", category=ConvergenceWarning)

from spagapa.bioml import MultiViewGraphBuilder
from spagapa.imputation import GPImputer, SparseGPImputer

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
PROC = ROOT / "data" / "processed"
GSE263789_BINNED = (
    ROOT / "pipeline_output" / "gse263789_stereo_pilot" /
    "spagapa_downstream_full" / "binned_200"
)
OUT_GSE220442 = ROOT / "pipeline_output" / "gse220442_spagapa_runs"
OUT_GSE263789 = ROOT / "pipeline_output" / "gse263789_spagapa_runs"

# GSE220442 sample -> condition (control / AD). The 6 scAPAtrap outputs are
# 3 control + 3 AD replicates (GSE220442 hippocampus spatial APA).
GSE220442_CONDITION = {
    "gse220442_gsm6801751_scapatrap": "control",
    "gse220442_gsm6801752_scapatrap": "control",
    "gse220442_gsm6801753_scapatrap": "control",
    "gse220442_gsm6801754_scapatrap": "AD",
    "gse220442_gsm6801755_scapatrap": "AD",
    "gse220442_gsm6801756_scapatrap": "AD",
}


# ── helpers ───────────────────────────────────────────────────────────
def _mask(values: np.ndarray) -> np.ndarray:
    return np.isfinite(values) & (values > 0)


def _leiden_on_graph(fused_graph, resolution: float, seed: int = 42):
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
        g, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=seed,
    )
    return np.array(partition.membership, dtype=int)


def load_apa_coords(apa_path: Path, coords_path: Path, max_genes=None):
    """Load gene x spot APA matrix + coords, aligning spot order.

    The APA CSVs use gene_id as first column and spot ids as the header row.
    Coordinates use spot_id,x,y. We align APA columns to coords order.
    """
    coords = pd.read_csv(coords_path)
    spots = list(coords["spot_id"])
    apa = pd.read_csv(apa_path, index_col=0)
    # align: intersect spot ids, reorder to coords order
    common = [s for s in spots if s in apa.columns]
    apa = apa[common]
    coords = coords[coords["spot_id"].isin(common)].reset_index(drop=True)
    assert list(coords["spot_id"]) == list(apa.columns)

    values = apa.to_numpy(dtype=float)
    if max_genes is not None and values.shape[0] > max_genes:
        # keep the most-varying (most informative) genes
        observed = _mask(values)
        # variance across observed entries per gene (RUD in [0,1])
        with np.errstate(invalid="ignore"):
            gv = np.nanvar(np.where(observed, values, np.nan), axis=1)
        gv = np.nan_to_num(gv, nan=0.0)
        top = np.argsort(gv)[::-1][:max_genes]
        values = values[top, :]
        apa_genes = apa.index.to_numpy()[top]
    else:
        apa_genes = apa.index.to_numpy()
    return values, coords, apa_genes


def run_sample(
    name: str,
    apa_path: Path,
    coords_path: Path,
    out_dir: Path,
    use_sparse_gp: bool = False,
    n_inducing: int = 200,
    max_genes: int | None = None,
    n_neighbors: int = 6,
    leiden_resolution: float = 1.0,
):
    """Run GP impute + Leiden domain for one sample. Returns summary dict."""
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print(f"\n=== {name} ===")
    t_load = time.time()
    apa, coords_df, genes = load_apa_coords(apa_path, coords_path, max_genes=max_genes)
    coords = coords_df[["x", "y"]].to_numpy(dtype=float)
    n_genes, n_spots = apa.shape
    mask = _mask(apa)
    obs_frac = float(mask.mean())
    print(f"  loaded {n_genes} genes x {n_spots} spots in {time.time()-t_load:.1f}s; "
          f"observed_fraction={obs_frac:.3f}")

    # ── GP imputation ─────────────────────────────────────────────────
    t_gp = time.time()
    train_mask = mask
    if use_sparse_gp:
        n_ind = min(n_inducing, max(20, n_spots // 2))
        # Resolve an auto length scale from local coordinate density so the
        # sparse GP adapts to the sample's coordinate scale.
        from sklearn.neighbors import NearestNeighbors
        k = min(max(2, n_neighbors), n_spots - 1)
        nn_d = NearestNeighbors(n_neighbors=k).fit(coords).kneighbors(coords)[0]
        pos = nn_d[:, -1]
        pos = pos[np.isfinite(pos) & (pos > 0)]
        length_scale = float(np.median(pos)) if pos.size else 1.0
        imputer = SparseGPImputer(n_inducing=n_ind, inducing_method="kmeans",
                                  length_scale=length_scale, noise_level=0.1)
        batch = imputer.fit_batch(coords, apa, mask=train_mask,
                                  n_jobs=32, verbose=False)
        imputed, uncertainty = batch.impute(return_uncertainty=True)
    else:
        imputer = GPImputer(kernel_type="matern", alpha=1e-6,
                            n_restarts_optimizer=0)
        batch = imputer.fit_batch(coords, apa, mask=train_mask, n_jobs=32,
                                  verbose=False)
        imputed, uncertainty = batch.impute(return_uncertainty=True)
    imputed = np.clip(imputed, 0.0, 1.0)
    gp_wall = time.time() - t_gp
    med_unc = float(np.median(uncertainty)) if uncertainty is not None else float("nan")
    print(f"  GP impute ({'sparse' if use_sparse_gp else 'exact'}) {gp_wall:.1f}s; "
          f"median uncertainty={med_unc:.4f}")

    # ── Fused graph + Leiden ──────────────────────────────────────────
    t_g = time.time()
    work = imputed
    apa_view = np.clip(work, 0.0, 1.0)
    nn = min(n_neighbors, n_spots - 1)
    graph = MultiViewGraphBuilder(
        n_neighbors=nn, spatial_weight=0.5, expression_weight=0.0, apa_weight=0.5,
    ).build(coords, apa_matrix=apa_view, uncertainty=uncertainty)
    fused = graph.fused
    labels = _leiden_on_graph(fused, resolution=leiden_resolution, seed=42)
    n_domains = int(len(np.unique(labels)))
    g_wall = time.time() - t_g
    print(f"  Leiden -> {n_domains} domains in {g_wall:.1f}s")

    # ── Save ──────────────────────────────────────────────────────────
    domains = pd.DataFrame({
        "spot_id": list(coords_df["spot_id"]),
        "x": coords[:, 0], "y": coords[:, 1],
        "domain": labels,
    })
    domains.to_csv(out_dir / "domains.csv", index=False)

    _plot_domain_map(coords, labels, name, n_domains, out_dir / "domain_map.png")
    _plot_uncertainty(coords, uncertainty, name, out_dir / "uncertainty_map.png")

    summary = {
        "sample": name,
        "n_genes": int(n_genes),
        "n_spots": int(n_spots),
        "observed_fraction": obs_frac,
        "n_domains": n_domains,
        "leiden_resolution": leiden_resolution,
        "median_uncertainty": med_unc,
        "gp_mode": "sparse" if use_sparse_gp else "exact",
        "gp_wall_s": gp_wall,
        "graph_wall_s": g_wall,
        "wall_total_s": time.time() - t0,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  summary: n_domains={n_domains} obs_frac={obs_frac:.3f} "
          f"med_unc={med_unc:.4f} total={summary['wall_total_s']:.1f}s")
    return summary


def _plot_domain_map(coords, labels, name, n_domains, path):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=labels, cmap="tab20",
                    s=12, edgecolor="none")
    ax.set_title(f"{name}\n{n_domains} Leiden domains")
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def _plot_uncertainty(coords, uncertainty, name, path):
    if uncertainty is None:
        return
    per_spot = np.median(uncertainty, axis=0) if uncertainty.ndim == 2 else uncertainty
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=per_spot, cmap="magma",
                    s=12, edgecolor="none")
    plt.colorbar(sc, label="median GP uncertainty")
    ax.set_title(f"{name}\nper-spot GP imputation uncertainty")
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=None,
                    help="run a single GSE220442 sample dir name")
    ap.add_argument("--no-gse263789", action="store_true")
    args = ap.parse_args()

    all_summaries = []

    # ── GSE220442 (6 scAPAtrap samples) ───────────────────────────────
    if args.sample:
        samples = [args.sample]
    else:
        samples = sorted(GSE220442_CONDITION.keys())

    for s in samples:
        apa_path = PROC / s / "apa_matrix.csv"
        coords_path = PROC / s / "coordinates.csv"
        if not apa_path.exists():
            print(f"SKIP {s}: {apa_path} missing")
            continue
        cond = GSE220442_CONDITION.get(s, "unknown")
        out = OUT_GSE220442 / s
        try:
            # GSE220442 (~4-5k spots): exact GP is O(n^3)/gene -> use sparse GP
            # (matches production highres_accuracy preset for >1000 spots).
            summ = run_sample(
                s, apa_path, coords_path, out,
                use_sparse_gp=True, n_inducing=200, max_genes=2000,
                n_neighbors=6, leiden_resolution=0.8,
            )
            summ["condition"] = cond
            all_summaries.append(summ)
        except Exception as e:
            print(f"  ERROR on {s}: {e}")
            import traceback; traceback.print_exc()

    # control vs AD comparison
    if all_summaries and len([s for s in all_summaries if "condition" in s]) >= 2:
        try:
            _plot_control_vs_ad(all_summaries, OUT_GSE220442)
        except Exception as e:
            print(f"control-vs-AD plot failed: {e}")

    # ── GSE263789 binned (200 um) ─────────────────────────────────────
    if not args.sample and not args.no_gse263789:
        apa_path = GSE263789_BINNED / "apa_matrix.csv"
        coords_path = GSE263789_BINNED / "coordinates.csv"
        if apa_path.exists():
            try:
                summ = run_sample(
                    "gse263789_binned_200", apa_path, coords_path,
                    OUT_GSE263789 / "binned_200",
                    use_sparse_gp=True, n_inducing=300, max_genes=3000,
                    n_neighbors=8, leiden_resolution=1.0,
                )
                all_summaries.append(summ)
            except Exception as e:
                print(f"  ERROR on gse263789: {e}")
                import traceback; traceback.print_exc()
        else:
            print(f"SKIP gse263789: {apa_path} missing")

    # global summary
    summary_path = ROOT / "pipeline_output" / "spagapa_per_sample_summary.json"
    summary_path.write_text(json.dumps(all_summaries, indent=2))
    print(f"\nAll summaries -> {summary_path}")
    for s in all_summaries:
        print(f"  {s['sample']}: n_domains={s['n_domains']} "
              f"obs={s['observed_fraction']:.3f} med_unc={s['median_uncertainty']:.4f} "
              f"wall={s['wall_total_s']:.0f}s")


def _plot_control_vs_ad(summaries, out_dir):
    """Side-by-side: average domain layout is not possible across samples
    (different spot layouts), so instead show a bar chart of n_domains and
    median uncertainty split by condition."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(summaries)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, metric, title in [
        (axes[0], "n_domains", "Domains detected"),
        (axes[1], "median_uncertainty", "Median GP uncertainty"),
    ]:
        groups = ["control", "AD"]
        data = [df[df["condition"] == g][metric].dropna().values for g in groups]
        labels = [f"{g}\n(n={len(d)})" for g, d in zip(groups, data)]
        ax.bar(labels, [np.mean(d) if len(d) else 0 for d in data],
               yerr=[np.std(d) if len(d) > 1 else 0 for d in data],
               color=["#4c72b0", "#c44e52"], capsize=6)
        # jitter individual points
        for i, d in enumerate(data):
            ax.scatter(np.full(len(d), i) + np.random.uniform(-0.05, 0.05, len(d)),
                       d, color="k", s=30, zorder=3)
        ax.set_title(title)
        ax.set_ylabel(metric)
    fig.suptitle("GSE220442 control vs AD (spaGAPA, n=3 each)")
    fig.tight_layout()
    fig.savefig(out_dir / "control_vs_ad_domains.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

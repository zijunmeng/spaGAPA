#!/usr/bin/env python3
"""Step 3-5: AD vs WT Stereo-seq differential APA + per-condition spaGAPA.

Steps
-----
(3) Build gene-level distal-usage APA index for AD and WT binned data, using
    the canonical median-split proximal/distal build_gene_index (same as the
    Visium differential + head-to-head benchmark).  Peaks are annotated to
    mm10 genes via the 10x mouse GTF; genes present in BOTH conditions are
    intersected.
(4) Differential APA (AD vs WT): pool AD binned spots vs WT binned spots,
    DifferentialAPAAnalyzer(method='t-test', min_spots_per_group=10), FDR.
    -> differential_apa_results.csv + volcano.png
(5) Per-condition spaGAPA: for AD and WT separately, SparseGPImputer
    (highres_fast: n_inducing=min(150, n_spots//200), length_scale=NN*5) +
    Leiden on the fused spatial+APA graph.
    -> {ad,wt}_domain_map.png + {ad,wt}_uncertainty_map.png

Inputs
------
AD  : pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200_raw/
       apa_matrix.csv (site x spot, RAW counts -- rebuilt from apa_site_counts.csv.gz
       because the original binned_200 stores usage fractions which can't feed
       build_gene_index), apa_sites.csv (no gene_name -> annotated here)
WT  : pipeline_output/gse263789_wt_control/binned_200/
       apa_matrix.csv (site x spot, raw counts), coordinates.csv,
       apa_sites.csv (no gene_name -> annotated here)

Output dir: pipeline_output/gse263789_ad_vs_wt_differential/
"""
from __future__ import annotations

import bisect
import gzip
import json
import os
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
AD_DIR = ROOT / "pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200_raw"
WT_DIR = ROOT / "pipeline_output/gse263789_wt_control/binned_200"
OUT = ROOT / "pipeline_output/gse263789_ad_vs_wt_differential"
OUT.mkdir(parents=True, exist_ok=True)

MOUSE_GTF = "/s1/SHARE/00_ref_genecode/refdata-gex-mm10-2020-A/genes/genes.gtf"

MIN_PARENT = 5          # build_gene_index threshold (matches Visium diff)
MIN_OBS_SPOTS = 10      # drop genes observed in < 10 spots


# ── GTF gene annotation (reused from run_scapatrap_spaceranger_mouse.py) ────
@dataclass
class GeneRecord:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    width: int


def parse_gtf_attrs(attr: str) -> dict:
    parsed = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def read_genes_from_gtf(gtf_path):
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    genes = {}
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_name = attrs.get("gene_name", gene_id)
            start = int(fields[3]); end = int(fields[4])
            rec = GeneRecord(fields[0], start, end, fields[6],
                             gene_id, gene_name, end - start + 1)
            genes.setdefault((rec.chrom, rec.strand), []).append(rec)
    for records in genes.values():
        records.sort(key=lambda r: r.start)
    return genes


def assign_genes(peaks: pd.DataFrame, genes: dict) -> pd.DataFrame:
    """Assign each peak to the smallest gene whose body contains its coord."""
    peaks = peaks.copy()
    if "coord" not in peaks.columns:
        peaks["coord"] = np.where(peaks["strand"].astype(str) == "+",
                                  peaks["end"], peaks["start"])
    rows = []
    for (chrom, strand), sub in peaks.groupby(["chr", "strand"], dropna=False):
        records = genes.get((str(chrom), str(strand)), [])
        starts = [r.start for r in records]
        active = []
        cursor = 0
        sub = sub.sort_values("coord")
        for idx, row in sub.iterrows():
            coord = int(row["coord"])
            upto = bisect.bisect_right(starts, coord, lo=cursor)
            active.extend(records[cursor:upto])
            cursor = upto
            active = [r for r in active if r.end >= coord]
            hits = [r for r in active if r.start <= coord <= r.end]
            best = min(hits, key=lambda r: r.width) if hits else None
            rows.append({"index": idx,
                         "gene_id": best.gene_id if best else np.nan,
                         "gene_name": best.gene_name if best else np.nan})
    anno = pd.DataFrame(rows).set_index("index")
    peaks[["gene_id", "gene_name"]] = anno.reindex(peaks.index)[["gene_id", "gene_name"]]
    return peaks


# ── canonical gene-level APA index (median-split proximal/distal) ───────────
def build_gene_index(counts: pd.DataFrame, sites: pd.DataFrame, min_parent: int):
    """Gene x spot distal-usage ratio. NaN where parent total < min_parent.

    Mirrors benchmark_stapaminer_headtohead.build_gene_index and
    analyze_gse220442_differential_apa.build_gene_index.
    """
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float),
                                 sites["coord"].astype(float))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    spot_cols = counts.columns
    vals = counts.values  # peak x spot
    gene_rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        rows = [peak_idx[r["site_id"]] for _, r in sub.iterrows()
                if r["site_id"] in peak_idx]
        if len(rows) < 2:
            continue
        oriented = np.array([
            float(sub.loc[(sub["site_id"] == counts.index[r]), "oriented"].iloc[0])
            for r in rows
        ])
        rows = np.array(rows)
        order = np.argsort(oriented)
        rows = rows[order]
        mid = len(rows) // 2
        if mid < 1:
            mid = 1
        prox_rows = rows[:mid]
        dist_rows = rows[mid:]
        prox_sum = vals[prox_rows].sum(axis=0)
        dist_sum = vals[dist_rows].sum(axis=0)
        total = prox_sum + dist_sum
        ratio = np.where(total >= min_parent,
                         dist_sum / np.where(total == 0, 1, total), np.nan)
        gene_rows[gene] = ratio
    if not gene_rows:
        raise ValueError("No multi-site genes with usable peaks")
    index = pd.DataFrame(gene_rows, index=spot_cols).T  # gene x spot
    index.index.name = "gene"
    return index


def load_binned(binned_dir: Path, label: str):
    """Load site x spot matrix + coords + sites; annotate peaks with GTF genes."""
    t = time.time()
    # matrix: site x spot (AD uses usage fraction, WT uses raw counts -- both fine
    # for build_gene_index since it sums per peak then takes a ratio).
    matrix_path = binned_dir / "apa_matrix.csv"
    # AD stores apa_sites as .csv.gz, WT as .csv; handle both.
    sites_path = binned_dir / "apa_sites.csv"
    if not sites_path.exists():
        sites_path = binned_dir / "apa_sites.csv.gz"
    coords_path = binned_dir / "coordinates.csv"

    counts = pd.read_csv(matrix_path, index_col=0)
    counts.columns = counts.columns.astype(str)
    opener = lambda p: gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)
    with opener(sites_path) as fh:
        sites = pd.read_csv(fh)
    sites["site_id"] = sites["site_id"].astype(str)
    sites = sites.set_index("site_id").loc[counts.index].reset_index()
    print(f"  [{label}] loaded matrix {counts.shape} sites {sites.shape} "
          f"in {time.time()-t:.1f}s", flush=True)
    return counts, sites, coords_path


def annotate_and_index(counts, sites, genes, label):
    """Annotate peaks with gene_name, then build gene x spot APA index."""
    t = time.time()
    # ensure required columns
    for col in ("chr", "strand", "coord"):
        if col not in sites.columns:
            raise ValueError(f"sites missing column {col}")
    before = sites["gene_name"].notna().sum() if "gene_name" in sites.columns else 0
    if "gene_name" not in sites.columns or before == 0:
        sites = assign_genes(sites, genes)
    after = sites["gene_name"].notna().sum()
    print(f"  [{label}] gene annotation: {before} -> {after} peaks annotated "
          f"in {time.time()-t:.1f}s", flush=True)

    t2 = time.time()
    idx = build_gene_index(counts, sites, MIN_PARENT)
    obs = np.isfinite(idx.values).sum(1)
    idx = idx.loc[idx.index[obs >= MIN_OBS_SPOTS]]
    print(f"  [{label}] gene index {idx.shape} (>= {MIN_OBS_SPOTS} obs spots) "
          f"in {time.time()-t2:.1f}s; obs frac="
          f"{np.isfinite(idx.values).mean():.3f}", flush=True)
    return idx, sites


# ── per-condition spaGAPA (GP impute + Leiden fused graph) ──────────────────
def run_per_condition_spagapa(name, gene_idx, coords_df, out_dir,
                              n_neighbors=6, leiden_resolution=1.0,
                              max_genes=1500):
    """SparseGP impute (highres_fast) + Leiden on fused spatial+APA graph.

    Mirrors run_spagapa_per_sample.run_sample but with the highres_fast
    preset: n_inducing = min(150, n_spots//200); length_scale = NN_dist * 5.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.neighbors import NearestNeighbors
    from spagapa.imputation import SparseGPImputer
    from spagapa.bioml import MultiViewGraphBuilder
    import igraph as ig
    import leidenalg

    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # align spots to coords order
    coords_df = coords_df.copy()
    coords_df["spot_id"] = coords_df["spot_id"].astype(str)
    common = [s for s in coords_df["spot_id"] if s in gene_idx.columns]
    gene_idx = gene_idx[common]
    coords_df = coords_df[coords_df["spot_id"].isin(common)].reset_index(drop=True)
    assert list(coords_df["spot_id"]) == list(gene_idx.columns)

    apa = gene_idx.to_numpy(dtype=float)
    coords = coords_df[["x", "y"]].to_numpy(dtype=float)
    n_genes_full, n_spots = apa.shape

    # keep the most-varying genes (cap for runtime)
    if max_genes is not None and apa.shape[0] > max_genes:
        observed = np.isfinite(apa) & (apa > 0)
        with np.errstate(invalid="ignore"):
            gv = np.nanvar(np.where(observed, apa, np.nan), axis=1)
        gv = np.nan_to_num(gv, nan=0.0)
        top = np.argsort(gv)[::-1][:max_genes]
        apa = apa[top, :]
        gene_names = gene_idx.index.to_numpy()[top]
        print(f"  [{name}] capped to top-{max_genes} varying genes "
              f"(of {n_genes_full})", flush=True)
    else:
        gene_names = gene_idx.index.to_numpy()

    mask = np.isfinite(apa) & (apa > 0)
    obs_frac = float(mask.mean())
    print(f"  [{name}] {apa.shape[0]} genes x {n_spots} spots "
          f"obs_frac={obs_frac:.3f}", flush=True)

    # ── GP impute (highres_fast) ──
    t_gp = time.time()
    n_ind = min(150, max(20, n_spots // 200))
    k = min(max(2, n_neighbors), n_spots - 1)
    nn_d = NearestNeighbors(n_neighbors=k).fit(coords).kneighbors(coords)[0]
    pos = nn_d[:, -1]
    pos = pos[np.isfinite(pos) & (pos > 0)]
    length_scale = float(np.median(pos)) * 5.0 if pos.size else 5.0
    print(f"  [{name}] SparseGP n_inducing={n_ind} "
          f"length_scale={length_scale:.1f} (=NN*5)", flush=True)
    imputer = SparseGPImputer(n_inducing=n_ind, inducing_method="kmeans",
                              length_scale=length_scale, noise_level=0.1)
    batch = imputer.fit_batch(coords, apa, mask=mask, n_jobs=16, verbose=False)
    imputed, uncertainty = batch.impute(return_uncertainty=True)
    imputed = np.clip(imputed, 0.0, 1.0)
    med_unc = float(np.median(uncertainty))
    print(f"  [{name}] GP impute {time.time()-t_gp:.0f}s med_unc={med_unc:.4f}",
          flush=True)

    # ── fused graph + Leiden ──
    t_g = time.time()
    apa_view = np.clip(imputed, 0.0, 1.0)
    nn = min(n_neighbors, n_spots - 1)
    graph = MultiViewGraphBuilder(
        n_neighbors=nn, spatial_weight=0.5, expression_weight=0.0, apa_weight=0.5,
    ).build(coords, apa_matrix=apa_view, uncertainty=uncertainty)
    fused = graph.fused
    coo = fused.tocoo()
    g_ig = ig.Graph(n=fused.shape[0],
                    edges=list(zip(coo.row.tolist(), coo.col.tolist())),
                    directed=False)
    g_ig.es["weight"] = coo.data.tolist()
    part = leidenalg.find_partition(
        g_ig, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=leiden_resolution, seed=42,
    )
    labels = np.array(part.membership, dtype=int)
    n_domains = int(len(np.unique(labels)))
    print(f"  [{name}] Leiden -> {n_domains} domains in {time.time()-t_g:.0f}s",
          flush=True)

    # ── plots ──
    domains_df = pd.DataFrame({
        "spot_id": list(coords_df["spot_id"]),
        "x": coords[:, 0], "y": coords[:, 1], "domain": labels,
    })
    domains_df.to_csv(out_dir / f"{name}_domains.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=labels, cmap="tab20",
                    s=8, edgecolor="none")
    ax.set_title(f"{name.upper()} Stereo-seq\n{n_domains} Leiden domains")
    ax.set_aspect("equal"); ax.set_xlabel("x (bin)"); ax.set_ylabel("y (bin)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{name}_domain_map.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    per_spot = np.median(uncertainty, axis=0)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=per_spot, cmap="magma",
                    s=8, edgecolor="none")
    plt.colorbar(sc, label="median GP uncertainty")
    ax.set_title(f"{name.upper()} Stereo-seq\nper-spot GP imputation uncertainty")
    ax.set_aspect("equal"); ax.set_xlabel("x (bin)"); ax.set_ylabel("y (bin)")
    fig.tight_layout()
    fig.savefig(out_dir / f"{name}_uncertainty_map.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "condition": name,
        "n_genes_input": int(n_genes_full),
        "n_genes_used": int(apa.shape[0]),
        "n_spots": int(n_spots),
        "observed_fraction": obs_frac,
        "n_domains": n_domains,
        "gp_n_inducing": int(n_ind),
        "gp_length_scale": round(length_scale, 2),
        "median_uncertainty": med_unc,
        "wall_s": round(time.time() - t0, 1),
    }
    (out_dir / f"{name}_spagapa_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  [{name}] done in {time.time()-t0:.0f}s", flush=True)
    return summary


def main():
    t0 = time.time()
    print("=== AD vs WT Stereo-seq differential APA + per-condition spaGAPA ===",
          flush=True)

    print("\n[annot] loading mouse GTF ...", flush=True)
    tg = time.time()
    genes = read_genes_from_gtf(MOUSE_GTF)
    n_genes_ref = sum(len(v) for v in genes.values())
    print(f"  GTF: {n_genes_ref} gene records in {time.time()-tg:.0f}s", flush=True)

    print("\n[step3] build gene-level APA index (AD, WT) ...", flush=True)
    ad_counts, ad_sites, ad_coords_path = load_binned(AD_DIR, "AD")
    wt_counts, wt_sites, wt_coords_path = load_binned(WT_DIR, "WT")

    ad_idx, ad_sites_ann = annotate_and_index(ad_counts, ad_sites, genes, "AD")
    wt_idx, wt_sites_ann = annotate_and_index(wt_counts, wt_sites, genes, "WT")

    # intersect genes present in both
    common_genes = sorted(set(ad_idx.index) & set(wt_idx.index))
    print(f"\n  common genes (AD & WT): {len(common_genes)}", flush=True)
    ad_idx_c = ad_idx.loc[common_genes]
    wt_idx_c = wt_idx.loc[common_genes]

    ad_idx_c.to_csv(OUT / "ad_gene_level_index.csv")
    wt_idx_c.to_csv(OUT / "wt_gene_level_index.csv")

    print("\n[step4] differential APA (AD vs WT, pooled spots, t-test) ...",
          flush=True)
    from spagapa.analysis import DifferentialAPAAnalyzer

    pooled = pd.concat([ad_idx_c, wt_idx_c], axis=1)
    cond = pd.concat([
        pd.Series("AD", index=ad_idx_c.columns),
        pd.Series("WT", index=wt_idx_c.columns),
    ])
    g_ad = np.where(cond.values == "AD")[0]
    g_wt = np.where(cond.values == "WT")[0]
    print(f"  pooled {pooled.shape} | AD spots={len(g_ad)} WT spots={len(g_wt)}",
          flush=True)

    ana = DifferentialAPAAnalyzer(method="t-test", min_spots_per_group=10)
    res = ana.test_differential_apa(pooled.values, g_ad, g_wt,
                                    gene_names=common_genes)
    res = ana.adjust_pvalues(res, method="fdr_bh")
    res = res.rename(columns={"mean_group1": "mean_AD", "mean_group2": "mean_WT"})
    res["delta_AD_minus_WT"] = res["mean_AD"] - res["mean_WT"]

    sig = res[(res["padj"] < 0.05) & (res["delta_AD_minus_WT"].abs() > 0.05)]
    res.to_csv(OUT / "differential_apa_results.csv", index=False)

    print(f"\n  genes tested: {len(res)}", flush=True)
    print(f"  padj<0.05 & |delta|>0.05: {len(sig)} "
          f"(proximal-shift AD<WT: {(sig['delta_AD_minus_WT']<0).sum()}, "
          f"distal-shift AD>WT: {(sig['delta_AD_minus_WT']>0).sum()})", flush=True)
    print(f"\n  top 20 differential APA genes:", flush=True)
    top = sig.reindex(sig["delta_AD_minus_WT"].abs()
                      .sort_values(ascending=False).index).head(20)
    if len(top):
        print(top[["gene", "mean_AD", "mean_WT", "delta_AD_minus_WT",
                   "log2fc", "padj"]].to_string(index=False), flush=True)
    else:
        # fall back: report top by padj regardless of effect threshold
        top = res.reindex(res["padj"].sort_values().index).head(20)
        print(top[["gene", "mean_AD", "mean_WT", "delta_AD_minus_WT",
                   "log2fc", "padj"]].to_string(index=False), flush=True)

    # volcano
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 6))
        ns = (res["padj"] >= 0.05) | (res["delta_AD_minus_WT"].abs() <= 0.05)
        ax.scatter(res.loc[ns, "delta_AD_minus_WT"],
                   -np.log10(res.loc[ns, "padj"].clip(lower=1e-300)),
                   s=6, c="#bdc3c7", alpha=0.5, label="ns")
        ax.scatter(sig["delta_AD_minus_WT"],
                   -np.log10(sig["padj"].clip(lower=1e-300)),
                   s=12, c="#c0392b", alpha=0.8, label="padj<0.05 & |Δ|>0.05")
        ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
        ax.axvline(0, ls="--", c="grey", lw=0.8)
        ax.axvline(-0.05, ls=":", c="grey", lw=0.6)
        ax.axvline(0.05, ls=":", c="grey", lw=0.6)
        ax.set_xlabel("Δ distal-usage (AD − WT)")
        ax.set_ylabel("-log10(padj)")
        ax.set_title("GSE263789 differential APA: AD vs WT Stereo-seq (bin200)",
                     fontweight="bold")
        ax.legend()
        for _, r in top.head(12).iterrows():
            ax.annotate(r["gene"], (r["delta_AD_minus_WT"],
                                    -np.log10(max(r["padj"], 1e-300))),
                        fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / "volcano.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  -> {OUT/'volcano.png'}", flush=True)
    except Exception as e:
        print(f"  (volcano skipped: {e})", flush=True)

    print("\n[step5] per-condition spaGAPA (AD) ...", flush=True)
    ad_coords = pd.read_csv(ad_coords_path)
    ad_summ = run_per_condition_spagapa("ad", ad_idx, ad_coords, OUT)

    print("\n[step5] per-condition spaGAPA (WT) ...", flush=True)
    wt_coords = pd.read_csv(wt_coords_path)
    wt_summ = run_per_condition_spagapa("wt", wt_idx, wt_coords, OUT)

    overall = {
        "n_genes_tested": int(len(res)),
        "n_significant": int(len(sig)),
        "proximal_shift_AD_lt_WT": int((sig["delta_AD_minus_WT"] < 0).sum()),
        "distal_shift_AD_gt_WT": int((sig["delta_AD_minus_WT"] > 0).sum()),
        "n_common_genes": int(len(common_genes)),
        "ad": ad_summ,
        "wt": wt_summ,
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "summary.json").write_text(json.dumps(overall, indent=2))
    print(f"\n=== ALL DONE {time.time()-t0:.0f}s ===", flush=True)
    print(json.dumps(overall, indent=2), flush=True)


if __name__ == "__main__":
    main()

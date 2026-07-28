#!/usr/bin/env python3
"""Expert-reviewer must-add result #3 (Fig 7F): Stereo-seq binning consistency.

Memory-efficient rewrite: representative APA sites are chosen first from the
existing bin200 matrix, then ONLY those sites are re-binned at bin50 / bin100 /
bin200 directly from the raw subcellular counts. This avoids ever holding a
full-resolution matrix (the original bin50 dense matrix was ~14 GB and its
npz compression was the bottleneck).

Protocol
--------
1. From the existing bin200 matrix, pick 8 representative APA sites = those
   with the highest spatial variance (truly spatially patterned), observed in
   many bins.
2. Stream the raw apa_site_counts.csv.gz ONCE; keep only rows whose peak_id is
   one of the representative sites.
3. For each bin size in {50, 100, 200}: aggregate the representative sites'
   counts onto the bin grid, convert to per-bin usage fraction (site_count /
   total peak count within bin), and align across the three resolutions on the
   nested (shared) coarse cells.
4. For each site and each pair of resolutions: Pearson r of the APA usage
   profile over shared bins. Also Moran's I at each resolution.
5. Plot: (a) spatial maps of up to 3 sites at the three bin sizes; (b) Moran's I
   vs bin size; (c) cross-bin Pearson r box plot.

Outputs (pipeline_output/stereo_binning_consistency/)
    binning_correlation.csv  gene, bin_size1, bin_size2, pearson_r,
                             morans_i_bin{50,100,200}, n_bins
    figure.png               publication quality, 3 panels, 300 DPI
    summary.json             cross-bin correlation, consistency assessment
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
RAW = Path(REPO) / "pipeline_output/gse263789_stereo_pilot/gsm8199179_full/scapatrap_raw"
BIN200 = Path(REPO) / "pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200"
OUT = Path(REPO) / "pipeline_output/stereo_binning_consistency"

BIN_SIZES = [50, 100, 200]
N_REP_GENES = 8
MORAN_K = 8
SEED = 42


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ===========================================================================
# Select representative sites from existing bin200 matrix
# ===========================================================================
def select_representative_sites():
    apa = pd.read_csv(BIN200 / "apa_matrix.csv", index_col=0)
    v = apa.values.astype(float)
    obs = np.isfinite(v).sum(axis=1)
    cand = np.where(obs >= 100)[0]
    var_sp = np.array([
        np.var(v[i][np.isfinite(v[i])]) if np.isfinite(v[i]).any() else 0.0
        for i in cand
    ])
    order = cand[np.argsort(-var_sp)]
    rep = order[:N_REP_GENES]
    site_ids = [str(x) for x in apa.index[rep]]
    log(f"  representative sites ({N_REP_GENES}): {site_ids}")
    return site_ids, apa, rep


# ===========================================================================
# Stream raw counts once, keep only representative-site rows.
# Build per-bin usage profiles for ALL three bin sizes from the same filtered
# rows. For each bin size we need both the site counts AND the column totals
# (over ALL peaks, not just representative) to compute usage fraction.
# So we DO need all-peak counts for the column-sum denominator. We compute the
# column sum at each bin size from the FULL stream, and the site rows from the
# filtered subset.
# ===========================================================================
def stream_and_aggregate(site_ids, bin_sizes):
    """One pass over apa_site_counts.csv.gz.

    Returns for each bin size b:
        site_usage[b] : dict site_id -> 1-D array (per-bin usage fraction)
        bin_centers[b]: (bx_centers, by_centers) arrays in column order
        col_keys[b]   : 1-D int64 array = bx*1e9+by per column (column order)
    The denominator (col_sums) is computed over ALL peaks.
    """
    site_set = set(site_ids)
    # peak -> row index map (for representative sites only)
    site_to_row = {s: i for i, s in enumerate(site_ids)}
    n_sites = len(site_ids)

    # accumulators per bin size
    # we aggregate (peak, bin) -> count. For representative sites we keep a
    # dict keyed by (site_row, bin_key) -> count. For col sums we keep a dict
    # keyed by bin_key -> total.  Using python dicts first, then densify.
    # But 97M rows -> too slow with pure python. Instead: vectorize per chunk.

    # Plan: read full csv, derive bin keys for all three sizes at once, then
    # use np.add.at style aggregation. To keep memory bounded, process in
    # chunks of rows.
    # For col_sums (over all peaks) we need every row's count -> full stream.
    # For site counts we only need representative-site rows.

    # We'll do it in chunks of 5M rows.
    log(f"  Streaming apa_site_counts.csv.gz in chunks (keeping rep-site rows + all-peak col sums) ...")
    t0 = time.time()

    # We accumulate per bin size:
    #   site_counts[b]: dict (site_row, col_code) -> float
    #   col_sum[b]    : dict col_code -> float
    # col_code = bx * BIG + by  (BIG chosen > max grid extent)
    site_counts = {b: {} for b in bin_sizes}
    col_sums = {b: {} for b in bin_sizes}

    cols = ["peak_id", "spot_id", "count"]
    reader = pd.read_csv(
        RAW / "apa_site_counts.csv.gz",
        dtype={"peak_id": "string", "spot_id": "string", "count": np.int32},
        engine="c", chunksize=5_000_000,
    )
    n_rows = 0
    n_rep_rows = 0
    for chunk in reader:
        n_rows += len(chunk)
        # parse spot coords once
        sp = chunk["spot_id"].str.split("_", expand=True)
        cx = sp[0].astype(np.int64).to_numpy()
        cy = sp[1].astype(np.int64).to_numpy()
        counts = chunk["count"].astype(np.float64).to_numpy()
        pids = chunk["peak_id"].astype(str).to_numpy()
        is_rep = pd.Series(pids).isin(site_set).to_numpy()

        for b in bin_sizes:
            bx = cx // b
            by = cy // b
            code = bx * np.int64(10**9) + by
            # column sums over ALL peaks
            # use np.bincount on unique codes -> need mapping. Use pandas groupby.
            # Faster: use np.add.at on a dict via python loop is slow; use
            # pandas groupby on this chunk's codes, then merge into dict.
            tmp = pd.DataFrame({"c": code, "v": counts})
            cs = tmp.groupby("c")["v"].sum()
            for k, val in cs.items():
                col_sums[b][k] = col_sums[b].get(k, 0.0) + float(val)
            # representative site rows
            if is_rep.any():
                rep_codes = code[is_rep]
                rep_counts = counts[is_rep]
                rep_pids = pids[is_rep]
                rep_rows = np.array([site_to_row[p] for p in rep_pids])
                # aggregate (site_row, code) within chunk first to limit dict ops
                # combine site_row into a single big key: site_row * BIG2 + code
                # BIG2 = 10**13 (> max col_code ~ 10**9, and site_row < 8)
                key2 = rep_rows.astype(np.int64) * np.int64(10**13) + rep_codes
                tmp2 = pd.DataFrame({"k": key2, "v": rep_counts})
                agg2 = tmp2.groupby("k")["v"].sum()
                for k, val in agg2.items():
                    site_counts[b][k] = site_counts[b].get(k, 0.0) + float(val)
                n_rep_rows += int(is_rep.sum())
        if n_rows % 20_000_000 < 5_000_000:
            log(f"    ... {n_rows/1e6:.0f}M rows  ({time.time()-t0:.0f}s)  rep_rows={n_rep_rows}")
    log(f"  stream done: {n_rows/1e6:.1f}M rows, rep_rows={n_rep_rows}  ({time.time()-t0:.0f}s)")

    # densify per bin size: build column index (sorted), then per-site usage
    result = {}
    for b in bin_sizes:
        cs = col_sums[b]
        col_keys = np.array(sorted(cs.keys()), dtype=np.int64)
        col_vals = np.array([cs[k] for k in col_keys], dtype=np.float64)
        # decode centers
        bx = (col_keys // np.int64(10**9)).astype(np.int64)
        by = (col_keys % np.int64(10**9)).astype(np.int64)
        bx_centers = (bx * b + b // 2).astype(np.int64)
        by_centers = (by * b + b // 2).astype(np.int64)
        col_pos = {k: j for j, k in enumerate(col_keys.tolist())}
        n_cols = len(col_keys)

        site_usage = np.full((n_sites, n_cols), np.nan, dtype=np.float32)
        sc = site_counts[b]
        # group site_counts keys by site_row
        per_site = {i: {} for i in range(n_sites)}
        for k, val in sc.items():
            sr = int(k // np.int64(10**13))
            cc = int(k % np.int64(10**13))
            per_site[sr][cc] = val
        for sr in range(n_sites):
            d = per_site[sr]
            if not d:
                continue
            for cc, val in d.items():
                j = col_pos.get(cc)
                if j is not None and col_vals[j] > 0:
                    site_usage[sr, j] = val / col_vals[j]
        result[b] = {
            "site_usage": site_usage,
            "bx_centers": bx_centers,
            "by_centers": by_centers,
            "col_keys": col_keys,
            "col_vals": col_vals,
        }
        log(f"    bin{b}: n_bins={n_cols}, finite site-cells="
            f"{int(np.isfinite(site_usage).sum())}/{site_usage.size}")
    return result


# ===========================================================================
# Moran's I
# ===========================================================================
def morans_i(values, bx, by, k=MORAN_K):
    mask = np.isfinite(values)
    if mask.sum() < k + 5:
        return float("nan")
    v = values[mask].astype(float)
    coords = np.stack([bx[mask], by[mask]], axis=1).astype(float)
    if v.std() == 0:
        return float("nan")
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    k_use = min(k, len(coords) - 1)
    _, idx = tree.query(coords, k=k_use + 1)
    nbr = idx[:, 1:]
    n = len(v)
    z = v - v.mean()
    # vectorized cross term
    cross = (z[:, None] * z[nbr]).sum()
    wsum = float(n * k_use)
    S2 = float((z ** 2).sum())
    if wsum == 0 or S2 == 0:
        return float("nan")
    return float((n / wsum) * (cross / S2))


def coarse_index_from_center(center, size):
    """Convert a bin-CENTER coordinate to that bin's INDEX at resolution `size`.

    center = idx*size + size//2  ->  idx = (center - size//2) // size.
    """
    return (np.asarray(center, dtype=np.int64) - size // 2) // size


def aggregate_to_coarse(site_usage_fine, bx_fine, by_fine, fine_size, target):
    """Mean of finite fine-bin usages within each coarse (bx_idx, by_idx) cell.

    Uses bin INDICES (not centers) so the returned coarse cells are directly
    comparable to the native target-resolution cells (also converted to indices).

    Returns (coarse_usage, coarse_bx_idx, coarse_by_idx).
    """
    finite = np.isfinite(site_usage_fine)
    if finite.sum() < 2:
        return None, None, None
    # fine-bin index from center, then coarse index = fine_idx * fine_size // target
    fine_bx_idx = coarse_index_from_center(bx_fine[finite], fine_size)
    fine_by_idx = coarse_index_from_center(by_fine[finite], fine_size)
    cbx = (fine_bx_idx * fine_size) // target
    cby = (fine_by_idx * fine_size) // target
    vals = site_usage_fine[finite]
    key = cbx.astype(np.int64) * np.int64(10**9) + cby.astype(np.int64)
    tmp = pd.DataFrame({"k": key, "v": vals})
    agg = tmp.groupby("k")["v"].mean()
    keys = agg.index.to_numpy()
    means = agg.to_numpy()
    cbx_u = (keys // np.int64(10**9)).astype(np.int64)
    cby_u = (keys % np.int64(10**9)).astype(np.int64)
    return means, cbx_u, cby_u


# ===========================================================================
# Main
# ===========================================================================
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    site_ids, apa200, rep_idx200 = select_representative_sites()

    cache = OUT / "binned_aggregated.pkl"
    if cache.exists():
        import pickle
        log(f"  loading cached aggregated data {cache}")
        with open(cache, "rb") as fh:
            site_ids, binned = pickle.load(fh)
    else:
        binned = stream_and_aggregate(site_ids, BIN_SIZES)
        import pickle
        with open(cache, "wb") as fh:
            pickle.dump((site_ids, binned), fh, protocol=4)
        log(f"  cached aggregated data -> {cache}")

    # ---- cross-bin correlation: bin50 (native) aggregated up vs bin100/bin200 native ----
    cor_rows = []
    moran_records = {s: {} for s in site_ids}
    for si, sid in enumerate(site_ids):
        for b in BIN_SIZES:
            u = binned[b]["site_usage"][si]
            moran_records[sid][f"bin{b}"] = morans_i(
                u, binned[b]["bx_centers"], binned[b]["by_centers"])

    bx50 = binned[50]["bx_centers"]
    by50 = binned[50]["by_centers"]
    for si, sid in enumerate(site_ids):
        u50 = binned[50]["site_usage"][si]
        for target in (100, 200):
            agg, cbx_u, cby_u = aggregate_to_coarse(u50, bx50, by50, 50, target)
            if agg is None:
                continue
            unat = binned[target]["site_usage"][si]
            nat_bx = binned[target]["bx_centers"]
            nat_by = binned[target]["by_centers"]
            # native target-resolution cells -> coarse INDEX (== their own index
            # at resolution `target`), matching the aggregated coarse indices.
            nat_bx_idx = coarse_index_from_center(nat_bx, target)
            nat_by_idx = coarse_index_from_center(nat_by, target)
            key_nat = nat_bx_idx.astype(np.int64) * np.int64(10**9) + nat_by_idx.astype(np.int64)
            nat_lookup = {}
            finite_nat = np.isfinite(unat)
            for kk, val in zip(key_nat[finite_nat], unat[finite_nat]):
                nat_lookup[int(kk)] = float(val)
            xs, ys = [], []
            for j, key in enumerate(cbx_u * np.int64(10**9) + cby_u):
                kint = int(key)
                if kint in nat_lookup:
                    xs.append(float(agg[j]))
                    ys.append(nat_lookup[kint])
            xs = np.array(xs); ys = np.array(ys)
            if len(xs) >= 5 and xs.std() > 0 and ys.std() > 0:
                from scipy.stats import pearsonr
                r, _ = pearsonr(xs, ys)
                r = float(r)
            else:
                r = float("nan")
            cor_rows.append(dict(
                gene=sid, bin_size1=50, bin_size2=target,
                pearson_r=r, n_bins=int(len(xs)),
                morans_i_bin50=moran_records[sid]["bin50"],
                morans_i_bin100=moran_records[sid]["bin100"],
                morans_i_bin200=moran_records[sid]["bin200"]))

    df_cor = pd.DataFrame(cor_rows)
    df_cor.to_csv(OUT / "binning_correlation.csv", index=False)
    log(f"  wrote binning_correlation.csv ({len(df_cor)} rows)")

    valid_r = df_cor["pearson_r"].dropna().values
    mean_mi = {b: float(np.nanmean([moran_records[s][f"bin{b}"] for s in site_ids]))
               for b in BIN_SIZES}
    summary = {
        "sample": "GSE263789 GSM8199179 (mouse AD brain, Stereo-seq, scAPAtrap)",
        "bin_sizes": BIN_SIZES,
        "n_representative_sites": int(len(site_ids)),
        "representative_sites": site_ids,
        "selection_criterion": "top sites by spatial variance of usage across bins (bin200)",
        "n_pairs_tested": int(len(df_cor)),
        "mean_cross_bin_pearson_r": float(np.mean(valid_r)) if len(valid_r) else float("nan"),
        "median_cross_bin_pearson_r": float(np.median(valid_r)) if len(valid_r) else float("nan"),
        "min_cross_bin_pearson_r": float(np.min(valid_r)) if len(valid_r) else float("nan"),
        "max_cross_bin_pearson_r": float(np.max(valid_r)) if len(valid_r) else float("nan"),
        "frac_pairs_r_above_0.7": float(np.mean(valid_r > 0.7)) if len(valid_r) else float("nan"),
        "mean_morans_i": {f"bin{b}": mean_mi[b] for b in BIN_SIZES},
        "n_bins_per_resolution": {f"bin{b}": int(len(binned[b]["col_keys"])) for b in BIN_SIZES},
        "consistency_assessment": (
            "STRONG: APA spatial patterns are consistent across binning "
            "resolutions (mean cross-bin r > 0.7)." if np.nanmean(valid_r) > 0.7
            else ("MODERATE: cross-bin consistency is partial." if np.nanmean(valid_r) > 0.4
                  else "WEAK: results are sensitive to binning resolution.")),
    }
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    _plot(df_cor, moran_records, site_ids, binned, summary, OUT / "figure.png")
    log(f"  wrote figure.png and summary.json")
    log("\n=== DONE ===")


def _plot(df_cor, moran_records, site_ids, binned, summary, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("apa", ["#f7f7f7", "#4575b4", "#0d2a4a"])
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 4, hspace=0.55, wspace=0.35)

    show_sites = list(range(min(3, len(site_ids))))
    for row_i, si in enumerate(show_sites):
        sid = site_ids[si]
        for col_j, b in enumerate(BIN_SIZES):
            ax = fig.add_subplot(gs[row_i, col_j])
            u = binned[b]["site_usage"][si]
            bx = binned[b]["bx_centers"]
            by = binned[b]["by_centers"]
            finite = np.isfinite(u)
            sc = ax.scatter(bx[finite], by[finite], c=u[finite], s=1.0,
                            cmap=cmap, edgecolors="none", rasterized=True)
            ax.set_title(f"{sid}  bin{b}", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_aspect("equal")
            if col_j == 0:
                ax.set_ylabel(f"site {row_i+1}", fontsize=9)
            if row_i == 0 and col_j == 0:
                cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label("usage", fontsize=7)
                cbar.ax.tick_params(labelsize=6)

    ax_mi = fig.add_subplot(gs[0, 3])
    mi_arr = np.array([[moran_records[s][f"bin{b}"] for b in BIN_SIZES]
                       for s in site_ids])
    for k in range(mi_arr.shape[0]):
        ax_mi.plot(BIN_SIZES, mi_arr[k], "o-", color="#999999", alpha=0.6, ms=4)
    ax_mi.plot(BIN_SIZES, np.nanmean(mi_arr, axis=0), "o-",
               color="#D55E00", lw=2.5, ms=8, label="mean")
    ax_mi.set_xlabel("Bin size", fontsize=9)
    ax_mi.set_ylabel("Moran's I", fontsize=9)
    ax_mi.set_title("Spatial autocorrelation\nvs bin size", fontsize=10)
    ax_mi.legend(fontsize=8)
    ax_mi.tick_params(labelsize=8)
    ax_mi.grid(alpha=0.25)

    ax_r = fig.add_subplot(gs[1, 3])
    sub100 = df_cor[df_cor["bin_size2"] == 100]["pearson_r"].dropna().values
    sub200 = df_cor[df_cor["bin_size2"] == 200]["pearson_r"].dropna().values
    pos = np.array([1, 2])
    parts = [sub100, sub200]
    bp = ax_r.boxplot(parts, positions=pos, widths=0.5,
                      showfliers=True, patch_artist=True)
    for patch, col in zip(bp["boxes"], ["#0072B2", "#009E73"]):
        patch.set_facecolor(col); patch.set_alpha(0.5)
    # overlay individual points
    for j, arr in enumerate(parts):
        x = np.full_like(arr, pos[j], dtype=float) + np.random.default_rng(0).uniform(-0.05, 0.05, len(arr))
        ax_r.plot(x, arr, "ko", ms=4, alpha=0.7)
    ax_r.axhline(0.7, color="#D55E00", ls="--", lw=1, label="r=0.7")
    ax_r.set_xticks(pos)
    ax_r.set_xticklabels(["bin50→100", "bin50→200"], fontsize=8)
    ax_r.set_ylabel("Pearson r (cross-bin)", fontsize=9)
    ax_r.set_title("Cross-bin APA\nprofile correlation", fontsize=10)
    ax_r.set_ylim(-0.2, 1.02)
    ax_r.legend(fontsize=8)
    ax_r.tick_params(labelsize=8)
    ax_r.grid(alpha=0.25, axis="y")

    ax_cap = fig.add_subplot(gs[2, :])
    ax_cap.axis("off")
    meanr = summary["mean_cross_bin_pearson_r"]
    frac = summary["frac_pairs_r_above_0.7"]
    cap = (f"Cross-bin consistency: mean Pearson r = {meanr:.3f}  "
           f"|  fraction of site-pairs with r > 0.7 = {frac:.0%}\n"
           f"Moran's I (mean): bin50={summary['mean_morans_i']['bin50']:.2f}, "
           f"bin100={summary['mean_morans_i']['bin100']:.2f}, "
           f"bin200={summary['mean_morans_i']['bin200']:.2f}\n"
           f"n_bins per resolution: bin50={summary['n_bins_per_resolution']['bin50']:,}, "
           f"bin100={summary['n_bins_per_resolution']['bin100']:,}, "
           f"bin200={summary['n_bins_per_resolution']['bin200']:,}\n"
           f"Verdict: {summary['consistency_assessment']}")
    ax_cap.text(0.02, 0.5, cap, fontsize=10, va="center", ha="left",
                family="monospace")

    fig.suptitle("Stereo-seq binning consistency (GSE263789 AD brain, scAPAtrap PAS)",
                 fontsize=12, y=0.995)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

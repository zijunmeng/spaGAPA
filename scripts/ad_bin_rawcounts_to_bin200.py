#!/usr/bin/env python3
"""Step 1b: Re-bin AD (GSM8199179) raw scAPAtrap counts -> peak x binned-spot
RAW COUNT matrix (comparable to the WT binning).

The existing binned_200/apa_matrix.csv for AD stores per-bin USAGE FRACTIONS
(site_count / col_sum), which cannot feed build_gene_index (it sums raw counts
and thresholds total >= min_parent).  This script rebuilds the AD matrix from
the raw long-format apa_site_counts.csv.gz using the SAME binning as WT, so
both conditions carry summed raw UMI counts per (peak, bin).

Input (pipeline_output/gse263789_stereo_pilot/gsm8199179_full/scapatrap_raw/):
  - apa_site_counts.csv.gz: peak_id, spot_id(Cx_Cy), count  (~23M rows)
  - peaks_meta.csv.gz: peakID, chr, start, end, strand, coord (21,455 peaks)

Output (pipeline_output/gse263789_stereo_pilot/.../binned_200_raw/):
  - apa_matrix.csv : peak x binned-spot (raw UMI counts)
  - coordinates.csv: spot_id, x (bin), y (bin)
  - apa_sites.csv  : site_id, peakID, chr, start, end, strand, coord
  - qc_summary.json
"""
from __future__ import annotations
import gzip
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd

RAW = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/"
           "pipeline_output/gse263789_stereo_pilot/gsm8199179_full/scapatrap_raw")
OUT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/"
           "pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200_raw")
OUT.mkdir(parents=True, exist_ok=True)
BIN = 200


def main():
    t0 = time.time()
    counts_gz = RAW / "apa_site_counts.csv.gz"
    meta_gz = RAW / "peaks_meta.csv.gz"

    print(f"[1/4] reading peaks_meta ({meta_gz.name}) ...", flush=True)
    with gzip.open(meta_gz, "rt") as fh:
        peaks = pd.read_csv(fh)
    if "peakID" not in peaks.columns:
        peaks = peaks.rename(columns={peaks.columns[0]: "peakID"})
    peaks["peakID"] = peaks["peakID"].astype(str)
    peaks["site_id"] = peaks["peakID"]
    peaks["coord"] = np.where(peaks["strand"].astype(str) == "+",
                              peaks["end"], peaks["start"])
    valid_peaks = set(peaks["peakID"])
    print(f"      {len(peaks)} curated peaks", flush=True)

    # Stream apa_site_counts.csv.gz (COMMA-separated).  Filter to curated peaks,
    # bin Cx_Cy, SUM per (peak, bin).  Read in chunks to bound memory.
    print(f"[2/4] streaming {counts_gz.name} (CSV) -> aggregated peak-bin dict "
          f"(bin={BIN}, peaks filtered to curated set) ...", flush=True)
    t_stream = time.time()
    agg = {}
    n_rows = 0
    n_kept = 0
    reader = pd.read_csv(
        counts_gz,
        usecols=["peak_id", "spot_id", "count"],
        dtype={"peak_id": "category", "spot_id": "string", "count": np.int32},
        chunksize=2_000_000,
        engine="c",
    )
    for chunk in reader:
        n_rows += len(chunk)
        # filter to curated peaks
        peak_str = chunk["peak_id"].astype(str)
        keep = peak_str.isin(valid_peaks)
        if not keep.any():
            continue
        n_kept += int(keep.sum())
        sub = chunk.loc[keep].copy()
        sub_peaks = peak_str.loc[keep].to_numpy()
        # vectorized Cx_Cy parse via pandas .str
        spot_str = sub["spot_id"].astype(str)
        split = spot_str.str.split("_", expand=True)
        cx = pd.to_numeric(split[0], errors="coerce").astype(np.int64) // BIN
        cy = pd.to_numeric(split[1], errors="coerce").astype(np.int64) // BIN
        cnts = sub["count"].astype(np.int64).to_numpy()
        peaks_c = sub_peaks
        bx_arr = cx.to_numpy()
        by_arr = cy.to_numpy()
        for pk, b1, b2, cnt in zip(peaks_c, bx_arr, by_arr, cnts):
            spot = f"{int(b1)}_{int(b2)}"
            key = (pk, spot)
            agg[key] = agg.get(key, 0) + int(cnt)
        if n_rows % 10_000_000 < 2_000_000:
            print(f"      {n_rows:,} rows read, {n_kept:,} kept, "
                  f"{len(agg):,} peak-bin keys in {time.time()-t_stream:.0f}s",
                  flush=True)
    print(f"      {n_rows:,} rows -> {n_kept:,} kept -> "
          f"{len(agg):,} unique (peak,bin) keys in {time.time()-t_stream:.0f}s",
          flush=True)

    print("[3/4] building triplets + pivoting to peak x spot ...", flush=True)
    t_piv = time.time()
    peaks_arr = np.fromiter((k[0] for k in agg), dtype=object, count=len(agg))
    spots_arr = np.fromiter((k[1] for k in agg), dtype=object, count=len(agg))
    counts_arr = np.fromiter(agg.values(), dtype=np.int64, count=len(agg))
    del agg
    trips = pd.DataFrame({"peak": peaks_arr, "spot_id": spots_arr,
                          "count": counts_arr})
    del peaks_arr, spots_arr, counts_arr
    print(f"      {len(trips):,} (peak,spot) triplets; "
          f"{trips['peak'].nunique()} peaks, {trips['spot_id'].nunique()} bins",
          flush=True)

    matrix = trips.pivot_table(
        index="peak", columns="spot_id", values="count",
        aggfunc="sum", fill_value=0,
    ).astype(np.int64)
    # reindex to full curated peak set (rows with no counts -> all zero)
    matrix = matrix.reindex(index=peaks["peakID"].tolist(), fill_value=0)
    matrix.index.name = "site_id"
    print(f"      matrix {matrix.shape[0]} peaks x {matrix.shape[1]} bins, "
          f"nnz={(matrix.values!=0).sum():,}, pivot {time.time()-t_piv:.0f}s",
          flush=True)

    print("[4/4] writing outputs ...", flush=True)
    spots = matrix.columns.tolist()
    xs, ys = [], []
    for s in spots:
        a, b = s.split("_")
        xs.append(int(a)); ys.append(int(b))
    coords = pd.DataFrame({"spot_id": spots, "x": xs, "y": ys})

    matrix.to_csv(OUT / "apa_matrix.csv")
    coords.to_csv(OUT / "coordinates.csv", index=False)

    sites_out = peaks[["site_id", "peakID", "chr", "start", "end",
                       "strand", "coord"]]
    sites_out.to_csv(OUT / "apa_sites.csv", index=False)

    summary = {
        "dataset": "gse263789_ad_raw_binned",
        "platform": "Stereo-seq",
        "sample": "GSM8199179_AD",
        "bin_size": BIN,
        "n_peaks_curated": int(len(peaks)),
        "n_peaks_with_counts": int((matrix.values != 0).any(axis=1).sum()),
        "n_bins": int(matrix.shape[1]),
        "n_triples_collapsed": int(len(trips)),
        "n_input_rows": int(n_rows),
        "n_rows_kept_curated": int(n_kept),
        "matrix_nnz": int((matrix.values != 0).sum()),
        "matrix_sparsity": float((matrix.values == 0).mean()),
        "value_type": "raw_umi_count_per_bin",
        "apa_ready": True,
        "files": {
            "apa_matrix": "apa_matrix.csv",
            "coordinates": "coordinates.csv",
            "apa_sites": "apa_sites.csv",
        },
        "notes": [
            "Re-binned from raw apa_site_counts.csv.gz (NOT the usage-fraction matrix).",
            "Stereo-seq DNB spots binned on a 200x200 grid; spot_id = 'xbin_ybin'.",
            "Peaks filtered to the 21,455 curated peaks_meta set.",
            "Values = summed raw UMI counts per (peak, bin).",
            "Comparable to WT binned_200 (same raw-count convention).",
        ],
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "qc_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n=== DONE {time.time()-t0:.0f}s ===")
    print(f"  n_peaks={matrix.shape[0]} n_bins={matrix.shape[1]} "
          f"nnz={summary['matrix_nnz']:,}")
    print(f"  -> {OUT/'apa_matrix.csv'}")
    print(f"  -> {OUT/'coordinates.csv'}")
    print(f"  -> {OUT/'apa_sites.csv'}")
    print(f"  -> {OUT/'qc_summary.json'}")


if __name__ == "__main__":
    main()

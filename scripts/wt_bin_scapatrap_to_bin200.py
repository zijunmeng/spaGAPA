#!/usr/bin/env python3
"""Step 1+2: Parse WT (GSM8199181) scAPAtrap counts -> peak x binned-spot matrix.

Input (pipeline_output/gse263789_wt_control/scapatrap_raw/):
  - counts.tsv.gz: umi_tools `gene\\tcell\\tcount` (peak_id, Cx_Cy barcode, count).
    53.7M rows, 93,451 peaks, cell tag encodes Stereo-seq DNB position "Cx_Cy".
  - peaks.saf: peak_id, chr, start, end, strand (93,451 peaks; matches counts).

Bin size 200 (matching the AD GSM8199179 binned_200 pilot):
  x_bin = floor(Cx / 200), y_bin = floor(Cy / 200)
  spot_id = f"{x_bin}_{y_bin}"   (same convention as AD binned matrix)
  Aggregate counts per (peak, bin) by SUM.

Output (pipeline_output/gse263789_wt_control/binned_200/):
  - apa_matrix.csv : peak x binned-spot (raw UMI counts, integer)
  - coordinates.csv: spot_id, x (x_bin), y (y_bin)
  - apa_sites.csv  : site_id, peakID, chr, start, end, strand, coord
  - qc_summary.json

Memory: streamed with awk pre-aggregation + pandas pivot of the reduced
(peak, spot, count) triplets.  ~53.7M input rows collapse to peak x spot.
"""
from __future__ import annotations
import gzip
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd

RAW = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control/scapatrap_raw")
OUT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_wt_control/binned_200")
OUT.mkdir(parents=True, exist_ok=True)
BIN = 200


def main():
    t0 = time.time()
    counts_gz = RAW / "counts.tsv.gz"
    peaks_saf = RAW / "peaks.saf"

    print(f"[1/4] reading peaks.saf ({peaks_saf.name}) ...", flush=True)
    peaks = pd.read_csv(
        peaks_saf, sep="\t", header=None,
        names=["peakID", "chr", "start", "end", "strand"],
    )
    peaks["site_id"] = peaks["peakID"].astype(str)
    peaks["coord"] = np.where(
        peaks["strand"].astype(str) == "+", peaks["end"], peaks["start"]
    )
    print(f"      {len(peaks)} peaks", flush=True)

    # Stream counts.tsv.gz in Python: bin Cx_Cy on the fly and SUM per
    # (peak, bin) using a dict.  This collapses the 53.7M-row file into at most
    # (n_peaks * n_bins) entries.  Memory for the dict is bounded (~hundreds
    # of MB for ~10M unique peak-bin keys).
    print(f"[2/4] streaming {counts_gz.name} -> aggregated peak-bin dict "
          f"(bin={BIN}) ...", flush=True)
    t_stream = time.time()
    agg: dict[tuple[str, str], np.int64] = {}
    n_rows = 0
    with gzip.open(counts_gz, "rt") as fh:
        header = fh.readline()  # 'gene\tcell\tcount'
        for line in fh:
            n_rows += 1
            # line: peak \t Cx_Cy \t count
            tab1 = line.find("\t")
            tab2 = line.find("\t", tab1 + 1)
            peak = line[:tab1]
            cell = line[tab1 + 1:tab2]
            cnt = int(line[tab2 + 1:])
            # parse Cx_Cy -> bin
            us = cell.find("_")
            cx = int(cell[:us])
            cy = int(cell[us + 1:])
            spot = f"{cx // BIN}_{cy // BIN}"
            key = (peak, spot)
            agg[key] = agg.get(key, 0) + cnt
            if n_rows % 5_000_000 == 0:
                print(f"      {n_rows:,} rows, {len(agg):,} peak-bin keys "
                      f"in {time.time()-t_stream:.0f}s", flush=True)
    print(f"      {n_rows:,} rows -> {len(agg):,} unique (peak,bin) keys "
          f"in {time.time()-t_stream:.0f}s", flush=True)

    print("[3/4] building triplets DataFrame + pivoting to peak x spot ...",
          flush=True)
    t_piv = time.time()
    peaks_arr = np.fromiter((k[0] for k in agg), dtype=object, count=len(agg))
    spots_arr = np.fromiter((k[1] for k in agg), dtype=object, count=len(agg))
    counts_arr = np.fromiter(agg.values(), dtype=np.int64, count=len(agg))
    del agg
    trips = pd.DataFrame({"peak": peaks_arr, "spot_id": spots_arr,
                          "count": counts_arr})
    print(f"      {len(trips):,} (peak,spot) triplets; "
          f"{trips['peak'].nunique():,} peaks, {trips['spot_id'].nunique():,} bins",
          flush=True)

    # Restrict to peaks present in peaks.saf (they should all be, but be safe).
    valid_peaks = set(peaks["peakID"].astype(str))
    trips = trips[trips["peak"].isin(valid_peaks)]

    matrix = trips.pivot_table(
        index="peak", columns="spot_id", values="count",
        aggfunc="sum", fill_value=0,
    ).astype(np.int64)
    matrix.index.name = "site_id"
    print(f"      matrix {matrix.shape[0]} peaks x {matrix.shape[1]} bins, "
          f"nnz={(matrix.values!=0).sum():,}, "
          f"pivot {time.time()-t_piv:.0f}s", flush=True)

    print("[4/4] writing outputs ...", flush=True)
    # coordinates.csv : spot_id, x (bin), y (bin) in coordinate order
    spots = matrix.columns.tolist()
    xs, ys = [], []
    for s in spots:
        a, b = s.split("_")
        xs.append(int(a)); ys.append(int(b))
    coords = pd.DataFrame({"spot_id": spots, "x": xs, "y": ys})

    # apa_matrix.csv : peak x spot  (matches AD layout: site_id, spot...)
    matrix.to_csv(OUT / "apa_matrix.csv")
    coords.to_csv(OUT / "coordinates.csv", index=False)

    # apa_sites.csv : peak annotation table
    sites_out = peaks[["site_id", "peakID", "chr", "start", "end", "strand", "coord"]]
    sites_out.to_csv(OUT / "apa_sites.csv", index=False)

    summary = {
        "dataset": "gse263789_wt_control_binned",
        "platform": "Stereo-seq",
        "sample": "GSM8199181_WT",
        "bin_size": BIN,
        "n_peaks": int(matrix.shape[0]),
        "n_bins": int(matrix.shape[1]),
        "n_triples_collapsed": int(len(trips)),
        "n_input_rows": int(n_rows),
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
            "Stereo-seq WT DNB spots binned on a 200x200 grid.",
            "Cell tag 'Cx_Cy' parsed: x_bin=floor(Cx/200), y_bin=floor(Cy/200).",
            "spot_id = f'{x_bin}_{y_bin}' (matches AD GSM8199179 binned_200 convention).",
            "Values are summed raw UMI counts per (peak, bin).",
            "Source: scapatrap_raw/counts.tsv.gz (umi_tools count output).",
        ],
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "qc_summary.json").write_text(json.dumps(summary, indent=2))

    # no intermediate file to clean up (aggregation done in-memory)

    print(f"\n=== DONE {time.time()-t0:.0f}s ===")
    print(f"  n_peaks={matrix.shape[0]} n_bins={matrix.shape[1]} "
          f"nnz={summary['matrix_nnz']:,}")
    print(f"  -> {OUT/'apa_matrix.csv'}")
    print(f"  -> {OUT/'coordinates.csv'}")
    print(f"  -> {OUT/'apa_sites.csv'}")
    print(f"  -> {OUT/'qc_summary.json'}")


if __name__ == "__main__":
    main()

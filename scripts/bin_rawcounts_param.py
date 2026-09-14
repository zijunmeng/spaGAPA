#!/usr/bin/env python3
"""Parameterized raw-count re-binner: scAPAtrap apa_site_counts.csv.gz -> peak x bin raw UMI matrix.

Same logic as ad_bin_rawcounts_to_bin200.py but takes CLI args so it can serve
E4 / 3M_E1 / 3M_E2 (the original is hardcoded for the E3 pilot).

Usage:
  python bin_rawcounts_param.py --raw-dir <scapatrap_raw> --out-dir <binned_200_raw> --bin 200
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
from scipy import sparse


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True, help="scapatrap_raw dir with apa_site_counts.csv.gz + peaks_meta.csv.gz")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--bin", type=int, default=200)
    ap.add_argument("--chunk", type=int, default=5_000_000)
    a = ap.parse_args()
    RAW, OUT, BIN = Path(a.raw_dir), Path(a.out_dir), a.bin
    OUT.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    # 1. peaks_meta
    print("[1/4] reading peaks_meta ...", flush=True)
    peaks = pd.read_csv(RAW / "peaks_meta.csv.gz")
    peaks = peaks[["peakID", "chr", "start", "end", "strand", "coord"]].copy()
    peaks["site_id"] = peaks["peakID"]
    peak_ids = peaks["peakID"].tolist()
    peak_pos = {p: i for i, p in enumerate(peak_ids)}
    print(f"      {len(peaks)} curated peaks")

    # 2. stream counts -> dict[(peak, bin)] -> sum
    print("[2/4] streaming apa_site_counts.csv.gz ...", flush=True)
    agg = {}
    n_rows = 0
    cols = pd.read_csv(RAW / "apa_site_counts.csv.gz", nrows=0).columns.tolist()
    # expected: peak_id, spot_id (Cx_Cy), count — tolerate naming variants
    with gzip.open(RAW / "apa_site_counts.csv.gz", "rt") as fh:
        for chunk in pd.read_csv(fh, chunksize=a.chunk, names=None, header=0):
            n_rows += len(chunk)
            pid = chunk[cols[0]].map(peak_pos)
            ok = pid.notna()
            if not ok.any():
                continue
            sp = chunk.loc[ok, cols[1]].astype(str)
            cnt = chunk.loc[ok, cols[2]].astype(float)
            # spot "Cx_Cy" -> bin coords
            xy = sp.str.split("_", expand=True).astype(int)
            bx, by = xy[0] // BIN, xy[1] // BIN
            keys = pd.Series(list(zip(pid[ok].astype(int), bx * 10_000_000 + by)))
            s = keys.map(cnt.values).groupby(keys).sum() if False else None
            df = pd.DataFrame({"k": bx.values * 10_000_000 + by.values,
                               "p": pid[ok].astype(int).values,
                               "c": cnt.values})
            g = df.groupby(["p", "k"])["c"].sum()
            for (p, k), v in g.items():
                agg[(p, k)] = agg.get((p, k), 0.0) + float(v)
            if n_rows % 20_000_000 < a.chunk:
                print(f"      {n_rows:,} rows read, {len(agg):,} keys in {time.time()-t0:.0f}s", flush=True)
    print(f"      total {n_rows:,} rows -> {len(agg):,} peak-bin keys")

    # 3. sparse matrix
    print("[3/4] building sparse matrix ...", flush=True)
    if not agg:
        raise SystemExit("no peak-bin counts aggregated")
    pk = np.fromiter((k[0] for k in agg), dtype=np.int64)
    bk = np.fromiter((k[1] for k in agg), dtype=np.int64)
    vv = np.fromiter(agg.values(), dtype=np.float64)
    bin_x = (bk // 10_000_000).astype(np.int64)
    bin_y = (bk % 10_000_000).astype(np.int64)
    uniq_bins, bin_code = np.unique(bk, return_inverse=True)
    nP, nB = len(peak_ids), len(uniq_bins)
    M = sparse.csr_matrix((vv, (pk, bin_code)), shape=(nP, nB))
    nnz = M.nnz
    print(f"      {nP} peaks x {nB} bins, nnz={nnz:,}")

    # 4. write outputs
    print("[4/4] writing outputs ...", flush=True)
    # derive bin labels from uniq_bins codes
    ux = (uniq_bins // 10_000_000).astype(int)
    uy = (uniq_bins % 10_000_000).astype(int)
    bin_labels = [f"{x}_{y}" for x, y in zip(ux, uy)]
    order = np.argsort([int(a_) * 1_000_000 + int(b_) for a_, b_ in zip(ux, uy)])
    bin_labels_sorted = [bin_labels[i] for i in order]
    M_sorted = M[:, order]

    # apa_matrix.csv (row-chunked write)
    with open(OUT / "apa_matrix.csv", "w") as fh:
        fh.write("site_id," + ",".join(bin_labels_sorted) + "\n")
        arr = M_sorted.toarray()
        for i, pid in enumerate(peak_ids):
            row = arr[i]
            fh.write(pid + "," + ",".join(f"{v:g}" if v else "0" for v in row) + "\n")
        del arr
    pd.DataFrame({"spot_id": bin_labels_sorted,
                  "x": [int(s.split("_")[0]) for s in bin_labels_sorted],
                  "y": [int(s.split("_")[1]) for s in bin_labels_sorted]}
                 ).to_csv(OUT / "coordinates.csv", index=False)
    sites = peaks[["site_id", "peakID", "chr", "start", "end", "strand", "coord"]]
    sites.to_csv(OUT / "apa_sites.csv", index=False)
    json.dump({"n_peaks": int(nP), "n_bins": int(nB), "nnz": int(nnz),
               "value_type": "raw_umi_counts", "bin": BIN,
               "source": str(RAW)}, open(OUT / "qc_summary.json", "w"), indent=2)
    print(f"[DONE] {time.time()-t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generic Stereo-seq scAPAtrap output -> spaGAPA bin-N aggregation.

Parameterized version of scripts/ad_bin_rawcounts_to_bin200.py (GSE263789 AD
pilot) — works for any scAPAtrap raw output with:
  <raw_dir>/apa_site_counts.csv.gz   (peak_id, spot_id "Cx_Cy", count)
  <raw_dir>/peaks_meta.csv.gz        (peakID, chr, start, end, strand)

Outputs to <out_dir>:
  apa_matrix.csv   site x binned-spot raw UMI counts
  coordinates.csv  spot_id, x(bin), y(bin)
  apa_sites.csv    site_id, peakID, chr, start, end, strand, coord
  qc_summary.json

Usage:
  python stereo_bin200.py --raw-dir RAW --out-dir OUT \
      --dataset gse293464 --sample GSM8882884 [--bin 200]
"""
import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bin", type=int, default=200)
    args = ap.parse_args()
    raw, out, bin_ = Path(args.raw_dir), Path(args.out_dir), args.bin
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    counts_gz, meta_gz = raw / "apa_site_counts.csv.gz", raw / "peaks_meta.csv.gz"
    for f in (counts_gz, meta_gz):
        if not f.is_file():
            raise SystemExit(f"missing input: {f}")

    print(f"[1/4] reading peaks_meta ...", flush=True)
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

    print(f"[2/4] streaming counts -> (peak,bin) aggregation, bin={bin_} ...", flush=True)
    t_stream = time.time()
    agg = {}
    n_rows = n_kept = 0
    reader = pd.read_csv(
        counts_gz, usecols=["peak_id", "spot_id", "count"],
        dtype={"peak_id": "category", "spot_id": "string", "count": np.int32},
        chunksize=2_000_000, engine="c",
    )
    for chunk in reader:
        n_rows += len(chunk)
        peak_str = chunk["peak_id"].astype(str)
        keep = peak_str.isin(valid_peaks)
        if not keep.any():
            continue
        n_kept += int(keep.sum())
        sub = chunk.loc[keep]
        split = sub["spot_id"].astype(str).str.split("_", expand=True)
        bx = (pd.to_numeric(split[0], errors="coerce") // bin_).to_numpy()
        by = (pd.to_numeric(split[1], errors="coerce") // bin_).to_numpy()
        cnts = sub["count"].astype(np.int64).to_numpy()
        pks = peak_str.loc[keep].to_numpy()
        for pk, b1, b2, cnt in zip(pks, bx, by, cnts):
            key = (pk, f"{int(b1)}_{int(b2)}")
            agg[key] = agg.get(key, 0) + int(cnt)
        if n_rows % 10_000_000 < 2_000_000:
            print(f"      {n_rows:,} rows, {n_kept:,} kept, {len(agg):,} keys "
                  f"in {time.time()-t_stream:.0f}s", flush=True)
    print(f"      {n_rows:,} rows -> {len(agg):,} (peak,bin) keys "
          f"in {time.time()-t_stream:.0f}s", flush=True)

    print("[3/4] pivoting ...", flush=True)
    t_piv = time.time()
    trips = pd.DataFrame({
        "peak": [k[0] for k in agg],
        "spot_id": [k[1] for k in agg],
        "count": list(agg.values()),
    })
    del agg
    matrix = trips.pivot_table(index="peak", columns="spot_id", values="count",
                               aggfunc="sum", fill_value=0).astype(np.int64)
    matrix = matrix.reindex(index=peaks["peakID"].tolist(), fill_value=0)
    matrix.index.name = "site_id"
    print(f"      {matrix.shape[0]} peaks x {matrix.shape[1]} bins, "
          f"nnz={(matrix.values != 0).sum():,}, {time.time()-t_piv:.0f}s", flush=True)

    print("[4/4] writing outputs ...", flush=True)
    spots = matrix.columns.tolist()
    coords = pd.DataFrame({
        "spot_id": spots,
        "x": [int(s.split("_")[0]) for s in spots],
        "y": [int(s.split("_")[1]) for s in spots],
    })
    matrix.to_csv(out / "apa_matrix.csv")
    coords.to_csv(out / "coordinates.csv", index=False)
    peaks[["site_id", "peakID", "chr", "start", "end", "strand", "coord"]] \
        .to_csv(out / "apa_sites.csv", index=False)

    summary = {
        "dataset": args.dataset,
        "platform": "Stereo-seq",
        "sample": args.sample,
        "bin_size": bin_,
        "n_peaks_curated": int(len(peaks)),
        "n_peaks_with_counts": int((matrix.values != 0).any(axis=1).sum()),
        "n_bins": int(matrix.shape[1]),
        "n_triples_collapsed": int(len(trips)),
        "n_input_rows": int(n_rows),
        "matrix_nnz": int((matrix.values != 0).sum()),
        "matrix_sparsity": float((matrix.values == 0).mean()),
        "value_type": "raw_umi_count_per_bin",
        "apa_ready": True,
        "wall_seconds": round(time.time() - t0, 1),
        "notes": [
            "Binned from scAPAtrap apa_site_counts.csv.gz (raw UMI counts).",
            f"DNB spots binned on a {bin_}x{bin_} grid; spot_id = 'xbin_ybin'.",
        ],
    }
    with open(out / "qc_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"[done] {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()

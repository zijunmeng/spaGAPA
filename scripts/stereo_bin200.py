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
    peak_to_row = {p: i for i, p in enumerate(peaks["peakID"])}
    rows_acc, cols_acc, vals_acc = [], [], []
    n_rows = n_kept = 0
    usecols = ["peak_id", "spot_id", "count"]
    peek = pd.read_csv(counts_gz, nrows=2)
    if "site" in peek.columns:
        usecols = ["site", "spot", "count"]
    reader = pd.read_csv(counts_gz, usecols=usecols,
                         dtype={"peak_id": str, "spot_id": str, "count": np.int64,
                                "site": str, "spot": str},
                         chunksize=5_000_000, engine="c")
    for chunk in reader:
        if "site" in chunk.columns:
            chunk = chunk.rename(columns={"site": "peak_id", "spot": "spot_id"})
        n_rows += len(chunk)
        rc = chunk["peak_id"].map(peak_to_row)
        keep = rc.notna()
        if not keep.any():
            continue
        n_kept += int(keep.sum())
        sub = chunk.loc[keep]
        r = rc.loc[keep].astype(np.int32).to_numpy()
        split = sub["spot_id"].str.split("_", expand=True)
        bx = (pd.to_numeric(split[0], errors="coerce") // bin_).to_numpy()
        by = (pd.to_numeric(split[1], errors="coerce") // bin_).to_numpy()
        ok = np.isfinite(bx) & np.isfinite(by)
        rows_acc.append(r[ok])
        cols_acc.append(by[ok].astype(np.int64) * 1_000_000 + bx[ok].astype(np.int64))
        vals_acc.append(sub["count"].to_numpy()[ok])
        if n_rows % 50_000_000 < 5_000_000:
            print(f"      {n_rows:,} rows, {n_kept:,} kept "
                  f"in {time.time()-t_stream:.0f}s", flush=True)
    rows = np.concatenate(rows_acc); cols = np.concatenate(cols_acc)
    vals = np.concatenate(vals_acc)
    del rows_acc, cols_acc, vals_acc
    print(f"      {n_rows:,} rows read in {time.time()-t_stream:.0f}s", flush=True)
    df = pd.DataFrame({"r": rows, "c": cols, "v": vals})
    del rows, cols, vals
    g = df.groupby(["r", "c"], sort=False, as_index=False)["v"].sum()
    del df
    agg = g  # columns r/c/v — downstream factorize block consumes equivalent
    print(f"      {len(g):,} unique (peak,bin) pairs", flush=True)
    print("[3/4] building sparse matrix (scipy) + usage fractions ...", flush=True)
    t_piv = time.time()
    import scipy.sparse as sp
    agg_r = agg["r"].to_numpy(); agg_c = agg["c"].to_numpy(); agg_v = agg["v"].to_numpy()
    del agg
    bx_arr = (agg_c % 1_000_000).astype(np.int64)
    by_arr = (agg_c // 1_000_000).astype(np.int64)
    del agg_c
    composite = bx_arr * 4_000_000 + by_arr   # rat y-bin < 4e6；解码 bx=c//4e6, by=c%4e6
    spot_codes, spot_uniq = pd.factorize(pd.Series(composite), sort=True)
    M = sp.coo_matrix((agg_v, (agg_r.astype(np.int64), spot_codes)),
                      shape=(len(peaks), len(spot_uniq))).tocsr()
    M.sum_duplicates()
    del agg_r, agg_v, spot_codes
    # 行对齐到完整 curated peak 集（agg_r 已是该索引，直接构造全形状）
    Mfull = sp.csr_matrix((len(peaks), M.shape[1]), dtype=np.int64)
    Mfull[:, :] = M
    del M
    n_bins = Mfull.shape[1]
    print(f"      {Mfull.shape[0]} peaks x {n_bins} bins, nnz={Mfull.nnz:,}, "
          f"{time.time()-t_piv:.0f}s", flush=True)

    print("[4/4] writing outputs (usage CSV, row-chunked) ...", flush=True)
    uniq_codes = np.asarray(spot_uniq, dtype=np.int64)
    spot_ids = [f"{int(c // 4_000_000)}_{int(c % 4_000_000)}" for c in uniq_codes]
    coords = pd.DataFrame({
        "spot_id": spot_ids,
        "x": [int(s.split("_")[0]) for s in spot_ids],
        "y": [int(s.split("_")[1]) for s in spot_ids],
    })
    coords.to_csv(out / "coordinates.csv", index=False)
    peaks[["site_id", "peakID", "chr", "start", "end", "strand", "coord"]] \
        .to_csv(out / "apa_sites.csv", index=False)
    col_sums = np.asarray(Mfull.sum(axis=0)).ravel()
    sp.csr_matrix(Mfull, dtype=np.float64)  # no-op type hint
    with open(out / "apa_matrix.csv", "w") as fh:
        fh.write("site_id," + ",".join(spot_ids) + "\n")
        site_ids = peaks["peakID"].tolist()
        inv = np.where(col_sums > 0, 1.0 / np.maximum(col_sums, 1), np.nan)
        chunk = 512
        for s in range(0, Mfull.shape[0], chunk):
            e = min(s + chunk, Mfull.shape[0])
            blk = np.asarray(Mfull[s:e].todense(), dtype=np.float64)
            blk = blk * inv[None, :]
            for i in range(e - s):
                fh.write(site_ids[s + i] + "," +
                         ",".join("nan" if np.isnan(v) else f"{v:.6g}" for v in blk[i]) + "\n")
            if (s // chunk) % 10 == 0:
                print(f"      wrote {e}/{Mfull.shape[0]} rows", flush=True)
    nnz = int(Mfull.nnz)
    n_with_counts = int((np.asarray(Mfull.sum(axis=1)).ravel() > 0).sum())

    summary = {
        "dataset": args.dataset,
        "platform": "Stereo-seq",
        "sample": args.sample,
        "bin_size": bin_,
        "n_peaks_curated": int(len(peaks)),
        "n_peaks_with_counts": int(n_with_counts),
        "n_bins": int(n_bins),
        "n_triples_collapsed": nnz,
        "n_input_rows": int(n_rows),
        "matrix_nnz": nnz,
        "matrix_sparsity": float(1.0 - nnz / max(1, Mfull.shape[0] * n_bins)),
        "value_type": "per_bin_site_usage_fraction",
        "apa_ready": True,
        "wall_seconds": round(time.time() - t0, 1),
        "notes": [
            "Binned from scAPAtrap apa_site_counts.csv.gz; usage = count/col_sum per bin (pilot convention).",
            f"DNB spots binned on a {bin_}x{bin_} grid; spot_id = 'xbin_ybin'.",
        ],
    }
    with open(out / "qc_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"[done] {time.time()-t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""HD bin 聚合：scAPAtrap apa_site_counts (peak, spot=CB, count) → peak × 8µm-bin raw 计数矩阵。

HD 空间条码 → 2µm 网格坐标：取自 spaceranger 的 spatial/tissue_positions 或
binned_outputs 元数据。兼容回退：若坐标文件缺失，从 CB 结构解析（HD 条码不含坐标，
必须用 spaceranger 映射）。

用法: hd_bin8um.py --raw-dir <scapatrap_raw> --out-dir <binned_8um> --um 8
"""
from __future__ import annotations
import argparse
import gzip
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--um", type=int, default=8)
    ap.add_argument("--sr-dir", default=None, help="spaceranger out (找坐标映射)")
    a = ap.parse_args()
    RAW, OUT = Path(a.raw_dir), Path(a.out_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    BIN_UM = a.um
    HD_BIN_UM = 2  # HD 原生 2µm 网格
    factor = BIN_UM // HD_BIN_UM  # 4x4

    sr = Path(a.sr_dir) if a.sr_dir else RAW.parent / "sr" / "outs"
    # 坐标来源优先级: spatial/tissue_positions.parquet > csv > binned_outputs src
    pos = None
    for cand in [sr / "spatial/tissue_positions.parquet",
                 sr / "spatial/tissue_positions_list.csv",
                 sr / "spatial/tissue_positions.csv"]:
        if cand.exists():
            pos = cand; break
    if pos is None:
        raise SystemExit(f"缺坐标文件（需 spaceranger outs/spatial）: {sr}")
    if str(pos).endswith(".parquet"):
        tpos = pd.read_parquet(pos)
    else:
        tpos = pd.read_csv(pos, header=None,
                           names=["barcode", "in_tissue", "array_row", "array_col",
                                  "pxl_row", "pxl_col"], index_col=0) \
              if (pos.stat().st_size and open(pos).readline().count(",") == 4) \
              else pd.read_csv(pos, index_col=0)
    # HD: array_row/array_col = 2µm 网格索引
    bc2xy = {b: (int(r.array_row) // factor, int(r.array_col) // factor)
             for b, r in tpos.iterrows()}

    peaks = pd.read_csv(RAW / "peaks_meta.csv.gz")
    prow = peaks["peakID"].tolist()
    pidx = {p: i for i, p in enumerate(prow)}

    agg = {}
    n = 0
    t0 = time.time()
    with gzip.open(RAW / "apa_site_counts.csv.gz", "rt") as fh:
        for chunk in pd.read_csv(fh, chunksize=5_000_000):
            n += len(chunk)
            xy = chunk["spot_id"].map(bc2xy)
            ok = xy.notna()
            if not ok.any():
                continue
            p = chunk.loc[ok, "peak_id"].map(pidx)
            ko = xy[ok]
            keys = list(zip(p[~p.isna()].astype(int),
                            [a_ * 10_000_000 + b_ for a_, b_ in ko[~p.isna()]]))
            for (pi, k), v in zip(keys, chunk.loc[ok, "count"].values[~p.isna().values]):
                agg[(pi, k)] = agg.get((pi, k), 0) + int(v)
            if n % 20_000_000 < 5_000_000:
                print(f"  {n:,} rows, {len(agg):,} keys, {time.time()-t0:.0f}s", flush=True)
    print(f"total {n:,} rows -> {len(agg):,} keys")

    from scipy import sparse
    pk = np.fromiter((k[0] for k in agg), dtype=np.int64)
    bk = np.fromiter((k[1] for k in agg), dtype=np.int64)
    vv = np.fromiter(agg.values(), dtype=np.float64)
    ub, code = np.unique(bk, return_inverse=True)
    M = sparse.csr_matrix((vv, (pk, code)), shape=(len(prow), len(ub)))
    ux, uy = ub // 10_000_000, ub % 10_000_000
    order = np.argsort(ux * 10_000_000 + uy)
    lab = [f"{x}_{y}" for x, y in zip(ux[order].tolist(), uy[order].tolist())]
    Ms = M[:, order]

    with open(OUT / "apa_matrix.csv", "w") as fh:
        fh.write("site_id," + ",".join(lab) + "\n")
        Mc = Ms.tocsr()
        for i, pid in enumerate(prow):
            row = np.zeros(len(lab))
            s_, e_ = Mc.indptr[i], Mc.indptr[i + 1]
            row[Mc.indices[s_:e_]] = Mc.data[s_:e_]
            fh.write(pid + "," + ",".join(f"{v:g}" if v else "0" for v in row) + "\n")
    pd.DataFrame({"spot_id": lab, "x": [int(s.split("_")[0]) for s in lab],
                  "y": [int(s.split("_")[1]) for s in lab]}
                 ).to_csv(OUT / "coordinates.csv", index=False)
    peaks.assign(site_id=peaks["peakID"]).to_csv(OUT / "apa_sites.csv", index=False)
    json.dump({"n_peaks": len(prow), "n_bins": len(lab), "nnz": int(Ms.nnz),
               "value_type": "raw_umi_counts", "bin_um": BIN_UM,
               "note": "HD 2um grid aggregated to 8um bins (4x4)"},
              open(OUT / "qc_summary.json", "w"), indent=2)
    print(f"[DONE] {time.time()-t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()

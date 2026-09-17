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

    # 坐标来源: BAM 的 (CB, sb:Z:s_002um_R_C) 标签对 → cb_to_position.tsv。
    # spaceranger 4.1.0 HD 的 tissue_positions.parquet barcode 列是合成 bin ID
    # （s_002um_...），不是序列，无法直接匹配 scAPAtrap 的 spot_id（=CB 序列）。
    cbpos = RAW.parent / "cb_to_position.tsv"
    if not cbpos.exists():
        bam = RAW.parent / "sr" / "outs" / "possorted_genome_bam.bam"
        if not bam.exists():
            raise SystemExit(f"缺 {cbpos} 且找不到 {bam}")
        print(f"[bin8um] cb_to_position.tsv 缺失，从 BAM 提取…")
        import subprocess
        subprocess.run(["/home/mengzijun/anaconda3/envs/spagapa/bin/python",
                        str(Path(__file__).parent / "hd_extract_cb_position.py"),
                        "--bam", str(bam), "--output", str(cbpos)], check=True)
    cpos = pd.read_csv(cbpos, sep="\t", dtype={"cb": str})
    # 2µm 网格 → 目标 bin（4x4 聚合）
    bc2xy = {b: (int(r) // factor, int(c) // factor)
             for b, r, c in zip(cpos["cb"], cpos["sb_row"], cpos["sb_col"])}
    print(f"[bin8um] CB→bin 映射: {len(bc2xy):,} barcodes")

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

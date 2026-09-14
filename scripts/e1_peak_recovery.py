#!/usr/bin/env python3
"""E1 peak-set recovery (memory-safe): re-curate from scAPAtrap long-format reduced counts.

(peak, cell) pairs verified unique -> per-peak n_cells = row count; total = sum(count).
Streams counts.tsv.gz.reduced (~75M rows), curates >=10 cells & >=10 counts,
writes peak x bin200 raw matrix (same layout as binned_200_raw).
"""
from __future__ import annotations
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/gse263789_expand/3m_e1")
RAW = ROOT / "scapatrap_raw"
OUT = ROOT / "binned_200_raw_recovered"
OUT.mkdir(parents=True, exist_ok=True)
BIN, MIN_CELLS, MIN_COUNT = 200, 10, 10


def main():
    t0 = time.time()
    saf = pd.read_csv(RAW / "peaks.saf.reduced", sep="\t",
                      names=["peakID", "chr", "start", "end", "strand"], header=None)
    pidx_all = {p: i for i, p in enumerate(saf["peakID"].tolist())}
    print(f"[1/4] {len(saf)} reduced peaks", flush=True)

    peak_n = {}   # peak -> n_cells (rows)
    peak_tot = {}  # peak -> total count
    agg = {}      # (peak_row_idx, bin_code) -> count
    n_rows = 0
    print("[2/4] streaming ...", flush=True)
    for chunk in pd.read_csv(RAW / "counts.tsv.gz.reduced", sep="\t", chunksize=5_000_000,
                             names=["peak", "cell", "count"], header=0,
                             dtype={"peak": "string", "cell": "string", "count": "float64"}):
        n_rows += len(chunk)
        # per-peak stats
        vc = chunk["peak"].value_counts()
        vs = chunk.groupby("peak")["count"].sum()
        for p, c in vc.items():
            peak_n[p] = peak_n.get(p, 0) + int(c)
        for p, s in vs.items():
            peak_tot[p] = peak_tot.get(p, 0.0) + float(s)
        # bin aggregation
        xy = chunk["cell"].str.split("_", expand=True)
        code = (xy[0].astype(int) // BIN).values * 10_000_000 + (xy[1].astype(int) // BIN).values
        df = pd.DataFrame({"p": chunk["peak"].map(pidx_all).values,
                           "k": code, "c": chunk["count"].values})
        g = df.groupby(["p", "k"], sort=False)["c"].sum()
        for (p, k), v in g.items():
            key = (int(p), int(k))
            agg[key] = agg.get(key, 0.0) + float(v)
        if n_rows % 20_000_000 == 0:
            print(f"      {n_rows:,} rows, {len(agg):,} keys, {time.time()-t0:.0f}s", flush=True)
    print(f"      total {n_rows:,} rows, {len(peak_n)} peaks, {len(agg):,} keys, {time.time()-t0:.0f}s", flush=True)

    keep = [p for p in peak_n if peak_n[p] >= MIN_CELLS and peak_tot.get(p, 0) >= MIN_COUNT]
    keep_idx = {pidx_all[p] for p in keep}
    print(f"[3/4] kept {len(keep)} / {len(peak_n)} peaks (>={MIN_CELLS} cells & >={MIN_COUNT} counts)", flush=True)

    keep_set = set(keep)
    saf_k = saf[saf["peakID"].isin(keep_set)].copy().reset_index(drop=True)
    saf_k["coord"] = np.where(saf_k["strand"] == "+", saf_k["end"], saf_k["start"])
    saf_k["site_id"] = saf_k["peakID"]
    prow = saf_k["peakID"].tolist()
    pidx_k = {p: i for i, p in enumerate(prow)}

    pk = np.fromiter((k[0] for k in agg), dtype=np.int64)
    bk = np.fromiter((k[1] for k in agg), dtype=np.int64)
    vv = np.fromiter(agg.values(), dtype=np.float64)
    # 向量化 remap: 旧行号 -> saf_k 行号（saf_k 保序子集）
    remap = np.full(len(saf), -1, dtype=np.int64)
    kept_old_idx = saf.index[saf["peakID"].isin(keep_set)].values
    remap[kept_old_idx] = np.arange(len(kept_old_idx))
    m = remap[pk] >= 0
    pk, bk, vv = remap[pk[m]], bk[m], vv[m]

    uniq_bins, bin_code = np.unique(bk, return_inverse=True)
    nP, nB = len(prow), len(uniq_bins)
    M = sparse.csr_matrix((vv, (pk, bin_code)), shape=(nP, nB))
    print(f"      matrix {nP} x {nB}, nnz={M.nnz:,}", flush=True)

    print("[4/4] writing ...", flush=True)
    ux = (uniq_bins // 10_000_000).astype(int)
    uy = (uniq_bins % 10_000_000).astype(int)
    order = np.argsort(ux * 1_000_000 + uy)
    lab = [f"{a}_{b}" for a, b in zip(ux[order].tolist(), uy[order].tolist())]
    M_s = M[:, order]
    with open(OUT / "apa_matrix.csv", "w") as fh:
        fh.write("site_id," + ",".join(lab) + "\n")
        M_c = M_s.tocsr()
        indptr, indices, data = M_c.indptr, M_c.indices, M_c.data
        for i, pid in enumerate(prow):
            row = np.zeros(nB)
            s_, e_ = indptr[i], indptr[i + 1]
            row[indices[s_:e_]] = data[s_:e_]
            fh.write(pid + "," + ",".join(f"{v:g}" if v else "0" for v in row) + "\n")
    pd.DataFrame({"spot_id": lab, "x": [int(s.split("_")[0]) for s in lab],
                  "y": [int(s.split("_")[1]) for s in lab]}
                 ).to_csv(OUT / "coordinates.csv", index=False)
    saf_k[["site_id", "peakID", "chr", "start", "end", "strand", "coord"]].to_csv(
        OUT / "apa_sites.csv", index=False)
    json.dump({"n_peaks": int(nP), "n_bins": int(nB), "nnz": int(M_s.nnz),
               "value_type": "raw_umi_counts", "bin": BIN,
               "curation": f">={MIN_CELLS} cells & >={MIN_COUNT} counts (uniform re-curation)",
               "source": "counts.tsv.gz.reduced (scAPAtrap final cut 133834->825 over-aggressive)"},
              open(OUT / "qc_summary.json", "w"), indent=2)
    print(f"[DONE] {time.time()-t0:.0f}s -> {OUT}  ({nP} peaks)", flush=True)


if __name__ == "__main__":
    main()

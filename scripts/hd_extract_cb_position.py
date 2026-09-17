#!/usr/bin/env python3
"""HD BAM → CB→2µm-bin 映射表（cb_to_position.tsv）。

spaceranger 4.1.0 HD BAM 的每条 read 携带:
  CB:Z:<31nt 校正条码>          ← scAPAtrap 的 spot_id 即此值
  sb:Z:s_002um_<row>_<col>-1    ← 空间 bin ID（2µm 网格坐标）
tissue_positions.parquet 的 barcode 列是合成 bin ID 而非序列，
因此 CB→坐标必须从 BAM 的 (CB, sb) 标签对提取。

用法: hd_extract_cb_position.py --bam <possorted_genome_bam.bam> --output <cb_to_position.tsv>
"""
from __future__ import annotations
import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

SAMTOOLS = "/home/mengzijun/anaconda3/envs/samtools/bin/samtools"
PAT = re.compile(rb"CB:Z:([ACGTNacgtn]+)\t[^\n]*?sb:Z:s_002um_(\d+)_(\d+)-1")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--threads", type=int, default=12)
    a = ap.parse_args()
    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    proc = subprocess.Popen(
        [SAMTOOLS, "view", "-@" + str(a.threads), a.bam],
        stdout=subprocess.PIPE, bufsize=0)
    seen: dict[bytes, tuple[int, int]] = {}
    n_reads = 0
    tail = b""
    assert proc.stdout is not None
    while True:
        chunk = proc.stdout.read(1 << 26)  # 64 MiB
        if not chunk:
            break
        data = tail + chunk
        nl = data.rfind(b"\n")
        if nl < 0:
            tail = data
            continue
        tail = data[nl + 1:]
        body = data[:nl + 1]
        n_reads += body.count(b"\n")
        for cb, r, c in PAT.findall(body):
            if cb not in seen:
                seen[cb] = (int(r), int(c))
    if tail:
        for cb, r, c in PAT.findall(tail):
            if cb not in seen:
                seen[cb] = (int(r), int(c))
    rc = proc.wait()
    if rc != 0:
        sys.exit(f"samtools view failed rc={rc}")

    with open(out, "w") as fh:
        fh.write("cb\tsb_row\tsb_col\n")
        for cb, (r, c) in seen.items():
            fh.write(f"{cb.decode()}\t{r}\t{c}\n")
    print(f"[cb_pos] reads={n_reads:,} unique_cb={len(seen):,} "
          f"wall={time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    main()

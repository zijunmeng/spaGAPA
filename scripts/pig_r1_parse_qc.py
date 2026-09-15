#!/usr/bin/env python3
"""无菌猪 KC 数据 R1 解析 + QC（无需 manifest 即可运行）。

结构（GeneMind 官方 + 实测验证）:
  R1[0:30]   空间条码
  R1[30:51]  固定序列 GTTCGCAACATGTCTGGCGTCATAGAATTC (实测锚: CAACATGTCTGGCGTCATAGA)
  R1[51:61]  UMI 10bp
  R1[61:~91] oligo-dT 残留（浮动）
  R1[91:150] mRNA 3' 末端（PAS 侧）

输出 QC: 条码质量/复杂度、UMI 复杂度、polyT 长度分布、捕获结构完整率、
        PAS 候选读段（T-连接点 + 下游 60bp）抽样导出。
"""
from __future__ import annotations
import gzip
import sys
import collections
import json
import re

FIXED_FULL = "GTTCGCAACATGTCTGGCGTCATAGAATTC"  # xlsx 官方
FIXED_ANCHOR = "CAACATGTCTGGCGTCATAGA"          # 实测稳定子串


def parse(r1_path: str, n_max: int = 2_000_000, sample_out: str | None = None):
    stats = collections.Counter()
    bc_count = collections.Counter()
    umi_count = collections.Counter()
    polyT_len = collections.Counter()
    fixed_offset = collections.Counter()
    pas_samples = []
    rng = 93
    out = open(sample_out, "w") if sample_out else None
    with gzip.open(r1_path, "rt") as f:
        n = 0
        h = s = q = None
        for line in f:
            n += 1
            if n % 4 == 2:
                s = line.strip()
            elif n % 4 == 0:
                q = line.strip()
                if n // 4 > n_max:
                    break
                stats["total"] += 1
                p = s.find(FIXED_ANCHOR)
                if p < 0:
                    stats["no_fixed"] += 1
                    continue
                fixed_offset[p - 4] += 1  # 完整 21mer 起点 ≈ p-4
                bc = s[:max(p - 4, 0)]
                if len(bc) < 20:
                    stats["short_bc"] += 1
                    continue
                stats["with_bc"] += 1
                bc_count[bc] += 1
                # UMI = 固定序列(21) 之后 10bp
                umi_start = (p - 4) + len(FIXED_FULL)
                umi = s[umi_start:umi_start + 10]
                if len(umi) == 10:
                    stats["with_umi"] += 1
                    umi_count[umi] += 1
                # polyT 残留 + PAS 连接
                rest = s[umi_start + 10:]
                m = re.match(r"(T{6,})", rest)
                if m:
                    stats["with_polyT"] += 1
                    polyT_len[min(len(m.group(1)), 40)] += 1
                    transcript = rest[m.end():]
                    if len(transcript) >= 25:
                        stats["pas_candidates"] += 1
                        if len(pas_samples) < 50000 and out:
                            out.write(f"{bc}\t{umi}\t{len(m.group(1))}\t{transcript}\n")
                else:
                    stats["no_polyT"] += 1
    top_off = fixed_offset.most_common(1)
    return dict(
        total=stats["total"],
        with_fixed=stats["total"] - stats["no_fixed"],
        with_bc=stats["with_bc"],
        with_umi=stats["with_umi"],
        with_polyT=stats["with_polyT"],
        pas_candidates=stats["pas_candidates"],
        n_unique_bc=len(bc_count),
        bc_top10_share=sum(c for _, c in bc_count.most_common(10)) / max(1, stats["with_bc"]),
        n_unique_umi=len(umi_count),
        fixed_offset_mode=top_off[0][0] if top_off else None,
        polyT_median_len=sorted(polyT_len.elements())[len(polyT_len) // 2] if polyT_len else None,
    )


if __name__ == "__main__":
    fq, out_json = sys.argv[1], sys.argv[2]
    sample = sys.argv[3] if len(sys.argv) > 3 else None
    res = parse(fq, sample_out=sample)
    json.dump(res, open(out_json, "w"), indent=2)
    print(json.dumps(res, indent=2))

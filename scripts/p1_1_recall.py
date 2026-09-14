#!/usr/bin/env python3
"""P1-1: E2a 补召回率 — PolyA_DB 组织可检测位点中 spaGAPA/scAPAtrap 找回多少。

公平分母定义: PolyA_DB 位点限制在【样本有表达证据的基因】内
（用样本的 apa_sites 基因注释 ∩ PolyA_DB 位点基因注释），±50bp 判命中。

对每个样本:
  recall = |DB sites (in detected genes) matched by our peaks| / |DB sites in detected genes|
  precision = 已有 hit_rate_50bp (74.7%)
  F1 = 2PR/(P+R)

输出: pipeline_output/polya_db_overlap/recall_by_sample.csv + 更新 summary
"""
from __future__ import annotations
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np
import pandas as pd

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT = ROOT / "pipeline_output/polya_db_overlap"

DB = {
    "hg38": ROOT / "data/external/polya_db_v4/human/HumanPas/hg38.PAS.max.tsv",
    "mm10": ROOT / "data/external/polya_db_v4/mouse/MousePas/mm10.PAS.max.tsv",
}

SAMPLES = {}
_ov = pd.read_csv(OUT / "overlap_by_sample.csv")
for _, r in _ov.iterrows():
    SAMPLES[r["sample"]] = dict(build=r["genome_build"], species=r["species"])
def load_db(build):
    df = pd.read_csv(DB[build], sep="\t")
    # PAS_ID = "chr:strand:pos"; GeneSymbol column
    parts = df["PAS_ID"].str.split(":", expand=True)
    df["_chr"] = parts[0].astype(str).str.replace("^chr", "", regex=True)
    df["_strand"] = parts[1]
    df["_pos"] = parts[2].astype(int)
    df["_sym"] = df["GeneSymbol"].astype(str)
    return df


def main():
    rows = []
    for sample, meta in SAMPLES.items():
        build = meta["build"]
        # find peaks_meta for the sample — reuse manifest from the original script
        import re
        src = (ROOT / "scripts/polya_db_overlap_validation.py").read_text()
        pat = re.compile(
            r'\(\s*"' + re.escape(sample) + r'"\s*,\s*"[^"]+"\s*,\s*"[^"]+"\s*,\s*"[^"]+"'
            r'\s*,\s*"[^"]+"\s*,\s*"' + re.escape(build) + r'"\s*,\s*f?"([^"]+)"')
        m = pat.search(src)
        if not m:
            continue
        peaks_path = ROOT / m.group(1).replace("{_S}", "raw_scapatrap/peaks_meta.csv.gz")
        if not peaks_path.exists():
            continue
        if build not in ("hg38", "mm10"):
            # GRCm39 样本需 liftOver（E2a 已对 precision 做）；recall 记为 deferred
            continue
        peaks = pd.read_csv(peaks_path)
        peaks["chr"] = peaks["chr"].astype(str).str.replace("^chr", "", regex=True)

        db = load_db(build)

        # 公平分母: DB 位点在样本检测到的基因内（搜索样本树的 apa_sites）
        detected_syms = set()
        base = peaks_path.parent.parent
        cands = list(base.rglob("apa_sites.csv"))[:1]
        if not cands and base.parent != ROOT:
            cands = list(base.parent.rglob("apa_sites.csv"))[:1]
        for p in cands:
            try:
                s = pd.read_csv(p)
                if "gene_name" in s.columns and s["gene_name"].notna().any():
                    detected_syms = set(s["gene_name"].dropna().astype(str)) - {""}
            except Exception:
                pass
            break
        if not detected_syms:
            # GTF 兜底: 直接注释峰
            import sys as _sys
            _sys.path.insert(0, str(ROOT / "scripts"))
            from gse263789_ad_vs_wt_analysis import read_genes_from_gtf, assign_genes
            gtf = "/s1/SHARE/01_software/SAW_refs/Homo_sapiens_index/genes/Homo_sapiens.GRCh38.93.saw.gtf" \
                if build == "hg38" else "/s1/SHARE/01_software/SAW_refs/Mus_musculus_index/genes/Mus_musculus.GRCm38.93.saw.gtf"
            try:
                genes = read_genes_from_gtf(gtf)
                ann = assign_genes(peaks.copy(), genes)
                detected_syms = set(ann["gene_name"].dropna().astype(str)) - {""}
            except Exception as e:
                print(f"  [warn] GTF annotation failed for {sample}: {e}")
        if detected_syms:
            db_mask = db["_sym"].isin(detected_syms)
        else:
            db_mask = db["_chr"].isin(set(peaks["chr"]))
        db_sub = db[db_mask]
        denom_total = len(db_sub)

        # 匹配: 对每个染色体, our peak coords vs db_sub positions, ±50bp
        matched = 0
        for chrom, sub in db_sub.groupby("_chr"):
            pk = peaks.loc[peaks["chr"] == chrom, "coord"].values.astype(int)
            if len(pk) == 0:
                continue
            pk.sort()
            pos = sub["_pos"].values
            idx = np.searchsorted(pk, pos)
            hit = np.zeros(len(pos), dtype=bool)
            for off in (-1, 0):
                ii = idx + off
                ok = (ii >= 0) & (ii < len(pk))
                hit[ok] |= np.abs(pos[ok] - pk[ii[ok]]) <= 50
            matched += int(hit.sum())

        recall = matched / denom_total if denom_total else np.nan
        prec_row = _ov[_ov["sample"] == sample].iloc[0]
        prec = prec_row["hit_rate_50bp"]
        f1 = 2 * prec * recall / (prec + recall) if (prec and recall and prec + recall) else np.nan
        rows.append(dict(sample=sample, build=build, n_db_denominator=denom_total,
                         n_db_matched=matched, recall_50bp=recall,
                         precision_50bp=prec, f1=f1,
                         denominator_basis="db sites in detected genes" if detected_syms else "db sites on detected chroms"))
        print(f"{sample}: denom={denom_total:,} matched={matched:,} recall={recall:.4f} prec={prec:.4f} F1={f1:.4f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "recall_by_sample.csv", index=False)
    summary = dict(
        n_samples=len(df),
        mean_recall=float(df["recall_50bp"].mean()),
        median_recall=float(df["recall_50bp"].median()),
        mean_f1=float(df["f1"].mean()),
        note="recall = fraction of PolyA_DB sites within sample-detected genes recovered within ±50bp; precision from E2a hit_rate_50bp",
    )
    json.dump(summary, open(OUT / "recall_summary.json", "w"), indent=2)
    print("\nmean recall:", summary["mean_recall"], " mean F1:", summary["mean_f1"])


if __name__ == "__main__":
    main()

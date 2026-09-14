#!/usr/bin/env python3
"""Phase 3 upgrade: AD 2v1 ({E3,E4} vs F5 WT) + 3M-vs-5M progression + Apoe genotyping.

Extends gse263789_ad_vs_wt_analysis.py (1v1) to:
  (A) AD 2v1: pool {E3 pilot + E4} binned spots vs F5 WT -> differential APA
      with doubled AD power; consistency check between E3 and E4 direction.
  (B) 3M-vs-5M progression (optional, auto-skips if 3M binned not present):
      {E1,E2} (3-month AD) vs {E3,E4} (5-month AD) -> early-vs-late APA shift.
  (C) Apoe/Trem2/Cst7 expression per condition from SAW gene matrices
      (microglia/plaque-response genotyping proxy; no plaque PNG needed).

Inputs
------
E3   pipeline_output/gse263789_stereo_pilot/spagapa_downstream_full/binned_200_raw/
E4   pipeline_output/gse263789_expand/ad18_e4/binned_200/
F5   pipeline_output/gse263789_wt_control/binned_200/
E1   pipeline_output/gse263789_expand/3m_e1/binned_200/      (optional)
E2   pipeline_output/gse263789_expand/3m_e2/binned_200/      (optional)

Output: pipeline_output/gse263789_ad_vs_wt_differential/phase3/
"""
from __future__ import annotations

import os, sys, json, gzip, warnings
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output")
OUT = ROOT / "gse263789_ad_vs_wt_differential" / "phase3"
OUT.mkdir(parents=True, exist_ok=True)

GTF = Path("/s1/SHARE/01_software/SAW_refs/Mus_musculus_index/genes/Mus_musculus.GRCm38.93.saw.gtf")

# reuse helpers from the 1v1 script
sys.path.insert(0, str(Path(__file__).parent))
from gse263789_ad_vs_wt_analysis import (
    read_genes_from_gtf, assign_genes, build_gene_index,
)

CONDITIONS = {
    "E3_5m_ad": ROOT / "gse263789_stereo_pilot/spagapa_downstream_full/binned_200_raw",
    "E4_5m_ad": ROOT / "gse263789_expand/ad18_e4/binned_200_raw_recovered",
    "F5_wt":    ROOT / "gse263789_wt_control/binned_200",
    "E1_3m_ad": ROOT / "gse263789_expand/3m_e1/binned_200",
    "E2_3m_ad": ROOT / "gse263789_expand/3m_e2/binned_200",
}

# plaque-response genotype markers (microglia DAM program)
GENOTYPE_GENES = ["Apoe", "Trem2", "Tyrobp", "Cst7", "Clec7a", "Lpl", "Csf1", "Axl"]


def load_condition(name: str, path: Path):
    """Load binned matrix (site x spot raw counts) + sites; annotate genes; gene-level index."""
    if not (path / "apa_matrix.csv").exists():
        print(f"[skip] {name}: {path}/apa_matrix.csv missing")
        return None
    counts = pd.read_csv(path / "apa_matrix.csv", index_col=0)
    sites = pd.read_csv(path / "apa_sites.csv")
    # normalize chr naming (scapatrap uses 'chr1', SAW GTF uses '1')
    if sites["chr"].astype(str).str.startswith("chr").any():
        sites["chr"] = sites["chr"].astype(str).str.replace("^chr", "", regex=True)
    if "gene_name" not in sites.columns or sites["gene_name"].isna().all():
        genes = read_genes_from_gtf(GTF)
        sites = assign_genes(sites, genes)
    gene_idx = build_gene_index(counts, sites, min_parent=10)  # gene x spot
    coords = pd.read_csv(path / "coordinates.csv", index_col=0) if (path / "coordinates.csv").exists() else None
    n_spots = counts.shape[1]
    print(f"[load] {name}: {counts.shape[0]} sites x {n_spots} spots -> {gene_idx.shape[0]} genes")
    return dict(name=name, gene_idx=gene_idx, coords=coords, n_spots=n_spots)


def pooled_diff(ad_conds, wt_conds, label, min_spots=10):
    """Pooled differential APA: AD group vs WT group on the gene-level distal index."""
    ad = pd.concat([c["gene_idx"] for c in ad_conds], axis=1)
    wt = pd.concat([c["gene_idx"] for c in wt_conds], axis=1)
    common = ad.index.intersection(wt.index)
    ad, wt = ad.loc[common], wt.loc[common]
    from scipy import stats
    rows = []
    for g in common:
        a, w = ad.loc[g].values.flatten(), wt.loc[g].values.flatten()
        a, w = a[~np.isnan(a)], w[~np.isnan(w)]
        if len(a) < min_spots or len(w) < min_spots:
            continue
        t, p = stats.mannwhitneyu(a, w, alternative="two-sided")
        rows.append(dict(gene=g, mean_ad=float(np.mean(a)), mean_wt=float(np.mean(w)),
                         delta=float(np.mean(a) - np.mean(w)),
                         p_value=float(p), n_ad=len(a), n_wt=len(w)))
    df = pd.DataFrame(rows)
    if len(df):
        # portable BH (skip NaN p-values; clip to [0,1] defensively)
        pv = np.asarray(df["p_value"].values, dtype=float)
        pv = np.clip(np.nan_to_num(pv, nan=1.0), 0.0, 1.0)
        m = len(df)
        order = np.argsort(pv)
        ranked = pv[order] * m / (np.arange(m) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adj = np.empty(m)
        adj[order] = np.clip(ranked, 0, 1)
        df["p_adj_BH"] = adj
    df.to_csv(OUT / f"{label}_differential_apa.csv", index=False)
    if len(df):
        sig = df[df["p_adj_BH"] < 0.05]
        print(f"[diff] {label}: {len(df)} genes tested, {len(sig)} FDR<0.05 "
              f"({(sig['delta']>0).sum()} distal-up / {(sig['delta']<0).sum()} proximal-up)")
    return df


def consistency_check(df_a, df_b, key="gene", val="delta"):
    """Direction consistency between two differential results (e.g. E3-only vs E4-only)."""
    m = df_a.merge(df_b, on=key, suffixes=("_a", "_b"))
    if not len(m):
        return None
    both_sig = m[(m["p_adj_BH_a"] < 0.05) & (m["p_adj_BH_b"] < 0.05)]
    consist = float((np.sign(both_sig[val + "_a"]) == np.sign(both_sig[val + "_b"])).mean()) \
        if len(both_sig) else None
    from scipy import stats as _st
    _ok = m[[val + "_a", val + "_b"]].dropna()
    rho = float(_st.spearmanr(_ok[val + "_a"], _ok[val + "_b"]).statistic) if len(_ok) > 2 else None
    return dict(n_overlap=int(len(m)), n_both_sig=int(len(both_sig)),
                direction_consistency=consist, delta_spearman=rho)


def genotype_markers(conds):
    """Per-condition expression of plaque-response genes from SAW matrices (if present)."""
    rows = []
    for c in conds:
        saw = None
        for cand in c.get("saw_paths", []):
            if cand.exists():
                saw = cand; break
        if saw is None:
            continue
        # SAW raw_barcode_gene_exp: long format x,y,geneID,MID,reads
        try:
            # SAW 格式: x y geneID(Ensembl) MIDIndex readCount (带 header)
            if "_ens2sym" not in genotype_markers.__dict__:
                from gse263789_ad_vs_wt_analysis import read_genes_from_gtf as _rg
                gmap = {}
                gtf = "/s1/SHARE/01_software/SAW_refs/Mus_musculus_index/genes/Mus_musculus.GRCm38.93.saw.gtf"
                for (chrom, strand), recs in _rg(gtf).items():
                    for r in recs:
                        gmap[r.gene_id] = r.gene_name
                genotype_markers._ens2sym = gmap
            gmap = genotype_markers._ens2sym
            chunks = pd.read_csv(saw, sep="\t", header=0,
                                 names=["x", "y", "gene", "umi", "reads"],
                                 chunksize=5_000_000, usecols=["gene", "reads"])
            tot = {}
            for ch in chunks:
                ch["gene"] = ch["gene"].map(gmap)
                s = ch.groupby("gene")["reads"].sum()
                for gn, v in s.items():
                    tot[gn] = tot.get(gn, 0) + int(v)
        except Exception as e:
            print(f"[genotype] {c['name']}: read fail {e}")
            continue
        # Ensembl->symbol naive: match gene names via GTF attrs if needed
        rows.append(dict(condition=c["name"],
                         **{g: tot.get(g, 0) for g in GENOTYPE_GENES}))
    if rows:
        pd.DataFrame(rows).to_csv(OUT / "genotype_markers.csv", index=False)
    return rows


def main():
    print("=== Phase 3: AD 2v1 + progression ===")
    genes_gtf = None

    conds = {}
    for name, path in CONDITIONS.items():
        c = load_condition(name, path)
        if c:
            conds[name] = c

    if "E3_5m_ad" in conds and "E4_5m_ad" in conds and "F5_wt" in conds:
        print("\n--- (A) AD 2v1: {E3,E4} vs F5 ---")
        df_2v1 = pooled_diff([conds["E3_5m_ad"], conds["E4_5m_ad"]], [conds["F5_wt"]], "ad2v1_vs_wt")

        print("\n--- consistency: E3-alone vs E4-alone ---")
        df_e3 = pooled_diff([conds["E3_5m_ad"]], [conds["F5_wt"]], "e3_vs_wt_check")
        df_e4 = pooled_diff([conds["E4_5m_ad"]], [conds["F5_wt"]], "e4_vs_wt_check")
        cc = consistency_check(df_e3, df_e4)
        if cc:
            print("consistency:", cc)
            json.dump(cc, open(OUT / "e3_e4_consistency.json", "w"), indent=2)
    else:
        print("!! need E3+E4+F5 for 2v1")

    if "E1_3m_ad" in conds and "E2_3m_ad" in conds and "E4_5m_ad" in conds:
        print("\n--- (B) progression: {E1,E2} 3m vs {E3,E4} 5m ---")
        pooled_diff([conds["E1_3m_ad"], conds["E2_3m_ad"]],
                    [conds["E3_5m_ad"], conds["E4_5m_ad"]], "3m_vs_5m_progression")
    else:
        print("[skip] 3M conditions not ready (E1/E2 binned missing)")

    print("\n--- (C) genotype markers ---")
    # attach SAW expression paths if we can find them
    saw_map = {
        "E3_5m_ad": ROOT / "gse263789_stereo_pilot",
        "E4_5m_ad": ROOT / "gse263789_expand/ad18_e4",
        "F5_wt": ROOT / "gse263789_wt_control",
        "E1_3m_ad": ROOT / "gse263789_expand/3m_e1",
        "E2_3m_ad": ROOT / "gse263789_expand/3m_e2",
    }
    for name, base in saw_map.items():
        if name in conds:
            hits = list(base.rglob("*raw_barcode_gene_exp*.txt"))[:1]
            conds[name]["saw_paths"] = hits
    genotype_markers(list(conds.values()))

    print(f"\ndone -> {OUT}")


if __name__ == "__main__":
    main()

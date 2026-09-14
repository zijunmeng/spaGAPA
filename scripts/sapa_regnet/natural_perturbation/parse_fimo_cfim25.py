#!/usr/bin/env python3
"""Parse NUDT21 (CFIm25) + neural-RBP FIMO outputs -> APA-region gene sets.

IMPORTANT (bug fix vs parse_fimo_rbp_genomewide.py): FIMO 5.5.9 auto-parses the
UCSC-style region suffix of the fasta headers (chr:start-end(strand)) and
reports start/stop ALREADY in genomic coordinates. The older parser added
utr3 offsets on top of them, corrupting coordinates (e.g. VAMP2 sites at
"1449") and forcing in_apa_region=True for every site. Here coordinates are
used directly, and APA-region membership is computed by genomic interval
overlap with Track-A's apa_region (proximal PAS -> gene TES) per gene.

Outputs (natural_perturbation/):
  cfim25_fimo/nudt21_sites.tsv             per-site, correct genomic coords
  cfim25_fimo/nudt21_apa_site_counts.csv   per-gene APA-region UUGUA counts
  _cache/rbp_apa_genes.csv                 per-RBP gene sets (15 neural + NUDT21)
"""
from pathlib import Path

import pandas as pd

PROJ = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT = PROJ / "pipeline_output/sapa_regnet/natural_perturbation"
CACHE = OUT / "_cache"
CACHE.mkdir(parents=True, exist_ok=True)
NEURAL_FIMO = PROJ / "pipeline_output/sapa_regnet/rbp_genomewide/fimo_neural_rbp/fimo.tsv"
CFIM25_FIMO = OUT / "cfim25_fimo/fimo.tsv"

ann = pd.read_csv(PROJ / "pipeline_output/sapa_regnet/apa_regulatory_annotation_rbp.csv",
                  low_memory=False,
                  usecols=["gene_name", "chr", "apa_region_start", "apa_region_end"])
ann = ann.dropna(subset=["apa_region_start"]).drop_duplicates(["gene_name", "chr"])


def load_fimo(path):
    f = pd.read_csv(path, sep="\t", comment="#",
                    usecols=["motif_id", "sequence_name", "start", "stop",
                             "strand", "score", "p-value"])
    f = f.rename(columns={"start": "gstart", "stop": "gend", "p-value": "pvalue"})
    parts = f.sequence_name.astype(str).str.extract(
        r"^(?P<gene_name>[^|]+)\|[^|]+\|(?P<chr>chr[^|]+)$")
    f = pd.concat([f, parts], axis=1).dropna(subset=["gene_name"])
    f["rbp"] = f.motif_id.str.split("_").str[0]
    return f


def apa_overlap(f):
    """Merge with Track-A apa_region intervals; flag genomic overlap."""
    m = f.merge(ann, on=["gene_name", "chr"], how="inner")
    m["in_apa_region"] = ((m.gstart <= m.apa_region_end) &
                          (m.gend >= m.apa_region_start))
    return m


def main():
    # ---- neural RBPs ----
    f = load_fimo(NEURAL_FIMO)
    print(f"neural FIMO hits: {len(f):,}")
    m = apa_overlap(f)
    print(f"  on annotated genes: {len(m):,}; in APA region: "
          f"{m.in_apa_region.sum():,} ({m.in_apa_region.mean():.1%})")
    f[["gene_name", "rbp", "chr", "gstart", "gend", "strand", "score", "pvalue"]].to_csv(
        CACHE / "neural_sites_genomic.tsv", sep="\t", index=False)
    apa = m[m.in_apa_region]
    sets = apa.groupby("rbp")["gene_name"].apply(set).to_dict()
    print("  per-RBP APA-region gene counts:",
          {k: len(v) for k, v in sorted(sets.items())})

    # ---- NUDT21 (CFIm25) ----
    fn = load_fimo(CFIM25_FIMO)
    # dedupe the two identical UUGUA motifs (same gene + same genomic window)
    fn = fn.drop_duplicates(subset=["gene_name", "chr", "gstart", "gend", "rbp"])
    mn = apa_overlap(fn)
    mn["rbp"] = "NUDT21"
    mn[["gene_name", "rbp", "chr", "gstart", "gend", "strand", "score",
        "pvalue", "in_apa_region"]].to_csv(OUT / "cfim25_fimo/nudt21_sites.tsv",
                                           sep="\t", index=False)
    print(f"NUDT21 hits: {len(fn):,}; on annotated genes: {len(mn):,}; "
          f"in APA region: {mn.in_apa_region.sum():,} ({mn.in_apa_region.mean():.1%})")
    cnt = (mn[mn.in_apa_region].groupby(["gene_name", "chr"]).size()
           .reset_index(name="n_nudt21_apa_sites"))
    cnt.to_csv(OUT / "cfim25_fimo/nudt21_apa_site_counts.csv", index=False)
    for k in (1, 2, 3):
        print(f"  genes with >= {k} APA-region UUGUA site(s): "
              f"{(cnt.n_nudt21_apa_sites >= k).sum():,}")
    sets["NUDT21"] = set(cnt.loc[cnt.n_nudt21_apa_sites >= 2, "gene_name"])

    rows = [(rbp, g) for rbp, genes in sets.items() for g in genes]
    pd.DataFrame(rows, columns=["rbp", "gene"]).to_csv(CACHE / "rbp_apa_genes.csv",
                                                       index=False)
    print(f"wrote {CACHE/'rbp_apa_genes.csv'} ({len(rows):,} pairs)")


if __name__ == "__main__":
    main()

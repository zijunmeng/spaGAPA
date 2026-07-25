#!/usr/bin/env python3
"""Parse the genome-wide neural-RBP FIMO output into APA-region annotation.

Fixes 2 bugs vs the original rescan parser:
  (1) this FIMO build names the end column `stop` (not `end`);
  (2) FIMO truncates sequence headers -> `GENE|gene_id|chr` (no :start-end(strand)),
      so genomic/strand context is recovered by joining utr3_per_gene (longest
      3'UTR per gene as representative) instead of from the header.

Output: apa_regulatory_annotation_rbp_neural.csv with RBP sites split into
APA-sensitive region (3' of proximal PAS) vs common, genome-wide, all chroms.
"""
import os, sys
import pandas as pd

PROJ = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
FIMO_TSV = f"{PROJ}/pipeline_output/sapa_regnet/rbp_genomewide/fimo_neural_rbp/fimo.tsv"
UTR3 = f"{PROJ}/pipeline_output/sapa_regnet/track_b_mirna/utr3_per_gene.tsv.gz"
# proximal PAS coord per gene (precomputed by Track A, 16079 genes):
RBP_TRACKA = f"{PROJ}/pipeline_output/sapa_regnet/apa_regulatory_annotation_rbp.csv"
APA_SITES = f"{PROJ}/data/processed/gse220442_gsm6801751_scapatrap/apa_sites.csv.gz"
OUT = f"{PROJ}/pipeline_output/sapa_regnet/apa_regulatory_annotation_rbp_neural.csv"
SITES_OUT = f"{PROJ}/pipeline_output/sapa_regnet/rbp_genomewide/neural_rbp_sites.tsv"

print("loading FIMO ...", flush=True)
fimo = pd.read_csv(FIMO_TSV, sep="\t", comment="#",
                   usecols=["motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value"])
fimo = fimo.rename(columns={"start": "fstart", "stop": "fstop", "p-value": "pvalue"})
print(f"  {len(fimo):,} FIMO hits", flush=True)

# parse sequence_name  GENE|gene_id|chr  (FIMO truncated the :start-end(strand))
parts = fimo["sequence_name"].astype(str).str.extract(
    r"^(?P<gene_name>[^|]+)\|(?P<gene_id>[^|]+)\|(?P<chr>chr[^|]+)$")
bad = parts["gene_name"].isna().sum()
print(f"  un-parseable sequence_name: {bad:,}", flush=True)
fimo = pd.concat([fimo, parts], axis=1).dropna(subset=["gene_name"])
fimo["rbp"] = fimo["motif_id"].astype(str).str.split("_").str[0]

# representative 3'UTR interval per gene = longest
print("loading utr3 intervals (longest per gene) ...", flush=True)
utr = pd.read_csv(UTR3, sep="\t")
utr["len"] = utr["utr3_end"] - utr["utr3_start"]
rep = utr.sort_values("len", ascending=False).drop_duplicates("gene_id")[
    ["gene_id", "strand", "utr3_start", "utr3_end"]].rename(columns={"strand": "gstrand"})

# proximal PAS coord per gene (strand-oriented, nearest TSS)
print("computing proximal PAS per gene from apa_sites ...", flush=True)
sites = pd.read_csv(APA_SITES)
def prox(g):
    g = g.dropna(subset=["coord"])
    if len(g) < 2:
        return None
    s = str(g["strand"].iloc[0])
    if s == "-":
        return float(g["coord"].max())  # nearest TSS = max coord on - strand
    return float(g["coord"].min())
pas = sites.groupby(["gene_name", "chr"]).apply(prox, include_groups=False).reset_index()
pas.columns = ["gene_name", "chr", "proximal_pas_coord"]
pas = pas.dropna(subset=["proximal_pas_coord"])

# merge: fimo -> representative utr3 (by gene_id) -> proximal PAS (by gene_name+chr)
fimo = fimo.merge(rep, on="gene_id", how="left")
fimo = fimo.merge(pas, on=["gene_name", "chr"], how="left")
fimo = fimo.dropna(subset=["utr3_start", "proximal_pas_coord"]).copy()
fimo["proximal_pas_coord"] = fimo["proximal_pas_coord"].astype(float)

# proximal PAS position within the SENSE 3'UTR sequence (1-based from 3'UTR 5' end)
plus = fimo["gstrand"] == "+"
prox_sense = pd.Series(0, index=fimo.index, dtype=float)
prox_sense[plus] = fimo.loc[plus, "proximal_pas_coord"] - fimo.loc[plus, "utr3_start"] + 1
prox_sense[~plus] = fimo.loc[~plus, "utr3_end"] - fimo.loc[~plus, "proximal_pas_coord"] + 1
fimo["prox_sense"] = prox_sense
# APA-sensitive region = sense position 3' of proximal PAS (fstart > prox_sense)
fimo["in_apa_region"] = fimo["fstart"] > fimo["prox_sense"]

print(f"  sites in APA region: {fimo['in_apa_region'].sum():,} / {len(fimo):,}", flush=True)

# per-site output (genomic coords via strand-aware conversion of fstart/fstop)
fp, fm = fimo["gstrand"] == "+", fimo["gstrand"] != "+"
gstart = pd.Series(0, index=fimo.index, dtype=int); gend = pd.Series(0, index=fimo.index, dtype=int)
gstart[fp] = (fimo.loc[fp, "utr3_start"] + fimo.loc[fp, "fstart"] - 1).astype(int)
gend[fp] = (fimo.loc[fp, "utr3_start"] + fimo.loc[fp, "fstop"] - 1).astype(int)
gstart[fm] = (fimo.loc[fm, "utr3_end"] - fimo.loc[fm, "fstop"] + 1).astype(int)
gend[fm] = (fimo.loc[fm, "utr3_end"] - fimo.loc[fm, "fstart"] + 1).astype(int)
fimo["genomic_start"] = gstart; fimo["genomic_end"] = gend
fimo[["gene_name", "rbp", "chr", "genomic_start", "genomic_end", "strand",
      "score", "pvalue", "in_apa_region"]].to_csv(SITES_OUT, sep="\t", index=False)
print(f"  wrote per-site {SITES_OUT}", flush=True)

# per-gene summary
apa_rbp = fimo[fimo["in_apa_region"]].groupby(["gene_name", "chr"])["rbp"].apply(
    lambda s: ",".join(sorted(set(s)))).reset_index().rename(columns={"rbp": "rbp_in_apa_region"})
g2 = fimo.assign(in_common=~fimo["in_apa_region"])
common_rbp = g2[g2["in_common"]].groupby(["gene_name", "chr"])["rbp"].apply(
    lambda s: ",".join(sorted(set(s)))).reset_index().rename(columns={"rbp": "rbp_in_common"})
allg = fimo[["gene_name", "chr"]].drop_duplicates()
agg = allg.merge(apa_rbp, on=["gene_name", "chr"], how="left").merge(
    common_rbp, on=["gene_name", "chr"], how="left").fillna("")
agg["n_rbp_apa"] = agg["rbp_in_apa_region"].apply(lambda s: len([x for x in str(s).split(",") if x]))
agg["n_rbp_common"] = agg["rbp_in_common"].apply(lambda s: len([x for x in str(s).split(",") if x]))
agg.to_csv(OUT, index=False)
print(f"\nDONE: {len(agg)} genes -> {OUT}", flush=True)
print(agg.head(8).to_string(index=False), flush=True)
# AD gene check
for g in ["APP", "SV2B", "YWHAZ", "ARPP19", "DCLK1", "TARDBP"]:
    r = agg[agg["gene_name"] == g]
    print(f"  {g}: APA-region RBPs = {str(r.iloc[0]['rbp_in_apa_region'])[:90]}" if len(r) else f"  {g}: not found", flush=True)

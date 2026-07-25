#!/usr/bin/env python3
"""Build a per-gene collapsed 3'UTR interval table from gencode v44 GTF.

For each protein_coding gene, take all UTR intervals that lie 3' of the
stop codon (on the transcript strand) across all transcripts, merge overlapping
intervals, and write one row per gene: gene_id, gene_name, chr, strand,
utr3_start, utr3_end (genome forward coords, start<=end).

Strategy: a UTR interval is a 3'UTR piece if, on the gene strand, it is
downstream of the gene's stop codons. Concretely:
  + strand: UTR.start > min(stop_codon.start)  -> 3' side
  - strand: UTR.end   < max(stop_codon.end)    -> 3' side
Union over transcripts, then merge.
"""
import gzip
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

GTF = Path("/s1/SHARE/00_ref_genecode/gencode.v44.annotation.gtf")
OUT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA/pipeline_output/sapa_regnet/track_b_mirna/utr3_per_gene.tsv.gz")

_GID = re.compile(r'gene_id "([^"]+)"')
_GNAME = re.compile(r'gene_name "([^"]+)"')
_GTYPE = re.compile(r'gene_type "([^"]+)"')


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def main():
    # per gene: strand, chr, gene_name, gene_type, stop positions (max end for -, min start for +), utr intervals
    g_strand = {}
    g_chr = {}
    g_name = {}
    g_type = {}
    stop_max_end = defaultdict(int)   # for - strand threshold
    stop_min_start = defaultdict(lambda: 10**12)  # for + strand threshold
    utr_by_gene = defaultdict(list)   # (start,end)

    with (gzip.open(GTF, "rt") if str(GTF).endswith(".gz") else open(GTF)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, src, feat, start, end, score, strand, frame, attr = f[:9]
            if feat not in ("UTR", "stop_codon"):
                continue
            mg = _GID.search(attr)
            if not mg:
                continue
            gid = mg.group(1).split(".")[0]
            g_strand[gid] = strand
            g_chr[gid] = chrom
            mn = _GNAME.search(attr)
            if mn:
                g_name[gid] = mn.group(1)
            mt = _GTYPE.search(attr)
            if mt:
                g_type[gid] = mt.group(1)
            s, e = int(start), int(end)
            if feat == "stop_codon":
                stop_max_end[gid] = max(stop_max_end[gid], e)
                stop_min_start[gid] = min(stop_min_start[gid], s)
            else:  # UTR
                utr_by_gene[gid].append((s, e))

    rows = []
    for gid, utrs in utr_by_gene.items():
        if g_type.get(gid) != "protein_coding":
            continue
        strand = g_strand[gid]
        if strand == "+":
            thr = stop_min_start.get(gid)
            if thr is None:
                continue
            utr3 = [(s, e) for (s, e) in utrs if s > thr]
        else:
            thr = stop_max_end.get(gid)
            if thr is None or thr == 0:
                continue
            utr3 = [(s, e) for (s, e) in utrs if e < thr]
        if not utr3:
            continue
        merged = merge_intervals(utr3)
        for (s, e) in merged:
            rows.append({
                "gene_id": gid,
                "gene_name": g_name.get(gid, ""),
                "chr": g_chr[gid],
                "strand": strand,
                "utr3_start": s,
                "utr3_end": e,
            })
    df = pd.DataFrame(rows)
    with gzip.open(OUT, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)
    n_genes = df["gene_id"].nunique()
    print(f"[info] {len(df)} 3'UTR intervals across {n_genes} protein_coding genes -> {OUT}")


if __name__ == "__main__":
    main()

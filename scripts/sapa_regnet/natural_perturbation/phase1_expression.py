#!/usr/bin/env python3
"""Phase 1: per-sample per-domain expression of 15 neural RBPs + NUDT21 (CFIm25).

Reads SAW raw expression long tables (x, y, geneID, MIDIndex, readCount) in
chunks, converts (x, y) -> bin200 spot id  f"{x//200}_{y//200}", aligns with
spagapa domain labels (domains.csv row-aligned to coordinates.csv), aggregates
to mean per-bin readCount for each (sample, domain, gene).

Outputs (natural_perturbation/):
  rbp_domain_expression.csv   sample, domain, gene, mean_expr, total_reads, n_bins
"""
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT = PROJ / "pipeline_output/sapa_regnet/natural_perturbation"
OUT.mkdir(parents=True, exist_ok=True)

GTF = "/s1/SHARE/01_software/SAW_refs/Homo_sapiens_index/genes/Homo_sapiens.GRCh38.93.saw.gtf"

SAMPLES = {
    "GSM8882884": "D02266B4",
    "GSM8882885": "D02266C2",
    "GSM8882886": "D02266D2",
    "GSM8882887": "D02266D4",
}
TARGETS = ["ELAVL1", "ELAVL2", "ELAVL3", "ELAVL4", "NOVA1", "NOVA2", "PTBP1",
           "FUS", "TARDBP", "HNRNPA1", "HNRNPK", "KHSRP", "FMR1", "TIA1",
           "TIAL1", "NUDT21"]


def gtf_target_ids():
    """symbol -> list of ENSG ids (GRCh38.93 GTF)."""
    want = set(TARGETS)
    pat = re.compile(r'gene_id "([^"]+)"; gene_version "[^"]+"; gene_name "([^"]+)"')
    m = {}
    with open(GTF) as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            hit = pat.search(line)
            if hit and hit.group(2) in want:
                m.setdefault(hit.group(2), []).append(hit.group(1))
    return m


def load_domains(sample):
    exp = next((PROJ / "pipeline_output/stereo_expansion_downstream").glob(
        f"gse293464_{sample}_*/spagapa_run"))
    coords = pd.read_csv(PROJ / f"pipeline_output/gse293464_retina/{sample}_binned/coordinates.csv")
    dom = pd.read_csv(exp / "domains.csv")
    assert len(coords) == len(dom), f"{sample}: coords {len(coords)} vs domains {len(dom)}"
    sd = pd.DataFrame({"spot_id": coords["spot_id"], "domain": dom["domain"].to_numpy()})
    sd["domain"] = sd["domain"].astype(int)
    return sd


def process_sample(sample, plate, sym2ids):
    ids = {i for lst in sym2ids.values() for i in lst}
    id2sym = {i: s for s, lst in sym2ids.items() for i in lst}
    expr_f = (PROJ / f"pipeline_output/gse293464_retina/{sample}_saw/outs/feature_expression/"
              f"{plate}_raw_barcode_gene_exp.txt")
    spots = load_domains(sample)
    spot2dom = dict(zip(spots.spot_id, spots.domain))
    dom_counts = spots.domain.value_counts().to_dict()

    acc = []  # (symbol, spot_id, readCount summed per bin)
    t0 = time.time()
    reader = pd.read_csv(expr_f, sep="\t", usecols=["x", "y", "geneID", "readCount"],
                         dtype={"x": np.int64, "y": np.int64, "geneID": str,
                                "readCount": np.int64},
                         chunksize=5_000_000, engine="c")
    n_rows = 0
    for chunk in reader:
        n_rows += len(chunk)
        chunk = chunk[chunk.geneID.isin(ids)]
        if not len(chunk):
            continue
        chunk["bin"] = (chunk.x // 200).astype(str) + "_" + (chunk.y // 200).astype(str)
        chunk["domain"] = chunk.bin.map(spot2dom)
        chunk = chunk.dropna(subset=["domain"])
        if not len(chunk):
            continue
        chunk["symbol"] = chunk.geneID.map(id2sym)
        g = chunk.groupby(["symbol", "domain", "bin"], as_index=False)["readCount"].sum()
        acc.append(g)
        print(f"  {sample}: {n_rows:,} rows read, {len(acc)} chunks kept "
              f"({time.time()-t0:.0f}s)", flush=True)
    df = pd.concat(acc, ignore_index=True)
    df = df.groupby(["symbol", "domain", "bin"], as_index=False)["readCount"].sum()
    # per (symbol, domain): total reads, mean per bin over ALL bins of domain
    df["domain"] = df.domain.astype(int)
    rows = []
    for (sym, dom), g in df.groupby(["symbol", "domain"]):
        nb = dom_counts[dom]
        rows.append((sample, dom, sym, g.readCount.sum() / nb, int(g.readCount.sum()), nb))
    out = pd.DataFrame(rows, columns=["sample", "domain", "gene", "mean_expr",
                                      "total_reads", "n_bins"])
    return out


def main():
    sym2ids = gtf_target_ids()
    missing = [t for t in TARGETS if t not in sym2ids]
    if missing:
        print("MISSING in GTF:", missing)
        sys.exit(1)
    print("symbol->ENSG:", {k: v for k, v in sym2ids.items()})
    res = []
    for sample, plate in SAMPLES.items():
        print(f"processing {sample} ({plate}) ...", flush=True)
        res.append(process_sample(sample, plate, sym2ids))
    allx = pd.concat(res, ignore_index=True)
    allx.to_csv(OUT / "rbp_domain_expression.csv", index=False)
    print(f"wrote {OUT/'rbp_domain_expression.csv'}  rows={len(allx)}")
    piv = allx.pivot_table(index=["sample", "domain"], columns="gene", values="mean_expr")
    print(piv.round(2).to_string())


if __name__ == "__main__":
    main()

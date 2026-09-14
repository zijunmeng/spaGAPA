#!/usr/bin/env python3
"""Phase 3 prep: per-sample per-domain distal-PAS usage per gene.

1. Map GSE293464 scAPAtrap peaks (apa_sites.csv) to genes via GENCODE v44
   gene bounds (body overlap, strand-aware, protein-coding preferred).
2. Per gene with >=2 peaks: proximal = peak coord nearest TSS, distal = nearest
   TES (strand-aware) — same convention as track_a_rbp_apa_annotation.py.
3. distal usage per spot = u_distal / (u_prox + u_distal) from apa_matrix.csv
   (per-spot column normalization cancels in the within-gene ratio).
4. Domain mean usage (>=5 covered spots), then per-RBP high/low domain contrast
   is applied later in the analysis script.

Outputs:
  _cache/peak_gene_map_<sample>.csv   site_id, gene, role(proximal/distal)
  _cache/domain_usage_<sample>.csv    gene, domain, distal_usage, n_spots
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT = PROJ / "pipeline_output/sapa_regnet/natural_perturbation"
CACHE = OUT / "_cache"
CACHE.mkdir(parents=True, exist_ok=True)
GENE_BOUNDS = PROJ / "pipeline_output/sapa_regnet/_gene_bounds_gencode_v44.tsv.gz"
SAMPLES = ["GSM8882884", "GSM8882885", "GSM8882886", "GSM8882887"]
MIN_SPOTS_PER_DOMAIN = 5

def load_bounds():
    return pd.read_csv(GENE_BOUNDS, sep="\t")


def map_peaks_to_genes(peaks, bounds):
    """Assign each peak to the gene whose body overlaps it most (strand aware)."""
    b = bounds[["gene_name", "gene_type", "chr", "gene_start", "gene_end", "strand"]]
    # interval overlap join via pyinterval-free approach: sort + searchsorted
    hits = []
    b_by_chr = {c: g for c, g in b.groupby("chr")}
    for chrom, pk in peaks.groupby("chr"):
        g = b_by_chr.get(chrom)
        if g is None:
            continue
        gs = g.gene_start.to_numpy()
        ge = g.gene_end.to_numpy()
        order = np.argsort(gs)
        gs_s, ge_s = gs[order], ge[order]
        pk_start = pk.start.to_numpy()
        pk_end = pk.end.to_numpy()
        # genes whose start <= pk_end
        lo = np.searchsorted(gs_s, pk_end, side="right")
        rows = []
        for i, (ps, pe) in enumerate(zip(pk_start, pk_end)):
            cand = np.arange(0, lo[i])
            if len(cand) == 0:
                rows.append(None)
                continue
            ov = np.minimum(ge_s[cand], pe) - np.maximum(gs_s[cand], ps)
            cand = cand[ov > 0]
            if len(cand) == 0:
                rows.append(None)
                continue
            sub = g.iloc[order[cand]]
            same = (sub.strand.to_numpy() == pk.strand.iloc[i])
            score = ov[cand] + same * 1e9 + (sub.gene_type == "protein_coding").to_numpy() * 1e6
            best = np.argmax(score)
            rows.append(sub.gene_name.iloc[best])
        pk = pk.assign(gene=rows)
        hits.append(pk)
    mapped = pd.concat(hits, ignore_index=True)
    return mapped


def assign_prox_dist(mapped):
    """For each gene with >=2 peaks define proximal (nearest TSS) and distal
    (nearest TES) peak."""
    recs = []
    for gene, g in mapped.dropna(subset=["gene"]).groupby("gene"):
        if len(g) < 2:
            continue
        strand = g.strand.mode().iloc[0]
        c = g.coord.to_numpy()
        site_ids = g.site_id.to_numpy()
        if strand == "+":
            prox_i, dist_i = np.argmin(c), np.argmax(c)
        else:
            prox_i, dist_i = np.argmax(c), np.argmin(c)
        if c[prox_i] == c[dist_i]:
            continue
        recs.append((gene, site_ids[prox_i], "proximal"))
        recs.append((gene, site_ids[dist_i], "distal"))
    return pd.DataFrame(recs, columns=["gene", "site_id", "role"])


def process_sample(sample):
    t0 = time.time()
    binned = PROJ / f"pipeline_output/gse293464_retina/{sample}_binned"
    peaks = pd.read_csv(binned / "apa_sites.csv", dtype={"site_id": str, "peakID": str})
    bounds = load_bounds()
    cache_f = CACHE / f"peak_gene_map_{sample}.csv"
    if cache_f.is_file():
        roles = pd.read_csv(cache_f, dtype={"site_id": str})
    else:
        mapped = map_peaks_to_genes(peaks, bounds)
        roles = assign_prox_dist(mapped)
        roles.to_csv(cache_f, index=False)
        mapped.to_csv(CACHE / f"peak_gene_map_full_{sample}.csv", index=False)
    print(f"[{sample}] genes with prox/dist pair: {roles.gene.nunique():,} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # domain per spot
    coords = pd.read_csv(binned / "coordinates.csv", dtype={"spot_id": str})
    dom_f = next((PROJ / "pipeline_output/stereo_expansion_downstream").glob(
        f"gse293464_{sample}_*/spagapa_run/domains.csv"))
    dom = pd.read_csv(dom_f)
    spot2dom = pd.Series(dom.domain.to_numpy(), index=coords.spot_id).to_dict()
    spot_dom = coords.spot_id.map(spot2dom).to_numpy()
    need = set(roles.site_id)
    mat = pd.read_csv(binned / "apa_matrix.csv", index_col=0)  # sites x spots
    mat = mat.loc[mat.index.intersection(need)]
    spots = mat.columns.to_numpy()
    spot_dom_arr = np.array([spot2dom.get(s, -1) for s in spots], dtype=int)
    M = mat.to_numpy(dtype=np.float32)

    dom_ids = sorted(set(spot_dom_arr[spot_dom_arr >= 0]))
    out_rows = []
    prox = roles[roles.role == "proximal"].set_index("gene").site_id.to_dict()
    dist = roles[roles.role == "distal"].set_index("gene").site_id.to_dict()
    idx = {s: i for i, s in enumerate(mat.index)}
    n_skip = 0
    n_zero_domain = 0
    for gene in prox:
        if gene not in dist or prox[gene] not in idx or dist[gene] not in idx:
            n_skip += 1
            continue
        d = M[idx[dist[gene]]]
        p = M[idx[prox[gene]]]
        tot = d + p
        frac = np.where(tot > 0, d / np.where(tot > 0, tot, 1), np.nan)
        n_ok_dom = 0
        for dm in dom_ids:
            m = spot_dom_arr == dm
            f = frac[m]
            k = np.isfinite(f)
            if k.sum() >= MIN_SPOTS_PER_DOMAIN:
                n_ok_dom += 1
                out_rows.append((sample, gene, int(dm), float(f[k].mean()), int(k.sum())))
        if n_ok_dom == 0:
            n_zero_domain += 1
    print(f"  DEBUG genes={len(prox)} skipped={n_skip} zero_domain={n_zero_domain} "
          f"out={len(out_rows)} dom_ids={len(dom_ids)} mat={mat.shape}", flush=True)
    res = pd.DataFrame(out_rows, columns=["sample", "gene", "domain",
                                          "distal_usage", "n_spots"])
    res.to_csv(CACHE / f"domain_usage_{sample}.csv", index=False)
    print(f"[{sample}] domain-usage rows: {len(res):,} genes: {res.gene.nunique():,} "
          f"({time.time()-t0:.0f}s)", flush=True)
    return res


def main():
    all_rows = []
    for s in SAMPLES:
        all_rows.append(process_sample(s))
    df = pd.concat(all_rows, ignore_index=True)
    df.to_csv(CACHE / "domain_usage_all.csv", index=False)
    print(f"wrote {CACHE/'domain_usage_all.csv'} rows={len(df):,} "
          f"genes={df.gene.nunique():,}")


if __name__ == "__main__":
    main()

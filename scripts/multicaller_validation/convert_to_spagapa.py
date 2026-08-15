#!/usr/bin/env python
"""Convert Sierra CountPeaks output to the spaGAPA processed-data format.

Inputs (pipeline_output/multicaller_validation/sierra/sierra_counts/):
  matrix.mtx.gz   sparse peak x barcode UMI counts (peaks in table order)
  barcodes.tsv.gz 3010 spot barcodes
  sitenames.tsv.gz peak ids 'Gene:Chr:start-end:Strand'

Plus sierra_peaks.txt for the fitted summit (Fit.max.pos / MaxPosition).

Outputs (same conventions as data/processed/gse183456_gsm6047774_scapatrap/):
  apa_sites.csv.gz   site_id, gene_id, gene_name, chr, start, end, strand, coord
  apa_site_counts.csv.gz  site x spot raw UMI counts (genes with >=2 sites only)
  apa_matrix.csv     site x spot usage fractions (NaN where parent < 5)
  qc_summary.json
Coordinates are reused verbatim from the scAPAtrap baseline (same Space
Ranger run, same 3010 spots).
"""
import os, sys, json, gzip
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd
import scipy.io as sio

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
OUT = os.path.join(ROOT, "pipeline_output/multicaller_validation/sierra")
COUNTS = os.path.join(OUT, "sierra_counts")
BASE = os.path.join(ROOT, "data/processed/gse183456_gsm6047774_scapatrap")
MIN_PARENT = 5


def main():
    # ---- sierra count matrix ----
    with gzip.open(os.path.join(COUNTS, "matrix.mtx.gz"), "rt") as f:
        mat = sio.mmread(f).tocsr()          # peaks x barcodes
    barcodes = [l.strip() for l in gzip.open(os.path.join(COUNTS, "barcodes.tsv.gz"), "rt")]
    sitenames = [l.strip() for l in gzip.open(os.path.join(COUNTS, "sitenames.tsv.gz"), "rt")]
    print(f"[conv] matrix {mat.shape}, {len(barcodes)} barcodes, {len(sitenames)} sitenames")

    counts = pd.DataFrame.sparse.from_spmatrix(mat, index=sitenames, columns=barcodes)
    counts = counts.sparse.to_dense()

    # ---- parse sitenames + join summit from sierra_peaks.txt ----
    peaks = pd.read_table(os.path.join(OUT, "sierra_peaks.txt"))
    peaks = peaks[peaks["Fit.start"].apply(lambda x: str(x).lstrip("-").isdigit())]
    for c in ("MaxPosition", "Fit.max.pos", "Fit.start", "Fit.end"):
        peaks[c] = pd.to_numeric(peaks[c], errors="coerce")
    peaks = peaks.dropna(subset=["Fit.start", "Fit.end"])
    peaks = peaks[peaks["Fit.start"] < peaks["Fit.end"]]

    site = pd.DataFrame({"sitename": counts.index})
    parts = site["sitename"].str.rsplit(":", n=3, expand=True)  # robust to ':' in gene symbols
    assert parts.shape[1] == 4, f"unexpected sitename format: {site['sitename'].iloc[:3].tolist()}"
    site["gene_name"] = parts[0]
    site["chr"] = parts[1]
    rng = parts[2].str.split("-", expand=True)
    site["start"] = rng[0].astype(int)
    site["end"] = rng[1].astype(int)
    site["strand"] = parts[3]

    # join the fitted summit (Fit.max.pos) on (gene, chr, start, end, strand)
    peaks["sitename"] = (peaks["Gene"] + ":" + peaks["Chr"] + ":" +
                         peaks["Fit.start"].astype(int).astype(str) + "-" +
                         peaks["Fit.end"].astype(int).astype(str) + ":" +
                         peaks["Strand"].astype(str))
    summ = peaks.drop_duplicates("sitename").set_index("sitename")
    # if multiple rows share the sitename key, take the max-coverage summit (first kept)
    site["coord"] = site["sitename"].map(summ["Fit.max.pos"])
    site["maxpos"] = site["sitename"].map(summ["MaxPosition"])
    site["coord"] = site["coord"].fillna(site["maxpos"])
    site["coord"] = site["coord"].fillna(((site["start"] + site["end"]) / 2).round().astype(int))
    n_missing = int(site["coord"].isna().sum())
    # Sierra strand convention is 1 / -1; normalise to the +/- used by the
    # baseline apa_sites table (needed for strand-aware overlap + oriented split).
    site["strand"] = site["strand"].astype(str).map({"1": "+", "-1": "-"}).fillna(site["strand"])

    # ---- gene ids: map gene_name -> ensembl id from baseline apa_sites ----
    base_sites = pd.read_csv(os.path.join(BASE, "apa_sites.csv.gz"))
    name2id = base_sites.dropna(subset=["gene_id"]).drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
    site["gene_id"] = site["gene_name"].map(name2id)

    # ---- site ids: sierra_0000... keep stable order ----
    site["site_id"] = ["sierra_%05d" % i for i in range(len(site))]

    # ---- restrict to genes with >= 2 peaks (baseline convention) ----
    ngene = site.groupby("gene_name")["site_id"].transform("count")
    keep = ngene >= 2
    site = site[keep].reset_index(drop=True)
    counts = counts.loc[site["sitename"]].reset_index(drop=True)
    counts.index = site["site_id"].values
    counts.columns.name = None

    # ---- usage matrix ----
    gene_of = site["gene_name"].values
    V = counts.values.astype(float)
    usage = np.full_like(V, np.nan)
    # vectorised per gene
    order = np.argsort(gene_of, kind="stable")
    sorted_genes = gene_of[order]
    bounds = np.flatnonzero(np.r_[True, sorted_genes[1:] != sorted_genes[:-1]])
    starts = np.r_[bounds, len(sorted_genes)]
    for i in range(len(bounds)):
        rows = order[starts[i]:starts[i + 1]]
        total = V[rows].sum(axis=0)
        usage[rows] = np.where(total >= MIN_PARENT, V[rows] / np.where(total == 0, 1, total), np.nan)
    usage_df = pd.DataFrame(usage, index=counts.index, columns=counts.columns)

    # ---- write ----
    sites_out = site[["site_id", "gene_id", "gene_name", "chr", "start", "end", "strand", "coord"]].copy()
    sites_out.insert(1, "sitename", site["sitename"])
    sites_out.to_csv(os.path.join(OUT, "apa_sites.csv.gz"), index=False, compression="gzip")
    counts.to_csv(os.path.join(OUT, "apa_site_counts.csv.gz"), compression="gzip")
    usage_df.to_csv(os.path.join(OUT, "apa_matrix.csv"))

    qc = {
        "dataset": "gse183456_gsm6047774_sierra",
        "caller": "Sierra 0.99.27 (Winnie09/Sierra)",
        "platform": "10x Visium",
        "n_spots": int(counts.shape[1]),
        "n_peaks_raw_all_genes": int(mat.shape[0]),
        "n_sites_ge2": int(counts.shape[0]),
        "n_genes_ge2": int(site["gene_name"].nunique()),
        "min_parent_count": MIN_PARENT,
        "summit_missing_filled_with_midpoint": n_missing,
        "files": {
            "apa_matrix": "apa_matrix.csv",
            "apa_sites": "apa_sites.csv.gz",
            "apa_site_counts": "apa_site_counts.csv.gz",
            "coordinates": "../../data/processed/gse183456_gsm6047774_scapatrap/coordinates.csv (reused)",
        },
    }
    with open(os.path.join(OUT, "qc_summary.json"), "w") as f:
        json.dump(qc, f, indent=2)
    print("[conv] wrote apa_sites.csv.gz / apa_site_counts.csv.gz / apa_matrix.csv / qc_summary.json")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()

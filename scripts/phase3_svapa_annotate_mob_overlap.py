#!/usr/bin/env python
"""Phase 3 follow-up #1: annotate SVAPA results + MOB SVAPA + competitor overlap.

Tasks
-----
1. Annotate existing peak-level SVAPA outputs (gse263789, gse220442) with
   gene names and aggregate to gene level.
2. Run SVAPA on the MOB gene-level APA matrix.
3. Overlap our MOB SVAPA hits with the competitor canonical MOB spatial-APA
   gene list (head-to-head validation).
4. Sanity report (overlap count, Pde1c recovery, etc.).

Environment (S91): OPENBLAS_NUM_THREADS=8, TMPDIR=/s3/mengzijun/tmp.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# project root = parent of scripts/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from spagapa.analysis.svapa import svapa  # noqa: E402

OUT = ROOT / "pipeline_output" / "svapa"
OUT.mkdir(parents=True, exist_ok=True)

GTF = Path("/s1/SHARE/00_ref_genecode/gencode.vM33.annotation.gtf")

COMP = ROOT / "pipeline_output" / "known_gene_validation" / "mob_spatial_apa_genes.csv"


# ---------------------------------------------------------------------------
# 1. Annotation helpers
# ---------------------------------------------------------------------------
def _parse_gene_name(attr: str) -> str:
    m = re.search(r'gene_name "([^"]+)"', attr)
    return m.group(1) if m else ""


def load_gtf_genes(gtf: Path = GTF) -> pd.DataFrame:
    """Parse GTF into a gene-level table: chr, start, end, strand, gene_id, gene_name."""
    rows = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
            name = _parse_gene_name(attr)
            if name:
                rows.append((chrom, start, end, strand, name))
    g = pd.DataFrame(rows, columns=["chr", "start", "end", "strand", "gene_name"])
    # collapse duplicate gene_name -> union interval (keep min start / max end)
    g = (
        g.sort_values(["chr", "gene_name", "start"])
        .groupby(["chr", "gene_name"], as_index=False)
        .agg(start=("start", "min"), end=("end", "max"),
             strand=("strand", "first"))
    )
    return g


def annotate_peaks_to_genes(sites: pd.DataFrame,
                            genes: pd.DataFrame,
                            flank: int = 5000) -> pd.DataFrame:
    """Assign each peak (chr, coord, strand) to a gene.

    Rule (in priority order):
      1. peak coord falls inside a gene body (same strand preferred)
      2. peak within `flank` bp upstream of a gene (3'UTR APA is downstream
         of gene body; scAPAtrap coord = end = polyA site, typically in/near
         3'UTR, so we also accept being within flank downstream)
      3. nearest gene overall (distance reported)
    Returns sites with added 'gene_name' (filled where empty) + 'dist'.
    """
    sites = sites.copy()
    # numeric coord
    sites["coord"] = pd.to_numeric(sites["coord"], errors="coerce")
    sites["start"] = pd.to_numeric(sites.get("start"), errors="coerce")
    sites["end"] = pd.to_numeric(sites.get("end"), errors="coerce")
    sites["gene_name_new"] = ""

    # index genes by chr for speed
    by_chr = {chrom: g.reset_index(drop=True) for chrom, g in genes.groupby("chr")}
    assigned = np.empty(len(sites), dtype=object)
    dists = np.full(len(sites), np.nan)

    for i, row in enumerate(sites.itertuples(index=False)):
        chrom = getattr(row, "chr", None)
        coord = getattr(row, "coord", np.nan)
        strand = getattr(row, "strand", None)
        if chrom not in by_chr or not np.isfinite(coord):
            continue
        g = by_chr[chrom]
        gs, ge, gname, gstrand = g["start"].values, g["end"].values, \
            g["gene_name"].values, g["strand"].values
        pos = int(coord)
        # 1. inside any gene body (prefer same strand)
        inside = (gs <= pos) & (ge >= pos)
        if inside.any():
            cand = np.where(inside)[0]
            if strand is not None and len(cand) > 1:
                same = cand[gstrand[cand] == strand]
                if len(same):
                    cand = same
            assigned[i] = gname[cand[0]]
            dists[i] = 0.0
            continue
        # 2. within flank upstream/downstream (polyA site near 3' end)
        # distance to gene body (0 if inside, else min gap)
        d_left = gs - pos  # >0 if peak is upstream of gene start
        d_right = pos - ge  # >0 if peak is downstream of gene end
        gap = np.where((d_left > 0), d_left,
                       np.where((d_right > 0), d_right, 0))
        within = np.where((gap >= 0) & (gap <= flank))[0]
        if len(within):
            # nearest among within-flank
            j = within[np.argmin(gap[within])]
            assigned[i] = gname[j]
            dists[i] = float(gap[j])
            continue
        # 3. nearest gene overall
        j = int(np.argmin(gap))
        assigned[i] = gname[j]
        dists[i] = float(gap[j])
    sites["gene_name_new"] = assigned
    sites["dist"] = dists
    return sites


def _read_matrix_index(apa_matrix_csv: Path) -> list[str]:
    """Read only the first column (site_id index) of a (possibly huge) csv."""
    import csv
    idx = []
    with open(apa_matrix_csv) as fh:
        reader = csv.reader(fh)
        next(reader, None)  # header
        for row in reader:
            idx.append(row[0])
    return idx


def map_svapa_peaks_to_gene(svapa_csv: Path, sites_csv: Path,
                            apa_matrix_csv: Path,
                            genes: pd.DataFrame | None) -> pd.DataFrame:
    """Join SVAPA peak ids (Gene_NNN) -> gene.

    SVAPA assigns ``Gene_{i}`` = the i-th row of the apa_matrix it was run on.
    The apa_matrix index name is ``site_id`` (peak ids like ``peak_27206``),
    and its row order is NOT the same as ``apa_sites`` row order, so we must
    load the apa_matrix index positionally, then look up the gene_name in
    ``apa_sites`` by ``site_id``.  If ``apa_sites.gene_name`` is empty
    (gse263789), fall back to GTF coordinate annotation.
    """
    res = pd.read_csv(svapa_csv)
    res["_pos"] = res["peak"].astype(str).str.extract(r"Gene_(\d+)").astype(int)

    mat_idx = _read_matrix_index(apa_matrix_csv)
    n_mat = len(mat_idx)
    bad = res["_pos"].max() >= n_mat
    if bad:
        raise RuntimeError(
            f"SVAPA positional index {res['_pos'].max()} >= apa_matrix rows "
            f"{n_mat} for {apa_matrix_csv}")
    res["site_id"] = [mat_idx[i] for i in res["_pos"].values]

    sites = pd.read_csv(sites_csv)
    keep = ["site_id", "chr", "start", "end", "strand", "coord", "gene_name"]
    keep = [c for c in keep if c in sites.columns]
    merged = res.merge(sites[keep], on="site_id", how="left")

    # fill missing gene_name via GTF coordinate annotation
    if merged["gene_name"].isna().any() or (merged["gene_name"].fillna("") == "").any():
        if genes is None:
            print(f"  [warn] missing gene_name in {svapa_csv.name} and no GTF")
        else:
            miss = merged["gene_name"].isna() | (merged["gene_name"].fillna("") == "")
            need = merged.loc[miss, ["chr", "start", "end", "strand", "coord"]].copy()
            if len(need):
                ann = annotate_peaks_to_genes(
                    need.rename(columns={"gene_name": "_drop"}), genes)
                merged.loc[miss, "gene_name"] = ann["gene_name_new"].values
                merged.loc[miss, "peak_dist"] = ann["dist"].values
    merged.drop(columns=["_pos"], inplace=True)
    return merged


def aggregate_to_gene(peaks_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate peak-level SVAPA -> gene-level (max Moran's I per gene)."""
    df = peaks_df.dropna(subset=["gene_name"])
    df = df[df["gene_name"].astype(str).str.strip() != ""]
    df = df.copy()
    df["morans_i"] = pd.to_numeric(df["morans_i"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    agg = (
        df.sort_values("morans_i", ascending=False)
        .groupby("gene_name", as_index=False)
        .agg(morans_i=("morans_i", "first"),
             pvalue=("pvalue", "first"),
             padj=("padj", "first"),
             n_peaks=("gene_name", "size"))
    )
    agg = agg.sort_values("morans_i", ascending=False,
                          na_position="last").reset_index(drop=True)
    return agg


# ---------------------------------------------------------------------------
# 2. MOB SVAPA
# ---------------------------------------------------------------------------
def run_mob_svapa() -> pd.DataFrame:
    mob_dir = ROOT / "data" / "processed" / "mob_st11"
    mat = pd.read_csv(mob_dir / "apa_matrix.csv", index_col=0)
    coords = pd.read_csv(mob_dir / "coordinates.csv", index_col="spot_id")
    # align columns (spots)
    common = [c for c in mat.columns if c in coords.index]
    mat = mat[common]
    coords = coords.loc[common]
    coord_arr = coords[["x", "y"]].values.astype(float)
    print(f"  MOB: {mat.shape[0]} genes x {mat.shape[1]} spots; "
          f"obs fraction={np.isfinite(mat.values.astype(float)).mean():.3f}")
    res = svapa(mat.astype(float), coord_arr, k=8, n_perm=200,
                seed=42, min_obs=10, include_geary=False)
    res = res.rename(columns={"morans_i_pvalue": "pvalue"})
    # ensure padj recomputed on pvalue
    try:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(res["pvalue"].fillna(1.0).values,
                                      method="fdr_bh")
        res["padj"] = padj
    except Exception:
        res["padj"] = res["pvalue"]
    res = res.sort_values(["padj", "morans_i"],
                          ascending=[True, False]).reset_index(drop=True)
    res.insert(0, "rank", np.arange(1, len(res) + 1))
    return res


# ---------------------------------------------------------------------------
# 3. Overlap with competitor list
# ---------------------------------------------------------------------------
def overlap_report(mob_res: pd.DataFrame, top_n: int = 50) -> str:
    comp = pd.read_csv(COMP)
    comp_genes = comp["gene"].astype(str).str.strip().tolist()
    n_comp = len(comp_genes)
    comp_set = set(comp_genes)

    # our top-N by padj then morans_i
    top = mob_res.head(top_n).copy()
    top_set = set(top["gene"].astype(str))

    # full overlap (any rank)
    mob_genes_all = set(mob_res["gene"].astype(str))
    recovered_all = [g for g in comp_genes if g in mob_genes_all]
    # top-N overlap
    recovered_top = []
    for g in comp_genes:
        if g in top_set:
            row = mob_res[mob_res["gene"].astype(str) == g].iloc[0]
            recovered_top.append((g, int(row["rank"]),
                                  float(row["morans_i"]),
                                  float(row["padj"])))

    md = []
    md.append("# MOB SVAPA vs Competitor Canonical Spatial-APA Genes\n")
    md.append(f"- Competitor canonical list: **{n_comp} genes** "
              f"(`mob_spatial_apa_genes.csv`; stAPAminer Ji 2023 + spvAPA Zhang 2025).\n")
    md.append(f"- Our MOB SVAPA list: **{len(mob_res)} genes** "
              f"(Moran's I, n_perm=200, BH-FDR).\n")
    md.append(f"- Head-to-head: top **{top_n}** our SVAPA genes vs competitor list.\n")
    md.append("")
    md.append(f"## Summary\n")
    md.append(f"- Competitor genes recovered in our **full** MOB SVAPA list: "
              f"**{len(recovered_all)}/{n_comp}** "
              f"({100*len(recovered_all)/n_comp:.1f}%).\n")
    md.append(f"- Competitor genes in our **top-{top_n}** SVAPA: "
              f"**{len(recovered_top)}/{n_comp}** "
              f"({100*len(recovered_top)/n_comp:.1f}%).\n")
    md.append("")
    # Pde1c sanity
    pde1c = mob_res[mob_res["gene"].astype(str) == "Pde1c"]
    if len(pde1c):
        pr = pde1c.iloc[0]
        md.append(f"## Sanity: Pde1c\n")
        md.append(f"- **Pde1c recovered** at rank {int(pr['rank'])}, "
                  f"Moran's I = {float(pr['morans_i']):.4f}, "
                  f"padj = {float(pr['padj']):.3g}.\n")
    else:
        md.append(f"## Sanity: Pde1c\n- **Pde1c NOT in MOB matrix/SVAPA output.**\n")
    md.append("")

    md.append(f"## Competitor genes in our top-{top_n} MOB SVAPA\n")
    md.append("| competitor gene | our rank | Moran's I | padj | evidence (competitor) |\n")
    md.append("|---|---:|---:|---:|---|\n")
    # map gene -> evidence
    ev = dict(zip(comp["gene"].astype(str), comp["evidence_strength"].astype(str)))
    if recovered_top:
        for g, rk, mi, pj in sorted(recovered_top, key=lambda x: x[1]):
            md.append(f"| {g} | {rk} | {mi:.4f} | {pj:.3g} | {ev.get(g,'')} |\n")
    else:
        md.append("| _(none in top-50; see full-list table below)_ | | | | |\n")
    md.append("")

    md.append(f"## Competitor genes recovered anywhere in our MOB SVAPA "
              f"({len(recovered_all)}/{n_comp})\n")
    md.append("| competitor gene | our rank | Moran's I | padj | n_obs | evidence |\n")
    md.append("|---|---:|---:|---:|---:|---|\n")
    for g in comp_genes:
        r = mob_res[mob_res["gene"].astype(str) == g]
        if len(r):
            rr = r.iloc[0]
            md.append(f"| {g} | {int(rr['rank'])} | {float(rr['morans_i']):.4f} "
                      f"| {float(rr['padj']):.3g} | {int(rr['n_obs'])} | "
                      f"{ev.get(g,'')} |\n")
        else:
            md.append(f"| {g} | — (not in matrix) | — | — | — | {ev.get(g,'')} |\n")
    md.append("")

    # competitor genes NOT in matrix at all
    not_in_mat = [g for g in comp_genes if g not in mob_genes_all]
    md.append(f"## Competitor genes absent from MOB apa_matrix ({len(not_in_mat)})\n")
    if not_in_mat:
        md.append(", ".join(not_in_mat) + "\n")
    md.append("")
    return "".join(md)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

    # ---- load GTF once (only needed for gse263789 coordinate annotation) ----
    print("[1/4] Loading GTF for peak->gene coordinate annotation ...")
    genes = load_gtf_genes(GTF)
    print(f"  GTF genes (collapsed by name): {len(genes)}")

    # ---- Task 1: annotate gse220442 ----
    print("[2/4] Annotating gse220442 SVAPA peaks ...")
    pdir_220442 = (ROOT / "data" / "processed" /
                   "gse220442_gsm6801751_scapatrap")
    sites_220442 = pdir_220442 / "apa_sites.csv.gz"
    mat_220442 = pdir_220442 / "apa_matrix.csv"
    pk_220442 = map_svapa_peaks_to_gene(OUT / "gse220442_svapa.csv",
                                        sites_220442, mat_220442, genes)
    gene_220442 = aggregate_to_gene(pk_220442)
    # write peak-level annotated + gene-level
    pk_220442.to_csv(OUT / "gse220442_svapa_annotated.csv", index=False)
    # gene-level summary table (reuse same filename stem per task spec, but
    # keep peak-level detail too): write gene-level as the *_annotated.csv and
    # a separate *_genelevel.csv for clarity.
    gene_220442.to_csv(OUT / "gse220442_svapa_genelevel.csv", index=False)
    print(f"  gse220442: {len(pk_220442)} peaks annotated; "
          f"{len(gene_220442)} unique genes. "
          f"Top genes: {', '.join(gene_220442['gene_name'].head(10).tolist())}")

    # ---- Task 1: annotate gse263789 ----
    print("[3/4] Annotating gse263789 SVAPA peaks (GTF coordinate fallback) ...")
    pdir_263789 = (ROOT / "pipeline_output" / "gse263789_stereo_pilot" /
                   "spagapa_downstream_full" / "binned_200")
    sites_263789 = pdir_263789 / "apa_sites.csv.gz"
    mat_263789 = pdir_263789 / "apa_matrix.csv"
    pk_263789 = map_svapa_peaks_to_gene(OUT / "gse263789_svapa.csv",
                                        sites_263789, mat_263789, genes)
    gene_263789 = aggregate_to_gene(pk_263789)
    pk_263789.to_csv(OUT / "gse263789_svapa_annotated.csv", index=False)
    gene_263789.to_csv(OUT / "gse263789_svapa_genelevel.csv", index=False)
    n_filled = (pk_263789["gene_name"].fillna("").astype(str).str.strip() != "").sum()
    print(f"  gse263789: {len(pk_263789)} peaks, {n_filled} with gene name; "
          f"{len(gene_263789)} unique genes. "
          f"Top genes: {', '.join(gene_263789['gene_name'].head(10).tolist())}")

    # ---- Task 2: MOB SVAPA ----
    print("[4/4] Running SVAPA on MOB gene-level APA matrix ...")
    mob_res = run_mob_svapa()
    mob_res.to_csv(OUT / "mob_svapa.csv", index=False)
    print(f"  MOB SVAPA: {len(mob_res)} genes. "
          f"Top: {', '.join(mob_res['gene'].head(10).tolist())}")

    # ---- Task 3 & 4: overlap report ----
    report = overlap_report(mob_res, top_n=50)
    (OUT / "mob_overlap_report.md").write_text(report)
    print("\n=== Overlap report written ===")
    # echo summary lines to stdout for the parent agent
    for line in report.splitlines():
        if line.startswith("- ") and ("recovered" in line or "top-" in line
                                      or "Pde1c" in line):
            print("  " + line)

    print("\nFiles written:")
    for f in ["gse220442_svapa_annotated.csv", "gse220442_svapa_genelevel.csv",
              "gse263789_svapa_annotated.csv", "gse263789_svapa_genelevel.csv",
              "mob_svapa.csv", "mob_overlap_report.md"]:
        p = OUT / f
        print(f"  {p}  ({p.stat().st_size if p.exists() else 'MISSING'} bytes)")


if __name__ == "__main__":
    main()

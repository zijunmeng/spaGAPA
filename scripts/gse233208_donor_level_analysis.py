#!/usr/bin/env python3
"""GSE233208 human Visium donor-level differential APA (external-review Phase 2).

Upgrade of the AD differential APA evidence from spot-level pseudoreplication
(GSE263789 mouse: 2 AD vs 1 WT) to donor-level statistics in an independent
human cohort with biological replicates:

  Primary cohort : 16 VisiumHuman L004 libraries = 16 independent donors
                   (Miyoshi/Morabito et al., Nat Genet 2024; GSE233208)
                   4 Control / 4 earlyAD / 4 AD / 4 AD_DS, PFC, fresh-frozen.
  Per sample     : scAPAtrap apa_site_counts (peak x spot counts) ->
                   gene annotation (dual GRCh38/GRCm39 GTF, GRCh38 only) ->
                   canonical median-split proximal/distal gene index ->
                   pseudobulk median distal-usage per gene (1 number per
                   sample -- the statistical unit is the DONOR, never the spot).
  Statistics     : Wilcoxon rank-sum across donors, Cliff's delta effect size,
                   bootstrap 95% CI, BH-FDR. Contrasts: AD vs Control (primary),
                   earlyAD vs Control, AD_DS vs Control, all-disease vs Control.
  Consistency    : direction concordance with GSE263789 mouse AD deltas
                   (mouse->human symbol capitalization) on shared genes.

Inputs : pipeline_output/gse233208_visium/<SRR>/{apa_site_counts.csv.gz,
         peaks_meta.csv.gz} for the 16 primary-cohort SRRs.
Outputs: pipeline_output/gse233208_donor_level/
         per_sample_gene_usage.csv, donor_level_differential_apa.csv,
         consistency_with_gse263789.csv, fig_donor_level_volcano.png,
         sample_qc.csv, summary.json
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("TMPDIR", "/s3/mengzijun/tmp")

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
SR_ROOT = ROOT / "pipeline_output/gse233208_visium"
OUT = ROOT / "pipeline_output/gse233208_donor_level"
OUT.mkdir(parents=True, exist_ok=True)

HUMAN_GTF = "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
MOUSE_REF = ROOT / "pipeline_output/gse263789_ad_vs_wt_differential/differential_apa_results.csv"

MIN_PARENT = 5      # canonical build_gene_index threshold (matches Visium diff)
MIN_OBS_SPOTS = 10  # gene usable in a sample if observed in >= this many spots
MIN_POOLED_READS = 20  # min pooled (prox+dist) reads for donor-level usage
N_BOOT = 2000
RNG = np.random.default_rng(20260917)


# ---------------------------------------------------------------- GTF genes
@dataclass
class GeneRecord:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    width: int


def parse_gtf_attrs(attr: str) -> dict:
    parsed = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def read_genes_from_gtf(gtf_path):
    """GRCh38 genes only from the dual-species GTF."""
    import bisect
    genes = {}
    n = 0
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            if not f[0].startswith("chr"):
                continue
            attrs = parse_gtf_attrs(f[8])
            gid, gname = attrs.get("gene_id"), attrs.get("gene_name", attrs.get("gene_id"))
            if not gid:
                continue
            start, end = int(f[3]), int(f[4])
            genes.setdefault((f[0], f[6]), []).append(GeneRecord(f[0], start, end, f[6], gid, gname, end - start + 1))
            n += 1
    for records in genes.values():
        records.sort(key=lambda r: r.start)
    print(f"[gtf] {n} human genes", flush=True)
    return genes


def assign_genes(peaks: pd.DataFrame, genes: dict) -> pd.DataFrame:
    """Assign each peak to the smallest gene body containing it (canonical)."""
    import bisect
    peaks = peaks.copy()
    peaks["coord"] = np.where(peaks["strand"].astype(str) == "-",
                              peaks["end"], peaks["start"])
    out = np.full(len(peaks), np.nan, dtype=object)
    # vectorized-by-contig sweep (same algorithm as gse263789 script)
    coord = peaks["coord"].to_numpy()
    for (chrom, strand), sub_pos in peaks.groupby(["chr", "strand"], sort=False).indices.items():
        records = genes.get((chrom, strand))
        if not records:
            continue
        starts = [r.start for r in records]
        active, cursor = [], 0
        for pos in sub_pos[np.argsort(coord[sub_pos])]:
            c = int(coord[pos])
            upto = bisect.bisect_right(starts, c, lo=cursor)
            active.extend(records[cursor:upto])
            cursor = upto
            active = [r for r in active if r.end >= c]
            hits = [r for r in active if r.start <= c <= r.end]
            if hits:
                best = min(hits, key=lambda r: r.width)
                out[pos] = best.gene_name
    peaks["gene_name"] = out
    return peaks


# ------------------------------------------------------- per-sample process
def load_sample(srr: str, genes: dict):
    """peak x spot sparse counts -> gene x spot distal-usage -> pseudobulk."""
    from scipy.sparse import coo_matrix
    t0 = time.time()
    cnt_path = SR_ROOT / srr / "apa_site_counts.csv.gz"
    meta_path = SR_ROOT / srr / "peaks_meta.csv.gz"
    counts = pd.read_csv(cnt_path)
    counts.columns = ["peak_id", "spot_id", "count"]
    peaks = pd.read_csv(meta_path)
    peaks = peaks.set_index("peakID").loc[  # align to unique peaks in counts
        pd.unique(counts["peak_id"])].reset_index()
    # dense pivot is safe: ~26k peaks x ~5k spots = 130M floats -> keep sparse
    pids = pd.unique(counts["peak_id"])
    pidx = {p: i for i, p in enumerate(pids)}
    sids = pd.unique(counts["spot_id"])
    sidx = {s: i for i, s in enumerate(sids)}
    pi = counts["peak_id"].map(pidx).to_numpy()
    si = counts["spot_id"].map(sidx).to_numpy()
    M = coo_matrix(
        (counts["count"].astype(float).to_numpy(), (pi, si)),
        shape=(len(pids), len(sids))).tocsr()

    peaks = assign_genes(peaks, genes)
    n_anno = peaks["gene_name"].notna().sum()
    print(f"  [{srr}] {len(pids)} peaks x {len(sids)} spots, "
          f"{n_anno} annotated ({time.time()-t0:.0f}s)", flush=True)

    # gene-level median-split index, sparse-friendly
    peaks["oriented"] = np.where(peaks["strand"].astype(str) == "-",
                                 -peaks["coord"].astype(float),
                                 peaks["coord"].astype(float))
    gene_rows = {}
    pooled_rows = {}
    sub_tab = peaks[["gene_name", "oriented"]].copy()
    sub_tab["row"] = np.arange(len(peaks))
    for gene, sub in sub_tab.dropna(subset=["gene_name"]).groupby("gene_name"):
        rows = sub["row"].to_numpy()
        if len(rows) < 2:
            continue
        order = rows[np.argsort(sub["oriented"].to_numpy())]
        mid = max(len(order) // 2, 1)
        prox = np.asarray(M[order[:mid]].sum(axis=0)).ravel()
        dist = np.asarray(M[order[mid:]].sum(axis=0)).ravel()
        tot = prox + dist
        usage = np.where(tot >= MIN_PARENT, dist / np.where(tot == 0, 1, tot), np.nan)
        gene_rows[gene] = usage
        # donor-level read-pooled usage (spots are technical subsamples):
        # sum proximal/distal reads over all spots -> single ratio per gene
        P, D = float(prox.sum()), float(dist.sum())
        pooled_rows[gene] = D / (P + D) if (P + D) >= MIN_POOLED_READS else np.nan

    gi = pd.DataFrame(gene_rows, index=sids).T  # gene x spot
    pooled = pd.Series(pooled_rows)             # gene -> pooled usage
    obs = np.isfinite(gi.values).sum(1)
    kept = gi.index[obs >= MIN_OBS_SPOTS]
    pseudo = pooled.reindex(pooled.index)       # all genes with pooled value
    qc = dict(srr=srr, n_peaks=len(pids), n_spots=len(sids),
              n_peaks_annotated=int(n_anno), n_genes_usable=int(pooled.notna().sum()),
              n_genes_spotmedian=len(kept),
              median_genes_per_spot=float(np.isfinite(gi.values).mean()))
    return pseudo, qc
def mouse_to_human_symbol(sym, human_set=None):
    """Quick ortholog mapping: mouse symbol -> human symbol.

    Rule: human orthologs of mouse Title-case symbols are the ALL-CAPS form
    (Cdk8->CDK8); symbols that keep mixed case in human (C11orf58) are tried
    Title-case. Validated against the GRCh38 gene set when provided; mouse-only
    ids (Riken/Gm/microRNA clusters/MHC) return None.
    """
    if not isinstance(sym, str) or not sym:
        return None
    if sym.endswith("Rik") or sym.startswith("Gm") or sym.startswith("BC0") \
            or sym.startswith("AI") or sym.startswith("AW") or sym.startswith("4930") \
            or sym.startswith("9130") or sym.startswith("Mir") or sym.startswith("Snrpn"):
        return None
    cands = []
    up = sym.upper()
    if up != sym:
        cands.append(up)
    cands.append(sym[0].upper() + sym[1:])  # Title-case form (C11orf58)
    if human_set is not None:
        for c in cands:
            if c in human_set:
                return c
        return None
    return cands[0]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)

def cliffs_delta(x, y):
    """Cliff's delta: P(x>y) - P(y>x)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return (gt - lt) / (len(x) * len(y))


# ------------------------------------------------------------------- main
def main():
    meta = pd.read_csv(OUT / "sample_metadata.csv")
    primary = meta[meta.in_primary_cohort == 1].copy()
    print(f"[main] primary cohort: {len(primary)} donors\n{primary.groupby('diagnosis').size()}", flush=True)

    genes = read_genes_from_gtf(HUMAN_GTF)
    human_set = {r.gene_name for recs in genes.values() for r in recs}

    series, qcs = {}, []
    for srr in primary["srr"]:
        if not (SR_ROOT / srr / "apa_site_counts.csv.gz").exists():
            print(f"  [{srr}] MISSING scAPAtrap output -- skipping", flush=True)
            continue
        s, qc = load_sample(srr, genes)
        series[srr] = s
        qcs.append(qc)
    if len(series) < 8:
        have_dx = primary.set_index("srr").loc[list(series), "diagnosis"]
        ok = ("AD" in set(have_dx)) and ("Control" in set(have_dx)) and len(series) >= 6
        if not ok:
            sys.exit(f"only {len(series)} samples processed; need >=6 covering AD+Control")
        print(f"[main] WARNING partial cohort n={len(series)}", flush=True)
    usage = pd.DataFrame(series)                       # gene x sample
    usage.index.name = "gene"
    usage.to_csv(OUT / "per_sample_gene_usage.csv")
    donor = primary.set_index("srr")["donor_id"].to_dict()
    dx = primary.set_index("srr")["diagnosis"].to_dict()
    donor_dx = {donor[s]: d for s, d in dx.items()}
    usage = usage.rename(columns=donor)                # gene x donor
    print(f"[main] usage matrix {usage.shape}; missingness "
          f"{usage.isna().mean().mean():.3f}", flush=True)

    contrasts = [
        ("AD_vs_Control", "AD", "Control"),
        ("earlyAD_vs_Control", "earlyAD", "Control"),
        ("AD_DS_vs_Control", "AD_DS", "Control"),
        ("all_disease_vs_Control", ["AD", "earlyAD", "AD_DS"], "Control"),
    ]
    cols = usage.columns
    results = {}
    for name, g1_dx, g0_dx in contrasts:
        g1 = [c for c in cols if donor_dx[c] in
              (([g1_dx] if isinstance(g1_dx, str) else g1_dx))]
        g0 = [c for c in cols if donor_dx[c] == g0_dx]
        if len(g1) < 2 or len(g0) < 2:
            print(f"[{name}] skipped (n_disease={len(g1)}, n_control={len(g0)})", flush=True)
            continue
        rows = []
        sub = usage[g1 + g0]
        for gene, r in sub.iterrows():
            x, y = r[g1].to_numpy(float), r[g0].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            x, y = x[ok], y[ok]
            if len(x) < 3 or len(y) < 3:
                continue
            if np.allclose(np.concatenate([x, y]), x[0]):
                p = 1.0
            else:
                p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
            rows.append(dict(
                gene=gene, n_disease=len(x), n_control=len(y),
                median_disease=float(np.median(x)), median_control=float(np.median(y)),
                delta_disease_minus_control=float(np.median(x) - np.median(y)),
                cliffs_delta=cliffs_delta(x, y), p_value=p))
        res = pd.DataFrame(rows).sort_values("p_value")
        res["fdr_BH"] = bh_fdr(res["p_value"].to_numpy())
        # bootstrap 95% CI of delta for the top signals (p < 0.05)
        lo, hi = [], []
        for _, row in res.iterrows():
            if row.p_value >= 0.05:
                lo.append(np.nan); hi.append(np.nan); continue
            x = usage.loc[row.gene, g1].dropna().to_numpy(float)
            y = usage.loc[row.gene, g0].dropna().to_numpy(float)
            boots = []
            for _ in range(N_BOOT):
                xb = x[RNG.integers(0, len(x), len(x))]
                yb = y[RNG.integers(0, len(y), len(y))]
                boots.append(np.median(xb) - np.median(yb))
            lo.append(float(np.percentile(boots, 2.5)))
            hi.append(float(np.percentile(boots, 97.5)))
        res["delta_ci95_lo"], res["delta_ci95_hi"] = lo, hi
        res.to_csv(OUT / f"donor_level_differential_apa_{name}.csv", index=False)
        results[name] = res
        print(f"[{name}] {len(res)} genes tested; "
              f"FDR<0.25: {(res.fdr_BH < 0.25).sum()}, "
              f"p<0.05: {(res.p_value < 0.05).sum()}", flush=True)

    # primary contrast gets the canonical filename too
    results["AD_vs_Control"].to_csv(OUT / "donor_level_differential_apa.csv", index=False)

    # ------------------------------------------------ consistency (mouse AD)
    cons_rows, mouse = None, None
    if MOUSE_REF.exists():
        mouse = pd.read_csv(MOUSE_REF)
        hum = results["AD_vs_Control"].set_index("gene")
        pairs = []
        for _, m in mouse.iterrows():
            h = mouse_to_human_symbol(m.gene, human_set)
            if h is None or h not in hum.index:
                continue
            hd = hum.loc[h]
            if not np.isfinite(hd.delta_disease_minus_control):
                continue
            pairs.append(dict(mouse_gene=m.gene, human_gene=h,
                              mouse_delta_AD_minus_WT=m.delta_AD_minus_WT,
                              human_delta_AD_minus_Control=hd.delta_disease_minus_control,
                              human_p=hd.p_value, human_fdr=hd.fdr_BH,
                              mouse_padj=m.padj))
        cons = pd.DataFrame(pairs)
        cons["same_direction"] = np.sign(cons.mouse_delta_AD_minus_WT) == \
            np.sign(cons.human_delta_AD_minus_Control)
        n_same = int(cons.same_direction.sum())
        binom_p = stats.binomtest(n_same, len(cons), 0.5).pvalue
        rho, rho_p = stats.spearmanr(cons.mouse_delta_AD_minus_WT,
                                     cons.human_delta_AD_minus_Control)
        # mouse spot-level-significant subset
        sig = cons[cons.mouse_padj < 0.05]
        n_same_sig = int(sig.same_direction.sum())
        binom_p_sig = stats.binomtest(n_same_sig, len(sig), 0.5).pvalue
        cons.to_csv(OUT / "consistency_with_gse263789.csv", index=False)
        cons_rows = dict(
            n_shared_genes=len(cons), n_same_direction=n_same,
            binom_p_direction=binom_p,
            spearman_rho=float(rho), spearman_p=float(rho_p),
            n_mouse_spotFDR=len(sig), n_same_direction_mouse_spotFDR=n_same_sig,
            binom_p_direction_mouse_spotFDR=binom_p_sig)
        print(f"[consistency] shared={len(cons)} same-dir={n_same} "
              f"(p={binom_p:.2g}); rho={rho:.3f} (p={rho_p:.2g}); "
              f"mouse-spotFDR subset: {n_same_sig}/{len(sig)}", flush=True)


    # --------------------------------------- published AD APA gene panel
    known_rows = None
    KNOWN = ROOT / "pipeline_output/known_gene_validation/ad_apa_genes.csv"
    if KNOWN.exists():
        known = pd.read_csv(KNOWN)
        ad_res = results["AD_vs_Control"].set_index("gene")
        mouse_syms = set()
        if MOUSE_REF.exists():
            mouse_syms = {mouse_to_human_symbol(g, human_set) for g in mouse.gene}
        ktab = []
        for gene in sorted(known.gene.unique()):
            if gene in ad_res.index:
                r = ad_res.loc[gene]
                ktab.append(dict(gene=gene,
                                 delta_AD_minus_Control=r.delta_disease_minus_control,
                                 p_value=r.p_value, fdr_BH=r.fdr_BH,
                                 present_in_mouse_ref=gene in mouse_syms))
            else:
                ktab.append(dict(gene=gene, delta_AD_minus_Control=np.nan,
                                 p_value=np.nan, fdr_BH=np.nan,
                                 present_in_mouse_ref=gene in mouse_syms))
        ktab = pd.DataFrame(ktab)
        ktab.to_csv(OUT / "known_ad_apa_genes_in_donor_level_results.csv", index=False)
        in_res = ktab.dropna(subset=["p_value"])
        hit = in_res[in_res.p_value < 0.05]
        # Fisher exact: p<0.05 enrichment within known panel vs all tested genes
        a = len(hit); b = len(in_res) - a
        c = int((ad_res.p_value < 0.05).sum()) - a
        d = len(ad_res) - a - b - c
        fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")[1]
        known_rows = dict(
            n_known_genes=len(ktab), n_tested_in_our_data=len(in_res),
            n_p_lt_05=len(hit), fisher_enrichment_p=float(fisher_p),
            hits=hit[["gene", "delta_AD_minus_Control", "p_value", "fdr_BH"]].to_dict("records"))
        print(f"[known-panel] {len(in_res)}/{len(ktab)} tested; "
              f"p<0.05: {len(hit)} (Fisher p={fisher_p:.2g})", flush=True)
    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.2))
    for ax, (name, res) in zip(axes, results.items()):
        res = res.dropna(subset=["delta_disease_minus_control"])
        x = res.delta_disease_minus_control
        y = -np.log10(np.clip(res.p_value, 1e-12, 1))
        sig = res.fdr_BH < 0.25
        ax.scatter(x[~sig], y[~sig], s=8, c="grey", alpha=.45, lw=0)
        up = sig & (x > 0)
        dn = sig & (x < 0)
        ax.scatter(x[up], y[up], s=14, c="#d43f3f", lw=0, label=f"distal-up FDR<0.25 ({up.sum()})")
        ax.scatter(x[dn], y[dn], s=14, c="#3f6fd4", lw=0, label=f"proximal-up FDR<0.25 ({dn.sum()})")
        for _, r in res[res.fdr_BH < 0.10].nsmallest(8, "fdr_BH").iterrows():
            ax.annotate(r.gene, (r.delta_disease_minus_control,
                                 -np.log10(max(r.p_value, 1e-12))),
                        fontsize=7, alpha=.9)
        ax.set_title(f"{name}\n(donor-level Wilcoxon)", fontsize=10)
        ax.set_xlabel("Δ median distal usage (disease − control)")
        ax.set_ylabel("−log10 p")
        ax.legend(fontsize=7, loc="upper left")
        ax.axhline(-np.log10(0.05), lw=.6, ls="--", c="k", alpha=.4)
        ax.axvline(0, lw=.6, c="k", alpha=.3)
    fig.suptitle("GSE233208 human PFC Visium — donor-level differential APA "
                 "(16 donors; unit = donor pseudobulk, not spot)", y=1.03, fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / "fig_donor_level_volcano.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------------------- summary
    def contrast_summary(name):
        res = results[name]
        return dict(
            n_genes_tested=int(len(res)),
            n_p_lt_05=int((res.p_value < 0.05).sum()),
            n_fdr_lt_25=int((res.fdr_BH < 0.25).sum()),
            n_fdr_lt_10=int((res.fdr_BH < 0.10).sum()),
            top_genes=res.nsmallest(10, "p_value")[
                ["gene", "delta_disease_minus_control", "cliffs_delta",
                 "p_value", "fdr_BH"]].to_dict("records"))

    summary = dict(
        experiment="gse233208_human_donor_level_differential_apa",
        date=time.strftime("%Y-%m-%d %H:%M:%S"),
        purpose=("External-review Phase 2: replace spot-level pseudoreplicated "
                 "AD APA evidence with donor-level tests in an independent "
                 "human cohort with biological replicates"),
        cohort=dict(source="GSE233208 (Miyoshi et al., Nat Genet 2024)",
                    platform="Visium fresh-frozen 55um, PFC",
                    n_donors=int(len(primary)),
                    diagnoses=primary.diagnosis.value_counts().to_dict(),
                    unit_of_analysis="donor (per-sample median spot distal-usage)"),
        methods=dict(
            pas_caller="scAPAtrap (TenX, peaks, cov.cutoff=10, min.cells=10)",
            gene_index="median-split proximal/distal, MIN_PARENT=5",
            pseudobulk="per-donor read-pooled distal usage: sum(distal)/sum(prox+dist) reads across spots (MIN_POOLED_READS=20); spots treated as technical subsamples",
            test="Wilcoxon rank-sum across donors (exact, two-sided)",
            effect_size="Cliff's delta + bootstrap 95% CI of median delta",
            multiple_testing="BH-FDR",
            note=("p-value floor for 4v4 exact Wilcoxon = 0.0286; "
                  "FDR<0.25 is the donor-level signal tier given n=4/group")),
        contrasts={k: contrast_summary(k) for k in results},
        consistency_with_gse263789_mouse_ad=cons_rows,
        published_ad_apa_gene_panel=known_rows,
    )
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print("[main] wrote summary.json + all outputs to", OUT, flush=True)


if __name__ == "__main__":
    main()

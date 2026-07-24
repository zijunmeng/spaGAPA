#!/usr/bin/env python3
"""Differential APA on GSE220442 (3 control vs 3 AD) using a UNIFIED peak set.

Fixes the v1 caveat (per-sample independent scAPAtrap peak calling -> AD1 called
~39k peaks vs 24-27k for the others, confounding distal-usage shifts). Instead
we:

  1. Build a UNIFIED consensus peak set across all 6 samples (bedtools-merge
     style: same strand, intervals whose start <= prev_end + MERGE_GAP get
     merged). Each consensus peak inherits a gene_name by majority vote of its
     constituent sample peaks.
  2. Re-quantify every sample's counts on the unified peaks: each sample's
     original peaks are mapped onto consensus peaks (1:many sample->consensus
     handled by summing; many:1 is impossible by construction since overlapping
     sample peaks collapse into one consensus peak -> their counts sum).
  3. Per sample, build the gene-level distal-usage index on the unified peaks
     (median-split proximal/distal peak groups per gene; NaN where parent total
     < MIN_PARENT; drop genes with < MIN_OBS_SPOTS observed spots).
  4. Differential APA (pooled spots) via spaGAPA's DifferentialAPAAnalyzer,
     running BOTH t-test and Wilcoxon with BH-FDR. Plus sample-level
     replication (3v3 direction agreement: min(AD) > max(control) or vice
     versa).
  5. Compare to v1 (per-sample-peak) gene set and magnitudes; check whether v1's
     extreme shifts (DDX24 0.18->0.97 etc.) survive on the unified peak set.
  6. (secondary) miRNA target gain/loss annotation for top differential genes,
     using the local miRanda-annotated DB. We proxy 3'UTR lengthening/shortening
     by the APA direction: distal-shift (delta>0) == distal PAS used more ==
     longer 3'UTR (more miRNA target sites exposed); proximal-shift == shorter.

Usage: python scripts/analyze_gse220442_diff_apa_unified.py
Output: pipeline_output/gse220442_differential_apa_unified/
"""
from __future__ import annotations
import json, time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
PROCESSED = PKG / "data/processed"
OUT = PKG / "pipeline_output/gse220442_differential_apa_unified"
OUT.mkdir(parents=True, exist_ok=True)

# control = GSM6801751/2/3 ; AD = GSM6801754/5/6
CONTROL = ["gsm6801751", "gsm6801752", "gsm6801753"]
AD = ["gsm6801754", "gsm6801755", "gsm6801756"]
ALL_SAMPLES = CONTROL + AD
COND_OF = {**{s: "control" for s in CONTROL}, **{s: "AD" for s in AD}}

MERGE_GAP = 50          # bp: merge adjacent peaks whose start <= prev_end + gap
MIN_PARENT = 5          # min total counts across a gene's peaks for a spot
MIN_OBS_SPOTS = 10      # drop genes observed in fewer spots in a sample
PADJ_THRESH = 0.05
DELTA_THRESH = 0.05

MIRNA_DB = Path("/s1/SHARE/apadata/MicroRNA/miRNA.mRNA.predict/human_miRNA/miRNA_results_annotated.csv")
V1_RESULTS = PKG / "pipeline_output/gse220442_differential_apa/differential_apa_results.csv"


# ---------------------------------------------------------------------------
# 1. UNIFIED PEAK SET
# ---------------------------------------------------------------------------
def load_sample_peaks(sid: str):
    """Return (sites_df, counts_df) for one sample.

    sites_df columns: site_id, gene_name, chr, start, end, strand, coord
    counts_df: peak x spot (index = site_id)
    """
    d = PROCESSED / f"gse220442_{sid}_scapatrap"
    sites = pd.read_csv(d / "apa_sites.csv.gz")
    counts = pd.read_csv(d / "apa_site_counts.csv.gz", index_col=0)
    counts.columns = counts.columns.astype(str)
    # prefix spot barcodes with sample id so pooled matrix has no collisions
    counts.columns = [f"{sid}:{c}" for c in counts.columns]
    sites = sites[sites["site_id"].isin(counts.index)].reset_index(drop=True)
    return sites, counts


def build_consensus_peaks(peaks_by_sample: dict) -> tuple[pd.DataFrame, dict]:
    """Merge all sample peaks (chr,start,end,strand) bedtools-style.

    Returns:
      consensus : DataFrame[consensus_id, chr, start, end, strand, gene_name,
                           n_samples, n_constituent]
      membership : {consensus_id: [(sid, site_id), ...]}  for traceability
    """
    all_peaks = []
    for sid, sites in peaks_by_sample.items():
        sub = sites[["site_id", "gene_name", "chr", "start", "end",
                     "strand"]].copy()
        sub["sid"] = sid
        all_peaks.append(sub)
    allp = pd.concat(all_peaks, ignore_index=True)
    n_raw = len(allp)

    rows = []
    membership = {}
    cid = 0
    # merge within each (chr, strand): sort by start, greedy merge when
    # start <= running_end + MERGE_GAP
    for (chrom, strand), grp in allp.groupby(["chr", "strand"]):
        grp = grp.sort_values("start").reset_index(drop=True)
        starts = grp["start"].values
        ends = grp["end"].values
        n = len(grp)
        i = 0
        while i < n:
            j = i
            run_end = ends[i]
            # extend while next start within MERGE_GAP of running end
            while j + 1 < n and starts[j + 1] <= run_end + MERGE_GAP:
                j += 1
                run_end = max(run_end, ends[j])
            block = grp.iloc[i:j + 1]
            # gene_name by majority vote (drop NaN just in case)
            genes = [g for g in block["gene_name"].tolist() if isinstance(g, str)]
            gene = Counter(genes).most_common(1)[0][0] if genes else np.nan
            sids = block["sid"].tolist()
            cid_str = f"cpeak_{cid}"
            rows.append({
                "consensus_id": cid_str,
                "chr": chrom,
                "start": int(starts[i]),
                "end": int(run_end),
                "strand": strand,
                "gene_name": gene,
                "n_samples": len(set(sids)),
                "n_constituent": len(block),
            })
            membership[cid_str] = list(zip(block["sid"].tolist(),
                                           block["site_id"].tolist()))
            cid += 1
            i = j + 1
    consensus = pd.DataFrame(rows)
    print(f"  merged {n_raw} raw peaks -> {len(consensus)} consensus peaks "
          f"({consensus['n_constituent'].sum()} constituents; "
          f"median {consensus['n_constituent'].median():.0f}/cpeak, "
          f"max {consensus['n_constituent'].max()})")
    multi_sample = (consensus["n_samples"] >= 2).sum()
    all6 = (consensus["n_samples"] == 6).sum()
    print(f"  consensus peaks present in >=2 samples: {multi_sample} ; "
          f"in all 6: {all6}")
    return consensus, membership


def requantify(peaks_by_sample, counts_by_sample, membership):
    """Build per-sample count vectors on consensus peaks.

    Each sample's original peaks map 1:1 to a consensus peak (by construction a
    sample peak belongs to exactly one consensus peak). Sample peaks that
    collapse into the same consensus peak sum. Consensus peaks a sample did not
    call get 0.
    """
    # invert membership: (sid, site_id) -> consensus_id
    sample_peak_to_cons = {}
    for cons_id, members in membership.items():
        for sid, site_id in members:
            sample_peak_to_cons[(sid, site_id)] = cons_id

    consensus_ids = list(membership.keys())
    requant = {sid: pd.DataFrame(0, index=consensus_ids,
                                 columns=counts_by_sample[sid].columns,
                                 dtype=np.float64)
               for sid in ALL_SAMPLES}

    for sid in ALL_SAMPLES:
        sites = peaks_by_sample[sid]
        counts = counts_by_sample[sid]
        # map each sample peak row -> consensus id
        site2cons = {row["site_id"]: sample_peak_to_cons.get((sid, row["site_id"]))
                     for _, row in sites.iterrows()}
        for site_id, cons_id in site2cons.items():
            if cons_id is None or site_id not in counts.index:
                continue
            requant[sid].loc[cons_id] = requant[sid].loc[cons_id].add(
                counts.loc[site_id].values, axis=0)
    return requant  # {sid: consensus_peak x spot DataFrame}


# ---------------------------------------------------------------------------
# 2. GENE-LEVEL DISTAL-USAGE INDEX ON UNIFIED PEAKS
# ---------------------------------------------------------------------------
def build_gene_index_unified(consensus: pd.DataFrame,
                             requant: dict) -> dict:
    """For each sample, compute gene x spot distal-usage on unified peaks.

    Mirrors v1 build_gene_index: peaks of each gene oriented by strand, sorted,
    median-split into proximal/distal groups; distal_sum/(distal+proximal);
    NaN where total < MIN_PARENT.
    """
    consensus = consensus.copy()
    consensus["oriented"] = np.where(
        consensus["strand"].astype(str) == "-",
        -consensus["end"].astype(float),    # use peak end as position
        consensus["end"].astype(float),
    )
    # keep only peaks with a gene_name
    consensus = consensus[consensus["gene_name"].notna()].reset_index(drop=True)

    indices = {}
    for sid in ALL_SAMPLES:
        counts = requant[sid]
        # restrict to peaks that belong to this gene map and are in counts index
        cons_in_counts = consensus[consensus["consensus_id"].isin(counts.index)]
        peak_idx = {p: i for i, p in enumerate(counts.index)}
        vals = counts.values
        rows = {}
        for gene, sub in cons_in_counts.groupby("gene_name"):
            ridx = [peak_idx[p] for p in sub["consensus_id"] if p in peak_idx]
            if len(ridx) < 2:
                continue
            oriented = sub.set_index("consensus_id").loc[
                [counts.index[r] for r in ridx], "oriented"].values
            ridx = np.array(ridx)[np.argsort(oriented)]
            mid = max(1, len(ridx) // 2)
            prox = vals[ridx[:mid]].sum(0)
            dist = vals[ridx[mid:]].sum(0)
            total = prox + dist
            rows[gene] = np.where(total >= MIN_PARENT,
                                  dist / np.where(total == 0, 1, total),
                                  np.nan)
        idx = pd.DataFrame(rows, index=counts.columns).T  # gene x spot
        obs = np.isfinite(idx.values).sum(1)
        idx = idx.loc[idx.index[obs >= MIN_OBS_SPOTS]]
        indices[sid] = idx
        print(f"    {sid}: {idx.shape[0]} genes ({MIN_OBS_SPOTS}+ obs spots)")
    return indices


# ---------------------------------------------------------------------------
# 3. DIFFERENTIAL APA (pooled + sample-level)
# ---------------------------------------------------------------------------
def pooled_test(sample_indices: dict, method: str, common_genes: list):
    from spagapa.analysis import DifferentialAPAAnalyzer
    mats, conds = [], []
    for sid in ALL_SAMPLES:
        m = sample_indices[sid].loc[common_genes]
        mats.append(m)
        conds.append(pd.Series(COND_OF[sid], index=m.columns))
    pooled = pd.concat(mats, axis=1)
    cond = pd.concat(conds)
    g1 = np.where(cond.values == "control")[0]
    g2 = np.where(cond.values == "AD")[0]
    ana = DifferentialAPAAnalyzer(method=method, min_spots_per_group=10)
    res = ana.test_differential_apa(pooled.values, g1, g2,
                                    gene_names=common_genes)
    res = ana.adjust_pvalues(res, method="fdr_bh")
    res = res.rename(columns={"mean_group1": "mean_control",
                              "mean_group2": "mean_AD"})
    res["delta_AD_minus_control"] = res["mean_AD"] - res["mean_control"]
    return res, len(g1), len(g2)


def sample_level(sample_indices: dict, common_genes: list) -> pd.DataFrame:
    samp_vals = {}
    for sid in ALL_SAMPLES:
        m = sample_indices[sid].loc[common_genes]
        grp = COND_OF[sid]
        for g in common_genes:
            v = m.loc[g].values
            v = v[np.isfinite(v)]
            samp_vals.setdefault(g, {"control": [], "AD": []})[grp].append(
                float(np.mean(v)) if len(v) else np.nan)
    rows = []
    for g in common_genes:
        c = samp_vals[g]["control"]
        a = samp_vals[g]["AD"]
        if any(np.isnan(c)) or any(np.isnan(a)):
            rows.append({"gene": g, "sample_delta": np.nan,
                         "direction_agree": np.nan})
            continue
        delta = float(np.mean(a) - np.mean(c))
        agree = int((min(a) > max(c)) or (max(a) < min(c)))
        rows.append({"gene": g, "sample_delta": delta,
                     "direction_agree": agree,
                     "control_vals": ";".join(f"{x:.3f}" for x in c),
                     "ad_vals": ";".join(f"{x:.3f}" for x in a)})
    return pd.DataFrame(rows).set_index("gene")


# ---------------------------------------------------------------------------
# 4. miRNA TARGET ANNOTATION (secondary)
# ---------------------------------------------------------------------------
def load_mirna_targets(top_genes: list):
    """Return {gene: set(mirna)} from local miRanda-annotated DB.

    DB row: miRNA, Target, Score, ..., Gene_Symbol (may contain 'G1;G2').
    Score >= 140 is the conventional miRanda high-confidence cutoff.
    """
    if not MIRNA_DB.exists():
        return None, None
    top_set = set(top_genes)
    # read only needed columns; Gene_Symbol may hold 'G1;G2'
    cols = ["miRNA", "Score", "Gene_Symbol"]
    db = pd.read_csv(MIRNA_DB, usecols=cols, low_memory=False)
    db["Score"] = pd.to_numeric(db["Score"], errors="coerce")
    db = db.dropna(subset=["Gene_Symbol", "Score"])
    db = db[db["Score"] >= 140]   # high-confidence miRanda cutoff
    # explode multi-gene rows
    db = db.assign(Gene_Symbol=db["Gene_Symbol"].astype(str).str.split(";"))
    db = db.explode("Gene_Symbol")
    db["Gene_Symbol"] = db["Gene_Symbol"].str.strip()
    db = db[db["Gene_Symbol"].isin(top_set)]
    targets = {}
    for gene, sub in db.groupby("Gene_Symbol"):
        targets[gene] = sorted(sub["miRNA"].unique().tolist())
    # background: total miRNAs per gene in DB (for the gene universe we tested)
    return targets, db


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=== GSE220442 unified-peak differential APA: 3 control vs 3 AD ===")

    # ---- 1. load all sample peaks ----
    print("\n[1] loading per-sample peaks/counts ...")
    t1 = time.time()
    peaks_by_sample, counts_by_sample = {}, {}
    for sid in ALL_SAMPLES:
        sites, counts = load_sample_peaks(sid)
        peaks_by_sample[sid] = sites
        counts_by_sample[sid] = counts
        print(f"  {sid} ({COND_OF[sid]}): {len(sites)} peaks x "
              f"{counts.shape[1]} spots")
    print(f"  load: {time.time()-t1:.1f}s")

    # ---- 2. build unified consensus peak set ----
    print("\n[2] building unified consensus peak set "
          f"(merge_gap={MERGE_GAP}bp, same strand) ...")
    t1 = time.time()
    consensus, membership = build_consensus_peaks(peaks_by_sample)
    print(f"  merge: {time.time()-t1:.1f}s")
    consensus.to_csv(OUT / "consensus_peaks.csv", index=False)
    summary_peak = {
        "n_raw_peaks_total": int(sum(len(p) for p in peaks_by_sample.values())),
        "n_consensus_peaks": int(len(consensus)),
        "n_per_sample_raw": {sid: int(len(p))
                             for sid, p in peaks_by_sample.items()},
        "merge_gap_bp": MERGE_GAP,
        "n_multi_sample": int((consensus["n_samples"] >= 2).sum()),
        "n_all6_samples": int((consensus["n_samples"] == 6).sum()),
        "n_genes_in_consensus": int(consensus["gene_name"].nunique()),
    }

    # ---- 3. re-quantify each sample on unified peaks ----
    print("\n[3] re-quantifying each sample on unified peaks ...")
    t1 = time.time()
    requant = requantify(peaks_by_sample, counts_by_sample, membership)
    # sanity: report total counts preserved per sample
    for sid in ALL_SAMPLES:
        orig = float(counts_by_sample[sid].values.sum())
        new = float(requant[sid].values.sum())
        print(f"  {sid}: orig counts={orig:.0f}  unified={new:.0f}  "
              f"ratio={new/orig:.4f}")
    print(f"  requantify: {time.time()-t1:.1f}s")

    # ---- 4. per-sample gene-level index on unified peaks ----
    print("\n[4] building gene-level distal-usage index per sample "
          "(unified peaks) ...")
    t1 = time.time()
    sample_indices = build_gene_index_unified(consensus, requant)
    print(f"  index build: {time.time()-t1:.1f}s")

    common_genes = set(sample_indices[ALL_SAMPLES[0]].index)
    for sid in ALL_SAMPLES[1:]:
        common_genes &= set(sample_indices[sid].index)
    common_genes = sorted(common_genes)
    print(f"  genes present in ALL 6 samples (unified): {len(common_genes)}")

    # ---- 5. differential APA: pooled, BOTH t-test and Wilcoxon ----
    print("\n[5a] pooled differential APA: t-test ...")
    t1 = time.time()
    res_t, nc, na = pooled_test(sample_indices, "t-test", common_genes)
    res_t = res_t.set_index("gene")
    print(f"  pooled matrix: {len(common_genes)} genes x {nc+na} spots "
          f"(control {nc} / AD {na})  [{time.time()-t1:.1f}s]")

    print("[5b] pooled differential APA: Wilcoxon ...")
    t1 = time.time()
    res_w, _, _ = pooled_test(sample_indices, "wilcoxon", common_genes)
    res_w = res_w.set_index("gene")
    print(f"  [{time.time()-t1:.1f}s]")

    print("[5c] sample-level replication ...")
    samp = sample_level(sample_indices, common_genes)

    # combine
    res = res_t.join(
        res_w[["pvalue", "padj"]].rename(columns={"pvalue": "pvalue_wilcox",
                                                  "padj": "padj_wilcox"}))
    res = res.rename(columns={"pvalue": "pvalue_ttest",
                              "padj": "padj_ttest"})
    res = res.join(samp)

    sig_t = res[(res["padj_ttest"] < PADJ_THRESH)
                & (res["delta_AD_minus_control"].abs() > DELTA_THRESH)]
    sig_w = res[(res["padj_wilcox"] < PADJ_THRESH)
                & (res["delta_AD_minus_control"].abs() > DELTA_THRESH)]
    sig_both = res.loc[sig_t.index.intersection(sig_w.index)]
    replicated_t = sig_t[sig_t["direction_agree"] == 1]
    replicated_w = sig_w[sig_w["direction_agree"] == 1]
    replicated_both = sig_both[sig_both["direction_agree"] == 1]

    res_out = res.reset_index()
    res_out.to_csv(OUT / "results.csv", index=False)
    print(f"\n  -> results written ({len(res)} genes)")

    print("\n=== differential APA results (unified peaks) ===")
    print(f"  genes tested (pooled):           {len(res)}")
    print(f"  sig t-test  (padj<.05 & |d|>.05): {len(sig_t)} "
          f"(prox {(sig_t['delta_AD_minus_control']<0).sum()}, "
          f"dist {(sig_t['delta_AD_minus_control']>0).sum()})")
    print(f"  sig Wilcoxon(padj<.05 & |d|>.05): {len(sig_w)} "
          f"(prox {(sig_w['delta_AD_minus_control']<0).sum()}, "
          f"dist {(sig_w['delta_AD_minus_control']>0).sum()})")
    print(f"  sig BOTH methods:                {len(sig_both)}")
    print(f"  + replicated (t-test, dir agree): {len(replicated_t)}")
    print(f"  + replicated (Wilcoxon):          {len(replicated_w)}")
    print(f"  + replicated (BOTH):              {len(replicated_both)}")

    print("\n  top 15 replicated (BOTH) AD-dysregulated APA genes:")
    if len(replicated_both):
        show = replicated_both.reindex(
            replicated_both["delta_AD_minus_control"].abs()
            .sort_values(ascending=False).index).head(15)
        print(show[["mean_control", "mean_AD", "delta_AD_minus_control",
                    "padj_ttest", "padj_wilcox", "sample_delta",
                    "direction_agree", "control_vals", "ad_vals"]].to_string())

    # ---- 6. comparison to v1 ----
    print("\n[6] comparison to v1 (per-sample-peak) ...")
    v1_top = ["DDX24", "TPD52", "ABHD12", "RPL14", "BAIAP2", "NCDN"]
    v1_survive = {}
    if V1_RESULTS.exists():
        v1 = pd.read_csv(V1_RESULTS)
        v1_map = {}
        for _, r in v1.iterrows():
            v1_map[r["gene"]] = r
        print(f"  v1 tested {len(v1)} genes, unified tests {len(res)}")
        print(f"  v1 sig (pooled): "
              f"{(v1['padj']<PADJ_THRESH).sum() if 'padj' in v1 else 'n/a'}, "
              f"v1 replicated: "
              f"{(v1['direction_agree']==1).sum() if 'direction_agree' in v1 else 'n/a'}")
        print("\n  do v1's top genes survive on the unified peak set?")
        print(f"  {'gene':8} {'v1_ctrl':>8} {'v1_AD':>8} {'v1_d':>7} | "
              f"{'un_ctrl':>8} {'un_AD':>8} {'un_d':>7} {'padj_t':>9} "
              f"{'padj_w':>9} {'diragree':>8}  verdict")
        for g in v1_top:
            if g in v1_map and g in res.index:
                v1r = v1_map[g]
                ur = res.loc[g]
                # verdict: same direction & still significant + replicated?
                v1_dir = np.sign(v1r.get("delta_AD_minus_control", np.nan))
                un_dir = np.sign(ur["delta_AD_minus_control"])
                same_dir = (v1_dir == un_dir) and not np.isnan(v1_dir)
                still_sig = (ur["padj_ttest"] < PADJ_THRESH) or \
                            (ur["padj_wilcox"] < PADJ_THRESH)
                still_rep = ur["direction_agree"] == 1
                if same_dir and still_rep:
                    verdict = "SURVIVES (replicated)"
                elif same_dir and still_sig:
                    verdict = "survives (sig, not rep)"
                elif same_dir:
                    verdict = "direction only"
                else:
                    verdict = "DROPS / flipped"
                v1_survive[g] = verdict
                print(f"  {g:8} {v1r.get('mean_control',float('nan')):8.2f} "
                      f"{v1r.get('mean_AD',float('nan')):8.2f} "
                      f"{v1r.get('delta_AD_minus_control',float('nan')):+7.2f} | "
                      f"{ur['mean_control']:8.2f} {ur['mean_AD']:8.2f} "
                      f"{ur['delta_AD_minus_control']:+7.2f} "
                      f"{ur['padj_ttest']:9.2e} {ur['padj_wilcox']:9.2e} "
                      f"{int(ur['direction_agree']) if not np.isnan(ur['direction_agree']) else -1:>8}  {verdict}")
            else:
                v1_survive[g] = "absent in unified"
                print(f"  {g:8}  -- absent from unified common genes")
    else:
        print(f"  (v1 results not found at {V1_RESULTS})")

    # ---- 7. miRNA target annotation (secondary) ----
    print("\n[7] miRNA target gain/loss annotation (secondary) ...")
    mirna_report = {"db_path": str(MIRNA_DB), "status": "skipped"}
    if MIRNA_DB.exists():
        top_genes_for_mirna = replicated_both.index.tolist()[:30] \
            if len(replicated_both) else \
            sig_both.index.tolist()[:30] if len(sig_both) else \
            sig_t.index.tolist()[:30]
        targets, db_sub = load_mirna_targets(top_genes_for_mirna)
        n_annotated = sum(1 for g in top_genes_for_mirna if g in (targets or {}))
        print(f"  DB: {MIRNA_DB.name} ({MIRNA_DB.stat().st_size//1_000_000}MB)")
        print(f"  annotated {n_annotated}/{len(top_genes_for_mirna)} top genes "
              f"with >=1 miRanda target (score>=140)")
        rows = []
        for g in top_genes_for_mirna:
            if g not in res.index:
                continue
            r = res.loc[g]
            d = r["delta_AD_minus_control"]
            # distal-shift (d>0) -> longer 3'UTR -> potential target GAIN
            # proximal-shift (d<0) -> shorter 3'UTR -> potential target LOSS
            shift = "distal-up (3'UTR lengthen)" if d > 0 else \
                    "proximal-up (3'UTR shorten)" if d < 0 else "none"
            gain_loss = "miRNA target GAIN (longer 3'UTR)" if d > 0 else \
                        "miRNA target LOSS (shorter 3'UTR)" if d < 0 else "-"
            mirs = targets.get(g, []) if targets else []
            rows.append({
                "gene": g,
                "delta_AD_minus_control": d,
                "direction": shift,
                "mirna_gain_loss_prediction": gain_loss,
                "n_mirna_targets": len(mirs),
                "top_mirna_targets": ";".join(mirs[:10]),
            })
        mirna_df = pd.DataFrame(rows)
        if len(mirna_df):
            mirna_df.to_csv(OUT / "mirna_target_annotation.csv", index=False)
            print(mirna_df[["gene", "delta_AD_minus_control", "direction",
                            "mirna_gain_loss_prediction",
                            "n_mirna_targets"]].head(20).to_string(index=False))
        mirna_report = {
            "db_path": str(MIRNA_DB),
            "status": "done",
            "n_genes_annotated": int(n_annotated),
            "n_top_genes": int(len(top_genes_for_mirna)),
        }
    else:
        print(f"  no miRNA DB at {MIRNA_DB} -> SKIPPED (remaining step)")

    # ---- volcano (t-test padj) ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 6))
        rep = res["direction_agree"] == 1
        sig = (res["padj_ttest"] < PADJ_THRESH) & \
              (res["delta_AD_minus_control"].abs() > DELTA_THRESH)
        ax.scatter(res.loc[~sig, "delta_AD_minus_control"],
                   -np.log10(res.loc[~sig, "padj_ttest"].clip(lower=1e-300)),
                   s=6, c="#bdc3c7", alpha=0.5, label="ns")
        ax.scatter(res.loc[sig & ~rep, "delta_AD_minus_control"],
                   -np.log10(res.loc[sig & ~rep, "padj_ttest"]
                             .clip(lower=1e-300)),
                   s=10, c="#e67e22", alpha=0.7, label="pooled sig (t-test)")
        ax.scatter(res.loc[sig & rep, "delta_AD_minus_control"],
                   -np.log10(res.loc[sig & rep, "padj_ttest"]
                             .clip(lower=1e-300)),
                   s=16, c="#c0392b", alpha=0.9,
                   label="pooled + replicated (3v3 agree)")
        ax.axhline(-np.log10(PADJ_THRESH), ls="--", c="grey", lw=0.8)
        ax.axvline(0, ls="--", c="grey", lw=0.8)
        ax.set_xlabel("delta distal-usage (AD - control)")
        ax.set_ylabel("-log10(padj_ttest)")
        ax.set_title("GSE220442 differential APA: AD vs control (3v3)\n"
                     "UNIFIED consensus peak set", fontweight="bold")
        ax.legend()
        top_lab = replicated_both if len(replicated_both) else sig_both
        if len(top_lab):
            for g, r in top_lab.reindex(
                    top_lab["delta_AD_minus_control"].abs()
                    .sort_values(ascending=False).index).head(10).iterrows():
                ax.annotate(g, (r["delta_AD_minus_control"],
                                -np.log10(max(r["padj_ttest"], 1e-300))),
                            fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / "volcano.png", dpi=150, bbox_inches="tight")
        print(f"\n  -> volcano.png")
    except Exception as e:
        print(f"  (volcano skipped: {e})")

    # ---- summary ----
    summary = {
        **summary_peak,
        "n_genes_common_all6": int(len(common_genes)),
        "n_spots_pooled": {"control": int(nc), "AD": int(na)},
        "n_sig_ttest": int(len(sig_t)),
        "n_sig_wilcox": int(len(sig_w)),
        "n_sig_both": int(len(sig_both)),
        "n_replicated_ttest": int(len(replicated_t)),
        "n_replicated_wilcox": int(len(replicated_w)),
        "n_replicated_both": int(len(replicated_both)),
        "proximal_shift_ttest": int((sig_t["delta_AD_minus_control"] < 0).sum()),
        "distal_shift_ttest": int((sig_t["delta_AD_minus_control"] > 0).sum()),
        "v1_top_gene_survival": v1_survive,
        "mirna_annotation": mirna_report,
        "params": {"merge_gap_bp": MERGE_GAP, "min_parent": MIN_PARENT,
                   "min_obs_spots": MIN_OBS_SPOTS,
                   "padj_thresh": PADJ_THRESH, "delta_thresh": DELTA_THRESH},
        "wall_s": round(time.time() - t0, 1),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n=== DONE ===")
    print(f"  -> {OUT}/  (results.csv, summary.json, consensus_peaks.csv, "
          f"volcano.png, mirna_target_annotation.csv)")
    print(f"  wall: {summary['wall_s']}s")


if __name__ == "__main__":
    main()

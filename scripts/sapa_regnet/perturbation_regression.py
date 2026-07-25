#!/usr/bin/env python3
"""Phase 1: regression-based interventional APA perturbation (CPU-only, go/no-go).

Replaces the hand-signed perturbation heuristic with a DATA-LEARNED model:
  log1p(expr_G, spot) = b0 + b_APA * distal_usage_G(spot) + domain fixed effects + eps
b_APA is learned per gene (the within-domain cis effect of G's own APA on G's
expression). Cell-type/domain confounding is controlled via domain fixed effects.

Then:
  - self-consistency: does b_APA sign match distal 3'UTR regulatory elements?
    (genes with more distal miRNAs -> repressive -> expect b_APA < 0)
  - this is the go/no-go: if APA predicts expression (significant b_APA,
    sign consistent with biology), the perturbation approach is viable.
"""
import os, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np
import pandas as pd

PROJ = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
SAMPLE = "gsm6801751"
D = f"{PROJ}/data/processed/gse220442_{SAMPLE}_scapatrap"
OUT = f"{PROJ}/pipeline_output/sapa_regnet/perturbation_regression"
os.makedirs(OUT, exist_ok=True)
MIN_SPOTS = 80  # per gene, with both APA and expression observed


def build_gene_index(counts, sites, min_parent=5):
    sites = sites.copy()
    sites["oriented"] = np.where(sites["strand"].astype(str) == "-",
                                 -sites["coord"].astype(float), sites["coord"].astype(float))
    peak_idx = {p: i for i, p in enumerate(counts.index)}
    vals = counts.values
    rows = {}
    for gene, sub in sites.groupby("gene_name"):
        sub = sub.dropna(subset=["gene_name"])
        ridx = [peak_idx[r["site_id"]] for _, r in sub.iterrows() if r["site_id"] in peak_idx]
        if len(ridx) < 2:
            continue
        oriented = np.array([float(sub.loc[(sub["site_id"] == counts.index[r]), "oriented"].iloc[0]) for r in ridx])
        ridx = np.array(ridx)[np.argsort(oriented)]
        mid = max(1, len(ridx) // 2)
        prox = vals[ridx[:mid]].sum(0); dist = vals[ridx[mid:]].sum(0)
        total = prox + dist
        rows[gene] = np.where(total >= min_parent, dist / np.where(total == 0, 1, total), np.nan)
    return pd.DataFrame(rows, index=counts.columns).T  # gene x spot


def main():
    t0 = time.time()
    print("loading APA + expression + domain ...", flush=True)
    counts = pd.read_csv(f"{D}/apa_site_counts.csv.gz", index_col=0); counts.columns = counts.columns.astype(str)
    sites = pd.read_csv(f"{D}/apa_sites.csv.gz")
    expr = pd.read_csv(f"{D}/expression_matrix.csv", index_col=0); expr.columns = expr.columns.astype(str)
    dom = pd.read_csv(f"{PROJ}/pipeline_output/gse220442_spagapa_runs/gse220442_{SAMPLE}_scapatrap/domains.csv")
    dom.columns = ["spot_id", "domain"] if list(dom.columns)[:2] == list(dom.columns) else dom.columns
    dom["spot_id"] = dom["spot_id"].astype(str)

    print("  building gene-level distal-usage ...", flush=True)
    du = build_gene_index(counts, sites)
    # align spots
    spots = sorted(set(du.columns) & set(expr.columns) & set(dom["spot_id"]))
    du = du[spots]; expr = expr[spots]; dom = dom.set_index("spot_id").loc[spots]
    dom_codes = pd.Categorical(dom["domain"]).codes
    n_dom = len(set(dom_codes))
    print(f"  spots={len(spots)} genes_du={du.shape[0]} genes_expr={expr.shape[0]} domains={n_dom}", flush=True)

    # domain design matrix (dummies)
    Dmat = np.column_stack([np.asarray([1.0 if c == k else 0.0 for c in dom_codes]) for k in sorted(set(dom_codes))])
    # per-spot depth covariate (controls detection/library-size confound: high-read
    # spots have both higher expression AND more-observable APA -> spurious +beta)
    spot_depth = np.log1p(expr.sum(axis=0).reindex(spots).values.astype(float))

    common_genes = sorted(set(du.index) & set(expr.index))
    print(f"  genes in both du+expr: {len(common_genes)}", flush=True)

    try:
        import statsmodels.api as sm
        HAVE_SM = True
    except Exception:
        HAVE_SM = False
        print("  (statsmodels unavailable; using numpy OLS)", flush=True)

    rows = []
    for i, g in enumerate(common_genes):
        duv = du.loc[g].values.astype(float)
        ev = expr.loc[g].values.astype(float)
        obs = np.isfinite(duv) & (ev > 0)
        if obs.sum() < MIN_SPOTS:
            continue
        y = np.log1p(ev[obs])
        x_du = duv[obs]
        Dd = Dmat[obs]
        depth = spot_depth[obs]
        X = np.column_stack([np.ones(obs.sum()), x_du, Dd, depth])
        # OLS
        XtX = X.T @ X
        try:
            beta = np.linalg.solve(XtX, X.T @ y)
            resid = y - X @ beta
            n, p = X.shape
            if n > p:
                sigma2 = resid @ resid / (n - p)
                cov = sigma2 * np.linalg.inv(XtX)
                se = np.sqrt(np.diag(cov))
                b_du = beta[1]; se_du = se[1]
                tstat = b_du / se_du if se_du > 0 else 0.0
                from scipy import stats as ss
                pval = 2 * ss.t.sf(abs(tstat), n - p)
            else:
                b_du = beta[1]; pval = np.nan
            rows.append({"gene": g, "beta_du": float(b_du), "pvalue": float(pval), "n": int(obs.sum())})
        except Exception:
            continue
        if (i + 1) % 500 == 0:
            print(f"    ...{i+1}/{len(common_genes)} genes", flush=True)

    res = pd.DataFrame(rows)
    # FDR
    from scipy import stats as ss
    from statsmodels.stats.multitest import multipletests
    valid = res["pvalue"].notna()
    padj = np.full(len(res), np.nan)
    _, padj[valid.values], _, _ = multipletests(res.loc[valid, "pvalue"].values, method="fdr_bh")
    res["padj"] = padj
    res.to_csv(f"{OUT}/per_gene_beta.csv", index=False)

    sig = res[res["padj"] < 0.05].copy()
    print(f"\n=== Phase 1 regression (go/no-go) ===", flush=True)
    print(f"  genes tested: {len(res)}", flush=True)
    print(f"  significant (padj<0.05): {len(sig)} ({100*len(sig)/max(len(res),1):.1f}%)", flush=True)
    if len(sig):
        print(f"  beta_du sign: neg={int((sig['beta_du']<0).sum())} pos={int((sig['beta_du']>0).sum())}", flush=True)
        print(f"  median |beta_du| (sig): {sig['beta_du'].abs().median():.3f}", flush=True)

    # self-consistency: join miRNA annotation; genes w/ more distal miRNAs -> beta more negative?
    try:
        mir = pd.read_csv(f"{PROJ}/pipeline_output/sapa_regnet/apa_regulatory_annotation_mirna.csv")
        if "n_mirna_apa" not in mir.columns:
            mir["n_mirna_apa"] = mir.get("mirna_in_apa_region", "").astype(str).apply(
                lambda s: len([x for x in s.split(",") if x and x != "nan"]))
        m = res.merge(mir[["gene_name", "n_mirna_apa"]], left_on="gene", right_on="gene_name", how="left")
        m["n_mirna_apa"] = m["n_mirna_apa"].fillna(0)
        has_mirna = m[m["n_mirna_apa"] > 0]
        no_mirna = m[m["n_mirna_apa"] == 0]
        print(f"\n  self-consistency (b_APA sign vs distal miRNA):", flush=True)
        print(f"    genes w/ distal miRNAs (n={len(has_mirna)}): median beta_du={has_mirna['beta_du'].median():.3f}, %neg={100*(has_mirna['beta_du']<0).mean():.0f}%", flush=True)
        print(f"    genes w/o distal miRNAs (n={len(no_mirna)}): median beta_du={no_mirna['beta_du'].median():.3f}, %neg={100*(no_mirna['beta_du']<0).mean():.0f}%", flush=True)
        # correlation beta vs n_mirna (expect negative)
        sub = m.dropna(subset=["beta_du", "n_mirna_apa"])
        if len(sub) > 20:
            r, p = ss.pearsonr(sub["n_mirna_apa"], sub["beta_du"])
            print(f"    corr(n_mirna_apa, beta_du) = {r:.3f} (p={p:.2e}) -- expect NEGATIVE if repressive", flush=True)
        m.to_csv(f"{OUT}/per_gene_beta_with_mirna.csv", index=False)
    except Exception as e:
        print(f"  (miRNA self-consistency skipped: {e})", flush=True)

    print(f"\n  wall {time.time()-t0:.1f}s", flush=True)
    print(f"  -> {OUT}/per_gene_beta.csv", flush=True)


if __name__ == "__main__":
    main()

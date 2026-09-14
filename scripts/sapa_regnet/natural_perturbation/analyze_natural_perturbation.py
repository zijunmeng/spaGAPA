#!/usr/bin/env python3
"""sAPA-RegNet natural-perturbation validation — Phases 2-6.

Uses natural spatial variation of 15 neural RBP (+NUDT21/CFIm25) expression
across spatial domains in GSE293464 retina organoids as "natural partial
knockdowns", testing whether sAPA-RegNet APA-region RBP annotations predict
the DIRECTION of distal-PAS usage shifts.

Design
  - Per (sample, RBP): rank domains by mean expression; high = top-5, low =
    bottom-5 (samples have 15-19 domains). Screening per spec: within-sample
    max/min domain FC > 1.5 (>= 2 domains differing > 1.5x). Low-expression
    domains near zero are kept — a near-absent domain is a strong natural
    perturbation; the FC baseline guard was dropped after it removed nearly
    all contrasts. RBP skipped entirely if mean expression < 1 across all
    domains.
  - Per gene: delta = mean(distal_usage | high domains) - mean(... | low),
    per sample requiring >= 2 valid domains on each side; gene delta = mean
    over all samples where measurable (>= 1).
  - Experimental set: genes with FIMO motif sites in the APA region
    (proximal PAS -> TES, genomic interval overlap; FIMO coordinates used
    directly — fixes the offset bug in parse_fimo_rbp_genomewide.py).
  - Control set: measurable genes without that RBP's APA-region sites.
  - IMPORTANT confound handling: high-RBP domains are shared across RBPs
    (cell-type composition), so the background delta distal usage is globally
    shifted (consistency of random gene sets != 0.5). All inference is
    therefore exp-vs-ctrl: two-sided Mann-Whitney U, rank-biserial effect,
    delta_consistency = consistency_exp - consistency_ctrl with 1000x
    label-shuffle null, BH-FDR over RBP tests.

Literature direction priors (sign of distal-usage change when RBP HIGH):
  PTBP1 -1 (Yeom 2019, Zhu 2022) | NOVA1/2 +1 (Licatalosi 2012, Hwang 2017)
  ELAVL1 -1 (Deng 2023) | KHSRP +1 | FUS/TARDBP -1 (Rot 2017) | HNRNPA1 -1
  FMR1 +1 | TIA1/TIAL1 -1 | NUDT21/CFIm25 +1 (Masamha 2014)
  ELAVL2/3/4 -1 (family-inferred) | HNRNPK none (exploratory)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
OUT = PROJ / "pipeline_output/sapa_regnet/natural_perturbation"
CACHE = OUT / "_cache"

SAMPLES = ["GSM8882884", "GSM8882885", "GSM8882886", "GSM8882887"]
N_GROUP = 5
MIN_VALID_DOMAINS = 2
MIN_SAMPLES = 1
MIN_EXPR_ALLDOMAINS = 1.0
FC_SCREEN = 1.5
N_SHUFFLE = 1000
SEED = 42

PREDICTION = {
    "PTBP1": (-1, "Yeom et al. 2019; Zhu et al. 2022"),
    "NOVA1": (+1, "Licatalosi et al. 2012; Hwang et al. 2017"),
    "NOVA2": (+1, "Licatalosi et al. 2012; Hwang et al. 2017"),
    "ELAVL1": (-1, "Deng et al. 2023 (review)"),
    "KHSRP": (+1, "promotes distal cleavage (spec)"),
    "FUS": (-1, "Rot et al. 2017"),
    "TARDBP": (-1, "Rot et al. 2017"),
    "HNRNPA1": (-1, "competitive distal inhibition (spec)"),
    "FMR1": (+1, "distal PAS selection (spec)"),
    "TIA1": (-1, "proximal PAS (spec)"),
    "TIAL1": (-1, "proximal PAS (spec)"),
    "NUDT21": (+1, "Masamha et al. 2014 (CFIm25)"),
    "ELAVL2": (-1, "family-inferred from ELAVL1"),
    "ELAVL3": (-1, "family-inferred from ELAVL1"),
    "ELAVL4": (-1, "family-inferred from ELAVL1"),
    "HNRNPK": (0, "no prior (exploratory)"),
}
RBP_ORDER = ["PTBP1", "NOVA1", "NOVA2", "ELAVL1", "ELAVL2", "ELAVL3", "ELAVL4",
             "FUS", "TARDBP", "HNRNPA1", "HNRNPK", "KHSRP", "FMR1", "TIA1",
             "TIAL1", "NUDT21"]


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z*z / (2*n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z*z / (4*n*n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def main():
    rng = np.random.default_rng(SEED)

    expr = pd.read_csv(OUT / "rbp_domain_expression.csv")
    usage = pd.read_csv(CACHE / "domain_usage_all.csv")
    genesets = pd.read_csv(CACHE / "rbp_apa_genes.csv")
    gs = genesets.groupby("rbp")["gene"].apply(set).to_dict()
    rbps = [r for r in RBP_ORDER if r in gs]

    upiv = {}
    for s in SAMPLES:
        u = usage[usage["sample"] == s]
        upiv[s] = u.pivot(index="gene", columns="domain", values="distal_usage")

    # ------------- Phase 2: screening + high/low domain assignment -------------
    contrasts = {}
    screen_rows = []
    for rbp in rbps:
        e = expr[expr.gene == rbp]
        overall_mean = float(e.mean_expr.mean())
        if overall_mean < MIN_EXPR_ALLDOMAINS:
            screen_rows.append(dict(rbp=rbp, status="skipped_low_expression",
                                    overall_mean_expr=round(overall_mean, 3),
                                    n_samples_passing=0))
            continue
        n_pass = 0
        for s in SAMPLES:
            es = e[e["sample"] == s].set_index("domain").mean_expr
            es = es.fillna(0.0).sort_values(ascending=False)
            if len(es) < 2 * N_GROUP:
                continue
            fc = es.max() / es.min() if es.min() > 0 else np.inf
            if fc <= FC_SCREEN:
                continue
            high, low = list(es.index[:N_GROUP]), list(es.index[-N_GROUP:])
            contrasts[(s, rbp)] = dict(
                high=[int(d) for d in high], low=[int(d) for d in low],
                fc_maxmin=float(fc), high_mean=float(es.loc[high].mean()),
                low_mean=float(es.loc[low].mean()))
            n_pass += 1
        screen_rows.append(dict(rbp=rbp,
                                status="tested" if n_pass else "no_passing_sample",
                                overall_mean_expr=round(overall_mean, 3),
                                n_samples_passing=n_pass))
    screen_df = pd.DataFrame(screen_rows)
    print("Phase 2 screening:")
    print(screen_df.to_string(index=False))
    with open(CACHE / "contrasts.json", "w") as fh:
        json.dump({f"{s}|{r}": v for (s, r), v in contrasts.items()}, fh,
                  indent=1)

    def gene_deltas(rbp, min_samples=MIN_SAMPLES):
        per_gene = {}
        for s in SAMPLES:
            if (s, rbp) not in contrasts:
                continue
            up = upiv[s]
            hi_v = up.reindex(columns=contrasts[(s, rbp)]["high"])
            lo_v = up.reindex(columns=contrasts[(s, rbp)]["low"])
            d_s = hi_v.mean(axis=1, skipna=True) - lo_v.mean(axis=1, skipna=True)
            ok = ((hi_v.notna().sum(axis=1) >= MIN_VALID_DOMAINS) &
                  (lo_v.notna().sum(axis=1) >= MIN_VALID_DOMAINS))
            for gene, val in d_s[ok].items():
                per_gene.setdefault(gene, []).append(float(val))
        return {g: (float(np.mean(v)), len(v)) for g, v in per_gene.items()
                if len(v) >= min_samples}

    all_deltas = {rbp: gene_deltas(rbp) for rbp in rbps}

    # ------------- Phase 4/5: per-RBP statistics -------------
    results, per_gene_rows = [], []
    for rbp in rbps:
        deltas = all_deltas[rbp]
        if not deltas or not any((s, rbp) in contrasts for s in SAMPLES):
            continue
        genes = sorted(deltas)
        dvals = pd.Series({g: deltas[g][0] for g in genes})
        nsamp = pd.Series({g: deltas[g][1] for g in genes})
        site_genes = gs[rbp]
        exp = dvals[dvals.index.isin(site_genes)]
        ctrl = dvals[~dvals.index.isin(site_genes)]
        pred, source = PREDICTION.get(rbp, (0, "none"))

        row = dict(rbp=rbp, prediction=pred, prediction_source=source,
                   n_contrast_samples=sum(1 for s in SAMPLES
                                          if (s, rbp) in contrasts),
                   n_genes_exp=len(exp), n_genes_ctrl=len(ctrl),
                   mean_n_samples_per_gene=float(nsamp.mean()))

        def consistency(series):
            nz = series[series != 0]
            if len(nz) == 0:
                return np.nan, 0
            return float((np.sign(nz) == pred).mean()), len(nz)

        if pred != 0 and len(exp) and len(ctrl):
            c_exp, n_exp = consistency(exp)
            c_ctrl, n_ctrl = consistency(ctrl)
            row["direction_consistency"] = c_exp
            row["direction_consistency_ctrl"] = c_ctrl
            row["delta_consistency"] = c_exp - c_ctrl
            lo_ci, hi_ci = wilson_ci(int(round(c_exp * n_exp)), n_exp)
            row["consistency_ci95_low"], row["consistency_ci95_high"] = lo_ci, hi_ci
            row["consistency_binom_p"] = float(stats.binomtest(
                int(round(c_exp * n_exp)), n_exp, 0.5).pvalue)
            # label-shuffle null on delta_consistency (background-corrected)
            vals = dvals.to_numpy()
            labels = dvals.index.isin(site_genes)
            n_draw = int(labels.sum())
            null_dc = np.empty(N_SHUFFLE)
            for i in range(N_SHUFFLE):
                perm = rng.permutation(len(vals))
                d_exp = vals[perm[:n_draw]]
                d_ctrl = vals[perm[n_draw:]]
                ce = (np.sign(d_exp[d_exp != 0]) == pred).mean()
                cc = (np.sign(d_ctrl[d_ctrl != 0]) == pred).mean()
                null_dc[i] = ce - cc
            row["null_delta_consistency_mean"] = float(null_dc.mean())
            row["null_delta_consistency_sd"] = float(null_dc.std())
            row["shuffle_p_value"] = float(
                (null_dc >= row["delta_consistency"]).mean())
        if len(exp) >= 10 and len(ctrl) >= 10:
            u, p = stats.mannwhitneyu(exp, ctrl, alternative="two-sided")
            row["mwu_p_value"] = float(p)
            row["effect_size_rank_biserial"] = float(
                2 * u / (len(exp) * len(ctrl)) - 1)
        row["median_delta_exp"] = float(exp.median()) if len(exp) else np.nan
        row["median_delta_ctrl"] = float(ctrl.median()) if len(ctrl) else np.nan
        row["mean_delta_exp"] = float(exp.mean()) if len(exp) else np.nan
        row["mean_delta_ctrl"] = float(ctrl.mean()) if len(ctrl) else np.nan
        results.append(row)
        for g in genes:
            per_gene_rows.append((rbp, g, float(dvals[g]), int(nsamp[g]),
                                  bool(g in site_genes)))

    res = pd.DataFrame(results)
    res["fdr_bh"] = bh_fdr(res["mwu_p_value"].fillna(1.0).to_numpy())
    res["_o"] = res.rbp.map({r: i for i, r in enumerate(RBP_ORDER)})
    res = res.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    per_gene = pd.DataFrame(per_gene_rows, columns=["rbp", "gene", "delta",
                                                    "n_samples", "in_site_set"])
    res.to_csv(OUT / "perturbation_results.csv", index=False)
    per_gene.to_csv(OUT / "per_gene_deltas.csv", index=False)
    print("\nPhase 4 results:")
    cols = ["rbp", "prediction", "n_genes_exp", "n_genes_ctrl",
            "direction_consistency", "direction_consistency_ctrl",
            "delta_consistency", "shuffle_p_value", "mwu_p_value",
            "effect_size_rank_biserial", "fdr_bh", "median_delta_exp",
            "median_delta_ctrl"]
    print(res[cols].round(4).to_string(index=False))

    # ------------- cross-RBP control -------------
    cross = []
    for rbp_a in rbps:
        pred, _ = PREDICTION.get(rbp_a, (0, "none"))
        if pred == 0:
            continue
        for rbp_b in rbps:
            if rbp_a == rbp_b or not all_deltas.get(rbp_b):
                continue
            exp = pd.Series({g: v[0] for g, v in all_deltas[rbp_b].items()
                             if g in gs[rbp_a] and v[0] != 0})
            ctrl = pd.Series({g: v[0] for g, v in all_deltas[rbp_b].items()
                              if g not in gs[rbp_a] and v[0] != 0})
            if len(exp) < 20 or len(ctrl) < 20:
                continue
            ce = (np.sign(exp) == pred).mean()
            cc = (np.sign(ctrl) == pred).mean()
            cross.append((rbp_a, rbp_b, float(ce), float(cc), float(ce - cc)))
    cross_df = pd.DataFrame(cross, columns=["rbp_genes", "rbp_contrast",
                                            "consistency_exp",
                                            "consistency_ctrl",
                                            "delta_consistency"])
    cross_df.to_csv(CACHE / "cross_rbp_consistency.csv", index=False)
    diag = res.set_index("rbp")["delta_consistency"]
    offdiag = cross_df[cross_df.rbp_genes != cross_df.rbp_contrast]
    print("\ncross-RBP control: mean off-diagonal delta_consistency = "
          f"{offdiag.delta_consistency.mean():.4f} vs mean diagonal "
          f"{diag.mean():.4f}")

    # ------------- summary -------------
    def g(rbp, col, default=np.nan):
        v = res.loc[res.rbp == rbp, col]
        return float(v.iloc[0]) if len(v) else default

    # direction concordance vs significance (background-corrected):
    # consistency already measures agreement with the prior, so a CONCORDANT
    # differential requires delta_consistency > 0 (and effect sign == prior
    # sign) for every RBP regardless of prior direction.
    res["concordant"] = (res["delta_consistency"] > 0) & (
        np.sign(res["effect_size_rank_biserial"]) == np.sign(res["prediction"]))
    res["sig"] = (res["fdr_bh"] < 0.10) | (res["shuffle_p_value"] < 0.05)
    res["sig_fdr"] = res["fdr_bh"] < 0.10
    wp = res[res.prediction != 0]
    conc_sig = wp[wp.concordant & wp.sig]
    conc_fdr = wp[wp.concordant & wp.sig_fdr]
    anti_sig = wp[~wp.concordant & wp.sig]
    print("\nconcordant & significant:", conc_sig.rbp.tolist())
    print("concordant & FDR<0.10:", conc_fdr.rbp.tolist())
    print("ANTI-concordant & significant:", anti_sig.rbp.tolist())
    n_fdr = int((res["fdr_bh"] < 0.10).sum())
    n_shuf = int((res["shuffle_p_value"] < 0.05).sum())
    with_prior = res[res.prediction != 0]
    crit = {
        "PTBP1": {
            "consistency": g("PTBP1", "direction_consistency"),
            "consistency_ctrl": g("PTBP1", "direction_consistency_ctrl"),
            "delta_consistency": g("PTBP1", "delta_consistency"),
            "shuffle_p": g("PTBP1", "shuffle_p_value", 1.0),
            "mwu_p": g("PTBP1", "mwu_p_value", 1.0),
            "fdr": g("PTBP1", "fdr_bh", 1.0),
            "pass": bool(g("PTBP1", "direction_consistency", 0) > 0.60 and
                         min(g("PTBP1", "shuffle_p_value", 1.0),
                             g("PTBP1", "mwu_p_value", 1.0)) < 0.05)},
        "CFIm25_positive_control": {
            "consistency": g("NUDT21", "direction_consistency"),
            "consistency_ctrl": g("NUDT21", "direction_consistency_ctrl"),
            "delta_consistency": g("NUDT21", "delta_consistency"),
            "median_delta_exp": g("NUDT21", "median_delta_exp"),
            "median_delta_ctrl": g("NUDT21", "median_delta_ctrl"),
            "shuffle_p": g("NUDT21", "shuffle_p_value", 1.0),
            "mwu_p": g("NUDT21", "mwu_p_value", 1.0),
            "fdr": g("NUDT21", "fdr_bh", 1.0),
            "concordant_masamha_2014": bool(
                g("NUDT21", "delta_consistency", 0) > 0 and
                g("NUDT21", "median_delta_exp", 0) >
                g("NUDT21", "median_delta_ctrl", 0)),
            "pass": bool(g("NUDT21", "delta_consistency", 0) > 0 and
                         min(g("NUDT21", "shuffle_p_value", 1.0),
                             g("NUDT21", "mwu_p_value", 1.0)) < 0.10)},
        "negative_control": {
            "median_abs_ctrl_delta": float(
                res["median_delta_ctrl"].abs().median()),
            "null_delta_consistency_abs_mean": float(
                res["null_delta_consistency_mean"].abs().mean()),
            "note": "background consistency != 0.5 by construction (shared high-RBP domains); specificity comes from exp-minus-ctrl contrast whose shuffle null is centred at 0",
            "pass": bool(res["null_delta_consistency_mean"].abs().mean() < 0.02)},
        "overall": {
            "n_rbp_fdr_lt_0.10_any_direction": n_fdr,
            "rbps_fdr_lt_0.10_any_direction": res.loc[
                res["fdr_bh"] < 0.10, "rbp"].tolist(),
            "n_rbp_shuffle_p_lt_0.05": n_shuf,
            "n_concordant_significant": int(len(conc_sig)),
            "rbps_concordant_significant": conc_sig.rbp.tolist(),
            "n_concordant_fdr_lt_0.10": int(len(conc_fdr)),
            "rbps_concordant_fdr_lt_0.10": conc_fdr.rbp.tolist(),
            "n_anti_concordant_significant": int(len(anti_sig)),
            "rbps_anti_concordant_significant": anti_sig.rbp.tolist(),
            "note": "spec criterion '>=3 RBP FDR<0.10' evaluated on direction-concordant RBPs only; anti-concordant significant RBPs counted against the model",
            "pass": bool(len(conc_fdr) >= 3)},
    }
    # verdict tiers (honesty-first):
    #   validated: PTBP1 pass + CFIm25 pass + >=3 concordant FDR<0.10
    #   partial:   CFIm25 pass + >=2 concordant significant (incl >=1 FDR<0.10)
    #   insufficient: otherwise
    full = (crit["PTBP1"]["pass"] and
            crit["CFIm25_positive_control"]["pass"] and len(conc_fdr) >= 3)
    partial = (crit["CFIm25_positive_control"]["pass"] and
               len(conc_sig) >= 2 and len(conc_fdr) >= 1)
    if full:
        verdict = "validated_with_natural_perturbation"
    elif partial:
        verdict = "partial_validation_positive_control_supported"
    else:
        verdict = "natural_variation_insufficient"

    summary = {
        "title": "sAPA-RegNet natural-perturbation validation (limitation #8)",
        "date": "2026-09-14",
        "design": {
            "samples": SAMPLES,
            "n_domains_per_sample": {s: int(usage[usage["sample"] == s]
                                            .domain.nunique()) for s in SAMPLES},
            "note_domain_count": "spec assumed 5 domains/sample; actual Leiden domain counts are 15/15/19/19 — all domains used",
            "high_low_domain_groups": f"top-{N_GROUP} vs bottom-{N_GROUP} domains by mean RBP expression per sample",
            "screening": f"within-sample max/min domain FC > {FC_SCREEN}; RBP mean expr >= {MIN_EXPR_ALLDOMAINS} across all domains",
            "gene_delta": f"mean distal usage (high - low domains); >= {MIN_VALID_DOMAINS} valid domains per side; gene delta averaged over all measurable samples",
            "distal_usage": "u_distal/(u_prox+u_distal) per bin200 spot; domain mean requiring >=5 covered spots",
            "site_definition": "FIMO motif sites overlapping Track-A apa_region (proximal PAS->TES); FIMO 5.5.9 genomic coords used directly (fixes offset bug in parse_fimo_rbp_genomewide.py)",
            "nudt21_target": "genes with >=2 UUGUA elements (ATtRACT NUDT21_420/575, consensus p<=1e-3) in APA region",
            "statistics": "two-sided Mann-Whitney U exp vs ctrl; rank-biserial; direction consistency exp and no-site ctrl; delta_consistency with 1000x label-shuffle null; BH-FDR",
            "confound_handling": "high-RBP domains are shared across RBPs (cell-type composition) -> background delta-consistency of random sets is not 0.5; all tests are background-corrected (exp minus no-site ctrl)"},
        "screening": screen_df.to_dict("records"),
        "acceptance_criteria": crit,
        "key_numbers": {
            "n_rbp_tested": int(res.rbp.nunique()),
            "n_rbp_with_prior": int(len(with_prior)),
            "rbps_skipped_low_expression": screen_df.loc[
                screen_df.status == "skipped_low_expression", "rbp"].tolist(),
            "mean_consistency_exp": float(
                with_prior["direction_consistency"].mean()),
            "mean_consistency_ctrl": float(
                with_prior["direction_consistency_ctrl"].mean()),
            "mean_delta_consistency": float(
                with_prior["delta_consistency"].mean()),
            "n_prior_delta_consistency_gt_0": int(
                (with_prior["delta_consistency"] > 0).sum()),
            "rbps_concordant_significant": conc_sig.rbp.tolist(),
            "rbps_concordant_fdr_lt_0.10": conc_fdr.rbp.tolist(),
            "rbps_anti_concordant_significant": anti_sig.rbp.tolist(),
            "all_rank_biserial_positive": bool(
                (res["effect_size_rank_biserial"] > 0).all()),
            "effect_size_range": [float(res["effect_size_rank_biserial"].min()),
                                  float(res["effect_size_rank_biserial"].max())],
            "cross_rbp_mean_offdiag_delta_consistency":
                float(offdiag.delta_consistency.mean()),
            "cross_rbp_mean_diag_delta_consistency": float(diag.mean()),
            "interpretation_caveats": [
                "all effect sizes (exp vs no-site ctrl) are positive -> genes with ANY RBP motif in the APA region shift slightly distal-up in high-RBP domains; the concordant signals partly ride on this shared component",
                "cross-RBP control: own-RBP contrast is not stronger than other RBPs' contrasts (off-diagonal ~= diagonal) -> motif-set specificity is weak in this natural-variation design",
                "PTBP1, TIA1, TIAL1 show significant shifts OPPOSITE to their priors",
            ],
        },
        "per_rbp": json.loads(res.to_json(orient="records")),
        "verdict": verdict,
    }
    with open(OUT / "validation_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    # ------------- figure -------------
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # A: heatmap (sample,domain) x RBP, z-scored per sample
    ax = axes[0, 0]
    piv = expr.pivot_table(index=["sample", "domain"], columns="gene",
                           values="mean_expr")
    piv = piv.reindex(columns=[c for c in RBP_ORDER if c in piv.columns])
    z = piv.groupby(level="sample", group_keys=False).apply(
        lambda x: (x - x.mean()) / (x.std() + 1e-9)).sort_index()
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdYlBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(z.shape[1]), z.columns, rotation=90, fontsize=8)
    run = 0
    for s in SAMPLES[:-1]:
        n_d = sum(1 for ss, _ in z.index if ss == s)
        ax.axhline(run + n_d - 0.5, color="k", lw=1.5)
        run += n_d
    tickpos, seen = [], set()
    for i, (s, d) in enumerate(z.index):
        if s not in seen:
            tickpos.append(i)
            seen.add(s)
    ax.set_yticks(tickpos, SAMPLES, fontsize=9)
    ax.set_title("A  RBP expression across spatial domains (z within sample)")
    plt.colorbar(im, ax=ax, shrink=0.7, label="z-score")

    # B: direction consistency exp vs ctrl, with CI
    ax = axes[0, 1]
    r = res[res.prediction != 0].dropna(subset=["direction_consistency"])
    r = r.sort_values("delta_consistency", ascending=False)
    x = np.arange(len(r))
    w = 0.38
    conc_color = ["#2166ac" if v > 0 else "#b2182b"
                  for v in r.delta_consistency]
    ax.bar(x - w / 2, r.direction_consistency, width=w, color=conc_color,
           alpha=0.9, label="site genes (blue=prior-concordant)")
    ax.bar(x + w / 2, r.direction_consistency_ctrl, width=w, color="#bdbdbd",
           alpha=0.9, label="no-site genes (background)")
    err = np.vstack([r.direction_consistency - r.consistency_ci95_low,
                     r.consistency_ci95_high - r.direction_consistency])
    ax.errorbar(x - w / 2, r.direction_consistency, yerr=err, fmt="none",
                ecolor="k", capsize=3, lw=1)
    ax.axhline(0.5, ls="--", color="k", lw=1)
    ax.set_xticks(x, r.rbp, rotation=90)
    ax.set_ylabel("direction consistency")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("B  Direction consistency vs prior — site genes vs background")
    for i, (sp, fd) in enumerate(zip(r.shuffle_p_value, r.fdr_bh)):
        star = "**" if sp < 0.01 else ("*" if sp < 0.05 else
                                       ("+" if fd < 0.10 else ""))
        if star:
            ax.text(i, 1.0, star, ha="center", fontsize=10)

    # C: exp vs ctrl median deltas per RBP
    ax = axes[1, 0]
    r2 = res.dropna(subset=["median_delta_exp"]).copy()
    r2["_o"] = r2.rbp.map({v: i for i, v in enumerate(RBP_ORDER)})
    r2 = r2.sort_values("_o")
    xs = np.arange(len(r2))
    for i, row in enumerate(r2.itertuples()):
        ax.plot([i - 0.18, i + 0.18],
                [row.median_delta_exp, row.median_delta_ctrl],
                color="k", lw=1, alpha=0.6)
    ax.scatter(xs - 0.18, r2.median_delta_exp, s=45, color="#d6604d", zorder=3,
               label="site genes (exp)", marker="o")
    ax.scatter(xs + 0.18, r2.median_delta_ctrl, s=45, color="#4393c3", zorder=3,
               label="no-site genes (ctrl)", marker="s")
    ax.axhline(0, ls="--", color="gray", lw=1)
    ax.set_xticks(xs, r2.rbp, rotation=90)
    ax.set_ylabel("median delta distal usage (high - low RBP domain)")
    ax.legend(fontsize=8)
    ax.set_title("C  Site genes vs no-site genes (medians, paired per RBP)")

    # D: CFIm25 positive control
    ax = axes[1, 1]
    deltas_n = all_deltas.get("NUDT21", {})
    nd = gs.get("NUDT21", set())
    exp_n = pd.Series({k2: v[0] for k2, v in deltas_n.items() if k2 in nd})
    ctrl_n = pd.Series({k2: v[0] for k2, v in deltas_n.items()
                        if k2 not in nd})
    bp = ax.boxplot([exp_n.dropna().to_numpy(), ctrl_n.dropna().to_numpy()],
                    tick_labels=[f"NUDT21 UUGUA x>=2\n(n={len(exp_n)})",
                                 f"no sites\n(n={len(ctrl_n)})"],
                    showfliers=False, widths=0.5, patch_artist=True)
    for patch, c in zip(bp["boxes"], ["#d6604d", "#4393c3"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    ax.axhline(0, ls="--", color="gray", lw=1)
    med_e = float(exp_n.median())
    med_c = float(ctrl_n.median())
    ax.annotate("", xy=(1.28, med_e), xytext=(1.28, med_c),
                arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.text(1.33, (med_e + med_c) / 2,
            f"Δmedians = {med_e - med_c:+.4f}\nMasamha 2014:\nCFIm25 high -> distal up",
            fontsize=9, va="center")
    if (res.rbp == "NUDT21").any():
        rr = res[res.rbp == "NUDT21"].iloc[0]
        ax.set_title("D  CFIm25/NUDT21 positive control — consistency "
                     f"{rr.direction_consistency:.2f} vs bg "
                     f"{rr.direction_consistency_ctrl:.2f}, "
                     f"MWU p={rr.mwu_p_value:.2g}")
    ax.set_xlim(0.5, 2.3)
    ax.set_ylabel("delta distal usage")

    fig.tight_layout()
    fig.savefig(OUT / "fig_regnet_validation.png", dpi=200)
    print(f"\nwrote figure + tables + summary -> {OUT}")
    print("VERDICT:", verdict)
    print("DONE")


if __name__ == "__main__":
    main()

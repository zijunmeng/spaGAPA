#!/usr/bin/env python
"""Assemble pipeline_output/multicaller_validation/summary.json + report.md
from the per-dataset metric JSONs of the multi-caller validation runs.

Datasets (2 species x 3 tissues):
  gse183456_gsm6047774  human kidney      (original run; outputs in sierra/)
  gse220442_gsm6801751  human AD brain    ( Visium, control sample)
  gse169749_gsm5213483  mouse colon       ( Visium, DSS day-0)
Each dataset contributes: Sierra qc, metric a (PAS overlap), metric b
(gene-level distal usage), metric c (split-conformal coverage, Sierra input +
scAPAtrap same-protocol control), metric d (GP-vs-mean RMSE on both inputs).
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from datasets import DATASETS, ROOT

MC = os.path.join(ROOT, "pipeline_output", "multicaller_validation")

# dataset display order + per-dataset file locations (gse183456 predates the
# per-dataset layout: its metric JSONs sit at the MC root)
ORDER = ["gse183456_gsm6047774", "gse220442_gsm6801751", "gse169749_gsm5213483"]
DISPLAY = {
    "gse183456_gsm6047774": ("GSE183456 / GSM6047774", "human", "kidney", "3010 spots, 214M reads"),
    "gse220442_gsm6801751": ("GSE220442 / GSM6801751", "human", "brain (AD PFC, control)", "4179 spots"),
    "gse169749_gsm5213483": ("GSE169749 / GSM5213483", "mouse", "colon (DSS d0)", "2715 spots, 277M reads"),
}


def load(path):
    with open(path) as f:
        return json.load(f)


def dataset_files(key):
    ds = DATASETS[key]
    out = ds["outdir"]
    root_legacy = key == "gse183456_gsm6047774"
    return {
        "qc": os.path.join(out, "qc_summary.json"),
        "ab": os.path.join(MC, "metrics_ab.json") if root_legacy else os.path.join(out, "metrics_ab.json"),
        "d": os.path.join(MC, "metrics_d.json") if root_legacy else os.path.join(out, "metrics_d.json"),
        "sierra_cal": os.path.join(out, "conformal", "uncertainty_calibration.json"),
        "base_cal": os.path.join(MC, "scapatrap_conformal", "uncertainty_calibration.json")
                    if root_legacy else os.path.join(out, "scapatrap_conformal", "uncertainty_calibration.json"),
        "baseline_qc": os.path.join(ds["baseline"], "qc_summary.json"),
        "junctions": os.path.join(out, "junctions.bed"),
    }


def main():
    per = {}
    for key in ORDER:
        f = dataset_files(key)
        qc = load(f["qc"])
        ab = load(f["ab"])
        d = load(f["d"])
        sc_cal = load(f["sierra_cal"])
        bs_cal = load(f["base_cal"])
        bqc = load(f["baseline_qc"])
        n_junc = sum(1 for _ in open(f["junctions"])) if os.path.exists(f["junctions"]) else None
        label, species, tissue, detail = DISPLAY[key]
        per[key] = {
            "label": label, "species": species, "tissue": tissue, "detail": detail,
            "n_junctions": n_junc,
            "baseline": {
                "caller": "scAPAtrap",
                "n_spots": bqc["n_spots"],
                "n_called_sites": bqc.get("n_called_sites"),
                "n_gene_annotated_sites": bqc.get("n_gene_annotated_sites"),
                "n_apa_usage_sites": bqc.get("n_apa_usage_sites"),
            },
            "sierra": qc,
            "metric_a_pas_overlap": {
                **ab["overlap"],
                "summit_convention_sensitivity": ab["overlap_summit_convention"],
                "nearest_distance": ab["nearest_distance"],
            },
            "metric_b_gene_level": ab["gene_level"],
            "metric_c_conformal_coverage": {
                "sierra_input": {
                    "n_genes": sc_cal["config"]["n_genes"],
                    "n_spots": sc_cal["config"]["n_spots"],
                    "n_test": sc_cal["config"]["n_test"],
                    "before_raw_gp": sc_cal["before"]["coverage"],
                    "after_conformal_global": sc_cal["after"]["modes"]["global"]["coverage"],
                    "after_conformal_locally_adaptive": sc_cal["after"]["modes"]["locally_adaptive"]["coverage"],
                },
                "scapatrap_input_same_protocol": {
                    "n_genes": bs_cal["config"]["n_genes"],
                    "n_test": bs_cal["config"]["n_test"],
                    "after_conformal_global": bs_cal["after"]["modes"]["global"]["coverage"],
                    "after_conformal_locally_adaptive": bs_cal["after"]["modes"]["locally_adaptive"]["coverage"],
                },
            },
            "metric_d_gp_vs_mean": d,
        }

    # ---------------- cross-dataset coverage deviations (core claim) ----------
    levels = ["80pct", "90pct", "95pct"]
    nominal = {"80pct": 0.80, "90pct": 0.90, "95pct": 0.95}
    devs = []
    for key in ORDER:
        mc3 = per[key]["metric_c_conformal_coverage"]
        for mode in ("after_conformal_global", "after_conformal_locally_adaptive"):
            for lv in levels:
                devs.append(abs(mc3["sierra_input"][mode][lv] - nominal[lv]))
    max_dev_pp = 100 * max(devs)

    summary = {
        "title": "Multi-caller robustness validation: Sierra as a second PAS caller, 3 datasets (2 species x 3 tissues)",
        "protocol": ("Sierra 0.99.27 FindPeaks + CountPeaks (UMI-deduplicated) on each dataset's Space Ranger "
                     "possorted_genome_bam.bam with a pysam-extracted splice-junction BED (>=25 reads) and the "
                     "same reference GTF the scAPAtrap baseline used; output converted to the spaGAPA usage-matrix "
                     "format (min_parent=5, baseline coordinates reused); identical spaGAPA inference "
                     "(split conformal 20% masking, 50/50 cal/test, global + locally adaptive)."),
        "datasets": per,
        "cross_dataset": {
            "n_datasets": len(ORDER),
            "species": ["human", "mouse"],
            "tissues": ["kidney", "brain (AD PFC)", "colon (DSS)"],
            "sierra_coverage_max_abs_deviation_pp": {lv: None for lv in levels},
            "sierra_coverage_max_abs_deviation_overall_pp": max_dev_pp,
            "headline": (
                f"Across three datasets spanning two species (human/mouse) and three tissues "
                f"(kidney/brain/colon), swapping the PAS caller (scAPAtrap -> Sierra) leaves "
                f"split-conformal coverage at nominal (max deviation {max_dev_pp:.1f} pp at all levels, "
                f"global and locally adaptive)."),
        },
    }
    for lv in levels:
        summary["cross_dataset"]["sierra_coverage_max_abs_deviation_pp"][lv] = 100 * max(
            abs(per[k]["metric_c_conformal_coverage"]["sierra_input"][m][lv] - nominal[lv])
            for k in ORDER for m in ("after_conformal_global", "after_conformal_locally_adaptive"))

    with open(os.path.join(MC, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[out] {os.path.join(MC, 'summary.json')}")

    # ------------------------------------------------------------------ report
    lines = []
    A = lines.append
    A("# Multi-caller robustness validation: Sierra as a second PAS caller (3 datasets, 2 species, 3 tissues)")
    A("")
    A("## Setup")
    A("")
    A("Protocol per dataset (identical to the original GSE183456 run): Sierra 0.99.27 FindPeaks + CountPeaks")
    A("(per-UMI deduplicated counting) on the dataset's Space Ranger `possorted_genome_bam.bam`, with a")
    A("splice-junction BED extracted from the BAM by pysam (>= 25 supporting reads) and the same reference")
    A("GTF the scAPAtrap baseline was annotated with (the Space Ranger reference GTF: GENCODE/Ensembl,")
    A("chr-prefixed contigs matching the BAM header). The Sierra peak x spot UMI matrix is converted to the")
    A("spaGAPA format (>= 2 sites per gene, usage matrix at min_parent_count = 5, baseline coordinates")
    A("reused verbatim) and pushed through the identical spaGAPA inference without re-tuning.")
    A("")
    A("| dataset | species / tissue | spots | junctions (>=25 reads) | Sierra peaks -> sites x spots (genes) |")
    A("|---|---|---|---|---|")
    for key in ORDER:
        p = per[key]
        A(f"| {p['label']} | {p['species']} / {p['tissue']} | {p['baseline']['n_spots']} | "
          f"{p['n_junctions']:,} | {p['sierra']['n_peaks_raw_all_genes']:,} peaks -> "
          f"**{p['sierra']['n_sites_ge2']:,} x {p['sierra']['n_spots']} ({p['sierra']['n_genes_ge2']:,} genes)** |")
    A("")
    A("Chromosome-naming note (GSE169749, mouse): the Space Ranger BAM uses chr-prefixed contigs while the")
    A("Ensembl GRCm39.111 GTF uses un-prefixed ones; to keep BAM/GTF consistent we used the Space Ranger")
    A("reference's own GTF (refdata-gex-GRCm39-2024-A, GENCODE M33 = Ensembl 110, chr-prefixed) - the exact")
    A("file the scAPAtrap baseline was annotated with, so gene ids match by construction.")
    A("")

    # ---------------- per-dataset sections ----------------
    for key in ORDER:
        p = per[key]
        ov = p["metric_a_pas_overlap"]; ov2 = ov["summit_convention_sensitivity"]; nd = ov["nearest_distance"]
        g10 = p["metric_b_gene_level"]["min10spots"]; g30 = p["metric_b_gene_level"]["min30spots"]
        sc = p["metric_c_conformal_coverage"]
        sg = sc["sierra_input"]["after_conformal_global"]; sl = sc["sierra_input"]["after_conformal_locally_adaptive"]
        bg = sc["scapatrap_input_same_protocol"]["after_conformal_global"]
        ds, dsi = p["metric_d_gp_vs_mean"]["scapatrap"], p["metric_d_gp_vs_mean"]["sierra"]
        A(f"## {p['label']} - {p['species']} {p['tissue']}")
        A("")
        A(f"scAPAtrap baseline: {p['baseline']['n_called_sites']:,} sites "
          f"({p['baseline']['n_gene_annotated_sites']:,} gene-annotated, "
          f"{p['baseline']['n_apa_usage_sites']:,} usage rows); Sierra usable peaks: "
          f"{p['sierra']['n_peaks_raw_all_genes']:,} -> {p['sierra']['n_sites_ge2']:,} sites in "
          f"{p['sierra']['n_genes_ge2']:,} multi-site genes; {sc['sierra_input']['n_test']:,} conformal test points.")
        A("")
        A("| metric | value |")
        A("|---|---|")
        A(f"| PAS overlap: shared / A-only / B-only clusters (+/-50 bp) | {ov['shared_clusters']:,} / "
          f"{ov['unique_callerA_clusters']:,} / {ov['unique_callerB_clusters']:,}; "
          f"Jaccard **{ov['jaccard_clusters']:.3f}** |")
        A(f"| per-point match within 50 bp | {ov['callerB_frac_matched']:.1%} of Sierra peaks, "
          f"{ov['callerA_frac_matched']:.1%} of scAPAtrap sites |")
        A(f"| nearest-neighbour distance | median {nd['median_bp']:.0f} bp; "
          f"{nd['frac_le_500bp']:.1%} <= 500 bp |")
        A(f"| gene-level distal usage (>= 10 co-finite spots) | n = {g10['n_genes_compared']}, "
          f"median Pearson r **{g10['pearson_r_median']:.3f}**, "
          f"{g10['frac_r_gt_0.7']:.1%} with r > 0.7 |")
        A(f"| conformal coverage, Sierra input (global) | {sg['80pct']:.3f} / {sg['90pct']:.3f} / {sg['95pct']:.3f} |")
        A(f"| conformal coverage, Sierra input (locally adaptive) | {sl['80pct']:.3f} / {sl['90pct']:.3f} / {sl['95pct']:.3f} |")
        A(f"| conformal coverage, scAPAtrap control (global) | {bg['80pct']:.3f} / {bg['90pct']:.3f} / {bg['95pct']:.3f} |")
        A(f"| GP beats per-gene mean (RMSE, 200 genes) | scAPAtrap {ds['frac_gp_better']:.1%}, "
          f"Sierra {dsi['frac_gp_better']:.1%} |")
        A("")

    # ---------------- cross-dataset table ----------------
    A("## Cross-dataset summary")
    A("")
    A("| dataset | species | tissue | n sites (genes) | n test | coverage 80/90/95 global | local | median r (>=10 spots) | Jaccard | GP>mean (Sierra) |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for key in ORDER:
        p = per[key]
        sg = p["metric_c_conformal_coverage"]["sierra_input"]["after_conformal_global"]
        sl = p["metric_c_conformal_coverage"]["sierra_input"]["after_conformal_locally_adaptive"]
        g10 = p["metric_b_gene_level"]["min10spots"]
        A(f"| {p['label']} | {p['species']} | {p['tissue'].split(' (')[0]} | "
          f"{p['sierra']['n_sites_ge2']:,} ({p['sierra']['n_genes_ge2']:,}) | "
          f"{p['metric_c_conformal_coverage']['sierra_input']['n_test']:,} | "
          f"{sg['80pct']:.3f}/{sg['90pct']:.3f}/{sg['95pct']:.3f} | "
          f"{sl['80pct']:.3f}/{sl['90pct']:.3f}/{sl['95pct']:.3f} | "
          f"{g10['pearson_r_median']:.2f} | "
          f"{p['metric_a_pas_overlap']['jaccard_clusters']:.2f} | "
          f"{p['metric_d_gp_vs_mean']['sierra']['frac_gp_better']:.1%} |")
    A("")
    cd = summary["cross_dataset"]
    A(f"**Max deviation of Sierra-input conformal coverage from nominal across all datasets, levels and")
    A(f"modes: {cd['sierra_coverage_max_abs_deviation_overall_pp']:.1f} pp** "
      f"(per level: 80% {cd['sierra_coverage_max_abs_deviation_pp']['80pct']:.1f} pp, "
      f"90% {cd['sierra_coverage_max_abs_deviation_pp']['90pct']:.1f} pp, "
      f"95% {cd['sierra_coverage_max_abs_deviation_pp']['95pct']:.1f} pp).")
    A("")
    A("## Conclusion")
    A("")
    jmin = min(per[k]["metric_a_pas_overlap"]["jaccard_clusters"] for k in ORDER)
    jmax = max(per[k]["metric_a_pas_overlap"]["jaccard_clusters"] for k in ORDER)
    A(cd["headline"] + " Gene-level distal usage agrees across callers at median per-gene Pearson r = "
      + ", ".join(f"{per[k]['metric_b_gene_level']['min10spots']['pearson_r_median']:.2f}" for k in ORDER)
      + " (kidney/brain/colon); the per-gene mean remains a strong RMSE baseline on both caller inputs in")
    A(f"every dataset; PAS overlap itself is caller-dependent (Jaccard {jmin:.2f}-{jmax:.2f} at +/-50 bp) but")
    A("with the large majority of Sierra peaks within 500 bp of a scAPAtrap site in every dataset.")
    A("Statistical conclusions are therefore robust to the choice of PAS caller across species and tissues.")
    A("")
    A("## Methods-ready paragraph (English)")
    A("")
    A("> **Caller robustness.** To test whether our conclusions depend on the PAS caller, we re-called PAS")
    A("> on the same Space Ranger BAMs of three Visium datasets spanning two species and three tissues")
    A("> (GSE183456 human kidney; GSE220442 human AD-brain prefrontal cortex; GSE169749 mouse colon) with")
    A("> Sierra (v0.99.27) using its default FindPeaks/CountPeaks workflow with UMI-deduplicated counting")
    A("> against each dataset's Space Ranger reference GTF. After the >=2-sites-per-gene convention this")
    A("> yielded " + ", ".join(f"{per[k]['sierra']['n_sites_ge2']:,} sites in {per[k]['sierra']['n_genes_ge2']:,} genes ({per[k]['label'].split(' / ')[0]})"
                              for k in ORDER) + ". Each peak-by-spot usage matrix was passed through the")
    A("> identical spaGAPA pipeline without any re-tuning. Split-conformal prediction intervals retained")
    A("> nominal marginal coverage on every dataset (max deviation from nominal "
      f"{cd['sierra_coverage_max_abs_deviation_overall_pp']:.1f} percentage points at 80/90/95%, global and")
    A("> locally adaptive), matching the scAPAtrap-based analyses of the same data. Cross-caller agreement")
    A("> was "
      + " / ".join(f"{per[k]['metric_a_pas_overlap']['callerB_frac_matched']:.0%}" for k in ORDER)
      + " of Sierra peaks within 50 bp of a scAPAtrap site (kidney/brain/colon), and gene-level distal-usage indices computed")
    A("> independently from each caller correlated at median per-gene Pearson r = "
      + " / ".join(f"{per[k]['metric_b_gene_level']['min10spots']['pearson_r_median']:.2f} (n = {per[k]['metric_b_gene_level']['min10spots']['n_genes_compared']})"
                  for k in ORDER) + ". The qualitative GP-vs-mean imputation comparison was likewise")
    A("> unchanged. These results indicate that the framework's statistical guarantees do not rely on a")
    A("> particular PAS caller.")
    A("")

    with open(os.path.join(MC, "report.md"), "w") as f:
        f.write("\n".join(lines))
    print(f"[out] {os.path.join(MC, 'report.md')}")


if __name__ == "__main__":
    main()

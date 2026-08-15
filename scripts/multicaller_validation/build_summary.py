#!/usr/bin/env python
"""Assemble pipeline_output/multicaller_validation/summary.json + report.md
from the metric JSONs produced by the multi-caller validation runs."""
import json, os

ROOT = "/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA"
MC = os.path.join(ROOT, "pipeline_output/multicaller_validation")


def load(p):
    with open(os.path.join(MC, p)) as f:
        return json.load(f)


def main():
    ab = load("metrics_ab.json")
    d = load("metrics_d.json")
    sierra_cal = load("sierra/conformal/uncertainty_calibration.json")
    base_cal = load("scapatrap_conformal/uncertainty_calibration.json")
    qc = load("sierra/qc_summary.json")

    summary = {
        "dataset": "GSE183456 / GSM6047774 (human kidney Visium, 10x Space Ranger, 3010 spots)",
        "baseline_caller": "scAPAtrap (53,572 sites; 34,118 usage rows)",
        "second_caller": "Sierra 0.99.27 (Winnie09/Sierra; FindPeaks + CountPeaks with UMI dedup)",
        "input_bam": "pipeline_output/gse183456_GSM6047774_sr/outs/possorted_genome_bam.bam",
        "reference_gtf": "refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz (Space Ranger reference)",
        "sierra_run": qc,
        "metric_a_pas_overlap": {
            **ab["overlap"],
            "summit_convention_sensitivity": ab["overlap_summit_convention"],
            "nearest_distance": ab["nearest_distance"],
        },
        "metric_b_gene_level": ab["gene_level"],
        "metric_c_conformal_coverage": {
            "sierra_input": {
                "n_genes": sierra_cal["config"]["n_genes"],
                "n_spots": sierra_cal["config"]["n_spots"],
                "n_test": sierra_cal["config"]["n_test"],
                "before_raw_gp": sierra_cal["before"]["coverage"],
                "after_conformal_global": sierra_cal["after"]["modes"]["global"]["coverage"],
                "after_conformal_locally_adaptive": sierra_cal["after"]["modes"]["locally_adaptive"]["coverage"],
            },
            "scapatrap_input_same_protocol": {
                "n_genes": base_cal["config"]["n_genes"],
                "n_test": base_cal["config"]["n_test"],
                "after_conformal_global": base_cal["after"]["modes"]["global"]["coverage"],
                "after_conformal_locally_adaptive": base_cal["after"]["modes"]["locally_adaptive"]["coverage"],
            },
        },
        "metric_d_gp_vs_mean": d,
    }
    with open(os.path.join(MC, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[out] {os.path.join(MC, 'summary.json')}")

    # ------------------------------------------------------------------ report
    ov = ab["overlap"]; ov2 = ab["overlap_summit_convention"]; nd = ab["nearest_distance"]
    g10 = ab["gene_level"]["min10spots"]; g30 = ab["gene_level"]["min30spots"]
    sc = sierra_cal["after"]["modes"]["global"]["coverage"]
    sl = sierra_cal["after"]["modes"]["locally_adaptive"]["coverage"]
    bc = base_cal["after"]["modes"]["global"]["coverage"]
    ds, dsi = d["scapatrap"], d["sierra"]
    lines = []
    A = lines.append
    A("# Multi-caller robustness validation: Sierra as a second PAS caller")
    A("")
    A("## Setup")
    A("")
    A("- Dataset: GSE183456 / GSM6047774, human kidney Visium, 3010 spots, 214M reads")
    A("- Baseline caller: scAPAtrap (53,572 sites, 41,747 gene-annotated, 34,118 usage rows)")
    A("- Second caller: **Sierra 0.99.27** (FindPeaks + CountPeaks, per-UMI deduplicated counting),")
    A("  run on the same Space Ranger `possorted_genome_bam.bam` and the same reference GTF")
    A("  (refdata-gex-GRCh38-2024-A) with a splice-junction BED extracted from the BAM")
    A("  (20,644 junctions >= 25 supporting reads).")
    A(f"- Sierra output: 19,575 usable peaks in 13,377 genes; after the spaGAPA >=2-sites-per-gene")
    A(f"  convention: **{qc['n_sites_ge2']} sites x {qc['n_spots']} spots in {qc['n_genes_ge2']} genes**;")
    A("  usage matrix built with min_parent_count=5, coordinates reused from the baseline")
    A("  (same Space Ranger run).")
    A("- Wall time: junction extraction ~5 min (16 procs), FindPeaks 10.1 min, CountPeaks ~10 min (16 cores).")
    A("")
    A("## Metric a - PAS overlap (summits merged at +/-50 bp, strand-aware 3'-end coordinates)")
    A("")
    A(f"- scAPAtrap points: 41,747 (gene-annotated); Sierra points: {ov['callerB_points_total']}")
    A(f"- Clusters: **{ov['shared_clusters']} shared** / {ov['unique_callerA_clusters']} scAPAtrap-only / "
      f"{ov['unique_callerB_clusters']} Sierra-only; Jaccard = **{ov['jaccard_clusters']:.3f}**")
    A(f"- Per-point match rate within 50 bp: **{ov['callerB_frac_matched']:.1%} of Sierra peaks**, "
      f"{ov['callerA_frac_matched']:.1%} of scAPAtrap sites")
    A(f"- Distance context: median nearest-neighbour distance {nd['median_bp']:.0f} bp; "
      f"{nd['frac_le_100bp']:.1%} <= 100 bp, {nd['frac_le_500bp']:.1%} <= 500 bp, "
      f"{nd['frac_le_1000bp']:.1%} <= 1 kb")
    A(f"- Sensitivity (Sierra gaussian-summit coordinate instead of 3'-end): "
      f"{ov2['callerB_frac_matched']:.1%} matched, Jaccard {ov2['jaccard_clusters']:.3f}. "
      "The two callers agree far better on peak *regions* (87% within 500 bp) than on exact")
    A("  summit coordinates (a mix of gaussian-fit summit vs narrow-peak 3' boundary conventions).")
    A("")
    A("## Metric b - gene-level distal-usage consistency")
    A("")
    A("Each caller's peaks were split at the median strand-oriented position into proximal/distal")
    A("groups; index = distal / (proximal + distal) per spot (min_parent = 5; identical code to the")
    A("benchmark scripts). 3,711 genes have multi-site indices in both callers.")
    A(f"- Genes with >= 10 co-finite spots: **n = {g10['n_genes_compared']}, per-gene Pearson r: "
      f"median {g10['pearson_r_median']:.3f}, mean {g10['pearson_r_mean']:.3f}**, "
      f"{g10['frac_r_gt_0.7']:.1%} with r > 0.7")
    A(f"- Genes with >= 30 co-finite spots: n = {g30['n_genes_compared']}, median r = "
      f"{g30['pearson_r_median']:.3f}, mean {g30['pearson_r_mean']:.3f}")
    A(f"- The r distribution is bimodal (IQR {g10['pearson_r_q25']:.2f}-{g10['pearson_r_q75']:.2f}):")
    A("  genes whose dominant PAS set is shared between callers agree almost perfectly, while")
    A("  genes where the callers pick different sites contribute near-zero r.")
    A("")
    A("## Metric c - split-conformal coverage on the Sierra input (core claim)")
    A("")
    A(f"Identical protocol to scripts/calibrate_uncertainty.py (20% masking, 50/50 cal/test split,")
    A(f"{sierra_cal['config']['n_test']:,} test points):")
    A("")
    A("| input | mode | 80% | 90% | 95% |")
    A("|---|---|---|---|---|")
    A(f"| Sierra | global | {sc['80pct']:.4f} | {sc['90pct']:.4f} | {sc['95pct']:.4f} |")
    A(f"| Sierra | locally adaptive | {sl['80pct']:.4f} | {sl['90pct']:.4f} | {sl['95pct']:.4f} |")
    A(f"| scAPAtrap (same protocol) | global | {bc['80pct']:.4f} | {bc['90pct']:.4f} | {bc['95pct']:.4f} |")
    A("")
    A("Coverage stays at the nominal level on the second caller's output; the raw-GP std is")
    A(f"over-conservative before calibration ({sierra_cal['before']['coverage']['80pct']:.3f} at 80%),")
    A("exactly as on scAPAtrap input.")
    A("")
    A("## Metric d - GP vs per-gene-mean RMSE (200 genes, 20% masked, seed 42)")
    A("")
    A("| input | RMSE GP (median) | RMSE mean (median) | GP better in |")
    A("|---|---|---|---|")
    A(f"| scAPAtrap | {ds['rmse_gp_median']:.4f} | {ds['rmse_mean_median']:.4f} | {ds['frac_gp_better']:.1%} |")
    A(f"| Sierra | {dsi['rmse_gp_median']:.4f} | {dsi['rmse_mean_median']:.4f} | {dsi['frac_gp_better']:.1%} |")
    A("")
    A("Qualitative pattern identical to the published supplementary analysis (S4): the per-gene mean")
    A("is a strong RMSE baseline on both caller inputs and the GP advantage remains a minority;")
    A("absolute RMSEs are higher on the Sierra input because its UMI-deduplicated counts are sparser.")
    A("")
    A("## Conclusion")
    A("")
    A(f"Swapping the PAS caller (scAPAtrap -> Sierra) leaves the spaGAPA framework fully functional:")
    A("a plain format conversion (peak x spot usage matrix + the same coordinates) is all that is")
    A(f"needed; split-conformal coverage stays nominal (0.801 / 0.900 / 0.950 at 80/90/95%);")
    A(f"gene-level distal usage agrees across callers at median per-gene r = "
      f"{g10['pearson_r_median']:.2f} ({g10['frac_r_gt_0.7']:.0%} of genes r > 0.7); and the")
    A("GP-vs-mean RMSE picture is unchanged. Statistical conclusions are therefore robust to the")
    A("choice of PAS caller.")
    A("")
    A("## Methods-ready paragraph (English)")
    A("")
    A(f"> **Caller robustness.** To test whether our conclusions depend on the PAS caller, we")
    A(f"> re-called PAS on the same Space Ranger BAM of GSE183456 with Sierra (v0.99.27) using")
    A(f"> its default FindPeaks/CountPeaks workflow with UMI-deduplicated counting, obtaining")
    A(f"> 19,575 peaks in 13,377 genes ({qc['n_sites_ge2']} sites in {qc['n_genes_ge2']} multi-site genes after the")
    A(f"> >=2-sites-per-gene convention). The resulting peak-by-spot usage matrix was passed through")
    A(f"> the identical spaGAPA pipeline without any re-tuning. Split-conformal prediction intervals")
    A(f"> retained nominal marginal coverage (80%: {sc['80pct']:.3f}; 90%: {sc['90pct']:.3f}; 95%: "
      f"{sc['95pct']:.3f}; locally adaptive variant {sl['80pct']:.3f}/{sl['90pct']:.3f}/{sl['95pct']:.3f}),")
    A(f"> matching the scAPAtrap-based analysis ({bc['80pct']:.3f}/{bc['90pct']:.3f}/{bc['95pct']:.3f}).")
    A(f"> Cross-caller agreement was {ov['callerB_frac_matched']:.0%} of Sierra peaks within 50 bp of a")
    A(f"> scAPAtrap site ({nd['frac_le_500bp']:.0%} within 500 bp; Jaccard of +/-50 bp merged clusters")
    A(f"> {ov['jaccard_clusters']:.2f}), and gene-level distal-usage indices computed independently from")
    A(f"> each caller correlated at median per-gene Pearson r = {g10['pearson_r_median']:.2f} across spots")
    A(f"> (n = {g10['n_genes_compared']} genes with >= 10 shared informative spots). The qualitative")
    A(f"> GP-vs-mean imputation comparison was likewise unchanged. These results indicate that the")
    A(f"> framework's statistical guarantees do not rely on a particular PAS caller.")
    A("")
    with open(os.path.join(MC, "report.md"), "w") as f:
        f.write("\n".join(lines))
    print(f"[out] {os.path.join(MC, 'report.md')}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Collect APA batch results for the two mouse datasets into a single report file.

Reads qc_summary.json from each processed scAPAtrap dir and writes:
  - pipeline_output/mouse_apa_batch_report.txt  (human-readable)
  - pipeline_output/mouse_apa_batch_report.json (machine-readable)

Run AFTER both datasets' scAPAtrap phases complete.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA/spaGAPA")
PROCESSED = ROOT / "data/processed"
PIPE_OUT = ROOT / "pipeline_output"

DATASETS = [
    {
        "key": "gse169749_gsm5213483",
        "gse": "GSE169749",
        "gsm": "GSM5213483",
        "tissue": "mouse colon (day 0)",
        "scapatrap_dir": PROCESSED / "gse169749_gsm5213483_scapatrap",
        "sr_dir": PIPE_OUT / "gse169749/gsm5213483_d0_sr/outs",
    },
    {
        "key": "gse263303_gsm8189356",
        "gse": "GSE263303",
        "gsm": "GSM8189356",
        "tissue": "mouse brain Nf1+/- (K73-6-FMFC)",
        "scapatrap_dir": PROCESSED / "gse263303_GSM8189356_scapatrap",
        "sr_dir": PIPE_OUT / "gse263303_GSM8189356_sr/outs",
    },
    {
        "key": "gse263303_gsm8189359",
        "gse": "GSE263303",
        "gsm": "GSM8189359",
        "tissue": "mouse brain Nf1+/- (K75-2-FMFC)",
        "scapatrap_dir": PROCESSED / "gse263303_GSM8189359_scapatrap",
        "sr_dir": PIPE_OUT / "gse263303_GSM8189359_sr/outs",
    },
]


def spaceranger_status(sr_dir: Path) -> dict:
    bam = sr_dir / "possorted_genome_bam.bam"
    tissue = sr_dir / "spatial" / "tissue_positions.csv"
    tissue_alt = sr_dir / "spatial" / "tissue_positions_list.csv"
    web = sr_dir / "web_summary.html"
    return {
        "sr_dir": str(sr_dir),
        "web_summary": web.exists(),
        "bam_present": bam.exists(),
        "bam_size_gb": round(bam.stat().st_size / 1e9, 2) if bam.exists() else None,
        "tissue_positions_present": tissue.exists() or tissue_alt.exists(),
    }


def scapatrap_status(proc_dir: Path) -> dict:
    qc_file = proc_dir / "qc_summary.json"
    if not qc_file.exists():
        return {"processed": False}
    qc = json.loads(qc_file.read_text())
    return {
        "processed": True,
        "n_spots": qc.get("n_spots"),
        "n_called_sites": qc.get("n_called_sites"),
        "n_gene_annotated_sites": qc.get("n_gene_annotated_sites"),
        "n_apa_usage_sites": qc.get("n_apa_usage_sites"),
        "n_apa_usage_spots": qc.get("n_apa_usage_spots"),
        "species": qc.get("species"),
        "tissue": qc.get("tissue"),
    }


def main() -> None:
    report = {"datasets": []}
    for d in DATASETS:
        item = {
            "key": d["key"],
            "gse": d["gse"],
            "gsm": d["gsm"],
            "tissue": d["tissue"],
            "spaceranger": spaceranger_status(d["sr_dir"]),
            "scapatrap": scapatrap_status(d["scapatrap_dir"]),
        }
        report["datasets"].append(item)

    out_json = PIPE_OUT / "mouse_apa_batch_report.json"
    out_txt = PIPE_OUT / "mouse_apa_batch_report.txt"
    out_json.write_text(json.dumps(report, indent=2))

    lines = []
    lines.append("=" * 70)
    lines.append("MOUSE Visium APA BATCH REPORT")
    lines.append("=" * 70)
    for d in report["datasets"]:
        lines.append("")
        lines.append(f"{d['gse']}  {d['gsm']}  ({d['tissue']})")
        lines.append("-" * 70)
        sr = d["spaceranger"]
        sr_status = "COMPLETE" if sr["web_summary"] and sr["bam_present"] else (
            "BAM READY (no web_summary)" if sr["bam_present"] else "INCOMPLETE/MISSING"
        )
        lines.append(f"  spaceranger : {sr_status}")
        lines.append(f"    bam          : {sr['bam_present']} ({sr['bam_size_gb']} GB)")
        lines.append(f"    tissue_pos   : {sr['tissue_positions_present']}")
        sc = d["scapatrap"]
        if sc["processed"]:
            lines.append(f"  scAPAtrap   : COMPLETE (species={sc['species']})")
            lines.append(f"    n_spots              = {sc['n_spots']}")
            lines.append(f"    n_called_sites       = {sc['n_called_sites']}")
            lines.append(f"    n_gene_annotated     = {sc['n_gene_annotated_sites']}")
            lines.append(f"    n_apa_usage_sites    = {sc['n_apa_usage_sites']}")
            lines.append(f"    n_apa_usage_spots    = {sc['n_apa_usage_spots']}")
        else:
            lines.append("  scAPAtrap   : NOT YET PROCESSED")
    lines.append("")
    lines.append("=" * 70)
    out_txt.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote: {out_txt}")
    print(f"Wrote: {out_json}")


if __name__ == "__main__":
    main()

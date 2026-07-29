#!/usr/bin/env python3
"""Run scAPAtrap on GSE237183 Visium Space Ranger BAMs (wrapper).

This is a thin dataset-specific driver that wraps the existing, reusable
launcher ``scripts/run_scapatrap_spaceranger.py``. It does NOT modify
scAPAtrap or the shared launcher; it only:

* maps a GSE237183 GSM -> its Space Ranger BAM / spatial dir,
* fills in the spaGAPA metadata (species/tissue/source-label),
* resolves server-specific env paths per CLAUDE.md (S90/S91/S97/S98),
* invokes ``run_scapatrap_spaceranger.py`` once per GSM (sequentially),

Outputs land at:
    spaGAPA/pipeline_output/gse237183_<GSM>_scapatrap/      (raw + stage + logs)
    spaGAPA/data/processed/gse237183_<gsm>_scapatrap/       (apa_matrix.csv etc.)

GSM7596588 is permanently SKIPPED (fiducial alignment failure).

Usage:
    python run_scapatrap_gse237183.py --gsm GSM7596590
    python run_scapatrap_gse237183.py --gsm GSM7596590 GSM7596601
    python run_scapatrap_gse237183.py --all            # all 18 valid samples
    python run_scapatrap_gse237183.py --gsm GSM7596590 --dry-run
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PACKAGE_ROOT / "scripts"
SR_COUNTS_DIR = PACKAGE_ROOT / "pipeline_output/gse237183_spaceranger_counts"
LAUNCHER = SCRIPTS_DIR / "run_scapatrap_spaceranger.py"
GTF = "/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"

# ---------------------------------------------------------------------------
# Sample metadata (GSE237183 — glioma Visium cohort)
# GSM7596588 is SKIPPED (fiducial alignment failure; partial data deleted).
# All samples: species=human, tissue=glioma.
# ---------------------------------------------------------------------------
ALL_GSMS = [
    "GSM7596587", "GSM7596589", "GSM7596590", "GSM7596591", "GSM7596592",
    "GSM7596593", "GSM7596594", "GSM7596595", "GSM7596596", "GSM7596597",
    "GSM7596598", "GSM7596599", "GSM7596600", "GSM7596601", "GSM7596602",
    "GSM7596603", "GSM7596604", "GSM7596605",
]
SKIPPED_GSMS = {"GSM7596588"}
SPECIES = "human"
TISSUE = "glioma"


def sr_dir(gsm: str) -> Path:
    return SR_COUNTS_DIR / f"gse237183_{gsm}_sr"


def bam_path(gsm: str) -> Path:
    return sr_dir(gsm) / "outs/possorted_genome_bam.bam"


def spatial_dir(gsm: str) -> Path:
    return sr_dir(gsm) / "outs/spatial"


def out_root(gsm: str) -> Path:
    return PACKAGE_ROOT / f"pipeline_output/gse237183_{gsm}_scapatrap"


def processed_dir(gsm: str) -> Path:
    return PACKAGE_ROOT / f"data/processed/gse237183_{gsm.lower()}_scapatrap"


def dataset_name(gsm: str) -> str:
    return f"gse237183_{gsm.lower()}_scapatrap"


def source_label(gsm: str) -> str:
    return f"GSE237183 {gsm} Space Ranger BAM + scAPAtrap"


# ---------------------------------------------------------------------------
# Server detection (CLAUDE.md)
# ---------------------------------------------------------------------------
def server_env() -> dict[str, str]:
    import socket
    host = socket.gethostname().split(".")[0].upper()
    base = {
        "S90": dict(RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript",
                    TOOL_PREFIX="/s1/mengzijun/anaconda3/envs/samtools/bin",
                    TMPDIR="/s2/mengzijun/tmp"),
        "S91": dict(RSCRIPT=str(Path.home() / "anaconda3/envs/r442/bin/Rscript"),
                    TOOL_PREFIX=str(Path.home() / "anaconda3/envs/samtools/bin"),
                    TMPDIR="/s3/mengzijun/tmp"),
        "S97": dict(RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript",
                    TOOL_PREFIX="/s1/mengzijun/anaconda3/envs/samtools/bin",
                    TMPDIR="/s972/mengzijun/tmp"),
        "S98": dict(RSCRIPT="/s1/mengzijun/anaconda3/envs/r442/bin/Rscript",
                    TOOL_PREFIX="/s1/mengzijun/anaconda3/envs/samtools/bin",
                    TMPDIR="/s982/mengzijun/tmp"),
    }
    if host not in base:
        raise SystemExit(f"错误：未知服务器 {host}，请检查环境配置 (expected S90/S91/S97/S98).")
    env = base[host]
    env["R_LIBS"] = "/s1/SHARE/01_software/R_442_SeuratV5/library"
    env["HOST"] = host
    return env


def validate_inputs(gsm: str) -> None:
    if gsm in SKIPPED_GSMS:
        raise SystemExit(f"GSM {gsm} is SKIPPED (fiducial alignment failure). Do not use.")
    bam = bam_path(gsm)
    spatial = spatial_dir(gsm)
    if not bam.exists():
        raise SystemExit(f"BAM not found for {gsm}: {bam}")
    if not (bam.with_suffix(bam.suffix + ".bai")).exists():
        # scAPAtrap reindexes if needed, but warn loudly
        print(f"WARNING: .bai index missing for {gsm}: {bam}.bai", file=sys.stderr)
    if not spatial.exists():
        raise SystemExit(f"spatial dir not found for {gsm}: {spatial}")
    tps = [spatial / "tissue_positions.csv", spatial / "tissue_positions_list.csv"]
    if not any(p.exists() for p in tps):
        raise SystemExit(f"tissue_positions file not found for {gsm} in {spatial}")


def run_one(gsm: str, env: dict, *, threads: int, force: bool, dry_run: bool,
            postprocess_only: bool) -> int:
    """Invoke run_scapatrap_spaceranger.py for a single GSM. Returns exit code."""
    validate_inputs(gsm)

    cmd = [
        sys.executable, str(LAUNCHER),
        "--bam", str(bam_path(gsm)),
        "--spatial-dir", str(spatial_dir(gsm)),
        "--gtf", GTF,
        "--output-root", str(out_root(gsm)),
        "--processed-dir", str(processed_dir(gsm)),
        "--dataset-name", dataset_name(gsm),
        "--source-label", source_label(gsm),
        "--species", SPECIES,
        "--tissue", TISSUE,
        "--rscript", env["RSCRIPT"],
        "--r-lib", env["R_LIBS"],
        "--tool-prefix", env["TOOL_PREFIX"],
        "--threads", str(threads),
    ]
    if force:
        cmd.append("--force")
    if dry_run:
        cmd.append("--dry-run")
    if postprocess_only:
        cmd.append("--postprocess-only")

    print(f"\n{'='*70}\n[{env['HOST']}] scAPAtrap -> {gsm}\n{'='*70}", flush=True)
    print("CMD: " + " ".join(cmd), flush=True)
    result = subprocess.run(cmd, cwd=str(PACKAGE_ROOT))
    return result.returncode


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--gsm", nargs="+", help="One or more GSE237183 GSM accessions.")
    g.add_argument("--all", action="store_true", help="Process all 18 valid GSE237183 samples.")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--force", action="store_true", help="Overwrite existing raw/processed dirs.")
    ap.add_argument("--dry-run", action="store_true", help="Validate + write manifest, do not run R.")
    ap.add_argument("--postprocess-only", action="store_true",
                    help="Skip scAPAtrap; export processed files from existing raw_scapatrap.")
    args = ap.parse_args()

    env = server_env()
    os.makedirs(env["TMPDIR"], exist_ok=True)

    gsms = ALL_GSMS if args.all else args.gsm
    # de-duplicate, preserve order
    seen = set()
    gsms = [g for g in gsms if not (g in seen or seen.add(g))]

    for g in gsms:
        if g in SKIPPED_GSMS:
            print(f"SKIP {g} (skipped sample)", flush=True)
            continue
        rc = run_one(g, env, threads=args.threads, force=args.force,
                     dry_run=args.dry_run, postprocess_only=args.postprocess_only)
        if rc != 0:
            print(f"ERROR: {g} failed (exit {rc}). Aborting sequential run.", file=sys.stderr)
            sys.exit(rc)
        print(f"DONE: {g} -> {processed_dir(g)}", flush=True)


if __name__ == "__main__":
    main()

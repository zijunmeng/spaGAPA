#!/usr/bin/env python3
"""Prepare and run parallel Space Ranger dry-runs for GSE237183."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/s1/SHARE/mengzijun/01_project/26_spaGAPA")
PACKAGE_ROOT = ROOT / "spaGAPA"
DATA_DIR = ROOT / "data/raw/gse237183"
FASTQ_DIR = DATA_DIR / "fastq_ena"
LOG_DIR = DATA_DIR / "logs"
RAW_TAR = DATA_DIR / "geo_supplementary/GSE237183_RAW.tar"
SPACERANGER = Path("/s2/mengzijun/pkg/spaceranger-4.1.0/spaceranger")
TRANSCRIPTOME = Path("/s1/SHARE/00_ref_genecode/refdata-gex-GRCh38-2024-A")
OUT_ROOT = PACKAGE_ROOT / "pipeline_output/gse237183_dryruns"


@dataclass(frozen=True)
class Sample:
    gsm: str
    title: str
    characteristics: str
    srx: str
    image_member: str
    runs: tuple[str, ...]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_samples() -> list[Sample]:
    sample_rows = read_tsv(LOG_DIR / "sample_manifest.tsv")
    run_rows = read_tsv(LOG_DIR / "ena_run_manifest.tsv")

    runs_by_srx: dict[str, list[str]] = {}
    for row in run_rows:
        runs_by_srx.setdefault(row["experiment_accession"], []).append(row["run_accession"])

    with tarfile.open(RAW_TAR) as tar:
        members = tar.getnames()
    image_by_gsm = {
        member.split("_", 1)[0]: member
        for member in members
        if member.endswith("_detected_tissue_image.jpg.gz")
    }

    samples: list[Sample] = []
    for row in sample_rows:
        gsm = row["gsm"]
        srx_values = [item for item in row["srx_accessions"].split(";") if item]
        if len(srx_values) != 1:
            raise ValueError(f"{gsm} should have exactly one SRX, got {srx_values}")
        srx = srx_values[0]
        runs = tuple(sorted(runs_by_srx.get(srx, [])))
        if not runs:
            raise ValueError(f"No runs found for {gsm}/{srx}")
        image_member = image_by_gsm.get(gsm)
        if not image_member:
            raise ValueError(f"No detected_tissue_image member found for {gsm}")
        samples.append(
            Sample(
                gsm=gsm,
                title=row["title"],
                characteristics=row["characteristics"],
                srx=srx,
                image_member=image_member,
                runs=runs,
            )
        )
    return samples


def ensure_inputs(sample: Sample) -> tuple[Path, Path]:
    sample_root = DATA_DIR / "spaceranger_inputs" / sample.gsm
    fastq_out = sample_root / "fastqs"
    image_out = sample_root / "detected_tissue_image.jpg"
    fastq_out.mkdir(parents=True, exist_ok=True)

    for lane_index, run in enumerate(sample.runs, start=1):
        for read_index in (1, 2):
            src = FASTQ_DIR / f"{run}_{read_index}.fastq.gz"
            if not src.exists():
                raise FileNotFoundError(src)
            dst = fastq_out / f"{sample.gsm}_S1_L{lane_index:03d}_R{read_index}_001.fastq.gz"
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(src, dst)

    if not image_out.exists() or image_out.stat().st_size == 0:
        with tarfile.open(RAW_TAR) as tar:
            extracted = tar.extractfile(sample.image_member)
            if extracted is None:
                raise FileNotFoundError(sample.image_member)
            with gzip.GzipFile(fileobj=extracted) as gz_in, image_out.open("wb") as out:
                shutil.copyfileobj(gz_in, out)
    return fastq_out, image_out


def run_dry(sample: Sample, threads: int, mem_gb: int, force: bool) -> dict[str, str]:
    fastq_dir, image = ensure_inputs(sample)
    sample_out = OUT_ROOT / sample.gsm
    work_dir = sample_out / "work"
    log_dir = sample_out / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    dry_id = f"gse237183_{sample.gsm}_sr"
    dry_pipestance = work_dir / dry_id
    dry_log = log_dir / "spaceranger_dry_run.log"
    metadata = log_dir / "sample_metadata.tsv"
    metadata.write_text(
        "field\tvalue\n"
        f"gsm\t{sample.gsm}\n"
        f"title\t{sample.title}\n"
        f"characteristics\t{sample.characteristics}\n"
        f"srx\t{sample.srx}\n"
        f"runs\t{','.join(sample.runs)}\n"
        f"fastq_dir\t{fastq_dir}\n"
        f"image\t{image}\n",
        encoding="utf-8",
    )
    if dry_pipestance.exists() and force:
        shutil.rmtree(dry_pipestance)

    cmd = [
        str(SPACERANGER),
        "count",
        "--id",
        dry_id,
        "--description",
        f"GSE237183 {sample.gsm} Visium dry-run",
        "--transcriptome",
        str(TRANSCRIPTOME),
        "--fastqs",
        str(fastq_dir),
        "--sample",
        sample.gsm,
        "--image",
        str(image),
        "--unknown-slide",
        "visium-1",
        "--create-bam",
        "true",
        "--localcores",
        str(threads),
        "--localmem",
        str(mem_gb),
        "--disable-cell-annotation",
        "--disable-ui",
        "--dry",
    ]

    start = time.time()
    with dry_log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\t" + " ".join(cmd) + "\n")
        handle.flush()
        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.time() - start
    return {
        "gsm": sample.gsm,
        "title": sample.title,
        "srx": sample.srx,
        "runs": ",".join(sample.runs),
        "returncode": str(proc.returncode),
        "status": "PASS" if proc.returncode == 0 else "FAIL",
        "elapsed_s": f"{elapsed:.2f}",
        "log": str(dry_log),
    }


def write_summary(rows: list[dict[str, str]]) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = OUT_ROOT / "dryrun_summary.tsv"
    fields = ["gsm", "title", "srx", "runs", "returncode", "status", "elapsed_s", "log"]
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Summary: {summary}")
    for row in rows:
        print(
            f"{row['status']}\t{row['gsm']}\treturncode={row['returncode']}\t"
            f"elapsed_s={row['elapsed_s']}\tlog={row['log']}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--mem-gb", type=int, default=24)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--samples", nargs="*", help="Optional GSM subset.")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()

    if not SPACERANGER.exists():
        raise FileNotFoundError(SPACERANGER)
    if not TRANSCRIPTOME.exists():
        raise FileNotFoundError(TRANSCRIPTOME)
    if not RAW_TAR.exists():
        raise FileNotFoundError(RAW_TAR)

    samples = load_samples()
    if args.samples:
        wanted = set(args.samples)
        samples = [sample for sample in samples if sample.gsm in wanted]
        missing = wanted - {sample.gsm for sample in samples}
        if missing:
            raise ValueError(f"Requested samples not found: {sorted(missing)}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Samples: {len(samples)}")
    print(f"Jobs: {args.jobs}")
    for sample in samples:
        ensure_inputs(sample)
    if args.prepare_only:
        print("Prepared inputs only.")
        return

    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_dry, sample, args.threads, args.mem_gb, args.force): sample
            for sample in samples
        }
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_summary(sorted(rows, key=lambda item: item["gsm"]))

    rows.sort(key=lambda item: item["gsm"])
    write_summary(rows)
    failed = [row for row in rows if row["status"] != "PASS"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

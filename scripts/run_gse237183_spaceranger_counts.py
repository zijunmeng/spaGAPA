#!/usr/bin/env python3
"""Run Space Ranger count for all or selected GSE237183 samples."""

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
OUT_ROOT = PACKAGE_ROOT / "pipeline_output/gse237183_spaceranger_counts"
SUMMARY = OUT_ROOT / "count_summary.tsv"


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
        srx_values = [item for item in row["srx_accessions"].split(";") if item]
        if len(srx_values) != 1:
            raise ValueError(f"{row['gsm']} should have one SRX, got {srx_values}")
        srx = srx_values[0]
        runs = tuple(sorted(runs_by_srx.get(srx, [])))
        if not runs:
            raise ValueError(f"No runs found for {row['gsm']}/{srx}")
        image_member = image_by_gsm.get(row["gsm"])
        if not image_member:
            raise ValueError(f"No detected_tissue_image member found for {row['gsm']}")
        samples.append(
            Sample(
                gsm=row["gsm"],
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
                if dst.is_symlink() and Path(os.readlink(dst)) == src:
                    continue
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


def output_status(pipestance: Path) -> str:
    if (pipestance / "outs/possorted_genome_bam.bam").exists():
        return "DONE"
    if (pipestance / "_log").exists() or pipestance.exists():
        return "EXISTS"
    return "NEW"


def run_count(
    sample: Sample,
    threads: int,
    mem_gb: int,
    force: bool,
    resume_existing: bool,
) -> dict[str, str]:
    fastq_dir, image = ensure_inputs(sample)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    log_root = OUT_ROOT / "_logs" / sample.gsm
    log_root.mkdir(parents=True, exist_ok=True)
    count_id = f"gse237183_{sample.gsm}_sr"
    pipestance = OUT_ROOT / count_id
    count_log = log_root / "spaceranger_count.log"
    metadata = log_root / "sample_metadata.tsv"
    metadata.write_text(
        "field\tvalue\n"
        f"gsm\t{sample.gsm}\n"
        f"title\t{sample.title}\n"
        f"characteristics\t{sample.characteristics}\n"
        f"srx\t{sample.srx}\n"
        f"runs\t{','.join(sample.runs)}\n"
        f"fastq_dir\t{fastq_dir}\n"
        f"image\t{image}\n"
        f"pipestance\t{pipestance}\n",
        encoding="utf-8",
    )

    status = output_status(pipestance)
    if status == "DONE" and not force:
        return {
            "gsm": sample.gsm,
            "title": sample.title,
            "srx": sample.srx,
            "runs": ",".join(sample.runs),
            "returncode": "0",
            "status": "SKIP_DONE",
            "elapsed_s": "0.00",
            "pipestance": str(pipestance),
            "log": str(count_log),
        }
    if status == "EXISTS" and not (force or resume_existing):
        return {
            "gsm": sample.gsm,
            "title": sample.title,
            "srx": sample.srx,
            "runs": ",".join(sample.runs),
            "returncode": "NA",
            "status": "SKIP_EXISTS",
            "elapsed_s": "0.00",
            "pipestance": str(pipestance),
            "log": str(count_log),
        }
    if pipestance.exists() and force:
        shutil.rmtree(pipestance)

    cmd = [
        str(SPACERANGER),
        "count",
        "--id",
        count_id,
        "--description",
        f"GSE237183 {sample.gsm} Visium count",
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
    ]

    start = time.time()
    with count_log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\t" + " ".join(cmd) + "\n")
        handle.flush()
        proc = subprocess.run(
            cmd,
            cwd=OUT_ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.time() - start
    final_status = "PASS" if proc.returncode == 0 else "FAIL"
    return {
        "gsm": sample.gsm,
        "title": sample.title,
        "srx": sample.srx,
        "runs": ",".join(sample.runs),
        "returncode": str(proc.returncode),
        "status": final_status,
        "elapsed_s": f"{elapsed:.2f}",
        "pipestance": str(pipestance),
        "log": str(count_log),
    }


def write_summary(rows: list[dict[str, str]]) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    fields = [
        "gsm",
        "title",
        "srx",
        "runs",
        "returncode",
        "status",
        "elapsed_s",
        "pipestance",
        "log",
    ]
    rows = sorted(rows, key=lambda item: item["gsm"])
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Summary: {SUMMARY}", flush=True)
    for row in rows:
        print(
            f"{row['status']}\t{row['gsm']}\treturncode={row['returncode']}\t"
            f"elapsed_s={row['elapsed_s']}\tpipestance={row['pipestance']}",
            flush=True,
        )


def check_only(samples: list[Sample]) -> None:
    print(f"samples={len(samples)}")
    for sample in samples:
        fastq_dir, image = ensure_inputs(sample)
        count_id = f"gse237183_{sample.gsm}_sr"
        pipestance = OUT_ROOT / count_id
        print(
            "\t".join(
                [
                    sample.gsm,
                    sample.title,
                    sample.srx,
                    ",".join(sample.runs),
                    str(fastq_dir),
                    str(image),
                    output_status(pipestance),
                    str(pipestance),
                ]
            )
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--mem-gb", type=int, default=128)
    parser.add_argument("--samples", nargs="*", help="Optional GSM subset.")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--resume-existing", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    if not SPACERANGER.exists():
        raise FileNotFoundError(SPACERANGER)
    if not TRANSCRIPTOME.exists():
        raise FileNotFoundError(TRANSCRIPTOME)
    if not RAW_TAR.exists():
        raise FileNotFoundError(RAW_TAR)
    if args.jobs < 1:
        raise ValueError("--jobs must be >= 1")

    samples = load_samples()
    if args.samples:
        wanted = set(args.samples)
        samples = [sample for sample in samples if sample.gsm in wanted]
        missing = wanted - {sample.gsm for sample in samples}
        if missing:
            raise ValueError(f"Requested samples not found: {sorted(missing)}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.check_only:
        check_only(samples)
        return

    print(f"Samples: {len(samples)}", flush=True)
    print(f"Jobs: {args.jobs}", flush=True)
    print(f"Threads per sample: {args.threads}", flush=True)
    print(f"Memory GB per sample: {args.mem_gb}", flush=True)
    print(f"Output root: {OUT_ROOT}", flush=True)

    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_count,
                sample,
                args.threads,
                args.mem_gb,
                args.force,
                args.resume_existing,
            ): sample
            for sample in samples
        }
        for future in as_completed(futures):
            rows.append(future.result())
            write_summary(rows)

    write_summary(rows)
    failed = [row for row in rows if row["status"] == "FAIL"]
    skipped_existing = [row for row in rows if row["status"] == "SKIP_EXISTS"]
    if failed or skipped_existing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

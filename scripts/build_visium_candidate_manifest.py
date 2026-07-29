#!/usr/bin/env python3
"""Build reproducible GEO/SRA/ENA manifests for candidate Visium datasets."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


def fetch(url: str, timeout: int = 90) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = response.read()
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def values(block: str, key: str) -> list[str]:
    prefix = key + " = "
    return [
        line[len(prefix) :].strip()
        for line in block.splitlines()
        if line.startswith(prefix)
    ]


def ncbi_ftp_to_https(url: str) -> str:
    return url.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")


def query_ena_for_srx(srx: str) -> list[dict[str, str]]:
    fields = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "library_strategy",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]
    query = urllib.parse.urlencode(
        {
            "result": "read_run",
            "query": f'experiment_accession="{srx}"',
            "fields": ",".join(fields),
            "format": "tsv",
        }
    )
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + query
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            data = fetch(url)
            return list(csv.DictReader(io.StringIO(data), delimiter="\t"))
        except Exception as exc:  # pragma: no cover - network retry path
            last_error = exc
            time.sleep(5 * attempt)
    raise RuntimeError(f"ENA query failed for {srx}: {last_error}")


def build_manifest(accession: str, out_dir: Path) -> None:
    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "geo_supplementary").mkdir(parents=True, exist_ok=True)
    (out_dir / "fastq_ena").mkdir(parents=True, exist_ok=True)

    soft_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{accession[:6]}nnn/{accession}/soft/{accession}_family.soft.gz"
    )
    text = fetch(soft_url)
    (log_dir / f"{accession}_family.soft").write_text(text, encoding="utf-8")
    (log_dir / "soft_url.txt").write_text(soft_url + "\n", encoding="utf-8")

    series = text.split("\n^SAMPLE", 1)[0]
    series_rows: dict[str, str | list[str]] = {
        "accession": accession,
        "soft_url": soft_url,
        "title": values(series, "!Series_title")[:1],
        "type": values(series, "!Series_type"),
        "summary": values(series, "!Series_summary"),
        "overall_design": values(series, "!Series_overall_design"),
        "relation": values(series, "!Series_relation"),
    }
    with (log_dir / "series_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["field", "value"])
        for key, value in series_rows.items():
            if isinstance(value, list):
                value = " | ".join(value)
            writer.writerow([key, value])

    supplementary_urls = [
        ncbi_ftp_to_https(item)
        for item in values(series, "!Series_supplementary_file")
        if item and item.lower() != "none"
    ]
    (log_dir / "geo_supplementary_urls.txt").write_text(
        "\n".join(supplementary_urls) + ("\n" if supplementary_urls else ""),
        encoding="utf-8",
    )

    samples = re.split(r"\n\^SAMPLE = ", text)[1:]
    sample_rows: list[list[str]] = []
    srx_accessions: list[str] = []
    for sample_block in samples:
        gsm = sample_block.split("\n", 1)[0].strip()
        title = " | ".join(values(sample_block, "!Sample_title"))
        source = " | ".join(values(sample_block, "!Sample_source_name_ch1"))
        organism = " | ".join(values(sample_block, "!Sample_organism_ch1"))
        characteristics = " | ".join(values(sample_block, "!Sample_characteristics_ch1"))
        relations = values(sample_block, "!Sample_relation")
        sample_srx = sorted(set(re.findall(r"SRX\d+", "\n".join(relations))))
        srx_accessions.extend(sample_srx)
        sample_rows.append(
            [gsm, title, organism, source, characteristics, ";".join(sample_srx)]
        )

    srx_accessions = sorted(set(srx_accessions))
    with (log_dir / "sample_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["gsm", "title", "organism", "source", "characteristics", "srx_accessions"]
        )
        writer.writerows(sample_rows)
    (log_dir / "srx_accessions.txt").write_text(
        "\n".join(srx_accessions) + ("\n" if srx_accessions else ""), encoding="utf-8"
    )

    fields = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "secondary_sample_accession",
        "library_strategy",
        "library_selection",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]
    run_rows: list[dict[str, str]] = []
    for srx in srx_accessions:
        run_rows.extend(query_ena_for_srx(srx))
    run_rows.sort(key=lambda row: (row.get("experiment_accession", ""), row.get("run_accession", "")))

    with (log_dir / "ena_run_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(run_rows)

    fastq_rows: list[list[str]] = []
    urls: list[str] = []
    total_bytes = 0
    for row in run_rows:
        fastq_field = row.get("fastq_ftp", "") or ""
        md5s = (row.get("fastq_md5", "") or "").split(";")
        sizes = (row.get("fastq_bytes", "") or "").split(";")
        for index, item in enumerate([entry for entry in fastq_field.split(";") if entry]):
            url = item if item.startswith(("http://", "https://", "ftp://")) else "https://" + item
            size = sizes[index] if index < len(sizes) else ""
            md5 = md5s[index] if index < len(md5s) else ""
            try:
                total_bytes += int(size)
            except ValueError:
                pass
            urls.append(url)
            fastq_rows.append(
                [
                    row.get("run_accession", ""),
                    row.get("experiment_accession", ""),
                    row.get("secondary_sample_accession", ""),
                    url,
                    Path(urllib.parse.urlparse(url).path).name,
                    size,
                    md5,
                ]
            )

    with (log_dir / "ena_fastq_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["run_accession", "experiment_accession", "gsm", "url", "filename", "bytes", "md5"])
        writer.writerows(fastq_rows)
    (log_dir / "ena_fastq_urls.txt").write_text(
        "\n".join(urls) + ("\n" if urls else ""), encoding="utf-8"
    )

    with (log_dir / "manifest_stats.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"accession\t{accession}\n")
        handle.write(f"samples\t{len(sample_rows)}\n")
        handle.write(f"srx\t{len(srx_accessions)}\n")
        handle.write(f"runs\t{len(run_rows)}\n")
        handle.write(f"fastq_files\t{len(urls)}\n")
        handle.write(f"fastq_total_bytes\t{total_bytes}\n")
        handle.write(f"fastq_total_gib\t{total_bytes / (1024**3):.3f}\n")
        handle.write(f"supplementary_files\t{len(supplementary_urls)}\n")

    print(
        f"{accession}: samples={len(sample_rows)} srx={len(srx_accessions)} "
        f"runs={len(run_rows)} fastq_files={len(urls)} "
        f"fastq_GiB={total_bytes / (1024**3):.2f} supp={len(supplementary_urls)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("accession")
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    build_manifest(args.accession, args.out_dir)


if __name__ == "__main__":
    main()

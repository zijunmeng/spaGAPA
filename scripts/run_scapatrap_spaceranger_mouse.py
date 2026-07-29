#!/usr/bin/env python
"""Run scAPAtrap on a Space Ranger BAM and export spaGAPA-ready files.

MOUSE-adapted fork of run_scapatrap_spaceranger.py:

* chromosome list is mouse (chr1-19 + chrX/chrY) by default, driven by
  --species (mouse|human);
* default GTF points at the 10x mouse mm10 reference;
* default species is "mouse".

This launcher is intentionally conservative:

* it uses barcode/UMI tags already present in a Space Ranger BAM;
* it writes an R script and logs for reproducibility;
* it does not delete BAM or scAPAtrap intermediate files;
* it exports site counts and a site-usage APA matrix derived from real PAS
  calling evidence.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_RSCRIPT = "~/anaconda3/envs/r442/bin/Rscript"
DEFAULT_R_LIB = "/s1/SHARE/01_software/R_442_SeuratV5/library"
DEFAULT_TOOL_PREFIX = "/home/mengzijun/anaconda3/envs/samtools/bin"
DEFAULT_BAM = (
    PACKAGE_ROOT
    / "pipeline_output/gse179572_GSM5420751_sr/outs/possorted_genome_bam.bam"
)
DEFAULT_SPATIAL = PACKAGE_ROOT / "pipeline_output/gse179572_GSM5420751_sr/outs/spatial"
DEFAULT_GTF = "/s1/SHARE/00_ref_genecode/refdata-gex-mm10-2020-A/genes/genes.gtf"
DEFAULT_OUTPUT_ROOT = PACKAGE_ROOT / "pipeline_output/gse179572_GSM5420751_scapatrap"
DEFAULT_PROCESSED_DIR = PACKAGE_ROOT / "data/processed/gse179572_gsm5420751_scapatrap"


def resolve(path: str | Path) -> Path:
    path = Path(path).expanduser()
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", default=str(DEFAULT_BAM))
    parser.add_argument("--spatial-dir", default=str(DEFAULT_SPATIAL))
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--processed-dir", default=str(DEFAULT_PROCESSED_DIR))
    parser.add_argument("--dataset-name", default="gse179572_gsm5420751_scapatrap")
    parser.add_argument(
        "--source-label",
        default="GSE179572 GSM5420751 Space Ranger BAM + scAPAtrap",
        help="Free-text source label written to metadata.csv.",
    )
    parser.add_argument("--species", default="mouse")
    parser.add_argument("--tissue", default="brain metastasis")
    parser.add_argument("--rscript", default=DEFAULT_RSCRIPT)
    parser.add_argument("--r-lib", default=DEFAULT_R_LIB)
    parser.add_argument("--tool-prefix", default=DEFAULT_TOOL_PREFIX)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--readlength", type=int, default=90)
    parser.add_argument("--cov-cutoff", type=int, default=10)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-parent-count", type=int, default=5)
    parser.add_argument("--tails-search", default="peaks", choices=["peaks", "genome", "no"])
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help="Skip scAPAtrap and export processed files from existing raw_scapatrap outputs.",
    )
    return parser.parse_args()


def ensure_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def tool_paths(tool_prefix: Path) -> dict[str, str]:
    tools = {
        "samtools": tool_prefix / "samtools",
        "umitools": tool_prefix / "umi_tools",
        "featureCounts": tool_prefix / "featureCounts",
        "star": tool_prefix / "STAR",
    }
    for name, path in tools.items():
        ensure_file(path, name)
    return {name: str(path) for name, path in tools.items()}


def human_chromosomes() -> list[str]:
    return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def mouse_chromosomes() -> list[str]:
    # Mouse: chr1-chr19, chrX, chrY (chrM is handled separately by scAPAtrap).
    return [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]


def chromosomes_for_species(species: str) -> list[str]:
    species = (species or "").strip().lower()
    if species in {"mouse", "mus_musculus", "mm", "mm10", "mm39", "grcm39", "grcm38"}:
        return mouse_chromosomes()
    return human_chromosomes()


def r_vector(values: Iterable[str]) -> str:
    quoted = ", ".join(json.dumps(str(value)) for value in values)
    return f"c({quoted})"


def write_r_script(args: argparse.Namespace, paths: dict[str, Path], tools: dict[str, str]) -> Path:
    rscript_path = paths["output_root"] / "run_scapatrap.R"
    chrs = r_vector(chromosomes_for_species(args.species))
    r_code = f"""
.libPaths(c({json.dumps(args.r_lib)}, .libPaths()))
suppressPackageStartupMessages(library(scAPAtrap))
suppressPackageStartupMessages(library(Matrix))

tools <- list(
  samtools = {json.dumps(tools["samtools"])},
  umitools = {json.dumps(tools["umitools"])},
  featureCounts = {json.dumps(tools["featureCounts"])},
  star = {json.dumps(tools["star"])}
)

trap.params <- setTrapParams(print=FALSE)
trap.params$TenX <- TRUE
trap.params$chrs <- {chrs}
trap.params$readlength <- {args.readlength}
trap.params$cov.cutoff <- {args.cov_cutoff}
trap.params$min.cells <- {args.min_cells}
trap.params$min.count <- {args.min_count}
trap.params$tails.search <- {json.dumps(args.tails_search)}
trap.params$thread <- {args.threads}

input_bam <- {json.dumps(str(paths["staged_bam"]))}
output_dir <- {json.dumps(str(paths["raw_dir"]))}
log_file <- {json.dumps(str(paths["internal_log"]))}

if (dir.exists(output_dir)) {{
  stop(paste("scAPAtrap output_dir already exists:", output_dir))
}}
if (file.exists(log_file)) {{
  stop(paste("scAPAtrap log_file already exists:", log_file))
}}

message("Running scAPAtrap")
message("input_bam=", input_bam)
message("output_dir=", output_dir)

scapatrap_rda <- scAPAtrap(
  tools = tools,
  trap.params = trap.params,
  inputBam = input_bam,
  outputDir = output_dir,
  logf = log_file,
  verbose = TRUE
)

load(scapatrap_rda)

write.csv(
  scAPAtrapData$peaks.meta,
  gzfile(file.path(output_dir, "peaks_meta.csv.gz")),
  quote = FALSE
)
write.csv(
  as.matrix(scAPAtrapData$peaks.count),
  gzfile(file.path(output_dir, "apa_site_counts.csv.gz")),
  quote = FALSE
)

qc <- list(
  scapatrap_rda = scapatrap_rda,
  n_sites = nrow(scAPAtrapData$peaks.meta),
  n_barcodes = ncol(scAPAtrapData$peaks.count),
  tails_search = trap.params$tails.search,
  readlength = trap.params$readlength,
  cov_cutoff = trap.params$cov.cutoff,
  min_cells = trap.params$min.cells,
  min_count = trap.params$min.count
)
writeLines(
  jsonlite::toJSON(qc, auto_unbox=TRUE, pretty=TRUE),
  file.path(output_dir, "scapatrap_qc.json")
)
message("scAPAtrap export complete")
"""
    rscript_path.write_text(r_code)
    return rscript_path


def load_10x_positions(spatial_dir: Path) -> pd.DataFrame:
    candidates = [spatial_dir / "tissue_positions.csv", spatial_dir / "tissue_positions_list.csv"]
    for path in candidates:
        if path.exists():
            break
    else:
        raise FileNotFoundError(f"No tissue_positions file found in {spatial_dir}")

    first_line = path.read_text().splitlines()[0]
    has_header = "barcode" in first_line.lower()
    if has_header:
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, header=None)
        df.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]

    rename = {
        "pxl_col_in_fullres": "x",
        "pxl_row_in_fullres": "y",
        "imagecol": "x",
        "imagerow": "y",
    }
    df = df.rename(columns=rename)
    if "barcode" not in df.columns:
        df = df.rename(columns={df.columns[0]: "barcode"})
    if "x" not in df.columns or "y" not in df.columns:
        if len(df.columns) >= 6:
            df = df.rename(columns={df.columns[5]: "x", df.columns[4]: "y"})
        else:
            raise ValueError(f"Cannot identify x/y columns in {path}")
    if "in_tissue" in df.columns:
        df = df[df["in_tissue"].astype(int) == 1].copy()
    return pd.DataFrame(
        {
            "spot_id": df["barcode"].astype(str),
            "spot_id_core": df["barcode"].astype(str).str.replace(r"-\d+$", "", regex=True),
            "x": pd.to_numeric(df["x"]),
            "y": pd.to_numeric(df["y"]),
        }
    )


@dataclass
class GeneRecord:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    width: int


def parse_gtf_attrs(attr: str) -> dict[str, str]:
    parsed = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def read_genes_from_gtf(gtf_path: Path) -> dict[tuple[str, str], list[GeneRecord]]:
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    genes: dict[tuple[str, str], list[GeneRecord]] = {}
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            gene_name = attrs.get("gene_name", gene_id)
            start = int(fields[3])
            end = int(fields[4])
            rec = GeneRecord(
                chrom=fields[0],
                start=start,
                end=end,
                strand=fields[6],
                gene_id=gene_id,
                gene_name=gene_name,
                width=end - start + 1,
            )
            genes.setdefault((rec.chrom, rec.strand), []).append(rec)
    for records in genes.values():
        records.sort(key=lambda rec: rec.start)
    return genes


def assign_genes(peaks: pd.DataFrame, genes: dict[tuple[str, str], list[GeneRecord]]) -> pd.DataFrame:
    peaks = peaks.copy()
    if "coord" not in peaks.columns:
        peaks["coord"] = np.where(peaks["strand"].astype(str) == "+", peaks["end"], peaks["start"])
    if "peakID" not in peaks.columns:
        peaks["peakID"] = peaks.index.astype(str)
    peaks["site_id"] = peaks["peakID"].astype(str)
    rows = []
    for (chrom, strand), sub in peaks.groupby(["chr", "strand"], dropna=False):
        records = genes.get((str(chrom), str(strand)), [])
        starts = [rec.start for rec in records]
        active: list[GeneRecord] = []
        cursor = 0
        sub = sub.sort_values("coord")
        for idx, row in sub.iterrows():
            coord = int(row["coord"])
            upto = bisect.bisect_right(starts, coord, lo=cursor)
            active.extend(records[cursor:upto])
            cursor = upto
            active = [rec for rec in active if rec.end >= coord]
            hits = [rec for rec in active if rec.start <= coord <= rec.end]
            best = min(hits, key=lambda rec: rec.width) if hits else None
            rows.append(
                {
                    "index": idx,
                    "gene_id": best.gene_id if best else np.nan,
                    "gene_name": best.gene_name if best else np.nan,
                }
            )
    anno = pd.DataFrame(rows).set_index("index")
    peaks[["gene_id", "gene_name"]] = anno.reindex(peaks.index)[["gene_id", "gene_name"]]
    return peaks


def build_site_usage(
    counts: pd.DataFrame,
    sites: pd.DataFrame,
    spot_ids: list[str],
    min_parent_count: int,
) -> pd.DataFrame:
    counts = counts.reindex(index=sites["site_id"], columns=spot_ids).fillna(0.0)
    usage = pd.DataFrame(index=counts.index, columns=counts.columns, dtype=float)
    grouped = sites.dropna(subset=["gene_id"]).groupby("gene_id")["site_id"].apply(list)
    valid_site_ids: list[str] = []
    for _, site_ids in grouped.items():
        site_ids = [site for site in site_ids if site in counts.index]
        if len(site_ids) < 2:
            continue
        parent_counts = counts.loc[site_ids].sum(axis=0)
        valid_spots = parent_counts >= min_parent_count
        values = counts.loc[site_ids].div(parent_counts.replace(0, np.nan), axis=1)
        values.loc[:, ~valid_spots] = np.nan
        usage.loc[site_ids] = values
        valid_site_ids.extend(site_ids)
    return usage.loc[valid_site_ids]


def postprocess(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    raw_dir = paths["raw_dir"]
    processed_dir = paths["processed_dir"]
    processed_dir.mkdir(parents=True, exist_ok=True)

    peaks = pd.read_csv(raw_dir / "peaks_meta.csv.gz", index_col=0)
    counts = pd.read_csv(raw_dir / "apa_site_counts.csv.gz", index_col=0)
    peaks.index = peaks.index.astype(str)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)

    coordinates = load_10x_positions(paths["spatial_dir"])
    coord_lookup = coordinates.copy()
    coord_lookup["spot_id_core"] = coord_lookup["spot_id_core"].astype(str)
    coord_lookup["spot_id"] = coord_lookup["spot_id"].astype(str)

    core_to_spot = {
        core: spot for core, spot in zip(coord_lookup["spot_id_core"], coord_lookup["spot_id"])
    }
    spot_ids = [core_to_spot[core] for core in counts.columns if core in core_to_spot]
    if not spot_ids:
        raise ValueError("No overlap between scAPAtrap barcodes and spatial coordinates")
    spot_cores = [spot.replace("-1", "") for spot in spot_ids]
    counts = counts.loc[:, spot_cores]
    counts.columns = spot_ids
    coordinates = coordinates.set_index("spot_id_core")
    coordinates = coordinates.loc[spot_cores].copy()
    coordinates["spot_id"] = spot_ids
    coordinates = coordinates.reset_index(drop=True)[["spot_id", "x", "y"]]

    genes = read_genes_from_gtf(paths["gtf"])
    sites = assign_genes(peaks, genes)
    sites = sites[sites["gene_id"].notna()].copy()
    if sites.empty:
        raise ValueError("No scAPAtrap sites overlap GTF gene intervals")
    sites["site_id"] = sites["site_id"].astype(str)
    sites = sites.drop_duplicates("site_id")
    counts = counts.reindex(index=sites["site_id"], columns=spot_ids).fillna(0)

    apa = build_site_usage(
        counts=counts,
        sites=sites,
        spot_ids=spot_ids,
        min_parent_count=args.min_parent_count,
    )
    if apa.empty:
        raise ValueError("No multi-site genes passed APA usage construction")

    site_cols = ["site_id", "peakID", "gene_id", "gene_name", "chr", "start", "end", "strand", "coord"]
    extra_cols = [col for col in sites.columns if col not in site_cols]
    sites[site_cols + extra_cols].to_csv(processed_dir / "apa_sites.csv.gz", index=False)
    counts.to_csv(processed_dir / "apa_site_counts.csv.gz")
    apa.to_csv(processed_dir / "apa_matrix.csv")
    coordinates.to_csv(processed_dir / "coordinates.csv", index=False)

    metadata = pd.DataFrame(
        {
            "spot_id": spot_ids,
            "spot_id_core": [spot.replace("-1", "") for spot in spot_ids],
            "dataset": args.dataset_name,
            "source": args.source_label,
        }
    )
    metadata.to_csv(processed_dir / "metadata.csv", index=False)

    qc = {
        "dataset": args.dataset_name,
        "platform": "10x Visium",
        "species": args.species,
        "tissue": args.tissue,
        "apa_source": "BAM_calling_scAPAtrap",
        "apa_ready": True,
        "expression_ready": False,
        "biological_labels_ready": False,
        "site_level_ready": True,
        "n_spots": len(spot_ids),
        "n_called_sites": int(len(peaks)),
        "n_gene_annotated_sites": int(len(sites)),
        "n_apa_usage_sites": int(apa.shape[0]),
        "n_apa_usage_spots": int(apa.shape[1]),
        "min_parent_count": args.min_parent_count,
        "files": {
            "apa_matrix": "apa_matrix.csv",
            "coordinates": "coordinates.csv",
            "metadata": "metadata.csv",
            "apa_sites": "apa_sites.csv.gz",
            "apa_site_counts": "apa_site_counts.csv.gz",
        },
        "notes": [
            "APA matrix rows are scAPAtrap PAS sites within genes that have at least two called sites.",
            "Values are site usage fractions: site_count / same-gene total site counts per spot.",
            "Spots below min_parent_count for a gene are set to NaN.",
        ],
    }
    (processed_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2))


def main() -> None:
    args = parse_args()
    paths = {
        "bam": resolve(args.bam),
        "spatial_dir": resolve(args.spatial_dir),
        "gtf": resolve(args.gtf),
        "output_root": resolve(args.output_root),
        "processed_dir": resolve(args.processed_dir),
    }
    paths["raw_dir"] = paths["output_root"] / "raw_scapatrap"
    paths["log_dir"] = paths["output_root"] / "logs"
    # scAPAtrap uses gsub(".bam", "", input) with regex semantics, so avoid
    # directory or stem names containing strings like "_bam".
    paths["stage_dir"] = paths["output_root"] / "stage"
    paths["staged_bam"] = paths["stage_dir"] / "spaceranger_input.bam"
    paths["internal_log"] = paths["log_dir"] / "scapatrap_internal.log"

    ensure_file(paths["bam"], "BAM")
    ensure_file(paths["gtf"], "GTF")
    if not paths["spatial_dir"].exists():
        raise FileNotFoundError(f"spatial-dir not found: {paths['spatial_dir']}")

    if args.postprocess_only:
        ensure_file(paths["raw_dir"] / "peaks_meta.csv.gz", "peaks_meta.csv.gz")
        ensure_file(paths["raw_dir"] / "apa_site_counts.csv.gz", "apa_site_counts.csv.gz")
        if args.force:
            shutil.rmtree(paths["processed_dir"], ignore_errors=True)
        postprocess(args, paths)
        print(f"Prepared dataset written to {paths['processed_dir']}")
        return

    if args.force:
        shutil.rmtree(paths["raw_dir"], ignore_errors=True)
        shutil.rmtree(paths["processed_dir"], ignore_errors=True)
        shutil.rmtree(paths["stage_dir"], ignore_errors=True)
        if paths["internal_log"].exists():
            paths["internal_log"].unlink()

    if paths["raw_dir"].exists():
        raise FileExistsError(f"raw scAPAtrap dir exists: {paths['raw_dir']}")
    paths["output_root"].mkdir(parents=True, exist_ok=True)
    paths["log_dir"].mkdir(parents=True, exist_ok=True)
    paths["stage_dir"].mkdir(parents=True, exist_ok=True)
    if not paths["staged_bam"].exists():
        paths["staged_bam"].symlink_to(paths["bam"])
    source_bai = Path(f"{paths['bam']}.bai")
    staged_bai = Path(f"{paths['staged_bam']}.bai")
    if source_bai.exists() and not staged_bai.exists():
        staged_bai.symlink_to(source_bai)

    tools = tool_paths(resolve(args.tool_prefix))
    rscript_path = write_r_script(args, paths, tools)
    manifest = {
        "dataset": args.dataset_name,
        "bam": str(paths["bam"]),
        "staged_bam": str(paths["staged_bam"]),
        "spatial_dir": str(paths["spatial_dir"]),
        "gtf": str(paths["gtf"]),
        "output_root": str(paths["output_root"]),
        "raw_dir": str(paths["raw_dir"]),
        "processed_dir": str(paths["processed_dir"]),
        "rscript": str(Path(args.rscript).expanduser()),
        "r_script": str(rscript_path),
        "tools": tools,
        "parameters": {
            "threads": args.threads,
            "readlength": args.readlength,
            "cov_cutoff": args.cov_cutoff,
            "min_cells": args.min_cells,
            "min_count": args.min_count,
            "min_parent_count": args.min_parent_count,
            "tails_search": args.tails_search,
        },
    }
    (paths["output_root"] / "run_manifest.json").write_text(json.dumps(manifest, indent=2))

    if args.dry_run:
        print(json.dumps(manifest, indent=2))
        return

    rscript = str(Path(args.rscript).expanduser())
    result = subprocess.run([rscript, str(rscript_path)], cwd=paths["output_root"], check=False)
    if result.returncode != 0:
        raise SystemExit(result.returncode)
    postprocess(args, paths)
    print(f"Prepared dataset written to {paths['processed_dir']}")


if __name__ == "__main__":
    main()

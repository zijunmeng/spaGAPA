#!/usr/bin/env python
"""Prepare the stAPAminer MOB example as a spaGAPA real-data benchmark.

This script uses the published stAPAminer package example files:

- APA.RDA: movAPA PACdataset with poly(A) site counts
- count.csv: gene expression matrix used by stAPAminer for KNN imputation
- position.txt: MOB spot coordinates and layer labels

It exports a standardized processed dataset:

- apa_matrix.csv: raw RUD APA index, genes x spots
- coordinates.csv: spot coordinates, spot_id/x/y
- metadata.csv: spot metadata with MOB layer labels
- expression_matrix.csv: gene expression, genes x spots
- stapaminer_rud_raw.csv: raw RUD from movAPA
- stapaminer_rud_imputed.csv: stAPAminer-compatible KNN-imputed RUD
- apa_site_counts.csv.gz: poly(A) site count matrix, sites x spots
- apa_sites.csv.gz: poly(A) site annotations
- qc_summary.json: dataset provenance and dimensions

The R code is generated at runtime so we can keep the source package untouched.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


DEFAULT_RSCRIPT = "~/anaconda3/envs/r442/bin/Rscript"
DEFAULT_R_LIB = "/s1/SHARE/01_software/R_442_SeuratV5/library"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stapaminer-dir",
        default="00_ref_packages/stAPAminer-main",
        help="Path to the unpacked stAPAminer source directory.",
    )
    parser.add_argument(
        "--output-dir",
        default="spaGAPA/data/processed/stapaminer_mob",
        help="Output directory for standardized spaGAPA files.",
    )
    parser.add_argument(
        "--rscript",
        default=DEFAULT_RSCRIPT,
        help="Rscript executable.",
    )
    parser.add_argument(
        "--r-lib",
        default=DEFAULT_R_LIB,
        help="R library path containing movAPA and stAPAminer dependencies.",
    )
    return parser.parse_args()


def run_r_export(
    stapaminer_dir: Path,
    output_dir: Path,
    rscript: str,
    r_lib: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    r_code = f"""
    .libPaths(c("{r_lib}", .libPaths()))
    suppressPackageStartupMessages(library(movAPA))

    stapaminer_dir <- normalizePath("{stapaminer_dir}")
    output_dir <- normalizePath("{output_dir}", mustWork=FALSE)
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

    load(file.path(stapaminer_dir, "inst", "extdata", "APA.RDA"))
    position <- read.table(
      file.path(stapaminer_dir, "inst", "extdata", "position.txt"),
      header=TRUE,
      check.names=FALSE
    )
    count <- read.csv(
      file.path(stapaminer_dir, "inst", "extdata", "count.csv"),
      row.names=1,
      check.names=FALSE
    )

    shared_spots <- Reduce(intersect, list(
      rownames(position),
      colnames(count),
      colnames(APA@counts)
    ))
    position <- position[shared_spots, , drop=FALSE]
    count <- count[, shared_spots, drop=FALSE]

    # stAPAminer's computeAPAIndex mutates APA@counts after subsetting.
    # With R 4.4 / movAPA 0.2.0 the slot requires AnyMatrix, so we use
    # as.matrix() while preserving the original method semantics.
    APA_subset <- APA
    APA_subset@counts <- as.matrix(APA_subset@counts[, shared_spots, drop=FALSE])
    rud_raw <- movAPAindex(
      APA_subset,
      method="RUD",
      choose2PA=NULL,
      RUD.includeNon3UTR=FALSE,
      clearPAT=0
    )
    rud_raw <- rud_raw[-nrow(rud_raw), , drop=FALSE]
    rud_raw <- as.matrix(rud_raw)

    imputeAPAIndex_compatible <- function(index, gene, k=10, init=TRUE) {{
      colNames <- colnames(index)[colnames(index) %in% colnames(gene)]
      index <- index[, colNames, drop=FALSE]
      gene <- gene[, colNames, drop=FALSE]
      if (init == TRUE) {{
        gene_some <- gene[rownames(index), , drop=FALSE]
        is_zero <- gene_some == 0
        index[is_zero & is.na(index)] <- 0
      }}
      scaleData <- as.data.frame(t(scale(gene)), stringsAsFactors=FALSE)
      scaleData[is.na(scaleData)] <- 0
      dist_mat <- as.matrix(dist(scaleData, method="euclidean"))
      k_eff <- min(k, nrow(dist_mat) - 1)
      weight <- data.frame()
      for (i in 1:nrow(dist_mat)) {{
        weight <- rbind(weight, order(dist_mat[i, ])[2:(k_eff + 1)])
      }}
      rownames(weight) <- rownames(dist_mat)
      colnames(weight) <- paste0("N", seq_len(k_eff))

      num <- 1
      while (sum(is.na(index)) > 0 && num <= 10) {{
        for (i in 1:ncol(index)) {{
          target <- index[, i]
          if (sum(is.na(target)) == 0) next
          fill_source <- index[, as.numeric(weight[i, ]), drop=FALSE]
          fill_mean <- apply(fill_source, 1, mean, na.rm=TRUE)
          target[is.na(target)] <- fill_mean[is.na(target)]
          index[, i] <- target
        }}
        num <- num + 1
      }}
      if (num > 10) {{
        index[is.na(index)] <- 0
      }}
      return(index)
    }}

    rud_imputed <- imputeAPAIndex_compatible(rud_raw, count, k=10, init=TRUE)
    rud_imputed <- as.matrix(rud_imputed)

    coords_out <- data.frame(
      spot_id=rownames(position),
      x=position$x,
      y=position$y,
      row.names=NULL,
      check.names=FALSE
    )
    meta_out <- data.frame(
      spot_id=rownames(position),
      layer=position$label,
      dataset="stAPAminer_MOB",
      source="stAPAminer inst/extdata",
      row.names=NULL,
      check.names=FALSE
    )

    write.csv(rud_raw, file.path(output_dir, "apa_matrix.csv"), quote=FALSE)
    write.csv(coords_out, file.path(output_dir, "coordinates.csv"), row.names=FALSE, quote=FALSE)
    write.csv(meta_out, file.path(output_dir, "metadata.csv"), row.names=FALSE, quote=FALSE)
    write.csv(count, file.path(output_dir, "expression_matrix.csv"), quote=FALSE)
    write.csv(rud_raw, file.path(output_dir, "stapaminer_rud_raw.csv"), quote=FALSE)
    write.csv(rud_imputed, file.path(output_dir, "stapaminer_rud_imputed.csv"), quote=FALSE)
    write.csv(as.matrix(APA_subset@counts), gzfile(file.path(output_dir, "apa_site_counts.csv.gz")), quote=FALSE)
    write.csv(as.data.frame(APA_subset@anno), gzfile(file.path(output_dir, "apa_sites.csv.gz")), quote=FALSE)

    qc <- list(
      dataset="stAPAminer_MOB",
      source_dir=stapaminer_dir,
      n_spots=nrow(position),
      n_expression_genes=nrow(count),
      n_polyA_sites=nrow(APA_subset@counts),
      n_rud_genes=nrow(rud_raw),
      raw_rud_na=sum(is.na(rud_raw)),
      raw_rud_missing_rate=sum(is.na(rud_raw)) / length(rud_raw),
      imputed_rud_na=sum(is.na(rud_imputed)),
      layers=as.list(table(position$label)),
      files=list(
        apa_matrix="apa_matrix.csv",
        coordinates="coordinates.csv",
        metadata="metadata.csv",
        expression_matrix="expression_matrix.csv",
        stapaminer_rud_raw="stapaminer_rud_raw.csv",
        stapaminer_rud_imputed="stapaminer_rud_imputed.csv",
        apa_site_counts="apa_site_counts.csv.gz",
        apa_sites="apa_sites.csv.gz"
      )
    )
    writeLines(jsonlite::toJSON(qc, auto_unbox=TRUE, pretty=TRUE),
               file.path(output_dir, "qc_summary.json"))

    cat("Exported stAPAminer MOB dataset to", output_dir, "\\n")
    cat("RUD genes:", nrow(rud_raw), "spots:", ncol(rud_raw), "\\n")
    cat("Raw RUD NA:", sum(is.na(rud_raw)), "\\n")
    cat("Imputed RUD NA:", sum(is.na(rud_imputed)), "\\n")
    """

    cmd = [str(Path(rscript).expanduser()), "-e", r_code]
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    run_r_export(
        stapaminer_dir=Path(args.stapaminer_dir).resolve(),
        output_dir=Path(args.output_dir).resolve(),
        rscript=args.rscript,
        r_lib=args.r_lib,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env Rscript
# Thin wrapper: run stAPAminer's imputeAPAIndex (expression-KNN) on a masked
# gene x spot APA index matrix and write the imputed matrix.
#
# Usage:
#   Rscript run_stapaminer_impute.R <masked_index.csv> <expression.csv> <out.csv> [k]
#
# masked_index.csv : gene x spot matrix, rownames=gene_symbol, colnames=spot,
#                    NA where masked/missing.
# expression.csv   : gene x spot COUNT matrix, rownames=gene_symbol (superset ok).
#                    Spots subset to those shared with index automatically.
# out.csv          : imputed gene x spot matrix (NA all filled).
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages({
  library(stAPAminer)
})

args <- commandArgs(trailingOnly = TRUE)
index_path <- args[[1]]
expr_path <- args[[2]]
out_path <- args[[3]]
k <- if (length(args) >= 4) as.integer(args[[4]]) else 10L

index <- as.matrix(read.csv(index_path, row.names = 1, check.names = FALSE))
expr <- as.matrix(read.csv(expr_path, row.names = 1, check.names = FALSE))
cat(sprintf("[stAPAminer] index %d x %d | expr %d x %d | k=%d\n",
            nrow(index), ncol(index), nrow(expr), ncol(expr), k))

# imputeAPAIndex subsets spots to the intersection and subsets expr rows to
# index rownames internally; pass the full expression matrix.
imp <- imputeAPAIndex(index, expr, k = k, init = TRUE)

write.csv(imp, out_path, quote = FALSE)
cat(sprintf("[stAPAminer] wrote imputed matrix %d x %d -> %s\n",
            nrow(imp), ncol(imp), out_path))

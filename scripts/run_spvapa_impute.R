#!/usr/bin/env Rscript
# Thin wrapper: run spvAPA's WNNImpute (multimodal weighted-nearest-neighbour,
# built on a Seurat WNN graph over RNA + APA assays) on a masked gene x spot
# APA index matrix and write the imputed matrix.
#
# Usage:
#   Rscript run_spvapa_impute.R <masked_index.csv> <expression.csv> <out.csv> [k]
#
# masked_index.csv : gene x spot matrix, rownames=gene_symbol, colnames=spot,
#                    NA where masked/missing.
# expression.csv   : gene x spot COUNT matrix, rownames=gene_symbol (superset ok).
#                    Spots subset to those shared with index automatically.
#                    NOTE: WNNImpute does gene[rownames(APA), ] for its `init`
#                    step, so every APA-row gene MUST be present in expression.
# out.csv          : imputed gene x spot matrix (NA all filled, k iterations).
#
# WNNImpute input contract (verified against installed spvAPA 0.1.0 source):
#   gene : gene x spot COUNT matrix (-> Seurat RNA assay, SCTransform + PCA)
#   APA  : gene x spot APA index matrix, NA where missing (-> Seurat PDUI assay)
#   k    : number of WNN neighbours to average (default 20)
#   init : TRUE -> set APA=0 where the matching gene has zero expression
#   is.weight=FALSE -> returns the imputed matrix (not a list).
# WNN derives neighbours from the multimodal (RNA+APA) PCA graph, NOT from
# spatial coordinates, so coords are not required here.
.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressPackageStartupMessages({
  library(spvAPA)
})

args <- commandArgs(trailingOnly = TRUE)
index_path <- args[[1]]
expr_path <- args[[2]]
out_path <- args[[3]]
k <- if (length(args) >= 4) as.integer(args[[4]]) else 15L

index <- as.matrix(read.csv(index_path, row.names = 1, check.names = FALSE))
expr <- as.matrix(read.csv(expr_path, row.names = 1, check.names = FALSE))
cat(sprintf("[spvAPA] index %d x %d | expr %d x %d | k=%d\n",
            nrow(index), ncol(index), nrow(expr), ncol(expr), k))

# shared spots (cols) in the same order in both matrices
shared <- intersect(colnames(index), colnames(expr))
if (length(shared) == 0) stop("no shared spots between index and expression")
index <- index[, shared, drop = FALSE]
expr <- expr[, shared, drop = FALSE]

# WNNImpute's init=TRUE step does gene[rownames(APA), ], so every index gene
# must exist in expression.  Subset expression rows to index genes that exist;
# drop index genes with no matching expression row (cannot be imputed by WNN).
missing_genes <- setdiff(rownames(index), rownames(expr))
if (length(missing_genes)) {
  cat(sprintf("[spvAPA] dropping %d index gene(s) absent from expression (init step requires them)\n",
              length(missing_genes)))
  keep <- rownames(index)[rownames(index) %in% rownames(expr)]
  index <- index[keep, , drop = FALSE]
}
expr_sub <- expr[rownames(index), , drop = FALSE]

# WNNImpute needs non-negative integer-like counts in the RNA assay for
# SCTransform; round to be safe (expression is already counts in our pipeline).
expr_sub[expr_sub < 0] <- 0

t0 <- proc.time()
imp <- WNNImpute(gene = expr_sub, APA = index, k = k, init = TRUE, is.weight = FALSE)
elapsed <- proc.time() - t0

# is.weight=FALSE -> a matrix; guard against accidental list returns.
if (is.list(imp)) imp <- imp[[1]]
imp <- as.matrix(imp)
write.csv(imp, out_path, quote = FALSE)
cat(sprintf("[spvAPA] wrote imputed matrix %d x %d (NA=%d) in %.1fs -> %s\n",
            nrow(imp), ncol(imp), sum(is.na(imp)), elapsed[["elapsed"]], out_path))

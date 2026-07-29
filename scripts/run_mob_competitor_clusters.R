#!/usr/bin/env Rscript
# Competitor domain-recovery on spvAPA's ST11 MOB data (Phase 3 Task Y).
#
# Runs the two published spatial-APA domain methods on the SAME 260-spot MOB
# data used by spaGAPA, and scores recovery vs the expert MOB layer labels
# (GCL/GL/MCL/ONL/OPL) with ARI/NMI:
#
#   * spvAPA::makeCluster      - Seurat Louvain on expression (ST11GEM) PCA
#   * stAPAminer::makeStCluster - Seurat clustering on expression (ST11GEM) PCA
#
# Both competitors cluster on gene EXPRESSION (their native mode); spaGAPA
# clusters on the fused spatial+APA graph. To make k comparable we scan a few
# Louvain resolutions and also report a fixed k=5 spectral clustering on the
# expression PCA embedding.
#
# Outputs append to pipeline_output/mob_domain_recovery/competitor_metrics.json
# (sibling to spaGAPA's spagapa_metrics.json).
#
# Run on S91:
#   source ~/anaconda3/etc/profile.d/conda.sh && conda activate r442
#   export R_LIBS=/s1/SHARE/01_software/R_442_SeuratV5/library
#   Rscript scripts/run_mob_competitor_clusters.R

.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressMessages({
    library(spvAPA)
    library(stAPAminer)
    library(Seurat)
    library(mclust)  # adjustedRandIndex
})

OUT <- "pipeline_output/mob_domain_recovery"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

data(ST11RUD, package = "spvAPA")
data(ST11GEM, package = "spvAPA")
data(ST11label, package = "spvAPA")

spots <- colnames(ST11RUD)
stopifnot(all(spots == rownames(ST11label)))
true_lab <- ST11label$label
gem <- as.matrix(ST11GEM[, spots])  # gene x spot raw counts

# adjustedRandIndex gives ARI; NMI via a small helper
nmi <- function(a, b) {
    tab <- table(a, b)
    p <- tab / sum(tab)
    pa <- rowSums(p); pb <- colSums(p)
    ha <- -sum(pa * log(pa + 1e-12))
    hb <- -sum(pb * log(pb + 1e-12))
    hab <- -sum(p * log(p + 1e-12))
    num <- 2 * (ha + hab)
    den <- ha + hb
    if (den <= 0) return(0)
    return(num / den)
}

results <- list(
    dataset = "spvAPA_ST11_MOB",
    n_spots = length(spots),
    method = "competitor_clusters",
    note = "Both competitors cluster on gene expression (native mode); spaGAPA uses spatial+APA graph.",
    spvAPA = list(),
    stAPAminer = list()
)

# ── Build a Seurat object from ST11GEM (spots x genes) ────────────────
# spvAPA/stAPAminer both feed a Seurat object; rownames = genes.
emat <- t(gem)  # spots x genes
seu <- CreateSeuratObject(counts = emat, meta.data = ST11label[spots, ])

# ── spvAPA::makeCluster across resolutions ────────────────────────────
cat("spvAPA makeCluster (expression-based Seurat clustering):\n")
spv_runs <- list()
for (res in c(0.1, 0.2, 0.3, 0.5, 0.8)) {
    set.seed(42)
    s2 <- spvAPA::makeCluster(seu, norm = TRUE, res = res, dims = 1:10, k.param = 20)
    cl <- as.integer(Idents(s2))
    ari <- mclust::adjustedRandIndex(true_lab, cl)
    n <- nmi(true_lab, cl)
    nd <- length(unique(cl))
    cat(sprintf("  res=%.1f -> %d clusters, ARI=%.3f, NMI=%.3f\n", res, nd, ari, n))
    spv_runs[[paste0("res_", res)]] <- list(n_domains = nd, ari = ari, nmi = n)
}
results$spvAPA$makeCluster <- spv_runs

# ── stAPAminer::makeStCluster + findLabels ───────────────────────────
cat("stAPAminer makeStCluster (expression-based):\n")
# createStAPAminerObject expects count (genes x spots) + metaData (spots x ...)
obj <- tryCatch(
    stAPAminer::createStAPAminerObject(count = gem, metaData = ST11label[spots, ]),
    error = function(e) { cat("  createStAPAminerObject failed:", conditionMessage(e), "\n"); NULL }
)
sta_runs <- list()
if (!is.null(obj)) {
    for (res in c(0.2, 0.3, 0.5, 0.8)) {
        set.seed(42)
        o2 <- tryCatch(stAPAminer::makeStCluster(obj, nfeatures = 2000,
                                                 dims = 1:10, resolution = res, k = 20),
                       error = function(e) { cat("  makeStCluster res=", res, "failed:", conditionMessage(e), "\n"); NULL })
        if (is.null(o2)) next
        cl <- as.integer(Idents(o2@seurat))
        ari <- mclust::adjustedRandIndex(true_lab, cl)
        n <- nmi(true_lab, cl)
        nd <- length(unique(cl))
        cat(sprintf("  res=%.1f -> %d clusters, ARI=%.3f, NMI=%.3f\n", res, nd, ari, n))
        sta_runs[[paste0("res_", res)]] <- list(n_domains = nd, ari = ari, nmi = n)
    }
}
results$stAPAminer$makeStCluster <- sta_runs

# ── write JSON ────────────────────────────────────────────────────────
jsonlite <- asNamespace("jsonlite")
out_path <- file.path(OUT, "competitor_metrics.json")
con <- file(out_path, "w")
writeLines(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE), con)
close(con)
cat("\nWrote", out_path, "\n")

#!/usr/bin/env Rscript
# Export spvAPA's MOB ST11RUD (gene x spot APA / RUD matrix) and ST11label
# (x, y, MOB-layer label) to CSV for the spaGAPA domain-recovery validation.
#
# Also exports ST11GEM (gene-expression matrix) so spaGAPA can build an
# expression-view graph, matching spvAPA/stAPAminer which cluster on
# expression.
#
# Output -> data/processed/mob_st11/
#   - apa_matrix.csv      : rows = gene, cols = spot, values = RUD in [0,1]
#                           (header row = spot ids, first col = gene id)
#   - coordinates.csv     : spot_id, x, y  (same spot order as apa_matrix cols)
#   - labels.csv          : spot_id, label  (MOB layer: GCL/GL/MCL/ONL/OPL)
#   - expression.csv      : rows = gene, cols = spot (raw counts)
#
# Run on S91:  source ~/anaconda3/etc/profile.d/conda.sh && conda activate r442
#              export R_LIBS=/s1/SHARE/01_software/R_442_SeuratV5/library
#              Rscript scripts/export_spvapa_st11.R

.libPaths(c("/s1/SHARE/01_software/R_442_SeuratV5/library", .libPaths()))
suppressMessages({
    library(spvAPA)
    library(stAPAminer)
})

OUT <- "data/processed/mob_st11"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Load the spvAPA MOB example data
data(ST11RUD, package = "spvAPA")
data(ST11label, package = "spvAPA")
data(ST11GEM, package = "spvAPA")

cat("ST11RUD:", dim(ST11RUD)[1], "genes x", dim(ST11RUD)[2], "spots\n")
cat("ST11label:", dim(ST11label)[1], "spots; layers:",
    paste(names(table(ST11label$label)), collapse = "/"), "\n")
cat("ST11GEM:", dim(ST11GEM)[1], "genes x", dim(ST11GEM)[2], "spots\n")

# Align labels to the APA-matrix column order (ST11RUD cols = spots)
spot_ids <- colnames(ST11RUD)
stopifnot(all(spot_ids %in% rownames(ST11label)))
labels_aligned <- ST11label[spot_ids, ]

stopifnot(all(spot_ids %in% colnames(ST11GEM)))
gem_aligned <- ST11GEM[, spot_ids]

# ── APA matrix: gene_id,spot1,spot2,... ───────────────────────────────
apa <- cbind(gene_id = rownames(ST11RUD), as.data.frame(ST11RUD))
write.csv(apa, file.path(OUT, "apa_matrix.csv"), row.names = FALSE)

# ── Coordinates (x, y) and labels ─────────────────────────────────────
coords <- data.frame(
    spot_id = spot_ids,
    x = labels_aligned$x,
    y = labels_aligned$y,
    stringsAsFactors = FALSE
)
write.csv(coords, file.path(OUT, "coordinates.csv"), row.names = FALSE)

labels <- data.frame(
    spot_id = spot_ids,
    label = labels_aligned$label,
    stringsAsFactors = FALSE
)
write.csv(labels, file.path(OUT, "labels.csv"), row.names = FALSE)

# ── Expression matrix (gene x spot, raw counts) ───────────────────────
gem <- cbind(gene_id = rownames(gem_aligned), as.data.frame(gem_aligned))
write.csv(gem, file.path(OUT, "expression.csv"), row.names = FALSE)

cat("\nExported to", OUT, ":\n")
cat("  apa_matrix.csv   (", nrow(apa) - 1, "genes x", ncol(apa) - 1, "spots)\n", sep = "")
cat("  coordinates.csv  (", nrow(coords), "spots)\n", sep = "")
cat("  labels.csv       layers:", paste(names(table(labels$label)), collapse = "/"), "\n")
cat("  expression.csv   (", nrow(gem) - 1, "genes)\n", sep = "")
cat("\nSpot-id alignment verified: APA cols == labels rows == expression cols.\n")

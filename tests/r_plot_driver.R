#!/usr/bin/env Rscript
# R reference plotSmoothers — uses pre-saved synthetic fixture.
suppressPackageStartupMessages({
  library(tradeSeq)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Reuse the manifest fixture
counts_path <- file.path(out_dir, "../../data/_sidecar_counts.csv")
pseudotime_path <- file.path(out_dir, "../../data/_sidecar_pseudotime.csv")
weights_path <- file.path(out_dir, "../../data/_sidecar_cellWeights.csv")

counts <- as.matrix(read.csv(counts_path, row.names=1, check.names=FALSE))
pseudotime <- as.matrix(read.csv(pseudotime_path, check.names=FALSE))
weights <- as.matrix(read.csv(weights_path, check.names=FALSE))
# tradeSeq requires row + col names
rownames(pseudotime) <- colnames(counts)
rownames(weights) <- colnames(counts)

cat("counts:", dim(counts), "  pseudotime:", dim(pseudotime), "  weights:", dim(weights), "\n")

set.seed(42)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = weights, nknots = 5)

# Pick a top gene
ar <- associationTest(sce)
top_gene <- rownames(ar)[order(ar$pvalue)[1]]
cat("top gene:", top_gene, "\n")

p <- plotSmoothers(sce, counts = counts, gene = top_gene)
ggsave(file.path(out_dir, "R_smoother.png"), p, width = 6, height = 4, dpi = 100)
write.csv(data.frame(gene = top_gene), file.path(out_dir, "top_gene.csv"), row.names = FALSE)
cat("R plot written\n")

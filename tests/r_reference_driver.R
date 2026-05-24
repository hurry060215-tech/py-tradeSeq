#!/usr/bin/env Rscript
# Reference runner — fits NB GAMs on tradeSeq's bundled `sds` fixture,
# runs associationTest + startVsEndTest, dumps everything to JSON.
#
# Usage: Rscript r_reference_driver.R <fixture.rds> <output.json>
# (fixture path is informational — we use tradeSeq's bundled data via data())

suppressMessages({
  library(tradeSeq)
  library(jsonlite)
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]
output_path  <- args[2]

# Load our pre-generated fixture: synthetic counts paired with the bundled `sds`
fixture <- readRDS(fixture_path)
counts <- fixture$counts
sds <- fixture$sds
cat("[ref] counts dims:", dim(counts), "\n")
cat("[ref] sds class:", class(sds), "\n")

# Sidecar: extract pseudotime + cellWeights from the S4 SlingshotDataSet so
# the Python candidate (which can't read S4) can rebuild the same input.
suppressMessages(library(slingshot))
sidecar_dir <- dirname(output_path)
write.csv(slingshot::slingPseudotime(sds, na = FALSE),
          file.path(sidecar_dir, "_sidecar_pseudotime.csv"), row.names = FALSE)
write.csv(slingshot::slingCurveWeights(sds),
          file.path(sidecar_dir, "_sidecar_cellWeights.csv"), row.names = FALSE)
write.csv(counts, file.path(sidecar_dir, "_sidecar_counts.csv"))
cat("[ref] sidecars written: pseudotime/cellWeights/counts\n")

set.seed(42)

# Fit GAMs — small (nknots=4) for speed
t0 <- proc.time()[["elapsed"]]
sce <- fitGAM(counts = counts, sds = sds, nknots = 4, verbose = FALSE)
t_fit <- proc.time()[["elapsed"]] - t0
cat("[ref] fitGAM:", t_fit, "s; output class:", class(sce), "\n")

# associationTest — does expression depend on pseudotime?
t0 <- proc.time()[["elapsed"]]
at <- associationTest(sce)
t_at <- proc.time()[["elapsed"]] - t0
cat("[ref] associationTest:", t_at, "s; output dims:", dim(at), "\n")
cat("[ref] head associationTest:\n"); print(head(at))

# startVsEndTest — does expression differ between trajectory endpoints?
t0 <- proc.time()[["elapsed"]]
svet <- startVsEndTest(sce)
t_svet <- proc.time()[["elapsed"]] - t0
cat("[ref] startVsEndTest:", t_svet, "s\n")

# diffEndTest — does expression differ between endpoints of different lineages?
t0 <- proc.time()[["elapsed"]]
det <- tryCatch(diffEndTest(sce),
                error = function(e) { cat("[ref] diffEndTest errored:", e$message, "\n"); NULL })
t_det <- proc.time()[["elapsed"]] - t0
cat("[ref] diffEndTest:", t_det, "s\n")

# Build flat JSON
gene_names <- rownames(counts)
at_pval  <- at$pvalue;  names(at_pval) <- rownames(at)
svet_pval <- svet$pvalue; names(svet_pval) <- rownames(svet)

# Top-50 genes by associationTest p-value
top50_at <- names(sort(at_pval, na.last = TRUE))[1:50]

# Same for startVsEndTest
top50_svet <- names(sort(svet_pval, na.last = TRUE))[1:50]

out <- list(
  gene_names = gene_names,
  n_genes = length(gene_names),
  n_cells = ncol(counts),
  n_lineages = ncol(slingshot::slingPseudotime(sds)),
  associationTest = list(
    pvalue = as.numeric(at_pval[gene_names]),
    waldStat = as.numeric(at$waldStat[match(gene_names, rownames(at))]),
    df = as.numeric(at$df[match(gene_names, rownames(at))]),
    top50 = top50_at
  ),
  startVsEndTest = list(
    pvalue = as.numeric(svet_pval[gene_names]),
    waldStat = as.numeric(svet$waldStat[match(gene_names, rownames(svet))]),
    df = as.numeric(svet$df[match(gene_names, rownames(svet))]),
    top50 = top50_svet
  ),
  timings = list(
    fitGAM = t_fit,
    associationTest = t_at,
    startVsEndTest = t_svet,
    diffEndTest = t_det
  )
)

write_json(out, output_path, auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE)
cat("[ref] wrote", output_path, "\n")

# Also save the SCE intermediate for any future debugging
saveRDS(sce, file.path(dirname(output_path), "sce_reference.rds"))
cat("[ref] saved sce_reference.rds for debugging\n")

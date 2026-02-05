#!/usr/bin/env Rscript
# ==============================================================================
# 04_bulk_rnaseq_gsea.R
# ==============================================================================
# Classic preranked GSEA analysis using fgsea.
#
# This script performs group-level pathway enrichment analysis using the
# Wald statistic from DESeq2 results as the ranking metric.
#
# Inputs:
#   - results/rnaseq/de/deseq2_results_*.tsv (DESeq2 results with stat column)
#   - references/gene_sets/*.gmt (pathway gene sets)
#
# Outputs (in results/rnaseq/gsea/):
#   - gsea_results.tsv (full GSEA results table)
#   - enrichment_*.pdf (enrichment plots for key pathways)
#
# Usage: Rscript scripts/04_bulk_rnaseq_gsea.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

cat("=== Preranked GSEA Analysis ===\n\n")

# ------------------------------------------------------------------------------
# Set up paths
# ------------------------------------------------------------------------------

repo_root <- getwd()

# Load config
config_path <- file.path(repo_root, "config", "config.yaml")
config <- read_yaml(config_path)

# Input files
de_dir <- file.path(repo_root, config$paths$results, "rnaseq", "de")
gmt_dir <- file.path(repo_root, "references", "gene_sets")

# Output directory
out_dir <- file.path(repo_root, config$paths$results, "rnaseq", "gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load DESeq2 results
# ------------------------------------------------------------------------------

cat("Loading DESeq2 results...\n")

# Use the first contrast
contrast_name <- config$contrasts$rnaseq[[1]]$name
de_path <- file.path(de_dir, paste0("deseq2_results_", contrast_name, ".tsv"))

if (!file.exists(de_path)) {
  stop("DESeq2 results not found: ", de_path, "\n  Run DESeq2 analysis first.")
}

de_results <- fread(de_path, sep = "\t")
cat("  Loaded", nrow(de_results), "genes from", contrast_name, "\n")

# ------------------------------------------------------------------------------
# Create ranked gene list
# ------------------------------------------------------------------------------

cat("Creating ranked gene list...\n")

# DESeq2 results now have gene_symbol as primary identifier
# Filter genes with valid symbols and stats
de_valid <- de_results[!is.na(gene_symbol) & gene_symbol != "" & !is.na(stat)]
cat("  Genes with valid symbols and stats:", nrow(de_valid), "\n")

# Create named vector: gene symbol -> Wald statistic
ranks <- setNames(de_valid$stat, de_valid$gene_symbol)
ranks <- sort(ranks, decreasing = TRUE)

cat("  Rank range:", round(min(ranks), 2), "to", round(max(ranks), 2), "\n\n")

# ------------------------------------------------------------------------------
# Load gene sets
# ------------------------------------------------------------------------------

cat("Loading gene sets...\n")

all_gene_sets <- list()

# Load all GMT files in the gene_sets directory
gmt_files <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)

for (gmt_file in gmt_files) {
  gmt_name <- basename(gmt_file)
  cat("  Loading:", gmt_name, "\n")

  # Read GMT file manually for fgsea format
  gmt_lines <- readLines(gmt_file)
  for (line in gmt_lines) {
    fields <- strsplit(line, "\t")[[1]]
    pathway_name <- fields[1]
    genes <- fields[3:length(fields)]  # Skip name and description
    genes <- genes[genes != ""]
    all_gene_sets[[pathway_name]] <- genes
  }
}

cat("  Total pathways:", length(all_gene_sets), "\n\n")

# ------------------------------------------------------------------------------
# Run fgsea
# ------------------------------------------------------------------------------

cat("Running fgsea...\n")

set.seed(42)

gsea_results <- fgsea(
  pathways = all_gene_sets,
  stats = ranks,
  minSize = 5,
  maxSize = 500,
  nPermSimple = 10000
)

# Sort by p-value
gsea_results <- gsea_results[order(pval)]

cat("  Pathways tested:", nrow(gsea_results), "\n")
cat("  Significant (padj < 0.05):", sum(gsea_results$padj < 0.05, na.rm = TRUE), "\n\n")

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

cat("Saving GSEA results...\n")

# Convert leadingEdge list to comma-separated string for TSV export
gsea_export <- copy(gsea_results)
gsea_export$leadingEdge <- sapply(gsea_export$leadingEdge, paste, collapse = ",")

results_path <- file.path(out_dir, "gsea_results.tsv")
fwrite(gsea_export, results_path, sep = "\t")
cat("  Saved:", results_path, "\n\n")

# ------------------------------------------------------------------------------
# Generate enrichment plots for key pathways
# ------------------------------------------------------------------------------

cat("Generating enrichment plots...\n")

# Key pathways to plot
key_pathways <- c(
  "INTERFERON_ALPHA_RESPONSE",
  "INTERFERON_GAMMA_RESPONSE",
  "INFLAMMATORY_RESPONSE",
  "ANTIGEN_PRESENTATION",
  "INTRINSIC_APOPTOSIS",
  "EXTRINSIC_APOPTOSIS",
  "VIRAL_MIMICRY_RESPONSE"
)

for (pathway in key_pathways) {
  if (pathway %in% names(all_gene_sets)) {
    # Get pathway stats
    pw_result <- gsea_results[pathway == pathway]

    if (nrow(pw_result) > 0) {
      # Create enrichment plot
      plot_path <- file.path(out_dir, paste0("enrichment_", tolower(pathway), ".pdf"))

      pdf(plot_path, width = 8, height = 5)
      p <- plotEnrichment(all_gene_sets[[pathway]], ranks) +
        labs(title = pathway,
             subtitle = sprintf("NES = %.2f, padj = %.2e",
                               pw_result$NES, pw_result$padj)) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))
      print(p)
      dev.off()

      cat("  Saved:", basename(plot_path), "\n")
    }
  } else {
    cat("  Pathway not found:", pathway, "\n")
  }
}

# ------------------------------------------------------------------------------
# Summary table of key pathways
# ------------------------------------------------------------------------------

cat("\nKey pathway results:\n")
key_results <- gsea_results[pathway %in% key_pathways]
if (nrow(key_results) > 0) {
  key_results <- key_results[, .(pathway, NES, pval, padj)]
  print(key_results)
}

cat("\n=== GSEA Analysis Complete ===\n")
cat("Output directory:", out_dir, "\n")

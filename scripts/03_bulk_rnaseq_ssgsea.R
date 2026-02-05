#!/usr/bin/env Rscript
# ==============================================================================
# 03_bulk_rnaseq_ssgsea.R
# ==============================================================================
# Single-sample GSEA (ssGSEA) analysis using GSVA.
#
# This script computes per-sample pathway activity scores for immune and
# viral mimicry gene sets using the ssGSEA method.
#
# Inputs:
#   - results/rnaseq/de/vst_matrix.tsv (VST-normalized expression from DESeq2)
#   - references/gene_sets/immune.gmt (immune pathway gene sets)
#   - references/gene_sets/viral_mimicry.gmt (viral mimicry gene sets)
#
# Outputs (in results/ssgsea/):
#   - ssgsea_scores.tsv (pathway scores per sample)
#   - ssgsea_plots.pdf (heatmap visualization)
#
# Usage: Rscript scripts/03_bulk_rnaseq_ssgsea.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(yaml)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

cat("=== ssGSEA Pathway Analysis ===\n\n")

# ------------------------------------------------------------------------------
# Set up paths
# ------------------------------------------------------------------------------

repo_root <- getwd()

# Input files
vst_path <- file.path(repo_root, "results", "rnaseq", "de", "vst_matrix.tsv")
immune_gmt_path <- file.path(repo_root, "references", "gene_sets", "immune.gmt")
viral_gmt_path <- file.path(repo_root, "references", "gene_sets", "viral_mimicry.gmt")
apoptosis_gmt_path <- file.path(repo_root, "references", "gene_sets", "apoptosis.gmt")

# Output directory
out_dir <- file.path(repo_root, "results", "ssgsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Output files
scores_path <- file.path(out_dir, "ssgsea_scores.tsv")
plots_path <- file.path(out_dir, "ssgsea_plots.pdf")
boxplots_path <- file.path(out_dir, "ssgsea_boxplots.pdf")

# Load config and metadata for condition information
config <- read_yaml(file.path(repo_root, "config", "config.yaml"))
meta_path <- file.path(repo_root, config$paths$metadata$rnaseq)
metadata <- fread(meta_path, sep = "\t")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

# ------------------------------------------------------------------------------
# Load VST-normalized expression matrix
# ------------------------------------------------------------------------------

cat("Loading VST-normalized expression matrix...\n")

if (!file.exists(vst_path)) {
  stop("VST matrix not found: ", vst_path, "\n  Run DESeq2 analysis first.")
}

vst_data <- fread(vst_path, sep = "\t")
cat("  Loaded expression data\n")

# First column is gene IDs
gene_col <- names(vst_data)[1]
genes <- vst_data[[gene_col]]

# Convert to matrix
expr_matrix <- as.matrix(vst_data[, -1, with = FALSE])
rownames(expr_matrix) <- genes

cat("  Genes:", nrow(expr_matrix), "\n")
cat("  Samples:", ncol(expr_matrix), "\n")
cat("  Gene IDs: symbols (from DESeq2)\n\n")

# ------------------------------------------------------------------------------
# Load gene sets from GMT files
# ------------------------------------------------------------------------------

cat("Loading gene sets...\n")

# Initialize empty list to hold all gene sets
all_gene_sets <- list()

# Load immune gene sets if file exists
if (file.exists(immune_gmt_path)) {
  cat("  Loading immune gene sets from:", immune_gmt_path, "\n")
  immune_gmt <- getGmt(immune_gmt_path)
  immune_list <- geneIds(immune_gmt)
  names(immune_list) <- names(immune_gmt)
  all_gene_sets <- c(all_gene_sets, immune_list)
  cat("    Loaded", length(immune_list), "immune gene sets\n")
} else {
  cat("  WARNING: Immune GMT file not found:", immune_gmt_path, "\n")
}

# Load viral mimicry gene sets if file exists
if (file.exists(viral_gmt_path)) {
  cat("  Loading viral mimicry gene sets from:", viral_gmt_path, "\n")
  viral_gmt <- getGmt(viral_gmt_path)
  viral_list <- geneIds(viral_gmt)
  names(viral_list) <- names(viral_gmt)
  all_gene_sets <- c(all_gene_sets, viral_list)
  cat("    Loaded", length(viral_list), "viral mimicry gene sets\n")
} else {
  cat("  WARNING: Viral mimicry GMT file not found:", viral_gmt_path, "\n")
}

# Load apoptosis gene sets if file exists
if (file.exists(apoptosis_gmt_path)) {
  cat("  Loading apoptosis gene sets from:", apoptosis_gmt_path, "\n")
  apoptosis_gmt <- getGmt(apoptosis_gmt_path)
  apoptosis_list <- geneIds(apoptosis_gmt)
  names(apoptosis_list) <- names(apoptosis_gmt)
  all_gene_sets <- c(all_gene_sets, apoptosis_list)
  cat("    Loaded", length(apoptosis_list), "apoptosis gene sets\n")
} else {
  cat("  WARNING: Apoptosis GMT file not found:", apoptosis_gmt_path, "\n")
}

# Check that we have gene sets to analyze
if (length(all_gene_sets) == 0) {
  stop("No gene sets loaded. Please provide GMT files in references/gene_sets/")
}

cat("  Total gene sets:", length(all_gene_sets), "\n\n")

# ------------------------------------------------------------------------------
# Check gene overlap between expression data and gene sets
# ------------------------------------------------------------------------------

cat("Checking gene overlap...\n")

# Get all genes in gene sets
all_geneset_genes <- unique(unlist(all_gene_sets))

# Find overlap with expression matrix
genes_in_expr <- intersect(all_geneset_genes, rownames(expr_matrix))
overlap_pct <- round(100 * length(genes_in_expr) / length(all_geneset_genes), 1)

cat("  Genes in gene sets:", length(all_geneset_genes), "\n")
cat("  Genes found in expression data:", length(genes_in_expr), "(", overlap_pct, "%)\n\n")

if (length(genes_in_expr) < 10) {
  stop("Too few genes overlap between gene sets and expression data.")
}

# ------------------------------------------------------------------------------
# Run ssGSEA using GSVA
# ------------------------------------------------------------------------------

cat("Running ssGSEA...\n")
cat("  This may take a few minutes...\n")

# Run ssGSEA
# method = "ssgsea" specifies single-sample GSEA
# kcdf = "Gaussian" is appropriate for continuous normalized data like VST
#
# Note: GSVA >= 1.46.0 changed the API. We handle both versions.
gsva_version <- packageVersion("GSVA")

if (gsva_version >= "1.46.0") {
  # New API: use ssgseaParam object
  ssgsea_param <- ssgseaParam(
    exprData = expr_matrix,
    geneSets = all_gene_sets,
    normalize = TRUE
  )
  ssgsea_results <- gsva(ssgsea_param, verbose = FALSE)
} else {
  # Old API: pass arguments directly
  ssgsea_results <- gsva(
    expr = expr_matrix,
    gset.idx.list = all_gene_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    verbose = FALSE
  )
}

cat("  ssGSEA complete\n")
cat("  Pathway scores computed:", nrow(ssgsea_results), "pathways x", ncol(ssgsea_results), "samples\n\n")

# ------------------------------------------------------------------------------
# Save ssGSEA scores
# ------------------------------------------------------------------------------

cat("Saving ssGSEA scores...\n")

# Convert to data frame with pathway names as first column
scores_df <- data.frame(
  pathway = rownames(ssgsea_results),
  ssgsea_results,
  check.names = FALSE
)

fwrite(scores_df, scores_path, sep = "\t")
cat("  Saved:", scores_path, "\n\n")

# ------------------------------------------------------------------------------
# Generate heatmap visualization
# ------------------------------------------------------------------------------

cat("Generating heatmap...\n")

# Set up color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Open PDF device
pdf(plots_path, width = 10, height = max(6, nrow(ssgsea_results) * 0.3))

# Create heatmap of ssGSEA scores
pheatmap(
  ssgsea_results,
  color = heatmap_colors,
  scale = "row",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "ssGSEA Pathway Scores\n(row-scaled)"
)

# Close PDF device
dev.off()

cat("  Saved:", plots_path, "\n\n")

# ------------------------------------------------------------------------------
# Generate boxplots by condition
# ------------------------------------------------------------------------------

cat("Generating boxplots by condition...\n")

# Reshape ssGSEA results for plotting
scores_long <- as.data.frame(t(ssgsea_results))
scores_long$sample_id <- rownames(scores_long)
scores_long <- merge(scores_long, metadata[, c("sample_id", "condition")], by = "sample_id")

# Melt to long format
scores_melt <- melt(as.data.table(scores_long),
                    id.vars = c("sample_id", "condition"),
                    variable.name = "pathway",
                    value.name = "score")

# Create faceted boxplot with all pathways
pdf(boxplots_path, width = 12, height = 8)

p <- ggplot(scores_melt, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 4) +
  labs(title = "ssGSEA Pathway Scores by Condition",
       x = "Condition", y = "Enrichment Score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 7))

# Add Wilcoxon p-values
conditions <- unique(scores_melt$condition)
if (length(conditions) == 2) {
  # Calculate p-values for each pathway
  pvals <- scores_melt[, .(pval = wilcox.test(score ~ condition)$p.value), by = pathway]
  pvals$label <- ifelse(pvals$pval < 0.001, "p < 0.001",
                        ifelse(pvals$pval < 0.01, sprintf("p = %.3f", pvals$pval),
                               sprintf("p = %.2f", pvals$pval)))

  # Add p-values as text annotation
  p <- p + geom_text(data = pvals,
                     aes(x = 1.5, y = Inf, label = label),
                     inherit.aes = FALSE,
                     vjust = 1.5, size = 2.5)
}

print(p)
dev.off()

cat("  Saved:", boxplots_path, "\n\n")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("=== ssGSEA Analysis Complete ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files created:\n")
cat("  - ssgsea_scores.tsv (pathway scores per sample)\n")
cat("  - ssgsea_plots.pdf (heatmap visualization)\n")
cat("  - ssgsea_boxplots.pdf (boxplots by condition)\n")

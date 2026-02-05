#!/usr/bin/env Rscript
# ==============================================================================
# 02_bulk_rnaseq_deseq2.R
# ==============================================================================
# Differential expression analysis of bulk RNA-seq data using DESeq2.
#
# Inputs:
#   - config/config.yaml (parameters)
#   - metadata/samples_rnaseq.tsv (sample information)
#   - results/rnaseq/counts/counts_raw.tsv (gene count matrix)
#
# Outputs (in results/rnaseq/de/):
#   - deseq2_results_{contrast}.tsv (full DE results per contrast)
#   - deseq2_results_{contrast}_sig.tsv (FDR < threshold, all significant)
#   - deseq2_results_{contrast}_sig_strong.tsv (FDR < threshold & |log2FC| > threshold)
#   - vst_matrix.tsv (variance-stabilized expression)
#   - pca_plot.pdf, volcano plots, MA plots, heatmaps
#   - session_info.txt
#
# Usage: Rscript scripts/01_bulk_rnaseq_deseq2.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(42)

repo_root <- getwd()

cat("=== Bulk RNA-seq Differential Expression Analysis ===\n\n")

# ------------------------------------------------------------------------------
# 1. Load configuration
# ------------------------------------------------------------------------------

cat("Loading configuration...\n")

config_path <- file.path(repo_root, "config", "config.yaml")
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
config <- read_yaml(config_path)

# Extract relevant parameters
de_params <- config$bulk_rnaseq$de
contrasts <- config$contrasts$rnaseq

# Output directory
out_dir <- file.path(repo_root, config$paths$results, "rnaseq", "de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("  FDR threshold:", de_params$fdr_threshold, "\n")
cat("  log2FC threshold (hypothesis testing):", de_params$lfc_threshold, "\n")
cat("  Note: Using lfcThreshold in DESeq2 tests H0: |log2FC| <= threshold\n")
cat("  Contrasts:", length(contrasts), "\n\n")

# ------------------------------------------------------------------------------
# 2. Load metadata
# ------------------------------------------------------------------------------

cat("Loading metadata...\n")

meta_path <- file.path(repo_root, config$paths$metadata$rnaseq)
if (!file.exists(meta_path)) {
  stop("Metadata file not found: ", meta_path)
}

metadata <- fread(meta_path, sep = "\t")
cat("  Samples in metadata:", nrow(metadata), "\n")

# Validate required columns
required_cols <- c("sample_id", "condition")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Convert to data.frame with sample_id as rownames
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

# Ensure condition is a factor
metadata$condition <- factor(metadata$condition)

cat("  Conditions:", paste(levels(metadata$condition), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# 3. Load count matrix
# ------------------------------------------------------------------------------

cat("Loading count matrix...\n")

counts_path <- file.path(repo_root, config$paths$results, "rnaseq", "counts", "counts_raw.tsv")
if (!file.exists(counts_path)) {
  stop("Count matrix not found: ", counts_path, "\n  Run alignment stage first.")
}

counts_raw <- fread(counts_path, sep = "\t")
cat("  Genes in count matrix:", nrow(counts_raw), "\n")

# First column should be gene IDs
gene_col <- names(counts_raw)[1]
genes <- counts_raw[[gene_col]]

# Convert to matrix (exclude gene ID column)
count_matrix <- as.matrix(counts_raw[, -1, with = FALSE])
rownames(count_matrix) <- genes

cat("  Samples in count matrix:", ncol(count_matrix), "\n\n")

# ------------------------------------------------------------------------------
# 4. Convert Ensembl IDs to gene symbols (if needed)
# ------------------------------------------------------------------------------

# Expected input: ENSEMBL IDs with version numbers (from preprocessing)
# Format: ENSG00000223972.5
cat("Converting Ensembl IDs to gene symbols...\n")

# Strip version suffixes (ENSG00000123456.7 -> ENSG00000123456)
ensembl_ids <- gsub("\\.\\d+$", "", rownames(count_matrix))

# Map Ensembl -> Symbol, Gene Name, and Gene Type (biotype)
symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_names <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                     column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
gene_types <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                     column = "GENETYPE", keytype = "ENSEMBL", multiVals = "first")

# Create gene annotation table
gene_info <- data.frame(
  gene_id = ensembl_ids,
  gene_symbol = as.character(symbols),
  gene_name = as.character(gene_names),
  gene_biotype = as.character(gene_types),
  stringsAsFactors = FALSE
)

# Count mappings
n_mapped <- sum(!is.na(gene_info$gene_symbol) & gene_info$gene_symbol != "")
cat("  Mapped", n_mapped, "of", nrow(gene_info), "genes to symbols\n")

# Filter to genes with valid symbols
valid_idx <- !is.na(gene_info$gene_symbol) & gene_info$gene_symbol != ""
count_matrix <- count_matrix[valid_idx, , drop = FALSE]
gene_info <- gene_info[valid_idx, ]

# Handle duplicate symbols: keep row with highest mean expression
mean_expr <- rowMeans(count_matrix)
gene_info$mean_expr <- mean_expr
gene_info <- gene_info[order(-gene_info$mean_expr), ]
count_matrix <- count_matrix[order(-mean_expr), ]

dup_symbols <- duplicated(gene_info$gene_symbol)
n_dups <- sum(dup_symbols)
if (n_dups > 0) {
  cat("  Removing", n_dups, "duplicate symbols (keeping highest expression)\n")
  count_matrix <- count_matrix[!dup_symbols, ]
  gene_info <- gene_info[!dup_symbols, ]
}

# Set rownames to gene symbols
rownames(count_matrix) <- gene_info$gene_symbol
rownames(gene_info) <- gene_info$gene_symbol
gene_info$mean_expr <- NULL

cat("  Final gene count:", nrow(count_matrix), "\n\n")

# ------------------------------------------------------------------------------
# 5. Match samples between metadata and counts
# ------------------------------------------------------------------------------

cat("Matching samples...\n")

# Find common samples
samples_meta <- rownames(metadata)
samples_counts <- colnames(count_matrix)
common_samples <- intersect(samples_meta, samples_counts)

if (length(common_samples) == 0) {
  stop("No matching samples between metadata and count matrix")
}

if (length(common_samples) < length(samples_meta)) {
  missing <- setdiff(samples_meta, samples_counts)
  cat("  WARNING: Samples in metadata but not in counts:", paste(missing, collapse = ", "), "\n")
}

if (length(common_samples) < length(samples_counts)) {
  extra <- setdiff(samples_counts, samples_meta)
  cat("  WARNING: Samples in counts but not in metadata:", paste(extra, collapse = ", "), "\n")
}

# Subset and order consistently
metadata <- metadata[common_samples, , drop = FALSE]
count_matrix <- count_matrix[, common_samples, drop = FALSE]

cat("  Using", length(common_samples), "matched samples\n\n")

# ------------------------------------------------------------------------------
# 6. Create DESeqDataSet
# ------------------------------------------------------------------------------

cat("Creating DESeqDataSet...\n")

# Build DESeqDataSet with design ~ condition
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

cat("  Initial genes:", nrow(dds), "\n")

# ------------------------------------------------------------------------------
# 7. Filter lowly expressed genes
# ------------------------------------------------------------------------------

cat("Filtering lowly expressed genes...\n")

min_counts <- de_params$min_counts

# Keep genes with at least min_counts total across all samples
keep <- rowSums(counts(dds)) >= min_counts
dds <- dds[keep, ]

cat("  Genes after filtering (>=", min_counts, "total counts):", nrow(dds), "\n\n")

# ------------------------------------------------------------------------------
# 8. Set reference level for contrasts
# ------------------------------------------------------------------------------

# For each contrast, the denominator is the reference
# Use the first contrast's denominator as the reference level
if (length(contrasts) > 0) {
  ref_level <- contrasts[[1]]$denominator
  if (ref_level %in% levels(dds$condition)) {
    dds$condition <- relevel(dds$condition, ref = ref_level)
    cat("Reference level set to:", ref_level, "\n\n")
  }
}

# ------------------------------------------------------------------------------
# 9. Run DESeq2
# ------------------------------------------------------------------------------

cat("Running DESeq2...\n")

# Use 'local' dispersion fit (better for n < 6)
# This adapts to gene-specific dispersion patterns better than parametric fit
dds <- DESeq(dds, fitType = 'local')

cat("  Dispersion estimation: complete (fitType='local' for n<6)\n")
cat("  Wald test: complete\n\n")

# Save DESeqDataSet object
dds_path <- file.path(out_dir, "dds.rds")
saveRDS(dds, dds_path)
cat("Saved DESeqDataSet to:", dds_path, "\n\n")

# ------------------------------------------------------------------------------
# 10. Extract and save results for each contrast
# ------------------------------------------------------------------------------

cat("Extracting results for each contrast...\n\n")

for (contrast in contrasts) {
  contrast_name <- contrast$name
  numerator <- contrast$numerator
  denominator <- contrast$denominator

  cat("  Contrast:", contrast_name, "(", numerator, "vs", denominator, ")\n")

  # Extract results using traditional approach optimized for small sample sizes
  # This tests H0: log2FC = 0 (any non-zero change)
  # Effect size threshold (lfc_threshold) will be used for filtering/prioritization, not hypothesis testing
  #
  # Special settings for n=3 replicates per group:
  # - independentFiltering=TRUE: ENABLE automatic filtering (matches raw script behavior)
  #   Paradoxically, this can INCREASE power by removing very noisy low-abundance genes
  #   which reduces the multiple testing burden on remaining genes
  # - cooksCutoff=FALSE: Disable Cook's distance outlier filtering
  #   (recommended by DESeq2 manual for n < 6, one outlier = 33% of data)
  res <- results(dds,
                 contrast = c("condition", numerator, denominator),
                 alpha = de_params$fdr_threshold,
                 independentFiltering = TRUE,
                 cooksCutoff = FALSE)

  cat("    NOTE: Using traditional testing approach (H0: log2FC = 0)\n")
  cat("    Independent filtering: ENABLED (matches raw script behavior)\n")
  cat("    Cook's distance filtering: DISABLED (recommended for n<6)\n")
  cat("    Effect size threshold |log2FC| >", de_params$lfc_threshold, "used for prioritization only\n")

  # Convert to data.frame - rownames are already gene symbols
  res_df <- as.data.frame(res)
  res_df$gene_symbol <- rownames(res_df)

  # Add gene_id, gene_name, and gene_biotype from gene_info table
  res_df$gene_id <- gene_info[res_df$gene_symbol, "gene_id"]
  res_df$gene_name <- gene_info[res_df$gene_symbol, "gene_name"]
  res_df$gene_biotype <- gene_info[res_df$gene_symbol, "gene_biotype"]

  # Reorder columns: gene info first, then statistics
  res_df <- res_df[, c("gene_symbol", "gene_id", "gene_name", "gene_biotype", "baseMean",
                        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

  # Sort by p-value
  res_df <- res_df[order(res_df$pvalue), ]

  # Count significant genes at different thresholds
  n_sig_fdr <- sum(res_df$padj < de_params$fdr_threshold, na.rm = TRUE)
  n_sig_strong <- sum(res_df$padj < de_params$fdr_threshold &
                      abs(res_df$log2FoldChange) > de_params$lfc_threshold, na.rm = TRUE)

  cat("    Genes tested:", nrow(res_df), "\n")
  cat("    Significant (FDR <", de_params$fdr_threshold, "):", n_sig_fdr, "\n")
  cat("    Strong effect (FDR <", de_params$fdr_threshold, "& |log2FC| >", de_params$lfc_threshold, "):", n_sig_strong, "\n")
  cat("    Note: Effect size threshold used for prioritization, not significance testing\n")

  # Report gene biotype distribution if available
  if (any(!is.na(res_df$gene_biotype))) {
    biotype_counts <- table(res_df$gene_biotype, useNA = "ifany")
    cat("    Gene biotype distribution (all tested genes):\n")
    for (bt in names(sort(biotype_counts, decreasing = TRUE))) {
      if (is.na(bt)) {
        cat("      Unknown:", biotype_counts[is.na(names(biotype_counts))], "\n")
      } else {
        cat("     ", bt, ":", biotype_counts[bt], "\n")
      }
    }
    # Report protein-coding percentage
    if ("protein-coding" %in% names(biotype_counts)) {
      pct_protein_coding <- round(100 * biotype_counts["protein-coding"] / sum(biotype_counts), 1)
      cat("    Protein-coding genes:", biotype_counts["protein-coding"],
          "(", pct_protein_coding, "%)\n")
    }
    cat("    Note: Consider filtering to protein_coding genes for cleaner DE analysis\n")
  }

  # Save full results
  full_path <- file.path(out_dir, paste0("deseq2_results_", contrast_name, ".tsv"))
  fwrite(res_df, full_path, sep = "\t")
  cat("    Saved:", basename(full_path), "\n")

  # Save filtered significant results (two files for different purposes)

  # File 1: All significant genes (FDR < threshold, any effect size)
  # Use for: Comprehensive gene lists, pathway analysis, GO enrichment
  res_sig_all <- res_df[!is.na(res_df$padj) &
                        res_df$padj < de_params$fdr_threshold, ]

  sig_all_path <- file.path(out_dir, paste0("deseq2_results_", contrast_name, "_sig.tsv"))
  fwrite(res_sig_all, sig_all_path, sep = "\t")
  cat("    Saved:", basename(sig_all_path), "(", nrow(res_sig_all), "genes - all significant)\n")

  # File 2: Strong effect genes (FDR < threshold AND |log2FC| > threshold)
  # Use for: Prioritized gene lists, focused analysis, validation targets
  res_sig_strong <- res_df[!is.na(res_df$padj) &
                            res_df$padj < de_params$fdr_threshold &
                            abs(res_df$log2FoldChange) > de_params$lfc_threshold, ]

  sig_strong_path <- file.path(out_dir, paste0("deseq2_results_", contrast_name, "_sig_strong.tsv"))
  fwrite(res_sig_strong, sig_strong_path, sep = "\t")
  cat("    Saved:", basename(sig_strong_path), "(", nrow(res_sig_strong), "genes - strong effect)\n")

  # --------------------------------------------------------------------------
  # Volcano plot
  # --------------------------------------------------------------------------

  cat("    Generating volcano plot...\n")

  # Prepare data for plotting
  # Note: Visual categories still use both FDR and log2FC thresholds for clarity
  # Even with lfcThreshold in results(), we separate Up/Down visually
  plot_df <- res_df
  plot_df$significant <- "Not Significant"
  plot_df$significant[!is.na(plot_df$padj) &
                      plot_df$padj < de_params$fdr_threshold &
                      plot_df$log2FoldChange > de_params$lfc_threshold] <- "Up"
  plot_df$significant[!is.na(plot_df$padj) &
                      plot_df$padj < de_params$fdr_threshold &
                      plot_df$log2FoldChange < -de_params$lfc_threshold] <- "Down"
  plot_df$significant <- factor(plot_df$significant, levels = c("Down", "Not Significant", "Up"))

  # Identify top 5 upregulated and top 5 downregulated genes for labeling
  sig_up <- plot_df[plot_df$significant == "Up" & !is.na(plot_df$padj), ]
  sig_down <- plot_df[plot_df$significant == "Down" & !is.na(plot_df$padj), ]

  top_up <- head(sig_up[order(sig_up$padj), ], 5)
  top_down <- head(sig_down[order(sig_down$padj), ], 5)
  top_genes <- rbind(top_up, top_down)

  # Use gene_symbol for labels if available, otherwise use gene ID
  top_genes$label <- ifelse(!is.na(top_genes$gene_symbol) & top_genes$gene_symbol != "",
                            top_genes$gene_symbol,
                            top_genes$gene)

  volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "blue", "Not Significant" = "grey", "Up" = "red")) +
    geom_vline(xintercept = c(-de_params$lfc_threshold, de_params$lfc_threshold),
               linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(de_params$fdr_threshold),
               linetype = "dashed", color = "black", alpha = 0.5) +
    geom_text_repel(data = top_genes,
                    aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    segment.color = "grey50",
                    show.legend = FALSE) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  volcano_path <- file.path(out_dir, paste0("volcano_", contrast_name, ".pdf"))
  ggsave(volcano_path, volcano_plot, width = 8, height = 6)
  cat("    Saved:", basename(volcano_path), "\n")

  # --------------------------------------------------------------------------
  # MA plot
  # --------------------------------------------------------------------------

  cat("    Generating MA plot...\n")

  ma_plot <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "blue", "Not Significant" = "grey", "Up" = "red")) +
    geom_hline(yintercept = c(-de_params$lfc_threshold, 0, de_params$lfc_threshold),
               linetype = c("dashed", "solid", "dashed"),
               color = c("black", "black", "black"),
               alpha = c(0.5, 1, 0.5)) +
    labs(
      title = paste("MA Plot:", contrast_name),
      x = "log10(Mean Expression + 1)",
      y = "log2 Fold Change",
      color = "Significance"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  ma_path <- file.path(out_dir, paste0("ma_plot_", contrast_name, ".pdf"))
  ggsave(ma_path, ma_plot, width = 8, height = 6)
  cat("    Saved:", basename(ma_path), "\n\n")
}

# ------------------------------------------------------------------------------
# 11. Generate normalized expression matrix (VST)
# ------------------------------------------------------------------------------

cat("Generating variance-stabilized expression matrix...\n")

# Use rlog for small sample sizes (n < 30), VST for larger datasets
n_samples <- ncol(dds)
if (n_samples < 30) {
  cat("  Using rlog transformation (n =", n_samples, "< 30 samples)\n")
  cat("  Note: rlog provides better variance stabilization for small sample sizes\n")
  vst <- rlog(dds, blind = FALSE)
} else {
  cat("  Using VST transformation (n =", n_samples, ">= 30 samples)\n")
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
}
vst_matrix <- assay(vst)

# Add gene column and save
vst_df <- data.frame(gene = rownames(vst_matrix), vst_matrix, check.names = FALSE)
vst_path <- file.path(out_dir, "vst_matrix.tsv")
fwrite(vst_df, vst_path, sep = "\t")

cat("  Saved:", vst_path, "\n\n")

# ------------------------------------------------------------------------------
# 12. QC: PCA plot
# ------------------------------------------------------------------------------

cat("Generating PCA plot...\n")

# Calculate PCA
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Create plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of VST-normalized expression") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Save
pca_path <- file.path(out_dir, "pca_plot.pdf")
ggsave(pca_path, pca_plot, width = 7, height = 5)
cat("  Saved:", pca_path, "\n\n")

# ------------------------------------------------------------------------------
# 13. QC: Sample distance heatmap
# ------------------------------------------------------------------------------

cat("Generating sample distance heatmap...\n")

# Calculate sample distances
sample_dists <- dist(t(assay(vst)))
sample_dist_matrix <- as.matrix(sample_dists)

# Annotation for heatmap
annotation_df <- data.frame(
  Condition = metadata$condition,
  row.names = rownames(metadata)
)

# Color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Create heatmap
heatmap_path <- file.path(out_dir, "sample_distance_heatmap.pdf")
pdf(heatmap_path, width = 8, height = 7)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = annotation_df,
  color = colors,
  main = "Sample-to-sample distances"
)
dev.off()

cat("  Saved:", heatmap_path, "\n\n")

# ------------------------------------------------------------------------------
# 14. Top 50 DE genes heatmap
# ------------------------------------------------------------------------------

cat("Generating top 50 DE genes heatmap...\n")

# Get top 50 DE genes by adjusted p-value from first contrast
res_all <- results(dds, contrast = c("condition", contrasts[[1]]$numerator, contrasts[[1]]$denominator))
top50 <- head(order(res_all$padj), 50)
top50_genes <- rownames(res_all)[top50]

# Extract VST values (rownames are already gene symbols)
top50_vst <- assay(vst)[top50_genes, ]

# Create heatmap
top50_path <- file.path(out_dir, "top50_de_heatmap.pdf")
pdf(top50_path, width = 8, height = 12)
pheatmap(top50_vst, scale = "row", annotation_col = annotation_df,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 8, main = "Top 50 Differentially Expressed Genes")
dev.off()

cat("  Saved:", top50_path, "\n\n")

# ------------------------------------------------------------------------------
# 15. Marker gene expression plots
# ------------------------------------------------------------------------------

cat("Generating marker gene expression plots...\n")

# Create markers output directory
markers_dir <- file.path(repo_root, config$paths$results, "rnaseq", "markers")
dir.create(markers_dir, recursive = TRUE, showWarnings = FALSE)

# Load marker genes from config
marker_lists <- config$marker_genes
all_markers <- unique(unlist(marker_lists))
cat("  Marker genes defined:", length(all_markers), "\n")

# Find which markers are present (rownames are already gene symbols)
markers_found <- all_markers[all_markers %in% rownames(vst_matrix)]
markers_missing <- setdiff(all_markers, markers_found)

if (length(markers_missing) > 0) {
  cat("  WARNING: Markers not found in data:", paste(markers_missing, collapse = ", "), "\n")
}
cat("  Markers found:", length(markers_found), "\n")

if (length(markers_found) > 0) {
  # Extract VST values for markers (rownames are already gene symbols)
  marker_vst <- vst_matrix[markers_found, , drop = FALSE]

  # --- Marker Heatmap ---
  marker_heatmap_path <- file.path(markers_dir, "marker_heatmap.pdf")
  pdf(marker_heatmap_path, width = 8, height = max(5, nrow(marker_vst) * 0.3))
  pheatmap(marker_vst, scale = "row", annotation_col = annotation_df,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize_row = 10, main = "Marker Gene Expression")
  dev.off()
  cat("  Saved:", marker_heatmap_path, "\n")

  # --- Marker Boxplots ---
  # Reshape for ggplot
  marker_long <- as.data.frame(t(marker_vst))
  marker_long$sample_id <- rownames(marker_long)
  marker_long$condition <- metadata[marker_long$sample_id, "condition"]

  marker_melt <- melt(as.data.table(marker_long),
                      id.vars = c("sample_id", "condition"),
                      variable.name = "gene",
                      value.name = "expression")

  # Create boxplots
  marker_boxplot_path <- file.path(markers_dir, "marker_boxplots.pdf")
  pdf(marker_boxplot_path, width = 10, height = 8)

  p <- ggplot(marker_melt, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    labs(title = "Marker Gene Expression by Condition",
         x = "Condition", y = "VST-normalized Expression") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  print(p)
  dev.off()
  cat("  Saved:", marker_boxplot_path, "\n")
}

cat("\n")

# ------------------------------------------------------------------------------
# 16. Save session info
# ------------------------------------------------------------------------------

cat("Saving session info...\n")

session_path <- file.path(out_dir, "session_info.txt")
sink(session_path)
cat("Bulk RNA-seq DESeq2 Analysis\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Parameters:\n")
cat("  FDR threshold:", de_params$fdr_threshold, "\n")
cat("  log2FC threshold (reporting):", de_params$lfc_threshold, "\n")
cat("  Min counts filter:", de_params$min_counts, "\n\n")
print(sessionInfo())
sink()

cat("  Saved:", session_path, "\n\n")

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------

cat("=== Analysis complete ===\n")
cat("Output directory:", out_dir, "\n")

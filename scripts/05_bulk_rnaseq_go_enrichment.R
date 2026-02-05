#!/usr/bin/env Rscript
# ==============================================================================
# 05_bulk_rnaseq_go_enrichment.R
# ==============================================================================
# Gene Ontology (GO) and KEGG pathway enrichment analysis of DE genes.
#
# Inputs:
#   - config/config.yaml (parameters)
#   - results/rnaseq/de/deseq2_results_{contrast}.tsv (DE results)
#
# Outputs (in results/rnaseq/go_enrichment/):
#   - go_enrichment_BP_{contrast}.tsv (Biological Process)
#   - go_enrichment_MF_{contrast}.tsv (Molecular Function)
#   - go_enrichment_CC_{contrast}.tsv (Cellular Component)
#   - kegg_enrichment_{contrast}.tsv (KEGG pathways)
#   - Various visualization PDFs (dot plots, bar plots, networks)
#
# Usage: Rscript scripts/05_bulk_rnaseq_go_enrichment.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
})

set.seed(42)

repo_root <- getwd()

cat("=== GO Term and KEGG Pathway Enrichment Analysis ===\n\n")

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

# Input/output directories
de_dir <- file.path(repo_root, config$paths$results, "rnaseq", "de")
out_dir <- file.path(repo_root, config$paths$results, "rnaseq", "go_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("  FDR threshold:", de_params$fdr_threshold, "\n")
cat("  Output directory:", out_dir, "\n\n")

# ------------------------------------------------------------------------------
# 2. Enrichment parameters
# ------------------------------------------------------------------------------

# GO/KEGG enrichment settings
PVALUE_CUTOFF <- 0.05        # P-value threshold for enrichment
QVALUE_CUTOFF <- 0.05        # FDR threshold for enrichment
MIN_GENE_SET_SIZE <- 10      # Minimum genes in a GO term
MAX_GENE_SET_SIZE <- 500     # Maximum genes in a GO term
TOP_N_PLOTS <- 20            # Number of terms to show in plots

cat("Enrichment parameters:\n")
cat("  P-value cutoff:", PVALUE_CUTOFF, "\n")
cat("  Q-value (FDR) cutoff:", QVALUE_CUTOFF, "\n")
cat("  Gene set size:", MIN_GENE_SET_SIZE, "-", MAX_GENE_SET_SIZE, "\n")
cat("  Top N for plots:", TOP_N_PLOTS, "\n\n")

# ------------------------------------------------------------------------------
# 3. Process each contrast
# ------------------------------------------------------------------------------

cat("Processing contrasts...\n\n")

for (contrast in contrasts) {
  contrast_name <- contrast$name

  cat("=== Contrast:", contrast_name, "===\n")

  # ----------------------------------------------------------------------------
  # Load DE results
  # ----------------------------------------------------------------------------

  cat("  Loading DE results...\n")

  de_file <- file.path(de_dir, paste0("deseq2_results_", contrast_name, ".tsv"))
  if (!file.exists(de_file)) {
    cat("    WARNING: DE results not found:", de_file, "\n")
    cat("    Skipping this contrast.\n\n")
    next
  }

  de_results <- fread(de_file, sep = "\t")
  cat("    Total genes:", nrow(de_results), "\n")

  # ----------------------------------------------------------------------------
  # Extract significant genes
  # ----------------------------------------------------------------------------

  cat("  Extracting significant genes...\n")

  # Filter to significant genes (FDR < threshold)
  sig_genes <- de_results[!is.na(padj) & padj < de_params$fdr_threshold]

  if (nrow(sig_genes) == 0) {
    cat("    WARNING: No significant genes found for", contrast_name, "\n")
    cat("    Skipping enrichment analysis.\n\n")
    next
  }

  cat("    Significant genes:", nrow(sig_genes), "\n")

  # Separate up- and down-regulated genes
  sig_up <- sig_genes[log2FoldChange > 0]
  sig_down <- sig_genes[log2FoldChange < 0]

  cat("    Up-regulated:", nrow(sig_up), "\n")
  cat("    Down-regulated:", nrow(sig_down), "\n")

  # ----------------------------------------------------------------------------
  # Convert gene symbols to Entrez IDs
  # ----------------------------------------------------------------------------

  cat("  Converting gene symbols to Entrez IDs...\n")

  # For enrichment, we need Entrez IDs
  # Try to get them from gene_id column (if Ensembl) or map from symbols

  # Check if gene_id column contains Ensembl IDs
  if (any(grepl("^ENSG", na.omit(sig_genes$gene_id)))) {
    # Map from Ensembl to Entrez
    ensembl_ids <- gsub("\\.\\d+$", "", sig_genes$gene_id)
    entrez_ids <- mapIds(org.Hs.eg.db, keys = ensembl_ids,
                         column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  } else {
    # Map from gene symbols to Entrez
    entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes$gene_symbol,
                         column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  }

  # Add Entrez IDs to significant genes
  sig_genes$entrez_id <- as.character(entrez_ids)

  # Filter to genes with valid Entrez IDs
  sig_genes_entrez <- sig_genes[!is.na(entrez_id)]

  n_mapped <- nrow(sig_genes_entrez)
  n_unmapped <- nrow(sig_genes) - n_mapped

  cat("    Mapped to Entrez:", n_mapped, "\n")
  if (n_unmapped > 0) {
    cat("    WARNING: Could not map", n_unmapped, "genes to Entrez IDs\n")
  }

  if (n_mapped < 5) {
    cat("    WARNING: Too few genes mapped for enrichment (n =", n_mapped, ")\n")
    cat("    Skipping enrichment analysis.\n\n")
    next
  }

  # Gene lists for enrichment
  all_genes_entrez <- sig_genes_entrez$entrez_id
  up_genes_entrez <- sig_genes_entrez[log2FoldChange > 0]$entrez_id
  down_genes_entrez <- sig_genes_entrez[log2FoldChange < 0]$entrez_id

  # Universe: all genes tested in DE analysis
  universe_symbols <- de_results[!is.na(padj)]$gene_symbol
  if (any(grepl("^ENSG", na.omit(de_results$gene_id)))) {
    universe_ensembl <- gsub("\\.\\d+$", "", de_results[!is.na(padj)]$gene_id)
    universe_entrez <- mapIds(org.Hs.eg.db, keys = universe_ensembl,
                              column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  } else {
    universe_entrez <- mapIds(org.Hs.eg.db, keys = universe_symbols,
                              column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  }
  universe_entrez <- unique(na.omit(as.character(universe_entrez)))

  cat("    Universe size:", length(universe_entrez), "genes\n")

  # ----------------------------------------------------------------------------
  # GO Enrichment: Biological Process
  # ----------------------------------------------------------------------------

  cat("  Running GO enrichment (Biological Process)...\n")

  tryCatch({
    ego_bp <- enrichGO(
      gene          = all_genes_entrez,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = PVALUE_CUTOFF,
      qvalueCutoff  = QVALUE_CUTOFF,
      minGSSize     = MIN_GENE_SET_SIZE,
      maxGSSize     = MAX_GENE_SET_SIZE,
      readable      = TRUE
    )

    if (nrow(ego_bp@result) > 0) {
      cat("    Significant BP terms:", nrow(ego_bp@result), "\n")

      # Save results
      bp_file <- file.path(out_dir, paste0("go_enrichment_BP_", contrast_name, ".tsv"))
      fwrite(as.data.frame(ego_bp), bp_file, sep = "\t")
      cat("    Saved:", basename(bp_file), "\n")

      # Dot plot
      if (nrow(ego_bp@result) > 0) {
        n_plot <- min(TOP_N_PLOTS, nrow(ego_bp@result))
        p <- dotplot(ego_bp, showCategory = n_plot) +
          ggtitle(paste("GO Biological Process -", contrast_name))

        plot_file <- file.path(out_dir, paste0("go_dotplot_BP_", contrast_name, ".pdf"))
        ggsave(plot_file, p, width = 10, height = max(6, n_plot * 0.3))
        cat("    Saved:", basename(plot_file), "\n")

        # Bar plot
        p_bar <- barplot(ego_bp, showCategory = n_plot) +
          ggtitle(paste("GO Biological Process -", contrast_name))

        plot_file_bar <- file.path(out_dir, paste0("go_barplot_BP_", contrast_name, ".pdf"))
        ggsave(plot_file_bar, p_bar, width = 10, height = max(6, n_plot * 0.3))
        cat("    Saved:", basename(plot_file_bar), "\n")
      }
    } else {
      cat("    No significant BP terms found.\n")
    }
  }, error = function(e) {
    cat("    ERROR in BP enrichment:", conditionMessage(e), "\n")
  })

  # ----------------------------------------------------------------------------
  # GO Enrichment: Molecular Function
  # ----------------------------------------------------------------------------

  cat("  Running GO enrichment (Molecular Function)...\n")

  tryCatch({
    ego_mf <- enrichGO(
      gene          = all_genes_entrez,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = PVALUE_CUTOFF,
      qvalueCutoff  = QVALUE_CUTOFF,
      minGSSize     = MIN_GENE_SET_SIZE,
      maxGSSize     = MAX_GENE_SET_SIZE,
      readable      = TRUE
    )

    if (nrow(ego_mf@result) > 0) {
      cat("    Significant MF terms:", nrow(ego_mf@result), "\n")

      # Save results
      mf_file <- file.path(out_dir, paste0("go_enrichment_MF_", contrast_name, ".tsv"))
      fwrite(as.data.frame(ego_mf), mf_file, sep = "\t")
      cat("    Saved:", basename(mf_file), "\n")

      # Dot plot
      if (nrow(ego_mf@result) > 0) {
        n_plot <- min(TOP_N_PLOTS, nrow(ego_mf@result))
        p <- dotplot(ego_mf, showCategory = n_plot) +
          ggtitle(paste("GO Molecular Function -", contrast_name))

        plot_file <- file.path(out_dir, paste0("go_dotplot_MF_", contrast_name, ".pdf"))
        ggsave(plot_file, p, width = 10, height = max(6, n_plot * 0.3))
        cat("    Saved:", basename(plot_file), "\n")
      }
    } else {
      cat("    No significant MF terms found.\n")
    }
  }, error = function(e) {
    cat("    ERROR in MF enrichment:", conditionMessage(e), "\n")
  })

  # ----------------------------------------------------------------------------
  # GO Enrichment: Cellular Component
  # ----------------------------------------------------------------------------

  cat("  Running GO enrichment (Cellular Component)...\n")

  tryCatch({
    ego_cc <- enrichGO(
      gene          = all_genes_entrez,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = PVALUE_CUTOFF,
      qvalueCutoff  = QVALUE_CUTOFF,
      minGSSize     = MIN_GENE_SET_SIZE,
      maxGSSize     = MAX_GENE_SET_SIZE,
      readable      = TRUE
    )

    if (nrow(ego_cc@result) > 0) {
      cat("    Significant CC terms:", nrow(ego_cc@result), "\n")

      # Save results
      cc_file <- file.path(out_dir, paste0("go_enrichment_CC_", contrast_name, ".tsv"))
      fwrite(as.data.frame(ego_cc), cc_file, sep = "\t")
      cat("    Saved:", basename(cc_file), "\n")

      # Dot plot
      if (nrow(ego_cc@result) > 0) {
        n_plot <- min(TOP_N_PLOTS, nrow(ego_cc@result))
        p <- dotplot(ego_cc, showCategory = n_plot) +
          ggtitle(paste("GO Cellular Component -", contrast_name))

        plot_file <- file.path(out_dir, paste0("go_dotplot_CC_", contrast_name, ".pdf"))
        ggsave(plot_file, p, width = 10, height = max(6, n_plot * 0.3))
        cat("    Saved:", basename(plot_file), "\n")
      }
    } else {
      cat("    No significant CC terms found.\n")
    }
  }, error = function(e) {
    cat("    ERROR in CC enrichment:", conditionMessage(e), "\n")
  })

  # ----------------------------------------------------------------------------
  # KEGG Pathway Enrichment
  # ----------------------------------------------------------------------------

  cat("  Running KEGG pathway enrichment...\n")

  tryCatch({
    kegg <- enrichKEGG(
      gene          = all_genes_entrez,
      universe      = universe_entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = PVALUE_CUTOFF,
      qvalueCutoff  = QVALUE_CUTOFF,
      minGSSize     = MIN_GENE_SET_SIZE,
      maxGSSize     = MAX_GENE_SET_SIZE
    )

    if (nrow(kegg@result) > 0) {
      # Convert Entrez IDs to gene symbols for readability
      kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

      cat("    Significant KEGG pathways:", nrow(kegg@result), "\n")

      # Save results
      kegg_file <- file.path(out_dir, paste0("kegg_enrichment_", contrast_name, ".tsv"))
      fwrite(as.data.frame(kegg), kegg_file, sep = "\t")
      cat("    Saved:", basename(kegg_file), "\n")

      # Dot plot
      if (nrow(kegg@result) > 0) {
        n_plot <- min(TOP_N_PLOTS, nrow(kegg@result))
        p <- dotplot(kegg, showCategory = n_plot) +
          ggtitle(paste("KEGG Pathways -", contrast_name))

        plot_file <- file.path(out_dir, paste0("kegg_dotplot_", contrast_name, ".pdf"))
        ggsave(plot_file, p, width = 10, height = max(6, n_plot * 0.3))
        cat("    Saved:", basename(plot_file), "\n")
      }
    } else {
      cat("    No significant KEGG pathways found.\n")
    }
  }, error = function(e) {
    cat("    ERROR in KEGG enrichment:", conditionMessage(e), "\n")
  })

  # ----------------------------------------------------------------------------
  # Enrichment map (if we have GO BP results with enough terms)
  # ----------------------------------------------------------------------------

  if (exists("ego_bp") && nrow(ego_bp@result) >= 5) {
    cat("  Generating enrichment map...\n")

    tryCatch({
      # Simplify GO terms to reduce redundancy
      ego_bp_simp <- simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

      if (nrow(ego_bp_simp@result) >= 3) {
        # Create pairwise similarity matrix
        ego_bp_simp <- pairwise_termsim(ego_bp_simp)

        # Enrichment map (network of related GO terms)
        p_emap <- emapplot(ego_bp_simp, showCategory = min(30, nrow(ego_bp_simp@result))) +
          ggtitle(paste("GO BP Enrichment Map -", contrast_name))

        plot_file <- file.path(out_dir, paste0("go_enrichmap_BP_", contrast_name, ".pdf"))
        ggsave(plot_file, p_emap, width = 12, height = 10)
        cat("    Saved:", basename(plot_file), "\n")
      }
    }, error = function(e) {
      cat("    WARNING: Could not generate enrichment map:", conditionMessage(e), "\n")
    })
  }

  # ----------------------------------------------------------------------------
  # Summary for this contrast
  # ----------------------------------------------------------------------------

  cat("  Completed enrichment for", contrast_name, "\n\n")
}

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("=== GO/KEGG Enrichment Analysis Complete ===\n")
cat("Output directory:", out_dir, "\n")
cat("\n")
cat("Output files per contrast:\n")
cat("  - go_enrichment_BP_{contrast}.tsv\n")
cat("  - go_enrichment_MF_{contrast}.tsv\n")
cat("  - go_enrichment_CC_{contrast}.tsv\n")
cat("  - kegg_enrichment_{contrast}.tsv\n")
cat("  - go_dotplot_BP_{contrast}.pdf\n")
cat("  - go_barplot_BP_{contrast}.pdf\n")
cat("  - go_dotplot_MF_{contrast}.pdf\n")
cat("  - go_dotplot_CC_{contrast}.pdf\n")
cat("  - kegg_dotplot_{contrast}.pdf\n")
cat("  - go_enrichmap_BP_{contrast}.pdf (if sufficient terms)\n")
cat("\n")

#!/usr/bin/env Rscript
# ==============================================================================
# Generate Synthetic Example Data for Demo Pipeline
# ==============================================================================
# Creates realistic RNA-seq count data with REAL ENSEMBL IDs
# ==============================================================================

set.seed(42)

library(org.Hs.eg.db)

cat("Generating example data with real ENSEMBL IDs...\n")

# Output path
out_dir <- "example_data"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Get real ENSEMBL IDs from database
# ------------------------------------------------------------------------------

# Get all ENSEMBL IDs
all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
cat("  Available ENSEMBL IDs:", length(all_ensembl), "\n")

# Sample 300 genes
set.seed(42)
selected_ensembl <- sample(all_ensembl, 300)

# Get gene symbols for these
gene_symbols <- mapIds(org.Hs.eg.db, keys = selected_ensembl,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Create mapping
gene_map <- data.frame(
  ensembl = selected_ensembl,
  symbol = as.character(gene_symbols),
  stringsAsFactors = FALSE
)

# Remove genes without symbols
gene_map <- gene_map[!is.na(gene_map$symbol), ]
n_genes <- nrow(gene_map)

cat("  Selected", n_genes, "genes with valid symbols\n")

# Define specific genes for differential expression
# Find ENSEMBL IDs for known cisplatin-responsive genes
known_genes <- c("IFIT1", "ISG15", "BAX", "CDKN1A", "GADD45A", "TP53I3",
                 "BBC3", "FAS", "CASP3", "MX1", "OAS1", "STAT1")

up_genes_idx <- which(gene_map$symbol %in% known_genes)
if (length(up_genes_idx) < 20) {
  # Add random genes to upregulated set
  remaining <- setdiff(1:n_genes, up_genes_idx)
  up_genes_idx <- c(up_genes_idx, sample(remaining, 20 - length(up_genes_idx)))
}

# Downregulated genes (random selection)
remaining <- setdiff(1:n_genes, up_genes_idx)
down_genes_idx <- sample(remaining, 10)

cat("  Upregulated genes:", length(up_genes_idx), "\n")
cat("  Downregulated genes:", length(down_genes_idx), "\n")

# ------------------------------------------------------------------------------
# Generate count matrix
# ------------------------------------------------------------------------------

samples <- c("C1", "C2", "C3", "V1", "V2", "V3")
conditions <- c(rep("cisplatin", 3), rep("vehicle", 3))

# Initialize count matrix
counts <- matrix(0, nrow = n_genes, ncol = 6)
rownames(counts) <- gene_map$ensembl
colnames(counts) <- samples

# Base expression levels (log-normal distributed)
base_means <- exp(rnorm(n_genes, mean = 6, sd = 2))
base_means <- pmax(base_means, 10)  # Minimum mean of 10

# Generate counts for each gene
for (i in 1:n_genes) {
  base <- base_means[i]

  # Determine fold change
  if (i %in% up_genes_idx) {
    # Upregulated in cisplatin: 2-4x higher
    fc <- runif(1, 2, 4)
    cis_mean <- base * fc
    veh_mean <- base
  } else if (i %in% down_genes_idx) {
    # Downregulated in cisplatin: 2-4x lower
    fc <- runif(1, 2, 4)
    cis_mean <- base
    veh_mean <- base * fc
  } else {
    # No change or small random variation
    cis_mean <- base * runif(1, 0.9, 1.1)
    veh_mean <- base * runif(1, 0.9, 1.1)
  }

  # Add biological variation
  cis_counts <- rnbinom(3, size = 10, mu = cis_mean * runif(3, 0.8, 1.2))
  veh_counts <- rnbinom(3, size = 10, mu = veh_mean * runif(3, 0.8, 1.2))

  counts[i, ] <- c(cis_counts, veh_counts)
}

# Ensure no negative counts
counts[counts < 0] <- 0

# ------------------------------------------------------------------------------
# Save count matrix with version numbers
# ------------------------------------------------------------------------------

# Add version numbers to match preprocessing output
# Use realistic version numbers (1-15)
gene_ids_with_version <- paste0(rownames(counts), ".", sample(1:15, n_genes, replace = TRUE))

counts_df <- data.frame(gene_id = gene_ids_with_version, counts, check.names = FALSE)

write.table(counts_df,
            file = file.path(out_dir, "counts_example.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("\nSaved:", file.path(out_dir, "counts_example.tsv"), "\n")
cat("  Genes:", nrow(counts), "\n")
cat("  Samples:", ncol(counts), "\n")

# ------------------------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------------------------

cat("\nSample total counts:\n")
print(colSums(counts))

cat("\nExample genes (first 10):\n")
example <- head(gene_map, 10)
example$ensembl_version <- head(gene_ids_with_version, 10)
print(example[, c("ensembl_version", "symbol")])

cat("\nDifferentially expressed genes include:\n")
up_symbols <- gene_map$symbol[up_genes_idx[1:min(10, length(up_genes_idx))]]
cat("  Upregulated:", paste(up_symbols, collapse = ", "), "\n")

cat("\nDone!\n")

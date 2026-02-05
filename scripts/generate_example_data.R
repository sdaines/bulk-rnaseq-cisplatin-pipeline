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
# Get real ENSEMBL IDs from database - prioritize marker/pathway genes
# ------------------------------------------------------------------------------

cat("  Collecting marker and pathway genes...\n")

# Define marker genes (from config) and key pathway genes
marker_genes <- c(
  # Immune markers
  "IFIT1", "ISG15", "MX1", "OAS1", "RSAD2", "PRF1", "NKG7",
  # Apoptosis markers
  "BAX", "BCL2", "CASP3", "TP53",
  # Housekeeping
  "GAPDH", "ACTB",
  # Additional response genes
  "CDKN1A", "GADD45A", "TP53I3", "BBC3", "FAS", "STAT1", "CXCL10",
  "IRF1", "IRF3", "IRF7", "DHX58", "DDX58", "IFIH1", "MAVS",
  "BID", "PUMA", "NOXA", "APAF1", "CYCS", "CASP9", "CASP8"
)

# Get ENSEMBL IDs for these priority genes
cat("    Mapping", length(marker_genes), "priority genes to ENSEMBL...\n")
priority_ensembl <- mapIds(org.Hs.eg.db, keys = marker_genes,
                           column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")

# Filter to valid mappings
priority_ensembl <- na.omit(priority_ensembl)
priority_symbols <- names(priority_ensembl)
cat("    Mapped", length(priority_ensembl), "priority genes\n")

# Get all available ENSEMBL IDs
all_ensembl <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

# Exclude priority genes from random sample pool
remaining_ensembl <- setdiff(all_ensembl, as.character(priority_ensembl))

# Sample additional genes to reach 300 total
n_additional <- 300 - length(priority_ensembl)
set.seed(42)
random_ensembl <- sample(remaining_ensembl, n_additional)

# Combine priority and random genes
selected_ensembl <- c(as.character(priority_ensembl), random_ensembl)

# Get gene symbols for all selected genes
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
cat("    Priority genes:", sum(gene_map$symbol %in% marker_genes), "\n")

# Define upregulated genes (prioritize immune/interferon response genes)
immune_genes <- c("IFIT1", "ISG15", "MX1", "OAS1", "RSAD2", "CXCL10",
                  "IRF1", "IRF3", "IRF7", "DHX58", "DDX58", "IFIH1", "MAVS", "STAT1")
apoptosis_genes <- c("BAX", "CASP3", "CASP9", "CASP8", "BID", "PUMA", "APAF1", "CYCS")

up_genes_idx <- which(gene_map$symbol %in% c(immune_genes, apoptosis_genes))

# Add random genes if needed to reach 20 upregulated
if (length(up_genes_idx) < 20) {
  remaining <- setdiff(1:n_genes, up_genes_idx)
  up_genes_idx <- c(up_genes_idx, sample(remaining, 20 - length(up_genes_idx)))
}

# Downregulated genes (random selection, excluding markers)
marker_idx <- which(gene_map$symbol %in% marker_genes)
remaining <- setdiff(setdiff(1:n_genes, up_genes_idx), marker_idx)
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

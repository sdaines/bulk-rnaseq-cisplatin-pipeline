# RUNBOOK.md

Step-by-step instructions for running the bulk RNA-seq pipeline.

---

## Prerequisites

### Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| STAR | 2.7+ | RNA-seq alignment |
| featureCounts | 2.0+ | Read quantification |
| samtools | 1.17+ | BAM processing |
| FastQC | 0.11+ | Quality control |
| R | 4.3+ | Statistical analysis |

### R Packages

```r
# Install required packages
install.packages(c("yaml", "data.table", "ggplot2", "pheatmap", "RColorBrewer"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GSVA", "GSEABase"))
```

---

## Setup

### 1. Clone Repository

```bash
git clone <repository-url>
cd bulk-rnaseq-pipeline
```

### 2. Download References

Download the reference genome and GTF annotation:

```bash
bash setup/download_references.sh
```

This downloads:
- GRCh38 primary assembly FASTA
- GENCODE v38 GTF annotation

### 3. Build STAR Index

Build the STAR genome index (requires ~32GB RAM):

```bash
bash setup/build_star_index.sh
```

### 4. Prepare Sample Metadata

Edit `metadata/samples_rnaseq.tsv` with your sample information:

| Column | Description |
|--------|-------------|
| sample_id | Unique sample identifier |
| fastq_r1 | R1 FASTQ filename |
| fastq_r2 | R2 FASTQ filename |
| condition | Experimental condition (e.g., vehicle, cisplatin) |
| batch | Batch identifier (optional) |

### 5. Place FASTQ Files

Copy your FASTQ files to `data/fastq/rnaseq/`.

### 6. Prepare Gene Sets (for ssGSEA)

Place GMT files in `references/gene_sets/`:
- `immune.gmt` — Immune-related pathway gene sets
- `viral_mimicry.gmt` — Viral mimicry / interferon response gene sets

---

## Running the Pipeline

### Full Pipeline

Run all stages (preprocessing, DESeq2, and ssGSEA):

```bash
bash scripts/run_pipeline.sh
```

### Individual Stages

Run only preprocessing (alignment and quantification):

```bash
bash scripts/run_pipeline.sh --stage preprocess
```

Run only DESeq2 analysis:

```bash
bash scripts/run_pipeline.sh --stage deseq2
```

Run only ssGSEA analysis:

```bash
bash scripts/run_pipeline.sh --stage ssgsea
```

### Dry Run

Preview commands without executing:

```bash
bash scripts/run_pipeline.sh --dry-run
```

---

## Pipeline Stages

### Stage 1: Preprocessing

**Script**: `scripts/01_bulk_rnaseq_preprocess.sh`

**What it does**:
1. Runs FastQC on raw FASTQ files
2. Aligns reads with STAR (2-pass mode)
3. Quantifies gene counts with featureCounts

**Outputs** (in `results/rnaseq/`):
- `fastqc/` — FastQC reports
- `bam/` — Aligned BAM files
- `counts/counts_raw.tsv` — Gene count matrix

### Stage 2: Differential Expression

**Script**: `scripts/02_bulk_rnaseq_deseq2.R`

**What it does**:
1. Loads count matrix and metadata
2. Filters lowly expressed genes
3. Runs DESeq2 differential expression
4. Generates QC plots (PCA, heatmap)

**Outputs** (in `results/rnaseq/de/`):
- `deseq2_results_cisplatin_vs_vehicle.tsv` — Full DE results
- `deseq2_results_cisplatin_vs_vehicle_sig.tsv` — Significant genes only
- `vst_matrix.tsv` — Normalized expression matrix
- `pca_plot.pdf` — PCA visualization
- `sample_distance_heatmap.pdf` — Sample clustering

### Stage 3: ssGSEA Pathway Analysis

**Script**: `scripts/03_bulk_rnaseq_ssgsea.R`

**What it does**:
Single-sample Gene Set Enrichment Analysis (ssGSEA) computes per-sample pathway activity scores. This allows you to quantify the activity of immune and viral mimicry pathways in each sample, independent of differential expression.

1. Loads VST-normalized expression from DESeq2
2. Loads gene sets from GMT files
3. Runs ssGSEA using GSVA
4. Generates heatmap visualization

**Outputs** (in `results/ssgsea/`):
- `ssgsea_scores.tsv` — Pathway activity scores per sample
- `ssgsea_plots.pdf` — Heatmap of pathway scores

---

## Output Summary

| File | Description |
|------|-------------|
| `results/rnaseq/counts/counts_raw.tsv` | Gene-level count matrix |
| `results/rnaseq/de/deseq2_results_*.tsv` | Differential expression results |
| `results/rnaseq/de/vst_matrix.tsv` | VST-normalized expression |
| `results/rnaseq/de/pca_plot.pdf` | Sample PCA plot |
| `results/ssgsea/ssgsea_scores.tsv` | ssGSEA pathway scores |
| `results/ssgsea/ssgsea_plots.pdf` | ssGSEA heatmap |
| `logs/` | Pipeline execution logs |

---

## Troubleshooting

### Common Issues

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| STAR index fails | Insufficient memory | Requires ~32GB RAM |
| featureCounts zero counts | Strandedness mismatch | Check library prep; update `strandedness` in config.yaml |
| DESeq2 errors | Missing samples | Ensure metadata matches count matrix columns |
| No BAM files created | Missing FASTQ files | Check FASTQ paths in metadata |
| ssGSEA no gene sets | GMT files missing | Place GMT files in `references/gene_sets/` |
| ssGSEA low overlap | Gene symbol mismatch | Ensure GMT uses same gene symbols as GTF |

### Log Files

Execution logs are saved to `logs/` with timestamps. Check these logs for error details when a stage fails.

---

## Configuration

All parameters are in `config/config.yaml`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `strandedness` | Library strandedness | unstranded |
| `star_threads` | STAR alignment threads | 8 |
| `min_counts` | Minimum gene counts for filtering | 10 |
| `fdr_threshold` | Significance threshold | 0.05 |
| `lfc_threshold` | Log2 fold change threshold | 1.0 |

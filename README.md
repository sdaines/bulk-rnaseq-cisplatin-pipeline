# Bulk RNA-seq Analysis Pipeline

A complete, reproducible pipeline for bulk RNA-seq differential expression and pathway analysis, optimized for small sample sizes (n=3-6 replicates).

[![Language](https://img.shields.io/badge/language-R%20%7C%20Bash-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## Overview

End-to-end analysis from FASTQ files to publication-ready results:

```
FASTQ → Alignment → Quantification → DESeq2 → Pathway Analysis → Plots
```

**Pipeline Stages:**
1. **Preprocessing**: Quality control, alignment (STAR), quantification (featureCounts)
2. **DESeq2**: Differential expression, PCA, volcano plots, heatmaps
3. **ssGSEA**: Single-sample pathway enrichment scores
4. **GSEA**: Gene set enrichment analysis
5. **GO/KEGG**: Pathway enrichment with clusterProfiler

---

## Quick Start (Demo)

Test the pipeline with synthetic data (no installation required except R/Bash):

```bash
# Clone repository
git clone https://github.com/your-username/bulk-rnaseq-pipeline.git
cd bulk-rnaseq-pipeline

# Run demo (uses 300 synthetic genes, completes in ~1 minute)
./scripts/run_pipeline.sh --demo
```

**Demo Output:**
- 42 differentially expressed genes
- Pathway enrichment analysis
- PCA, volcano plots, heatmaps
- All results in `results/`

---

## Features

- **Optimized for small n**: DESeq2 settings tuned for n=3-6 replicates
- **Complete workflow**: FASTQ → DE genes → pathway analysis
- **Custom pathways**: Apoptosis, interferon, viral mimicry gene sets
- **Reproducible**: Conda environment, version-controlled configs
- **Demo mode**: Test without real data

---

## Installation

### Option 1: Conda 

```bash
# Create environment from file
conda env create -f environment.yml
conda activate rnaseq

# Verify installation
Rscript -e "library(DESeq2); library(GSVA); library(fgsea)"
```

### Option 2: Manual Installation

**Required tools:**
- STAR (v2.7+)
- featureCounts/subread (v2.0+)
- samtools (v1.17+)
- FastQC, fastp, MultiQC
- R (v4.3+)

**R packages:**
```r
# CRAN packages
install.packages(c("yaml", "data.table", "ggplot2", "ggrepel", "pheatmap"))

# Bioconductor packages
BiocManager::install(c("DESeq2", "GSVA", "fgsea", "clusterProfiler",
                       "org.Hs.eg.db", "AnnotationDbi", "GSEABase"))
```

---

## Usage

### Full Pipeline

```bash
# 1. Download reference genome and build STAR index
./setup/setup_references.sh

# 2. Place FASTQ files in data/fastq/

# 3. Update metadata
nano metadata/samples_rnaseq.tsv

# 4. Run pipeline
./scripts/run_pipeline.sh
```

### Run Specific Stages

```bash
# Run only differential expression (if you already have counts)
./scripts/run_pipeline.sh --stage deseq2

# Run pathway analysis
./scripts/run_pipeline.sh --stage gsea
./scripts/run_pipeline.sh --stage go_enrichment
```

### Advanced Usage

```bash
# Use custom config
./scripts/run_pipeline.sh --config my_config.yaml

# Dry run (preview commands)
./scripts/run_pipeline.sh --dry-run

```

---

## Input Data

### Metadata Format

Edit `metadata/samples_rnaseq.tsv`:

| sample_id | fastq_r1 | fastq_r2 | condition | batch |
|-----------|----------|----------|-----------|-------|
| C1 | C1_R1.fastq.gz | C1_R2.fastq.gz | cisplatin | 1 |
| V1 | V1_R1.fastq.gz | V1_R2.fastq.gz | vehicle | 1 |

### Reference Files

Pipeline automatically downloads:
- GRCh38 primary assembly (GENCODE)
- GENCODE v38 comprehensive annotation

Or use your own by editing `config/config.yaml`

---

## Output Files

### DESeq2 Results (`results/rnaseq/de/`)

| File | Description |
|------|-------------|
| `deseq2_results_*.tsv` | Full differential expression results |
| `deseq2_results_*_sig.tsv` | Significant genes (FDR < 0.05) |
| `deseq2_results_*_sig_strong.tsv` | Strong effect (FDR < 0.05 & \|log2FC\| > 1) |
| `vst_matrix.tsv` | Variance-stabilized expression matrix |
| `volcano_*.pdf` | Volcano plot |
| `pca_plot.pdf` | PCA of samples |
| `top50_de_heatmap.pdf` | Heatmap of top DE genes |

### Pathway Analysis (`results/rnaseq/gsea/`, `results/ssgsea/`)

| File | Description |
|------|-------------|
| `gsea_results.tsv` | GSEA enrichment statistics |
| `ssgsea_scores.tsv` | Per-sample pathway scores |
| `enrichment_*.pdf` | GSEA enrichment plots |
| `go_enrichment_*.tsv` | GO term enrichment results |

---

## Configuration

Edit `config/config.yaml` to customize:

```yaml
bulk_rnaseq:
  alignment:
    strandedness: unstranded  # or "forward", "reverse"
    star_threads: 8

  de:
    min_counts: 5            # Pre-filter threshold
    fdr_threshold: 0.05      # Significance cutoff
    lfc_threshold: 1.0       # Fold-change for reporting

contrasts:
  rnaseq:
    - name: "treatment_vs_control"
      numerator: "treatment"
      denominator: "control"
```

---

## Example Outputs

See `docs/example_outputs/` for example figures from the demo:

<table>
<tr>
<td><img src="docs/example_outputs/volcano_example.pdf" width="200"/><br/><b>Volcano Plot</b></td>
<td><img src="docs/example_outputs/pca_example.pdf" width="200"/><br/><b>PCA Plot</b></td>
<td><img src="docs/example_outputs/heatmap_example.pdf" width="200"/><br/><b>DE Heatmap</b></td>
</tr>
</table>

*(Note: PDFs may not render on GitHub - download to view)*

---

## Key Methods

### DESeq2 Settings (Optimized for n=3)

- **Dispersion estimation**: `fitType='local'` 
- **Cook's distance**: Disabled (too harsh for n<6)
- **Pre-filtering**: Genes with <5 total counts removed

### Pathway Analysis

- **Custom gene sets**: 18 curated pathways (apoptosis, interferon, viral mimicry)
- **ssGSEA**: Sample-level pathway scores (GSVA package)
- **GSEA**: Preranked enrichment (fgsea package)
- **GO/KEGG**: Pathway enrichment (clusterProfiler)

See `docs/ASSUMPTIONS.md` for complete scientific rationale.

---

## Directory Structure

```
bulk-rnaseq-pipeline/
├── README.md                   # This file
├── CHANGELOG.md                # Version history
├── LICENSE                     # MIT license
├── environment.yml             # Conda environment
│
├── config/
│   ├── config.yaml            # Main configuration
│   └── config_example.yaml    # Demo configuration
│
├── scripts/
│   ├── run_pipeline.sh        # Main pipeline runner
│   ├── 01_bulk_rnaseq_preprocess.sh
│   ├── 02_bulk_rnaseq_deseq2.R
│   ├── 03_bulk_rnaseq_ssgsea.R
│   ├── 04_bulk_rnaseq_gsea.R
│   ├── 05_bulk_rnaseq_go_enrichment.R
│   ├── merge_individual_counts.R
│   └── generate_example_data.R
│
├── setup/
│   └── setup_references.sh    # Download genome/annotation
│
├── metadata/
│   └── samples_rnaseq.tsv     # Sample metadata template
│
├── references/
│   └── gene_sets/             # Custom pathway GMT files
│       ├── README.md          # Gene set documentation
│       ├── apoptosis.gmt
│       ├── immune.gmt
│       └── viral_mimicry.gmt
│
├── example_data/              # Synthetic demo data
│   ├── counts_example.tsv
│   └── metadata_example.tsv
│
├── docs/
│   ├── ASSUMPTIONS.md         # Scientific assumptions
│   ├── RUNBOOK.md             # Detailed usage guide
│   └── example_outputs/       # Example plots
│
├── data/                      # Your data (gitignored)
├── results/                   # Analysis outputs (gitignored)
└── logs/                      # Pipeline logs (gitignored)
```

---

## Troubleshooting

Check `logs/` directory for detailed error messages.

---

## Acknowledgments

This repository was organized and documented with assistance from [Claude Code](https://claude.ai/code) (Anthropic). The scientific analysis, pipeline development, and validation were performed by Saige Daines. Claude Code assisted with:
- Repository structure and organization
- Documentation formatting and clarity
- Code consolidation and standardization
- Creation of configuration files and examples

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **STAR**: Dobin et al. (2013) Bioinformatics. doi:10.1093/bioinformatics/bts635
- **DESeq2**: Love et al. (2014) Genome Biology. doi:10.1186/s13059-014-0550-8
- **fgsea**: Korotkevich et al. (2021) bioRxiv. doi:10.1101/060012
- **GSVA**: Hänzelmann et al. (2013) BMC Bioinformatics. doi:10.1186/1471-2105-14-7
- **clusterProfiler**: Wu et al. (2021) Innovation. doi:10.1016/j.xinn.2021.100141

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Documentation

- **Quick Start**: This README
- **Detailed Usage**: [docs/RUNBOOK.md](docs/RUNBOOK.md)
- **Scientific Methods**: [docs/ASSUMPTIONS.md](docs/ASSUMPTIONS.md)
- **Gene Set Provenance**: [references/gene_sets/README.md](references/gene_sets/README.md)

---

**Developed for analyzing cisplatin response in cancer cell lines, but adaptable to any 2-condition bulk RNA-seq experiment.**

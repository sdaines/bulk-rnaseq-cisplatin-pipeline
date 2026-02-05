# ASSUMPTIONS.md

Scientific and analytical assumptions for the bulk RNA-seq pipeline.

---

## 1. Biological System and Experimental Design

### 1.1 Study Objective

Characterize the transcriptional response of A2780 ovarian cancer cells to cisplatin treatment using bulk RNA-seq.

### 1.2 Model System

| System | Description |
|--------|-------------|
| A2780 cell line | Cisplatin-sensitive ovarian cancer cell line |

### 1.3 Experimental Conditions

| Parameter | Value | Notes |
|-----------|-------|-------|
| **Treatment** | Cisplatin vs. vehicle control | DMSO vehicle |
| **Biological replication** | 3 replicates per condition | Independent biological samples |
| **Treatment duration** | 24 hours ||
| **Cisplatin concentration** | 10 μM ||

### 1.4 Library Preparation

| Parameter | Assumed Value | Notes |
|-----------|---------------|-------|
| **Library type** | Poly-A selected ||
| **Strandedness** | Unstranded | Configured in `config.yaml` |
| **Read length** | ≥50 bp ||
| **Sequencing depth** | Target ≥20M reads/sample ||
| **Sequencing platform** | Illumina ||

---

## 2. Reference Genome and Annotation

### 2.1 Genome Build

All analyses use **GRCh38/hg38** as the reference genome.

- FASTA: GENCODE GRCh38 primary assembly
- GTF: GENCODE v38 comprehensive gene annotation

### 2.2 Annotation Consistency

- The same GTF annotation is used for alignment (STAR) and quantification (featureCounts).

### 2.3 Organism Databases

R/Bioconductor analyses use:

- `org.Hs.eg.db` for human gene identifier mapping

---

## 3. Alignment and Quantification

### 3.1 Alignment

- **Aligner**: STAR (2-pass mode) against GRCh38
- **Read type**: Paired-end Illumina short reads
- **Strandedness**: Unstranded (configurable in config.yaml)

### 3.2 Quantification

- **Method**: featureCounts at the gene level
- **Feature type**: Exon
- **Attribute**: gene_id
- **Quality filters**:
  - `-B`: Only count read pairs with both ends successfully aligned
  - `-C`: Exclude chimeric fragments (ends mapping to different chromosomes or strands)

### 3.3 Read Filtering

- **Quality threshold**: Q20 (99% base call accuracy)
- **Minimum read length**: 25 bp after trimming

### 3.4 QC Thresholds

Samples below these thresholds are flagged for review:

- **Minimum alignment rate**: 70% uniquely mapped
- **Minimum mapped reads**: 1 million
- **PCR duplication**: Marked for metrics but NOT removed
  - Duplicates are retained in count quantification
  - Rationale: For n=3 studies, keeping duplicates increases statistical power
  - Trade-off: More counts vs potential PCR bias
  - Duplication rates are reported for QC monitoring
  
---

## 4. Differential Expression

### 4.1 Normalization

- **Method**: DESeq2 with default size factor normalization
- **Design formula**: `~ condition`

### 4.2 Filtering

- **Pre-filtering**: Genes with fewer than 5 counts across all samples are excluded prior to normalization

- **Cook's distance filtering**: DISABLED (set to FALSE in results() call)
  - Recommended by DESeq2 manual for n<6
  - Outlier genes are not automatically removed from results

- **Dispersion estimation**: 'local' fit type
  - Uses local regression to fit dispersion-mean relationship

### 4.3 Statistical Thresholds

- **Significance threshold**: Adjusted p-value (Benjamini-Hochberg FDR) < 0.05
  - Hypothesis test: H0: log2FC = 0 vs H1: log2FC ≠ 0

- **Effect size threshold**: |log2 fold change| > 1.0 (2-fold change)
  - Used for **prioritization and filtering**, not significance testing
  - Genes are first tested for significance (FDR < 0.05)
  - Then filtered by effect size for focused analysis

- **Two output files generated**:
  - `*_sig.tsv`: All genes with FDR < 0.05 (used for pathway analysis, comprehensive lists)
  - `*_sig_strong.tsv`: Genes with FDR < 0.05 AND |log2FC| > 1.0 (used for prioritization, validation)

### 4.4 Variance Stabilization

- **Transformation method**:
  - `rlog()` for n < 30 samples (better variance stabilization for small sample sizes)
  - Used for visualization (PCA, heatmaps) only, not for differential expression testing

---

## 5. Reproducibility

### 5.1 Replication

- Minimum 3 biological replicates per condition

### 5.2 Configuration

- All parameters are stored in `config/config.yaml`

### 5.3 Software Environment

- R version: 4.3+
- Key packages: DESeq2, ggplot2, pheatmap, data.table, yaml
- Alignment tools: STAR 2.7+, featureCounts (subread 2.0+)

---

## 6. Outputs

### 6.1 Primary Outputs

| Output | Description |
|--------|-------------|
| `counts_raw.tsv` | Gene-level count matrix |
| `deseq2_results_*.tsv` | Full differential expression results |
| `deseq2_results_*_sig.tsv` | Filtered significant results |
| `vst_matrix.tsv` | Variance-stabilized expression matrix |
| `pca_plot.pdf` | PCA visualization |
| `sample_distance_heatmap.pdf` | Sample clustering heatmap |

---

## 7. Limitations and Considerations

### 7.1 Biological Limitations

- **Cell line model**: In vitro A2780 cell line models may not fully recapitulate:
  - Patient tumor heterogeneity
  - Tumor microenvironment interactions
  - In vivo drug metabolism and pharmacokinetics
- **Acute treatment**: Single timepoint analysis captures immediate response but not:
  - Time-course dynamics
  - Long-term adaptive responses
  - Cell fate decisions (survival vs death)

### 7.2 Technical Limitations

- **Bulk RNA-seq**: Measures average expression across cell population
  - Cannot resolve cell-to-cell heterogeneity
  - Cannot identify rare cell subpopulations
- **Poly-A selection**: Only captures polyadenylated RNAs
  - Misses non-coding RNAs without poly-A tails
  - May not capture all regulatory transcripts

### 7.3 Analytical Considerations

- **Sample size**: With n=3 per group, statistical power is limited
  - Can detect large effect sizes (>2-fold) reliably
  - May miss subtle but biologically relevant changes

### 7.4 Validation

- Results should be validated with orthogonal methods:
  - RT-qPCR for key genes
  - Western blot for protein-level validation
  - Functional assays for pathway activity

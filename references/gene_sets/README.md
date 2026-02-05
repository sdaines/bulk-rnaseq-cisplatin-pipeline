# Gene Set Provenance Documentation

This directory contains curated gene sets in GMT (Gene Matrix Transposed) format for pathway analysis of cisplatin-treated ovarian cancer cells.

---

## Overview

**Total Gene Sets:** 18 pathways across 3 biological themes
- **Immune response:** 6 pathways
- **Viral mimicry response:** 5 pathways
- **Apoptosis:** 7 pathways

**File Format:** GMT (Gene Matrix Transposed)
- Standard format for gene set enrichment analysis (GSEA)
- Each line represents one gene set: `PATHWAY_NAME\tDESCRIPTION\tGENE1\tGENE2\t...`
- Gene identifiers: HGNC gene symbols (human)

---

## File Contents

### 1. `immune.gmt` (6 pathways)

Gene sets related to immune and inflammatory responses expected in cisplatin-treated cells.

| Pathway Name | Genes | Source / Rationale |
|--------------|-------|-------------------|
| `INTERFERON_ALPHA_RESPONSE` | 21 | Type I IFN signaling; key to innate immune response |
| `INTERFERON_GAMMA_RESPONSE` | 20 | Type II IFN signaling; antiviral and antitumor immunity |
| `INFLAMMATORY_RESPONSE` | 20 | General inflammatory cytokine signaling |
| `ANTIGEN_PRESENTATION` | 25 | MHC-I/II presentation; immunogenic cell death marker |
| `T_CELL_ACTIVATION` | 18 | T cell co-stimulation and activation pathways |
| `CYTOTOXICITY` | 15 | NK cell and CTL effector molecules (perforin, granzymes) |

**Biological Context:**
- Cisplatin induces immunogenic cell death (ICD)
- DNA damage triggers cGAS-STING pathway → Type I IFN production
- Expected upregulation of interferon-stimulated genes (ISGs)

---

### 2. `viral_mimicry.gmt` (5 pathways)

Gene sets related to viral mimicry response - a hallmark of DNA-damaging chemotherapy.

| Pathway Name | Genes | Source / Rationale |
|--------------|-------|-------------------|
| `VIRAL_MIMICRY_RESPONSE` | 24 | Integrated viral mimicry signature from literature |
| `DSRNA_SENSING` | 16 | RIG-I, MDA5, PKR pathways; detect cytosolic dsRNA |
| `RETROTRANSPOSON_SILENCING` | 12 | Chromatin silencing of LINE-1, ALU, ERV elements |
| `TYPE_I_INTERFERON_PRODUCTION` | 18 | IRF3/IRF7-mediated IFN-α/β transcription |
| `CYTOSOLIC_DNA_SENSING` | 14 | cGAS-STING-TBK1-IRF3 axis; responds to DNA damage |

**Biological Context:**
- Cisplatin causes DNA damage → chromatin remodeling
- De-repression of endogenous retroelements (ERVs, LINE-1)
- Cytosolic dsRNA/DNA mimics viral infection → innate immune activation
- **Key Literature:**
  - Chiappinelli et al., *Cell* 2015 (5-aza induces viral mimicry)
  - Roulois et al., *Cell* 2015 (epigenetic therapy → dsRNA response)
  - Guler et al., *Immunity* 2017 (cisplatin activates STING)

---

### 3. `apoptosis.gmt` (7 pathways)

Gene sets related to programmed cell death pathways activated by cisplatin.

| Pathway Name | Genes | Source / Rationale |
|--------------|-------|-------------------|
| `EXTRINSIC_APOPTOSIS` | 16 | Death receptor pathway (FAS, TRAIL, TNF-R1) |
| `INTRINSIC_APOPTOSIS` | 22 | Mitochondrial pathway; BAX/BAK pore formation |
| `CASPASE_ACTIVATION` | 18 | Executioner caspases (CASP3, CASP7) and substrates |
| `TP53_DEPENDENT_APOPTOSIS` | 20 | p53 transcriptional targets (PUMA, NOXA, BAX) |
| `APOPTOSIS_VIA_TNF` | 14 | TNF superfamily-mediated cell death |
| `BCL2_FAMILY_APOPTOSIS` | 17 | BCL2, BCL-XL (anti-apoptotic) vs BAX, BAK (pro-apoptotic) |
| `DEATH_RECEPTOR_SIGNALING` | 15 | Upstream death receptor signaling cascades |

**Biological Context:**
- Cisplatin-DNA adducts → ATR/ATM kinase activation
- p53 stabilization → transcription of pro-apoptotic genes
- Expected upregulation of PUMA, NOXA, BAX in A2780 cells

---

## Curation Methodology

### Sources

These gene sets were curated from multiple trusted sources:

1. **MSigDB (Molecular Signatures Database)**
   - Version: v2023.1.Hs (accessed January 2025)
   - Collections used:
     - Hallmark gene sets (H)
     - Gene Ontology Biological Process (C5.BP)
     - Reactome pathways (C2.CP.REACTOME)
   - URL: https://www.gsea-msigdb.org/gsea/msigdb/

2. **Gene Ontology (GO) Consortium**
   - Version: 2024-09-08 release
   - Specific GO terms:
     - GO:0060337 (type I interferon signaling)
     - GO:0034340 (response to type I interferon)
     - GO:0071357 (cellular response to type I interferon)
     - GO:0006915 (apoptotic process)
     - GO:0006955 (immune response)
   - URL: http://geneontology.org/

3. **KEGG PATHWAY Database**
   - Human pathways:
     - hsa04210: Apoptosis
     - hsa04620: Toll-like receptor signaling
     - hsa04622: RIG-I-like receptor signaling
   - URL: https://www.genome.jp/kegg/pathway.html

4. **Literature-Derived Gene Lists**
   - Viral mimicry: Chiappinelli et al. (2015), Roulois et al. (2015)
   - Cisplatin response: Guler et al. (2017), Zhao et al. (2021)
   - Interferon-stimulated genes (ISGs): Interferome database v2.01

---

## Curation Process

### Step 1: Initial Gene List Compilation
- Downloaded relevant gene sets from MSigDB, GO, KEGG
- Extracted gene lists from key publications (Supplementary Tables)
- Filtered to human orthologs (Homo sapiens)

### Step 2: Gene Symbol Standardization
- Converted all gene identifiers to HGNC approved symbols
- Used org.Hs.eg.db (Bioconductor) for ID mapping
- Ensembl → HGNC symbol mapping (GRCh38.p13)
- Removed deprecated symbols, merged duplicates

### Step 3: Gene Set Size Filtering
- Minimum: 10 genes per pathway (avoid noisy small sets)
- Maximum: 30 genes per pathway (ensure specificity)

### Step 4: Redundancy Reduction
- Merged pathways with >70% gene overlap
- Retained most specific/interpretable pathway name

### Step 5: Biological Validation
- Cross-referenced with published cisplatin response data
- Verified expression in A2780 cells (baseline RNA-seq data)
- Confirmed relevance to ovarian cancer and DNA damage response

---

## Quality Metrics

| Metric | Value |
|--------|-------|
| Total unique genes | 247 |
| Genes per pathway (median) | 18 |
| Genes per pathway (range) | 12-25 |
| Overlapping genes between themes | 8.5% |
| Genes expressed in A2780 baseline | >90% |

**Note:** Gene sets are designed to have minimal overlap to maximize pathway specificity while maintaining biological coherence.

---

## Usage in This Pipeline

These gene sets are used in:

1. **Stage 3: ssGSEA** (`03_bulk_rnaseq_ssgsea.R`)
   - Single-sample gene set enrichment analysis
   - Calculates pathway activity scores for each sample
   - Independent of condition labels

2. **Stage 4: Preranked GSEA** (`04_bulk_rnaseq_gsea.R`)
   - Ranks genes by DESeq2 Wald statistic
   - Tests for pathway enrichment in cisplatin vs vehicle
   - Identifies directional pathway changes

3. **Stage 5: GO/KEGG Enrichment** (`05_bulk_rnaseq_go_enrichment.R`)
   - Uses standard GO/KEGG databases (not custom GMT files)
   - Complements custom pathway analysis with broader ontology

---

## Limitations and Considerations

### Biological Limitations
1. **Cell line specificity:** Gene sets optimized for A2780 cells; may not generalize to other cancer types
2. **Treatment specificity:** Focused on DNA damage/cisplatin response; not comprehensive cancer pathways
3. **Temporal dynamics:** Gene sets represent aggregate response, not time-course changes

### Technical Limitations
1. **Gene symbol dependency:** Relies on HGNC symbols; Ensembl IDs or Entrez may be preferred for some analyses
2. **Static annotation:** Based on 2024 genome annotation; retroelements and non-coding RNAs under-represented
3. **Size constraints:** 10-30 gene range excludes some valid small pathways (<10) or broad processes (>30)

### Statistical Considerations
1. **Pathway overlap:** Some genes appear in multiple pathways (e.g., STAT1 in both IFN-α and IFN-γ)
2. **Independence assumption:** ssGSEA and GSEA assume pathway independence (not strictly true)
3. **Multiple testing:** With 18 pathways, FDR correction is applied but power may be limited

---

## File Format Specification

**GMT Format (Gene Matrix Transposed):**
```
PATHWAY_NAME<TAB>DESCRIPTION<TAB>GENE1<TAB>GENE2<TAB>...<TAB>GENEN
```

**Example:**
```
INTERFERON_ALPHA_RESPONSE	Type I IFN signaling pathway	IFIT1	IFIT2	IFIT3	ISG15	MX1	OAS1	...
```
---

## Additional Resources

**Related Databases:**
- **Interferome:** Database of interferon-regulated genes (http://www.interferome.org/)
- **InnateDB:** Curated database of innate immunity genes (https://www.innatedb.com/)
- **MSigDB:** Comprehensive collection of gene sets (https://www.gsea-msigdb.org/)

---

**Last Updated:** 2025-02-04
**Curator:** Saige Daines
**Version:** 1.0
**Compatible with:** GRCh38/hg38 genome build

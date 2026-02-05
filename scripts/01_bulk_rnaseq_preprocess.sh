#!/bin/bash
# ==============================================================================
# 01_bulk_rnaseq_preprocess.sh
# ==============================================================================
# Preprocess bulk RNA-seq FASTQ files to produce a gene-level count matrix.
#
# Steps:
#   1. FastQC on raw FASTQ files
#   2. STAR alignment (2-pass mode)
#   3. featureCounts quantification
#
# Inputs:
#   - config/config.yaml (paths and parameters)
#   - metadata/samples_rnaseq.tsv (sample sheet)
#   - FASTQ files
#
# Outputs (in results/rnaseq/):
#   - fastqc/           FastQC reports
#   - bam/              Aligned BAM files
#   - counts/           Gene count matrix
#   - logs/             Processing logs
#
# Usage: bash scripts/01_bulk_rnaseq_preprocess.sh
# ==============================================================================

set -euo pipefail

echo "=== Bulk RNA-seq Preprocessing ==="
echo "Started: $(date)"
echo ""

# ------------------------------------------------------------------------------
# Helper: parse YAML config (simple grep-based, no dependencies)
# ------------------------------------------------------------------------------

get_config() {
    local key="$1"
    grep -E "^\s*${key}:" config/config.yaml | head -1 | sed 's/.*:\s*"\?\([^"]*\)"\?.*/\1/' | tr -d ' '
}

# ------------------------------------------------------------------------------
# Load configuration
# ------------------------------------------------------------------------------

echo "Loading configuration..."

# Check config exists
if [[ ! -f "config/config.yaml" ]]; then
    echo "ERROR: config/config.yaml not found"
    exit 1
fi

# Reference files
GENOME_FASTA=$(get_config "fasta")
GTF=$(get_config "gtf")
STAR_INDEX=$(get_config "star_index")

# Paths
FASTQ_DIR=$(get_config "fastq_rnaseq")
RESULTS_DIR=$(get_config "results")
# Parse metadata path (under paths.metadata.rnaseq)
METADATA=$(grep -E "^\s*rnaseq:\s*\"" config/config.yaml | head -1 | sed 's/.*:\s*"\([^"]*\)".*/\1/')

# Parameters
STRANDEDNESS=$(grep -E "^\s*strandedness:" config/config.yaml | head -1 | sed 's/.*:\s*"\?\([^"]*\)"\?.*/\1/' | tr -d ' ')
STAR_THREADS=$(grep -E "^\s*star_threads:" config/config.yaml | head -1 | sed 's/.*:\s*\([0-9]*\).*/\1/')

# QC thresholds
MIN_ALIGNMENT_RATE=$(grep -E "^\s*min_alignment_rate:" config/config.yaml | head -1 | sed 's/.*:\s*\([0-9.]*\).*/\1/')
MIN_MAPPED_READS=$(grep -E "^\s*min_mapped_reads:" config/config.yaml | head -1 | sed 's/.*:\s*\([0-9]*\).*/\1/')

# Set defaults if not found
STRANDEDNESS=${STRANDEDNESS:-"unstranded"}
STAR_THREADS=${STAR_THREADS:-4}
MIN_ALIGNMENT_RATE=${MIN_ALIGNMENT_RATE:-0.70}
MIN_MAPPED_READS=${MIN_MAPPED_READS:-1000000}

# Output directories
OUT_DIR="${RESULTS_DIR}/rnaseq"
FASTQC_DIR="${OUT_DIR}/fastqc"
TRIMMED_DIR="${OUT_DIR}/trimmed"
BAM_DIR="${OUT_DIR}/bam"
COUNTS_DIR="${OUT_DIR}/counts"
LOG_DIR="${OUT_DIR}/logs"

echo "  STAR index: ${STAR_INDEX}"
echo "  GTF: ${GTF}"
echo "  FASTQ dir: ${FASTQ_DIR}"
echo "  Strandedness: ${STRANDEDNESS}"
echo "  Threads: ${STAR_THREADS}"
echo ""

# ------------------------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------------------------

echo "Validating inputs..."

if [[ ! -d "${STAR_INDEX}" ]]; then
    echo "ERROR: STAR index not found: ${STAR_INDEX}"
    echo "  Run Stage 0 setup first."
    exit 1
fi

if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF file not found: ${GTF}"
    exit 1
fi

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: Metadata file not found: ${METADATA}"
    exit 1
fi

echo "  All inputs validated."
echo ""

# ------------------------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------------------------

mkdir -p "${FASTQC_DIR}"
mkdir -p "${TRIMMED_DIR}"
mkdir -p "${BAM_DIR}"
mkdir -p "${COUNTS_DIR}"
mkdir -p "${LOG_DIR}"

# ------------------------------------------------------------------------------
# Log tool versions
# ------------------------------------------------------------------------------

VERSION_LOG="${LOG_DIR}/tool_versions.txt"
echo "Recording tool versions..."

{
    echo "Tool Versions"
    echo "============="
    echo "Date: $(date)"
    echo ""
    echo "STAR:"
    STAR --version 2>&1 || echo "STAR not found"
    echo ""
    echo "featureCounts:"
    featureCounts -v 2>&1 || echo "featureCounts not found"
    echo ""
    echo "FastQC:"
    fastqc --version 2>&1 || echo "FastQC not found"
    echo ""
    echo "samtools:"
    samtools --version 2>&1 | head -2 || echo "samtools not found"
    echo ""
    echo "fastp:"
    fastp --version 2>&1 || echo "fastp not found"
    echo ""
    echo "multiqc:"
    multiqc --version 2>&1 || echo "multiqc not found"
} > "${VERSION_LOG}"

echo "  Saved: ${VERSION_LOG}"
echo ""

# ------------------------------------------------------------------------------
# Process each sample
# ------------------------------------------------------------------------------

echo "Processing samples..."
echo ""

# Read sample metadata (skip header)
SAMPLE_COUNT=0
BAM_FILES=""

while IFS=$'\t' read -r sample_id fastq_r1 fastq_r2 condition batch; do
    # Skip header
    if [[ "${sample_id}" == "sample_id" ]]; then
        continue
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo "--- Sample ${SAMPLE_COUNT}: ${sample_id} ---"

    # Construct full FASTQ paths
    R1="${FASTQ_DIR}/${fastq_r1}"
    R2="${FASTQ_DIR}/${fastq_r2}"

    # Check FASTQ files exist
    if [[ ! -f "${R1}" ]]; then
        echo "  WARNING: R1 not found: ${R1}"
        echo "  Skipping sample."
        continue
    fi

    # --------------------------------------------------------------------------
    # Step 1: FastQC on raw reads
    # --------------------------------------------------------------------------

    echo "  Running FastQC on raw reads..."
    FASTQC_LOG="${LOG_DIR}/${sample_id}_fastqc.log"

    if [[ -f "${R2}" ]]; then
        fastqc -o "${FASTQC_DIR}" -t 2 "${R1}" "${R2}" > "${FASTQC_LOG}" 2>&1
    else
        fastqc -o "${FASTQC_DIR}" "${R1}" > "${FASTQC_LOG}" 2>&1
    fi

    # --------------------------------------------------------------------------
    # Step 2: Adapter trimming with fastp
    # --------------------------------------------------------------------------

    echo "  Running fastp trimming..."
    FASTP_LOG="${LOG_DIR}/${sample_id}_fastp.log"
    FASTP_JSON="${TRIMMED_DIR}/${sample_id}_fastp.json"
    FASTP_HTML="${TRIMMED_DIR}/${sample_id}_fastp.html"

    # Output trimmed files
    TRIMMED_R1="${TRIMMED_DIR}/${sample_id}_R1_trimmed.fastq.gz"
    TRIMMED_R2="${TRIMMED_DIR}/${sample_id}_R2_trimmed.fastq.gz"

    if [[ -f "${R2}" ]]; then
        # Paired-end trimming
        fastp \
            -i "${R1}" \
            -I "${R2}" \
            -o "${TRIMMED_R1}" \
            -O "${TRIMMED_R2}" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 25 \
            --json "${FASTP_JSON}" \
            --html "${FASTP_HTML}" \
            > "${FASTP_LOG}" 2>&1
    else
        # Single-end trimming
        fastp \
            -i "${R1}" \
            -o "${TRIMMED_R1}" \
            --qualified_quality_phred 20 \
            --length_required 25 \
            --json "${FASTP_JSON}" \
            --html "${FASTP_HTML}" \
            > "${FASTP_LOG}" 2>&1
        TRIMMED_R2=""
    fi

    # --------------------------------------------------------------------------
    # Step 3: STAR alignment (using trimmed reads)
    # --------------------------------------------------------------------------

    echo "  Running STAR alignment..."
    STAR_OUT_PREFIX="${BAM_DIR}/${sample_id}."
    STAR_LOG="${LOG_DIR}/${sample_id}_star.log"

    # Build STAR command (using trimmed reads)
    STAR_CMD="STAR \
        --runThreadN ${STAR_THREADS} \
        --genomeDir ${STAR_INDEX} \
        --readFilesIn ${TRIMMED_R1}"

    # Add R2 if paired-end
    if [[ -n "${TRIMMED_R2}" ]] && [[ -f "${TRIMMED_R2}" ]]; then
        STAR_CMD="${STAR_CMD} ${TRIMMED_R2}"
    fi

    # Handle gzipped files (trimmed files are gzipped)
    STAR_CMD="${STAR_CMD} --readFilesCommand zcat"

    # Output settings
    STAR_CMD="${STAR_CMD} \
        --outFileNamePrefix ${STAR_OUT_PREFIX} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --twopassMode Basic"

    # Run STAR
    eval "${STAR_CMD}" > "${STAR_LOG}" 2>&1

    # Rename output BAM for consistency
    ALIGNED_BAM="${BAM_DIR}/${sample_id}.Aligned.sortedByCoord.out.bam"
    if [[ -f "${ALIGNED_BAM}" ]]; then
        # --------------------------------------------------------------------------
        # Mark and remove PCR duplicates
        # --------------------------------------------------------------------------
        echo "  Marking PCR duplicates..."
        # DUPLICATE REMOVAL DISABLED TO MATCH RAW SCRIPTS
        # Raw scripts did not remove PCR duplicates, keeping them retains higher counts
        # Trade-off: More statistical power vs potential PCR bias
        # For n=3 studies, the power gain often outweighs the bias risk

        # Mark duplicates for metrics but DO NOT remove them (no -r flag)
        DEDUP_METRICS="${LOG_DIR}/${sample_id}_dup_metrics.txt"
        samtools markdup -s "${ALIGNED_BAM}" /dev/null 2> "${DEDUP_METRICS}"

        # Use original aligned BAM (with duplicates) for quantification
        samtools index "${ALIGNED_BAM}"
        BAM_FILES="${BAM_FILES} ${ALIGNED_BAM}"
        echo "  Using aligned BAM (duplicates retained to match raw scripts)"

        # Still report duplication metrics for QC
        if [[ -f "${DEDUP_METRICS}" ]]; then
            echo "  Duplication metrics: ${DEDUP_METRICS}"
        fi
    else
        echo "  WARNING: BAM not created for ${sample_id}"
    fi

    echo ""

done < "${METADATA}"

echo "Processed ${SAMPLE_COUNT} samples."
echo ""

# ------------------------------------------------------------------------------
# QC: Check alignment metrics against thresholds
# ------------------------------------------------------------------------------

echo "Checking alignment QC metrics..."
QC_REPORT="${LOG_DIR}/qc_summary.tsv"

# Write QC header (added duplication_rate column)
echo -e "sample_id\ttotal_reads\tuniquely_mapped\tmapping_rate\tduplication_rate\tpassed_qc" > "${QC_REPORT}"

FLAGGED_SAMPLES=""

for bam in ${BAM_FILES}; do
    # Handle both .dedup.bam and .Aligned.sortedByCoord.out.bam names
    sample_name=$(basename "${bam}" .dedup.bam)
    sample_name=$(basename "${sample_name}" .Aligned.sortedByCoord.out.bam)
    star_final_log="${BAM_DIR}/${sample_name}.Log.final.out"

    if [[ -f "${star_final_log}" ]]; then
        # Extract metrics from STAR final log
        total_reads=$(grep "Number of input reads" "${star_final_log}" | awk -F'\t' '{print $2}')
        uniquely_mapped=$(grep "Uniquely mapped reads number" "${star_final_log}" | awk -F'\t' '{print $2}')
        mapping_pct=$(grep "Uniquely mapped reads %" "${star_final_log}" | awk -F'\t' '{print $2}' | tr -d '%')

        # Convert percentage to decimal for comparison
        mapping_rate=$(echo "scale=4; ${mapping_pct} / 100" | bc)

        # Extract duplication rate from samtools markdup metrics
        dup_metrics_file="${LOG_DIR}/${sample_name}_dup_metrics.txt"
        duplication_rate="NA"
        if [[ -f "${dup_metrics_file}" ]]; then
            # samtools markdup outputs: READ duplicates
            # Calculate duplication rate from metrics
            duplicates=$(grep -v "^#" "${dup_metrics_file}" | head -1 | awk '{print $3}')
            examined=$(grep -v "^#" "${dup_metrics_file}" | head -1 | awk '{print $2}')
            if [[ -n "${duplicates}" ]] && [[ -n "${examined}" ]] && [[ "${examined}" -gt 0 ]]; then
                duplication_rate=$(echo "scale=4; ${duplicates} / ${examined}" | bc)
            fi
        fi

        # Check against thresholds
        passed="PASS"
        if (( $(echo "${mapping_rate} < ${MIN_ALIGNMENT_RATE}" | bc -l) )); then
            passed="FAIL_ALIGNMENT_RATE"
            FLAGGED_SAMPLES="${FLAGGED_SAMPLES} ${sample_name}"
        fi
        if (( uniquely_mapped < MIN_MAPPED_READS )); then
            passed="FAIL_READ_COUNT"
            FLAGGED_SAMPLES="${FLAGGED_SAMPLES} ${sample_name}"
        fi

        echo -e "${sample_name}\t${total_reads}\t${uniquely_mapped}\t${mapping_rate}\t${duplication_rate}\t${passed}" >> "${QC_REPORT}"
    fi
done

echo "  QC summary: ${QC_REPORT}"

if [[ -n "${FLAGGED_SAMPLES}" ]]; then
    echo "  WARNING: The following samples failed QC thresholds:"
    echo "    ${FLAGGED_SAMPLES}"
    echo "  (min alignment rate: ${MIN_ALIGNMENT_RATE}, min reads: ${MIN_MAPPED_READS})"
else
    echo "  All samples passed QC thresholds."
fi
echo ""

# ------------------------------------------------------------------------------
# featureCounts quantification
# ------------------------------------------------------------------------------

echo "Running featureCounts..."

# Determine strand setting
case "${STRANDEDNESS}" in
    "unstranded")
        STRAND_OPT=0
        ;;
    "forward")
        STRAND_OPT=1
        ;;
    "reverse")
        STRAND_OPT=2
        ;;
    *)
        STRAND_OPT=0
        ;;
esac

COUNTS_RAW="${COUNTS_DIR}/counts_raw.txt"
COUNTS_LOG="${LOG_DIR}/featurecounts.log"

# Build file list (trim leading space)
BAM_LIST=$(echo "${BAM_FILES}" | sed 's/^ //')

if [[ -z "${BAM_LIST}" ]]; then
    echo "ERROR: No BAM files to quantify."
    exit 1
fi

# Run featureCounts
# -B: Only count read pairs with both ends successfully aligned
# -C: Do not count chimeric fragments (ends map to different chromosomes/strands)
# These flags match the raw script settings for quality filtering
featureCounts \
    -T "${STAR_THREADS}" \
    -a "${GTF}" \
    -o "${COUNTS_RAW}" \
    -t exon \
    -g gene_id \
    -s "${STRAND_OPT}" \
    -p \
    -B \
    -C \
    --countReadPairs \
    ${BAM_LIST} \
    > "${COUNTS_LOG}" 2>&1

echo "  Raw counts: ${COUNTS_RAW}"

# ------------------------------------------------------------------------------
# Format count matrix for DESeq2
# ------------------------------------------------------------------------------

echo "Formatting count matrix..."

COUNTS_TSV="${COUNTS_DIR}/counts_raw.tsv"

# featureCounts output format:
# Geneid Chr Start End Strand Length sample1 sample2 ...
# Extract gene ID and count columns, clean sample names

# Get header with cleaned sample names
head -2 "${COUNTS_RAW}" | tail -1 | \
    awk -F'\t' '{
        printf "gene_id";
        for (i=7; i<=NF; i++) {
            # Extract sample name from BAM path
            split($i, a, "/");
            sample = a[length(a)];
            gsub(/\.Aligned\.sortedByCoord\.out\.bam/, "", sample);
            printf "\t%s", sample;
        }
        printf "\n";
    }' > "${COUNTS_TSV}"

# Extract counts (skip first 2 header lines)
tail -n +3 "${COUNTS_RAW}" | \
    awk -F'\t' '{
        printf "%s", $1;
        for (i=7; i<=NF; i++) {
            printf "\t%s", $i;
        }
        printf "\n";
    }' >> "${COUNTS_TSV}"

echo "  Formatted counts: ${COUNTS_TSV}"
echo ""

# ------------------------------------------------------------------------------
# MultiQC: Aggregate QC reports
# ------------------------------------------------------------------------------

echo "Running MultiQC to aggregate QC reports..."

MULTIQC_DIR="${OUT_DIR}/multiqc"
mkdir -p "${MULTIQC_DIR}"

# Run MultiQC on FastQC, fastp, and STAR outputs
multiqc \
    "${FASTQC_DIR}" \
    "${TRIMMED_DIR}" \
    "${BAM_DIR}" \
    "${COUNTS_DIR}" \
    -o "${MULTIQC_DIR}" \
    -n "multiqc_report" \
    --force \
    > "${LOG_DIR}/multiqc.log" 2>&1 || echo "  WARNING: MultiQC encountered issues (check log)"

echo "  MultiQC report: ${MULTIQC_DIR}/multiqc_report.html"
echo ""

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

# Count genes and samples in output
N_GENES=$(tail -n +2 "${COUNTS_TSV}" | wc -l | tr -d ' ')
N_SAMPLES=$(head -1 "${COUNTS_TSV}" | awk -F'\t' '{print NF-1}')

echo "=== Preprocessing Complete ==="
echo "  Samples: ${N_SAMPLES}"
echo "  Genes: ${N_GENES}"
echo "  Count matrix: ${COUNTS_TSV}"
echo "  Trimmed reads: ${TRIMMED_DIR}/"
echo "  BAM files: ${BAM_DIR}/"
echo "  FastQC reports: ${FASTQC_DIR}/"
echo "  MultiQC report: ${MULTIQC_DIR}/"
echo "  Logs: ${LOG_DIR}/"
echo ""
echo "Finished: $(date)"

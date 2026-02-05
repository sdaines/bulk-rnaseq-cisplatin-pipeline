#!/bin/bash
# =============================================================================
# Build STAR Genome Index
# =============================================================================
# Builds the STAR genome index required for bulk RNA-seq alignment.
#
# Usage:
#   ./build_star_index.sh [THREADS]
#
# Arguments:
#   THREADS  - Number of threads for index building (default: 8)
#
# Prerequisites:
#   - STAR aligner installed (v2.7+)
#   - Reference files from download_references.sh
#   - ~32GB RAM for human genome
#   - ~30GB disk space for index
#
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

THREADS="${1:-8}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REF_DIR="${PROJECT_ROOT}/data/references"

# Input files (must match download_references.sh output)
GENOME_FASTA="${REF_DIR}/GRCh38.primary_assembly.fa"
GTF_FILE="${REF_DIR}/gencode.v38.annotation.gtf"

# Output directory (must match config.yaml)
STAR_INDEX_DIR="${REF_DIR}/star_index"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

main() {
    log_info "=============================================="
    log_info "STAR Genome Index Builder"
    log_info "=============================================="
    log_info "Threads: $THREADS"
    log_info "Output directory: $STAR_INDEX_DIR"
    log_info ""

    # Check if index already exists
    if [[ -d "$STAR_INDEX_DIR" ]] && [[ -f "${STAR_INDEX_DIR}/SA" ]]; then
        log_info "STAR index already exists at: $STAR_INDEX_DIR"
        log_info "To rebuild, remove the directory first:"
        log_info "  rm -rf $STAR_INDEX_DIR"
        exit 0
    fi

    # Check prerequisites
    log_info "Checking prerequisites..."

    if ! command -v STAR &> /dev/null; then
        log_error "STAR is required but not installed"
        log_error "Install with: conda install -c bioconda star"
        exit 1
    fi

    STAR_VERSION=$(STAR --version 2>&1 | head -1)
    log_info "STAR version: $STAR_VERSION"

    # Check input files
    if [[ ! -f "$GENOME_FASTA" ]]; then
        log_error "Genome FASTA not found: $GENOME_FASTA"
        log_error "Run download_references.sh first"
        exit 1
    fi

    if [[ ! -f "$GTF_FILE" ]]; then
        log_error "GTF annotation not found: $GTF_FILE"
        log_error "Run download_references.sh first"
        exit 1
    fi

    log_info "Input files:"
    log_info "  FASTA: $GENOME_FASTA"
    log_info "  GTF: $GTF_FILE"
    log_info ""

    # Create output directory
    log_info "Creating output directory..."
    mkdir -p "$STAR_INDEX_DIR"

    # Build STAR index
    log_info "=============================================="
    log_info "Building STAR index..."
    log_info "This may take 30-60 minutes and requires ~32GB RAM"
    log_info "=============================================="
    log_info ""

    STAR \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang 100

    # Verify output
    if [[ ! -f "${STAR_INDEX_DIR}/SA" ]]; then
        log_error "STAR index build failed - SA file not created"
        exit 1
    fi

    log_info ""
    log_info "=============================================="
    log_info "STAR INDEX BUILD COMPLETE"
    log_info "=============================================="
    log_info ""
    log_info "Index location: $STAR_INDEX_DIR"
    log_info ""
    log_info "Next step:"
    log_info "  Place FASTQ files in data/fastq/rnaseq/"
    log_info "  Then run: bash scripts/run_pipeline.sh"
    log_info ""
}

main "$@"

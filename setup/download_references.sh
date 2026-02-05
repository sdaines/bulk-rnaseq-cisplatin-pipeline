#!/bin/bash
# =============================================================================
# Download Reference Files (FASTA and GTF)
# =============================================================================
# Downloads the reference genome and annotation files needed for bulk RNA-seq.
#
# This script downloads:
#   - GRCh38 primary assembly FASTA from GENCODE
#   - GENCODE v38 GTF annotation
#
# Usage:
#   ./download_references.sh
#
# Prerequisites:
#   - wget or curl
#   - gunzip
#   - ~30GB disk space
#
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REF_DIR="${PROJECT_ROOT}/data/references"

# URLs from GENCODE
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"

# Output files
FASTA_GZ="${REF_DIR}/GRCh38.primary_assembly.fa.gz"
FASTA="${REF_DIR}/GRCh38.primary_assembly.fa"
GTF_GZ="${REF_DIR}/gencode.v38.annotation.gtf.gz"
GTF="${REF_DIR}/gencode.v38.annotation.gtf"

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

download_file() {
    local url="$1"
    local output="$2"
    local description="$3"

    if [[ -f "$output" ]]; then
        log_info "SKIP: $description already exists"
        return 0
    fi

    log_info "Downloading: $description"
    log_info "  URL: $url"
    log_info "  Output: $output"

    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "$output" "$url"
    elif command -v curl &> /dev/null; then
        curl -L -o "$output" "$url"
    else
        log_error "Neither wget nor curl available"
        exit 1
    fi

    if [[ ! -f "$output" ]]; then
        log_error "Download failed: $output"
        exit 1
    fi

    log_info "Downloaded: $description"
}

decompress_file() {
    local input="$1"
    local output="$2"
    local description="$3"

    if [[ -f "$output" ]]; then
        log_info "SKIP: $description already decompressed"
        return 0
    fi

    if [[ ! -f "$input" ]]; then
        log_error "Cannot decompress: $input does not exist"
        exit 1
    fi

    log_info "Decompressing: $description"
    gunzip -k "$input"
    log_info "Decompressed: $description"
}

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

main() {
    log_info "=============================================="
    log_info "Reference File Download Script"
    log_info "=============================================="
    log_info "Reference directory: $REF_DIR"
    log_info ""

    # Check prerequisites
    log_info "Checking prerequisites..."
    if ! command -v gunzip &> /dev/null; then
        log_error "gunzip is required but not installed"
        exit 1
    fi
    if ! command -v wget &> /dev/null && ! command -v curl &> /dev/null; then
        log_error "Either wget or curl is required"
        exit 1
    fi
    log_info "Prerequisites OK"
    log_info ""

    # Create reference directory
    log_info "Creating reference directory..."
    mkdir -p "$REF_DIR"

    # -------------------------------------------------------------------------
    # Download FASTA
    # -------------------------------------------------------------------------

    log_info "--- FASTA ---"
    download_file "$FASTA_URL" "$FASTA_GZ" "GRCh38 primary assembly FASTA"
    decompress_file "$FASTA_GZ" "$FASTA" "GRCh38 FASTA"
    log_info ""

    # -------------------------------------------------------------------------
    # Download GTF
    # -------------------------------------------------------------------------

    log_info "--- GTF ---"
    download_file "$GTF_URL" "$GTF_GZ" "GENCODE v38 GTF"
    decompress_file "$GTF_GZ" "$GTF" "GENCODE v38 GTF"
    log_info ""

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    log_info "=============================================="
    log_info "DOWNLOAD COMPLETE"
    log_info "=============================================="
    log_info ""
    log_info "Downloaded files:"
    log_info "  FASTA: $FASTA"
    log_info "  GTF:   $GTF"
    log_info ""
    log_info "Next step:"
    log_info "  Run: bash setup/build_star_index.sh"
    log_info ""
}

main "$@"

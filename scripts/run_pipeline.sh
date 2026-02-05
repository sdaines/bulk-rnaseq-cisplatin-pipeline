#!/usr/bin/env bash
# =============================================================================
# Bulk RNA-seq Pipeline Runner
# =============================================================================
# Main script that executes the bulk RNA-seq pipeline stages in order.
#
# Usage:
#   ./run_pipeline.sh [OPTIONS]
#
# Options:
#   --stage STAGE     Run only the specified stage (preprocess, deseq2, ssgsea, gsea)
#   --dry-run         Print commands without executing
#   --config FILE     Use alternate config file
#   --demo            Run in demo mode with example data (skips alignment)
#   --help            Show this help message
#
# Examples:
#   ./run_pipeline.sh                    # Run all stages
#   ./run_pipeline.sh --demo             # Run demo with example data
#   ./run_pipeline.sh --stage preprocess # Run only preprocessing
#   ./run_pipeline.sh --stage deseq2     # Run only DESeq2 analysis
#   ./run_pipeline.sh --stage ssgsea     # Run only ssGSEA analysis
#   ./run_pipeline.sh --dry-run          # Preview commands
#
# Environment variables:
#   DEMO_MODE=true    Enable demo mode (alternative to --demo flag)
#
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_FILE="${PROJECT_ROOT}/config/config.yaml"
LOG_DIR="${PROJECT_ROOT}/logs"

# Default settings
DRY_RUN=false
SINGLE_STAGE=""
DEMO_MODE="${DEMO_MODE:-false}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1" >&2
}

log_stage() {
    echo -e "${BLUE}[STAGE]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

show_help() {
    cat << EOF
Bulk RNA-seq Pipeline Runner

Usage:
  ./run_pipeline.sh [OPTIONS]

Options:
  --stage STAGE     Run only the specified stage
  --dry-run         Print commands without executing
  --config FILE     Use alternate config file (default: config/config.yaml)
  --demo            Run in demo mode with synthetic example data
  --help            Show this help message

Available Stages:
  preprocess        STAR alignment and featureCounts quantification
  deseq2            DESeq2 differential expression analysis
  ssgsea            ssGSEA pathway analysis
  gsea              Preranked GSEA pathway analysis
  go_enrichment     GO term and KEGG pathway enrichment analysis

Demo Mode:
  Use --demo or set DEMO_MODE=true to run the pipeline with synthetic
  example data. This skips alignment/counting and demonstrates the
  analysis stages (DESeq2, ssGSEA, GSEA, plots).

Examples:
  ./run_pipeline.sh                    # Run all stages (requires FASTQ data)
  ./run_pipeline.sh --demo             # Run demo with example data
  ./run_pipeline.sh --stage deseq2     # Run only DESeq2 analysis
  ./run_pipeline.sh --dry-run          # Preview commands without running
EOF
}

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

check_config() {
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "Config file not found: $CONFIG_FILE"
        exit 1
    fi
    log_info "Using config: $CONFIG_FILE"

    # Check if config has demo_mode: true
    if grep -q "demo_mode:\s*true" "$CONFIG_FILE" 2>/dev/null; then
        DEMO_MODE=true
    fi
}

setup_demo() {
    log_info "=========================================="
    log_info "DEMO MODE: Using synthetic example data"
    log_info "=========================================="

    # Check for example data
    local example_counts="${PROJECT_ROOT}/example_data/counts_example.tsv"
    local example_meta="${PROJECT_ROOT}/example_data/metadata_example.tsv"

    if [[ ! -f "$example_counts" ]]; then
        log_error "Example counts not found: $example_counts"
        log_error "Run: Rscript scripts/generate_example_data.R"
        exit 1
    fi

    if [[ ! -f "$example_meta" ]]; then
        log_error "Example metadata not found: $example_meta"
        exit 1
    fi

    # Create output directory for counts
    local counts_dir="${PROJECT_ROOT}/results/rnaseq/counts"
    mkdir -p "$counts_dir"

    # Copy example counts to expected location (simulating preprocessing output)
    # The example data uses gene symbols directly, so we copy as-is
    cp "$example_counts" "${counts_dir}/counts_raw.tsv"
    log_info "Copied example counts to: ${counts_dir}/counts_raw.tsv"

    log_info "Demo setup complete. Skipping preprocessing stage."
    echo ""
}

setup_directories() {
    log_info "Creating output directories..."
    mkdir -p "${LOG_DIR}"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/fastqc"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/trimmed"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/bam"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/counts"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/de"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/logs"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/multiqc"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/markers"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/gsea"
    mkdir -p "${PROJECT_ROOT}/results/rnaseq/go_enrichment"
    mkdir -p "${PROJECT_ROOT}/results/ssgsea"
    log_info "Directories created."
}

# -----------------------------------------------------------------------------
# Stage Execution
# -----------------------------------------------------------------------------

run_preprocess() {
    local script="${SCRIPT_DIR}/01_bulk_rnaseq_preprocess.sh"
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local log_file="${LOG_DIR}/preprocess_${timestamp}.log"

    log_stage "=========================================="
    log_stage "Stage 1: Preprocessing"
    log_stage "STAR alignment and featureCounts quantification"
    log_stage "=========================================="

    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        exit 1
    fi

    chmod +x "$script"

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would execute: bash $script"
        log_info "[DRY-RUN] Log would be written to: $log_file"
        return 0
    fi

    log_info "Log file: $log_file"

    local start_time=$(date +%s)

    if bash "$script" 2>&1 | tee "$log_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "Preprocessing completed successfully in ${duration}s"
    else
        log_error "Preprocessing failed. Check log: $log_file"
        exit 1
    fi
}

run_deseq2() {
    local script="${SCRIPT_DIR}/02_bulk_rnaseq_deseq2.R"
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local log_file="${LOG_DIR}/deseq2_${timestamp}.log"

    log_stage "=========================================="
    log_stage "Stage 2: Differential Expression"
    log_stage "DESeq2 analysis"
    log_stage "=========================================="

    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        exit 1
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would execute: Rscript $script"
        log_info "[DRY-RUN] Log would be written to: $log_file"
        return 0
    fi

    log_info "Log file: $log_file"

    local start_time=$(date +%s)

    if Rscript "$script" 2>&1 | tee "$log_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "DESeq2 analysis completed successfully in ${duration}s"
    else
        log_error "DESeq2 analysis failed. Check log: $log_file"
        exit 1
    fi
}

run_ssgsea() {
    local script="${SCRIPT_DIR}/03_bulk_rnaseq_ssgsea.R"
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local log_file="${LOG_DIR}/ssgsea_${timestamp}.log"

    log_stage "=========================================="
    log_stage "Stage 3: Pathway Analysis"
    log_stage "ssGSEA analysis"
    log_stage "=========================================="

    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        exit 1
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would execute: Rscript $script"
        log_info "[DRY-RUN] Log would be written to: $log_file"
        return 0
    fi

    log_info "Log file: $log_file"

    local start_time=$(date +%s)

    if Rscript "$script" 2>&1 | tee "$log_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "ssGSEA analysis completed successfully in ${duration}s"
    else
        log_error "ssGSEA analysis failed. Check log: $log_file"
        exit 1
    fi
}

run_gsea() {
    local script="${SCRIPT_DIR}/04_bulk_rnaseq_gsea.R"
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local log_file="${LOG_DIR}/gsea_${timestamp}.log"

    log_stage "=========================================="
    log_stage "Stage 4: Preranked GSEA"
    log_stage "fgsea pathway analysis"
    log_stage "=========================================="

    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        exit 1
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would execute: Rscript $script"
        log_info "[DRY-RUN] Log would be written to: $log_file"
        return 0
    fi

    log_info "Log file: $log_file"

    local start_time=$(date +%s)

    if Rscript "$script" 2>&1 | tee "$log_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "GSEA analysis completed successfully in ${duration}s"
    else
        log_error "GSEA analysis failed. Check log: $log_file"
        exit 1
    fi
}

run_go_enrichment() {
    local script="${SCRIPT_DIR}/05_bulk_rnaseq_go_enrichment.R"
    local timestamp=$(date '+%Y-%m-%d_%H-%M-%S')
    local log_file="${LOG_DIR}/go_enrichment_${timestamp}.log"

    log_stage "=========================================="
    log_stage "Stage 5: GO/KEGG Enrichment"
    log_stage "clusterProfiler enrichment analysis"
    log_stage "=========================================="

    if [[ ! -f "$script" ]]; then
        log_error "Script not found: $script"
        exit 1
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would execute: Rscript $script"
        log_info "[DRY-RUN] Log would be written to: $log_file"
        return 0
    fi

    log_info "Log file: $log_file"

    local start_time=$(date +%s)

    if Rscript "$script" 2>&1 | tee "$log_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "GO/KEGG enrichment analysis completed successfully in ${duration}s"
    else
        log_error "GO/KEGG enrichment analysis failed. Check log: $log_file"
        exit 1
    fi
}

# -----------------------------------------------------------------------------
# Argument Parsing
# -----------------------------------------------------------------------------

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --stage)
                SINGLE_STAGE="$2"
                if [[ "$SINGLE_STAGE" != "preprocess" && "$SINGLE_STAGE" != "deseq2" && "$SINGLE_STAGE" != "ssgsea" && "$SINGLE_STAGE" != "gsea" && "$SINGLE_STAGE" != "go_enrichment" ]]; then
                    log_error "Unknown stage: $SINGLE_STAGE"
                    log_error "Valid stages: preprocess, deseq2, ssgsea, gsea, go_enrichment"
                    exit 1
                fi
                shift 2
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --demo)
                DEMO_MODE=true
                CONFIG_FILE="${PROJECT_ROOT}/config/config_example.yaml"
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

main() {
    echo ""
    echo "=============================================="
    echo "  Bulk RNA-seq Pipeline"
    echo "  $(date)"
    echo "=============================================="
    echo ""

    parse_args "$@"
    check_config

    # Export variables for stage scripts
    export PROJECT_ROOT
    export CONFIG_FILE

    # Setup directories
    setup_directories

    # Setup demo mode if enabled
    if [[ "$DEMO_MODE" == "true" ]]; then
        setup_demo
    fi

    # Run stages
    local total_start=$(date +%s)

    if [[ -z "$SINGLE_STAGE" ]]; then
        # Run all stages (or skip preprocessing in demo mode)
        if [[ "$DEMO_MODE" == "true" ]]; then
            log_info "Running demo pipeline (skipping preprocessing)..."
            echo ""
        else
            log_info "Running full pipeline..."
            echo ""
            run_preprocess
            echo ""
        fi
        run_deseq2
        echo ""
        run_ssgsea
        echo ""
        run_gsea
        echo ""
        run_go_enrichment
    elif [[ "$SINGLE_STAGE" == "preprocess" ]]; then
        if [[ "$DEMO_MODE" == "true" ]]; then
            log_warn "Preprocessing skipped in demo mode (no FASTQ data)"
        else
            run_preprocess
        fi
    elif [[ "$SINGLE_STAGE" == "deseq2" ]]; then
        run_deseq2
    elif [[ "$SINGLE_STAGE" == "ssgsea" ]]; then
        run_ssgsea
    elif [[ "$SINGLE_STAGE" == "gsea" ]]; then
        run_gsea
    elif [[ "$SINGLE_STAGE" == "go_enrichment" ]]; then
        run_go_enrichment
    fi

    local total_end=$(date +%s)
    local total_duration=$((total_end - total_start))

    echo ""
    log_info "=========================================="
    log_info "Pipeline completed!"
    log_info "Total duration: ${total_duration}s"
    log_info "=========================================="
}

main "$@"

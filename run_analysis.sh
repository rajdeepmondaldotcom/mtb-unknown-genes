#!/bin/bash

# Mycobacterium Tuberculosis Homology Analysis Pipeline
# 
# This script executes the complete MTB homology analysis pipeline including
# data download, database setup, BLAST search, and homology analysis.
#
# Author: MTB Research Team
# License: MIT

set -euo pipefail  # Exit on error, undefined vars, pipe failures

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly DATA_DIR="$SCRIPT_DIR/data"
readonly MTB_PROTEOME_URL="https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000001584%29"

# Configuration
readonly MTB_PROTEOME_FILE="$DATA_DIR/mtb_h37rv_proteome.fasta"
readonly BLAST_RESULTS_FILE="$DATA_DIR/mtb_vs_swissprot.tsv"
readonly BLAST_EVALUE="1e-5"
readonly BLAST_THREADS="8"

# Colors for output formatting
readonly COLOR_RED='\033[0;31m'
readonly COLOR_GREEN='\033[0;32m'
readonly COLOR_YELLOW='\033[1;33m'
readonly COLOR_BLUE='\033[0;34m'
readonly COLOR_NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${COLOR_BLUE}[INFO]${COLOR_NC} $1"
}

log_success() {
    echo -e "${COLOR_GREEN}[SUCCESS]${COLOR_NC} $1"
}

log_warning() {
    echo -e "${COLOR_YELLOW}[WARNING]${COLOR_NC} $1"
}

log_error() {
    echo -e "${COLOR_RED}[ERROR]${COLOR_NC} $1" >&2
}

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Verify required tools are available
verify_dependencies() {
    log_info "Verifying required dependencies"
    
    local required_tools=("curl" "blastp" "makeblastdb" "update_blastdb.pl" "python3")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command_exists "$tool"; then
            missing_tools+=("$tool")
        fi
    done
    
    if [[ ${#missing_tools[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_error "Please run './setup_environment.sh' first"
        exit 1
    fi
    
    log_success "All required dependencies are available"
}

# Create and setup data directory
setup_data_directory() {
    log_info "Setting up data directory"
    
    if [[ ! -d "$DATA_DIR" ]]; then
        if mkdir -p "$DATA_DIR"; then
            log_success "Created data directory: $DATA_DIR"
        else
            log_error "Failed to create data directory: $DATA_DIR"
            exit 1
        fi
    else
        log_info "Data directory already exists: $DATA_DIR"
    fi
    
    cd "$DATA_DIR" || {
        log_error "Failed to change to data directory"
        exit 1
    }
}

# Download MTB H37Rv proteome from UniProt
download_mtb_proteome() {
    log_info "Downloading MTB H37Rv proteome"
    
    if [[ -f "mtb_h37rv_proteome.fasta" ]]; then
        log_info "MTB proteome already exists"
        return 0
    fi
    
    log_info "Downloading from UniProt..."
    if curl -f -o "mtb_h37rv_proteome.fasta" "$MTB_PROTEOME_URL"; then
        log_success "MTB proteome downloaded successfully"
        
        # Verify download
        if [[ -s "mtb_h37rv_proteome.fasta" ]]; then
            local sequence_count
            sequence_count=$(grep -c "^>" "mtb_h37rv_proteome.fasta" || echo "0")
            log_info "Downloaded $sequence_count protein sequences"
        else
            log_error "Downloaded file is empty"
            exit 1
        fi
    else
        log_error "Failed to download MTB proteome"
        log_error "Please check your internet connection and try again"
        exit 1
    fi
}

# Setup Swiss-Prot BLAST database
setup_swissprot_database() {
    log_info "Setting up Swiss-Prot BLAST database"
    
    # Check if database already exists
    if [[ -f "swissprot.pin" ]]; then
        log_info "Swiss-Prot database already exists"
        return 0
    fi
    
    log_info "Downloading Swiss-Prot database..."
    if update_blastdb.pl --decompress swissprot; then
        log_success "Swiss-Prot database setup complete"
        
        # Verify database
        if [[ -f "swissprot.pin" ]]; then
            local db_info
            db_info=$(blastdbcmd -db swissprot -info | head -3 || echo "Database info unavailable")
            log_info "Database information:"
            echo "$db_info"
        else
            log_error "Database setup failed - missing database files"
            exit 1
        fi
    else
        log_error "Failed to setup Swiss-Prot database"
        log_error "Please check your internet connection and try again"
        exit 1
    fi
}

# Run BLAST search
run_blast_search() {
    log_info "Running BLAST search"
    
    if [[ -f "mtb_vs_swissprot.tsv" ]]; then
        log_info "BLAST results already exist"
        return 0
    fi
    
    log_info "Executing BLASTP search (this may take some time)..."
    log_info "Parameters: E-value=$BLAST_EVALUE, Threads=$BLAST_THREADS"
    
    local blast_start_time
    blast_start_time=$(date +%s)
    
    if blastp \
        -query "mtb_h37rv_proteome.fasta" \
        -db "swissprot" \
        -out "mtb_vs_swissprot.tsv" \
        -evalue "$BLAST_EVALUE" \
        -num_threads "$BLAST_THREADS" \
        -outfmt "7 qseqid sseqid pident ppos length evalue bitscore"; then
        
        local blast_end_time
        blast_end_time=$(date +%s)
        local blast_duration=$((blast_end_time - blast_start_time))
        
        log_success "BLAST search completed successfully"
        log_info "BLAST search duration: ${blast_duration} seconds"
        
        # Verify results
        if [[ -s "mtb_vs_swissprot.tsv" ]]; then
            local hit_count
            hit_count=$(grep -c "^[^#]" "mtb_vs_swissprot.tsv" || echo "0")
            log_info "Generated $hit_count BLAST hits"
        else
            log_error "BLAST results file is empty"
            exit 1
        fi
    else
        log_error "BLAST search failed"
        exit 1
    fi
}

# Execute homology analysis
run_homology_analysis() {
    log_info "Running homology analysis"
    
    cd "$SCRIPT_DIR" || {
        log_error "Failed to change to script directory"
        exit 1
    }
    
    log_info "Executing MTB homology analyzer..."
    log_info "Features:"
    log_info "  - Raw file parsing for maximum performance"
    log_info "  - Smart caching system"
    log_info "  - Vectorized operations"
    log_info "  - Comprehensive error handling"
    echo
    
    if python3 mtb_homology_analyzer.py; then
        log_success "Homology analysis completed successfully"
    else
        log_error "Homology analysis failed"
        exit 1
    fi
}

# Display final results summary
display_results_summary() {
    log_info "Analysis Results Summary"
    echo "========================================"
    
    local results_5_percent="$DATA_DIR/mtb_homologs_top_5_percent.csv"
    local results_10_percent="$DATA_DIR/mtb_homologs_top_10_percent.csv"
    
    if [[ -f "$results_5_percent" && -f "$results_10_percent" ]]; then
        log_success "Results available:"
        
        # Count records in each file (excluding header)
        local count_5_percent
        local count_10_percent
        count_5_percent=$(tail -n +2 "$results_5_percent" | wc -l || echo "0")
        count_10_percent=$(tail -n +2 "$results_10_percent" | wc -l || echo "0")
        
        log_info "  - mtb_homologs_top_5_percent.csv ($count_5_percent matches)"
        log_info "  - mtb_homologs_top_10_percent.csv ($count_10_percent matches)"
        
        echo
        log_info "Complete annotations included:"
        log_info "  - MTB Rv numbers and gene names"
        log_info "  - Ortholog protein names and organisms"
        log_info "  - Comprehensive similarity metrics"
        
    else
        log_error "Expected result files not found"
        exit 1
    fi
}

# Display usage information
show_usage() {
    cat << EOF
Mycobacterium Tuberculosis Homology Analysis Pipeline

Usage: $0 [OPTIONS]

This script executes the complete MTB homology analysis pipeline:
  1. Downloads MTB H37Rv proteome from UniProt
  2. Sets up Swiss-Prot BLAST database
  3. Performs BLAST search
  4. Runs homology analysis

Options:
  -h, --help     Show this help message
  -v, --verbose  Enable verbose output

Requirements:
  - NCBI BLAST+ tools
  - Python 3.9+ with pandas and numpy
  - Internet connection for data download
  - Minimum 8GB RAM and 50GB disk space

EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_usage
                exit 0
                ;;
            -v|--verbose)
                set -x  # Enable verbose output
                shift
                ;;
            *)
                log_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
}

# Main execution function
main() {
    local start_time
    start_time=$(date +%s)
    
    echo "Mycobacterium Tuberculosis Homology Analysis Pipeline"
    echo "====================================================="
    echo
    
    # Parse command line arguments
    parse_arguments "$@"
    
    # Verify dependencies
    verify_dependencies
    
    # Setup data directory
    setup_data_directory
    
    # Execute pipeline steps
    download_mtb_proteome
    setup_swissprot_database
    run_blast_search
    
    # Return to script directory
    cd "$SCRIPT_DIR" || exit 1
    
    # Run analysis
    run_homology_analysis
    
    # Display results
    display_results_summary
    
    local end_time
    end_time=$(date +%s)
    local total_duration=$((end_time - start_time))
    
    echo
    log_success "ANALYSIS PIPELINE COMPLETE"
    log_info "Total execution time: ${total_duration} seconds"
    echo
    log_info "Next steps:"
    log_info "  - Review result files in $DATA_DIR"
    log_info "  - Analyze homology patterns in your favorite tool"
    log_info "  - Validate biological findings experimentally"
}

# Trap to handle script interruption
cleanup() {
    log_warning "Analysis interrupted by user"
    exit 130
}
trap cleanup INT

# Execute main function with all arguments
main "$@"
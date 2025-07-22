#!/bin/bash

# Mycobacterium Tuberculosis Analysis Environment Setup
# 
# This script sets up the development environment for MTB homology analysis
# including Python dependencies and NCBI BLAST+ tools.
#

set -euo pipefail  # Exit on error, undefined vars, pipe failures

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_NAME="mtb-analysis"
readonly PYTHON_VERSION="3.9"

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

# Detect operating system
detect_os() {
    case "${OSTYPE}" in
        darwin*) echo "macOS" ;;
        linux*) echo "Linux" ;;
        *) echo "Unknown" ;;
    esac
}

# Install BLAST+ on macOS using Homebrew
install_blast_macos() {
    log_info "Installing BLAST+ on macOS using Homebrew"
    
    if ! command_exists brew; then
        log_error "Homebrew is not installed"
        log_info "Please install Homebrew first:"
        log_info "/bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
        return 1
    fi
    
    if brew install blast; then
        log_success "BLAST+ installed successfully via Homebrew"
        return 0
    else
        log_error "Failed to install BLAST+ via Homebrew"
        return 1
    fi
}

# Install BLAST+ on Linux
install_blast_linux() {
    log_info "Installing BLAST+ on Linux"
    
    if command_exists apt-get; then
        log_info "Using apt-get package manager"
        sudo apt-get update && sudo apt-get install -y ncbi-blast+
    elif command_exists yum; then
        log_info "Using yum package manager"
        sudo yum install -y ncbi-blast+
    elif command_exists dnf; then
        log_info "Using dnf package manager"
        sudo dnf install -y ncbi-blast+
    else
        log_error "No supported package manager found"
        log_info "Please install BLAST+ manually from:"
        log_info "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
        return 1
    fi
}

# Install BLAST+ via conda
install_blast_conda() {
    log_info "Installing BLAST+ via conda"
    
    # Try different conda channels
    local channels=("bioconda" "conda-forge bioconda")
    
    for channel in "${channels[@]}"; do
        log_info "Trying channel: $channel"
        if conda install -c $channel blast -y; then
            log_success "BLAST+ installed successfully via conda"
            return 0
        fi
    done
    
    log_warning "Could not install BLAST+ via conda due to dependency conflicts"
    return 1
}

# Setup conda environment
setup_conda_environment() {
    log_info "Setting up conda environment: $PROJECT_NAME"
    
    # Create conda environment
    if conda create -n "$PROJECT_NAME" python="$PYTHON_VERSION" -y; then
        log_success "Conda environment created successfully"
    else
        log_error "Failed to create conda environment"
        return 1
    fi
    
    # Activate environment
    eval "$(conda shell.bash hook)"
    if conda activate "$PROJECT_NAME"; then
        log_success "Conda environment activated"
    else
        log_error "Failed to activate conda environment"
        return 1
    fi
    
    # Install BLAST+ in conda environment
    if ! install_blast_conda; then
        log_warning "Conda BLAST+ installation failed, trying system installation"
        install_blast_system
    fi
    
    # Install Python dependencies
    install_python_dependencies
}

# Install BLAST+ using system package manager
install_blast_system() {
    local os_type
    os_type=$(detect_os)
    
    case "$os_type" in
        "macOS")
            install_blast_macos
            ;;
        "Linux")
            install_blast_linux
            ;;
        *)
            log_error "Unsupported operating system: $os_type"
            log_info "Please install BLAST+ manually"
            return 1
            ;;
    esac
}

# Install Python dependencies
install_python_dependencies() {
    log_info "Installing Python dependencies"
    
    if [[ -f "requirements.txt" ]]; then
        if pip install -r requirements.txt; then
            log_success "Python dependencies installed successfully"
        else
            log_error "Failed to install Python dependencies"
            return 1
        fi
    else
        log_warning "requirements.txt not found, installing basic dependencies"
        if pip install pandas numpy; then
            log_success "Basic Python dependencies installed"
        else
            log_error "Failed to install basic Python dependencies"
            return 1
        fi
    fi
}

# Verify installations
verify_installations() {
    log_info "Verifying installations"
    
    local tools=("blastp" "makeblastdb" "update_blastdb.pl" "blastdbcmd")
    local missing_tools=()
    
    # Check BLAST+ tools
    for tool in "${tools[@]}"; do
        if command_exists "$tool"; then
            log_success "$tool is available"
        else
            log_warning "$tool is not available"
            missing_tools+=("$tool")
        fi
    done
    
    # Check Python packages
    if python3 -c "import pandas, numpy" 2>/dev/null; then
        log_success "Python dependencies (pandas, numpy) are available"
    else
        log_warning "Some Python dependencies may not be available"
    fi
    
    # Report BLAST+ version if available
    if command_exists blastp; then
        local blast_version
        blast_version=$(blastp -version 2>/dev/null | head -1 || echo "Version unknown")
        log_info "BLAST+ version: $blast_version"
    fi
    
    # Final status
    if [[ ${#missing_tools[@]} -eq 0 ]]; then
        log_success "All BLAST+ tools are available"
        return 0
    else
        log_warning "Missing BLAST+ tools: ${missing_tools[*]}"
        return 1
    fi
}

# Main execution function
main() {
    echo "Mycobacterium Tuberculosis Analysis - Environment Setup"
    echo "======================================================"
    echo
    
    log_info "Detected OS: $(detect_os)"
    
    if command_exists conda; then
        log_info "Conda found - setting up conda environment"
        setup_conda_environment
    else
        log_info "Conda not found - using system Python"
        install_blast_system
        install_python_dependencies
    fi
    
    echo
    verify_installations
    
    echo
    if verify_installations >/dev/null 2>&1; then
        log_success "Environment setup complete!"
        log_info "You can now run the analysis with: ./run_analysis.sh"
    else
        log_warning "Environment setup completed with some issues"
        log_info "Please check the warnings above and install missing components"
    fi
    
    if command_exists conda; then
        echo
        log_info "To use the conda environment in future sessions:"
        log_info "conda activate $PROJECT_NAME"
    fi
}

# Execute main function
main "$@" 
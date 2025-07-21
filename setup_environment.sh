#!/bin/bash

# Exit on error
set -e

echo "Setting up environment for Mycobacterium tuberculosis analysis..."

# Check if conda is installed
if command -v conda &> /dev/null; then
    echo "Setting up conda environment..."
    
    # Create a new conda environment
    conda create -n mtb-analysis python=3.9 -y
    
    # Activate the environment
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate mtb-analysis
    
    # Install BLAST+ via conda
    echo "Installing NCBI BLAST+..."
    conda install -c bioconda blast -y
    
    # Install Python dependencies
    echo "Installing Python dependencies..."
    pip install pandas biopython
    
    echo "Environment setup complete!"
    echo "To use this environment, run: conda activate mtb-analysis"
else
    echo "Conda not found. Installing dependencies using pip..."
    
    echo "Note: You'll need to install NCBI BLAST+ separately."
    echo "On macOS, you can use: brew install blast"
    echo "On Ubuntu/Debian, you can use: sudo apt-get install ncbi-blast+"
    
    # Install Python dependencies
    pip install pandas biopython
    
    echo "Python dependencies installed. Please ensure BLAST+ is installed before running the analysis."
fi

echo "Setup complete! Run './run_analysis.sh' to start the analysis." 
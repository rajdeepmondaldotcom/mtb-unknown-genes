#!/bin/bash

# Exit on error
set -e

echo "Starting Mycobacterium tuberculosis proteome analysis..."

# Create a directory for data
mkdir -p data
cd data

# Step 1: Download the M. tuberculosis H37Rv Proteome
echo "Downloading M. tuberculosis H37Rv proteome from UniProt..."
curl -o mtb_h37rv_proteome.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000001584%29"
echo "Proteome download complete."

# Step 2: Download the Swiss-Prot Database
# Check if the update_blastdb.pl script is available
if command -v update_blastdb.pl &> /dev/null; then
    echo "Downloading Swiss-Prot database using update_blastdb.pl..."
    update_blastdb.pl --decompress swissprot
else
    echo "Warning: update_blastdb.pl not found. Please install NCBI BLAST+ tools."
    echo "You can install it using: conda install -c bioconda blast"
    exit 1
fi

# Step 3: The Swiss-Prot database is already formatted from update_blastdb.pl
echo "Swiss-Prot database is ready for BLAST search..."

# Step 4: Execute the BLAST Search
echo "Running BLAST search (this may take some time)..."
blastp -query mtb_h37rv_proteome.fasta \
       -db swissprot \
       -out mtb_vs_swissprot.tsv \
       -evalue 1e-5 \
       -num_threads 8 \
       -outfmt "7 qseqid sseqid pident ppos length evalue bitscore"

cd ..

# Step 5: Run the Python script to filter and analyze results
echo "Processing BLAST results with Python..."
python3 filter_blast_results.py

echo "Analysis complete! Results are available in data/top_5_percent_homologs.csv" 
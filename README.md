# Mycobacterium Tuberculosis Homology Analyzer

A high-performance bioinformatics pipeline for analyzing protein sequence homology between Mycobacterium tuberculosis H37Rv proteome and the Swiss-Prot database. This tool identifies and annotates homologous proteins to provide insights into MTB protein function and evolutionary relationships.

## Overview

Mycobacterium tuberculosis remains one of the world's most significant pathogens. Understanding protein homology relationships enables researchers to:

- Identify functional orthologs across different organisms
- Predict protein functions based on well-annotated homologs
- Study evolutionary relationships and conservation patterns
- Guide drug discovery through identification of similar proteins in other organisms

This analyzer performs large-scale BLAST searches against Swiss-Prot and generates comprehensive, annotated results with optimized performance and guaranteed data completeness.

## Features

### Comprehensive Analysis
- Dual threshold analysis: Top 5% and 10% homology matches
- Complete annotations including:
  - MTB Rv numbers and gene names
  - Swiss-Prot protein names and source organisms
  - Detailed similarity metrics (identity, positives, E-value, bit score)
- Zero missing data guarantee through robust error handling

## System Requirements

### Hardware Requirements
- Python 3.9 or higher
- NCBI BLAST+ tools (blastp, makeblastdb, update_blastdb.pl)

### Software Dependencies
- pandas >= 2.0.0
- numpy >= 1.24.0

## Installation

### Quick Setup
```bash
git clone https://github.com/rajdeepmondaldotcom/mtb-unknown-genes.git
cd mtb-unknown-genes
./setup_environment.sh
```

### Manual Installation
```bash
# Create and activate conda environment
conda create -n mtb-analysis python=3.9 -y
conda activate mtb-analysis

# Install BLAST+ tools
conda install -c bioconda blast -y

# Install Python dependencies
pip install pandas numpy
```

## Usage

### Automated Pipeline
Execute the complete analysis pipeline:
```bash
./run_analysis.sh
```

### Manual Execution
For custom analysis or troubleshooting:

```bash
# Create data directory
mkdir -p data
cd data

# Download MTB proteome
curl -o mtb_h37rv_proteome.fasta \
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000001584%29"

# Setup Swiss-Prot database
update_blastdb.pl --decompress swissprot

# Run BLAST search
blastp -query mtb_h37rv_proteome.fasta \
       -db swissprot \
       -out mtb_vs_swissprot.tsv \
       -evalue 1e-5 \
       -num_threads 8 \
       -outfmt "7 qseqid sseqid pident ppos length evalue bitscore"

# Execute homology analysis
cd ..
python3 mtb_homology_analyzer.py
```

## Output Files

### Primary Results
- `data/mtb_homologs_top_5_percent.csv`: Top 5% homology matches
- `data/mtb_homologs_top_10_percent.csv`: Top 10% homology matches

### Result Schema
Each CSV contains the following columns:

| Column | Description |
|--------|-------------|
| mtb_uniprot_id | MTB protein UniProt identifier |
| mtb_rv_number | MTB gene Rv number |
| mtb_gene_name | MTB gene name |
| mtb_protein_description | MTB protein functional description |
| ortholog_swissprot_id | Swiss-Prot homolog identifier |
| ortholog_protein_name | Swiss-Prot protein name |
| ortholog_organism | Source organism of homolog |
| ortholog_full_description | Complete Swiss-Prot description |
| sequence_identity_percent | Percentage sequence identity |
| positive_substitutions_percent | Percentage positive substitutions |
| alignment_length | Length of sequence alignment |
| statistical_evalue | Statistical significance (E-value) |
| alignment_bitscore | BLAST alignment bit score |

### Cache Files
- `data/cache/mtb_genes.pkl`: Cached MTB gene information
- `data/cache/swissprot_annotations.pkl`: Cached Swiss-Prot annotations


## Scientific Context

### Background
Mycobacterium tuberculosis causes tuberculosis, affecting millions worldwide. Protein homology analysis enables:

- Functional annotation of hypothetical proteins
- Drug target identification through ortholog analysis
- Evolutionary studies of pathogenic mechanisms
- Comparative genomics across mycobacterial species

### Methodology
- BLASTP algorithm for protein sequence comparison
- E-value threshold: 1e-5 for statistical significance
- Swiss-Prot database: High-quality, manually annotated proteins
- Top-percentile filtering: Focus on most significant matches

## Acknowledgments

- NCBI BLAST+ for sequence alignment tools
- UniProt for high-quality protein databases
- Swiss-Prot for manually curated protein annotations
- MTB research community for proteome data and insights
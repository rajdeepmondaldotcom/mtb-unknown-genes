# Mycobacterium Tuberculosis Proteome Homology Analysis

This project performs a comprehensive homology analysis to identify the top 5% of homologous proteins for the complete *Mycobacterium tuberculosis* H37Rv proteome. The analysis uses BLAST (Basic Local Alignment Search Tool) to search against the Swiss-Prot database and applies statistical filtering to identify the most significant protein homologs.

## Project Overview

Tuberculosis remains one of the world's leading infectious disease killers. Understanding protein homology relationships in *M. tuberculosis* is crucial for:

- Functional annotation of hypothetical proteins
- Drug target identification and validation
- Understanding evolutionary relationships
- Comparative genomics studies
- Vaccine development research

This automated pipeline processes all ~4,000 proteins in the MTB proteome and identifies their most significant homologs based on statistical significance, sequence identity, and positive substitutions.

## Methodology

### 1. Data Sources
- **Query Dataset**: Complete *M. tuberculosis* H37Rv proteome (UP000001584) from UniProt
- **Target Database**: Swiss-Prot manually reviewed protein database from NCBI
- **Reference Strain**: H37Rv (the most widely studied laboratory strain)

### 2. Analysis Pipeline
1. **Environment Setup**: Install NCBI BLAST+ suite and Python dependencies
2. **Data Acquisition**: Download MTB proteome and Swiss-Prot database
3. **Homology Search**: Execute proteome-wide BLAST search using blastp
4. **Statistical Filtering**: Apply E-value threshold (1e-5) to remove spurious matches
5. **Result Ranking**: Sort hits using hierarchical criteria:
   - E-value (ascending - statistical significance)
   - Bit score (descending - alignment quality)
   - Percent identity (descending - sequence similarity)
   - Percent positives (descending - conservative substitutions)
6. **Top Selection**: Extract top 5% of results for each MTB protein

### 3. Quality Control
- Minimum E-value threshold of 1e-5 for statistical significance
- Multi-threaded processing for computational efficiency
- Comprehensive error handling and validation
- Hierarchical sorting ensures biologically meaningful results

## Installation and Setup

### Prerequisites
- Unix-like operating system (macOS, Linux)
- Internet connection for data download
- At least 4GB free disk space
- 8GB+ RAM recommended

### Quick Start
```bash
# Clone the repository
git clone https://github.com/rajdeepmondaldotcom/mtb-unknown-genes.git
cd mtb-unknown-genes

# Set up environment and dependencies
./setup_environment.sh

# Run the complete analysis
./run_analysis.sh
```

### Manual Installation
If you prefer manual setup:

```bash
# Install BLAST+ (choose one method)
# Via conda:
conda install -c bioconda blast

# Via homebrew (macOS):
brew install blast

# Via apt (Ubuntu/Debian):
sudo apt-get install ncbi-blast+

# Install Python dependencies
pip install pandas biopython
```

## Usage

### Automated Analysis
The complete analysis can be run with a single command:
```bash
./run_analysis.sh
```

This script will:
1. Create a `data/` directory
2. Download the MTB proteome (~1.8MB)
3. Download and setup Swiss-Prot database (~400MB compressed)
4. Execute BLAST search (may take 30-60 minutes depending on system)
5. Process and filter results
6. Generate final output file

### Manual Execution
For step-by-step execution or troubleshooting:

```bash
# Create data directory
mkdir -p data && cd data

# Download MTB proteome
curl -o mtb_h37rv_proteome.fasta "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000001584%29"

# Download Swiss-Prot database
update_blastdb.pl --decompress swissprot

# Run BLAST search
blastp -query mtb_h37rv_proteome.fasta \
       -db swissprot \
       -out mtb_vs_swissprot.tsv \
       -evalue 1e-5 \
       -num_threads 8 \
       -outfmt "7 qseqid sseqid pident ppos length evalue bitscore"

# Process results
cd .. && python3 filter_blast_results.py
```

## Output Files

### Primary Output
- **`data/top_5_percent_homologs.csv`**: Final filtered results containing the top 5% homologs for each MTB protein

### Intermediate Files
- **`data/mtb_h37rv_proteome.fasta`**: MTB proteome sequences
- **`data/mtb_vs_swissprot.tsv`**: Raw BLAST results
- **`data/swissprot.*`**: Swiss-Prot database files

### Output Format
The main results file contains the following columns:
- **qseqid**: Query sequence ID (MTB protein)
- **sseqid**: Subject sequence ID (homologous protein)
- **pident**: Percentage of identical matches
- **ppos**: Percentage of positive-scoring matches
- **length**: Alignment length
- **evalue**: Expect value (statistical significance)
- **bitscore**: Bit score (alignment quality)

## Results Summary

Based on the latest analysis:
- **Total MTB proteins analyzed**: 3,267 (with significant matches)
- **Total homologous proteins identified**: 20,021
- **Average homologs per MTB protein**: 6.1
- **Database searched**: Swiss-Prot (manually curated, high-quality)
- **Statistical threshold**: E-value ≤ 1e-5

## File Structure
```
mtb-unknown-genes/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── setup_environment.sh         # Environment setup script
├── run_analysis.sh             # Main analysis pipeline
├── filter_blast_results.py     # Results processing script
├── .gitignore                  # Git ignore patterns
└── data/                       # Analysis data and results
    ├── mtb_h37rv_proteome.fasta
    ├── mtb_vs_swissprot.tsv
    ├── top_5_percent_homologs.csv
    └── swissprot.*             # BLAST database files
```

## Technical Specifications

### System Requirements
- **CPU**: Multi-core processor recommended (8+ cores optimal)
- **Memory**: 8GB RAM minimum, 16GB recommended
- **Storage**: 4GB free space for databases and results
- **Network**: Stable internet connection for downloads

### Performance Notes
- BLAST search time: 30-60 minutes (system dependent)
- Database download: 5-10 minutes (network dependent)
- Results processing: <5 minutes
- Total runtime: 45-75 minutes

### Computational Parameters
- **E-value threshold**: 1e-5 (stringent statistical cutoff)
- **Output format**: Tabular with custom fields
- **Threading**: 8 cores (adjustable in script)
- **Database**: Swiss-Prot (manually reviewed entries only)

## Troubleshooting

### Common Issues
1. **BLAST+ installation fails**: Try alternative installation method (homebrew, apt, manual)
2. **Database download timeout**: Check internet connection, retry download
3. **Memory issues**: Reduce thread count in blastp command
4. **Permission errors**: Ensure scripts are executable (`chmod +x *.sh`)

### Error Messages
- **"File swissprot does not exist"**: Database download incomplete, re-run setup
- **"Command not found: blastp"**: BLAST+ not installed or not in PATH
- **"No module named pandas"**: Python dependencies not installed

## Citation and References

If you use this pipeline in your research, please cite:
- NCBI BLAST+ suite
- UniProt database
- Swiss-Prot database
- This analysis pipeline

## Contributing

Contributions are welcome. Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is provided as-is for research and educational purposes. Please respect the terms of use for:
- NCBI BLAST+ software
- UniProt database
- Swiss-Prot database

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review the error messages carefully
3. Ensure all dependencies are properly installed
4. Verify internet connectivity for downloads

## Version History

- **v1.0.0**: Initial release with complete analysis pipeline
- Automated environment setup
- Comprehensive error handling
- Statistical filtering and ranking 
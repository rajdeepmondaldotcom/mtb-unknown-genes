import pandas as pd
import os

# --- Configuration ---
DATA_DIR = 'data'
BLAST_OUTPUT_FILE = os.path.join(DATA_DIR, 'mtb_vs_swissprot.tsv')
FINAL_OUTPUT_FILE = os.path.join(DATA_DIR, 'top_5_percent_homologs.csv')

# Define the column names based on the -outfmt string used in the BLAST command
COLUMN_NAMES = [
    'qseqid', 'sseqid', 'pident', 'ppos', 'length', 'evalue', 'bitscore'
]

# --- Main Script ---
def process_blast_results():
    """
    Reads, sorts, and filters BLAST results to find the top 5% of hits for each query protein.
    """
    print(f"Reading BLAST output from '{BLAST_OUTPUT_FILE}'...")

    # Check if the input file exists and is not empty
    if not os.path.exists(BLAST_OUTPUT_FILE) or os.path.getsize(BLAST_OUTPUT_FILE) == 0:
        print(f"Error: The BLAST output file '{BLAST_OUTPUT_FILE}' was not found or is empty.")
        print("Please ensure the blastp command in Step 3 completed successfully.")
        return

    # Read the tabular BLAST output file, skipping the commented header lines
    try:
        df = pd.read_csv(BLAST_OUTPUT_FILE, sep='\t', comment='#', header=None, names=COLUMN_NAMES)
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{BLAST_OUTPUT_FILE}' contains no data to parse.")
        print("This might happen if the BLAST search yielded no results.")
        return

    print(f"Successfully loaded {len(df)} total hits.")

    # Group the results by the query protein ID ('qseqid')
    grouped = df.groupby('qseqid')
    
    top_hits_list = []

    print("Processing each gene to find top 5% of hits...")
    for query_id, group in grouped:
        # Sort the hits for each gene using a hierarchical approach:
        # 1. E-value (ascending - lower is better)
        # 2. Bit Score (descending - higher is better)
        # 3. Percent Identity (descending - higher is better)
        # 4. Percent Positives (descending - higher is better)
        # This ensures statistical significance is prioritized over simple identity scores.
        sorted_group = group.sort_values(
            by=['evalue', 'bitscore', 'pident', 'ppos'],
            ascending=[True, False, False, False]
        )
        
        # Calculate how many hits represent the top 5%
        num_hits = len(sorted_group)
        top_5_percent_count = max(1, int(num_hits * 0.05))  # Ensure at least one hit is kept
        
        # Select the top N hits and add them to our list
        top_hits_list.append(sorted_group.head(top_5_percent_count))

    # Combine the top hits from all groups into a single final DataFrame
    if top_hits_list:
        final_df = pd.concat(top_hits_list)
        
        # Save the final, shortlisted results to a CSV file
        final_df.to_csv(FINAL_OUTPUT_FILE, index=False)
        print(f"\nSuccess! The final shortlist has been saved to '{FINAL_OUTPUT_FILE}'.")
        print(f"The file contains {len(final_df)} top hits across {len(grouped)} genes.")
    else:
        print("No results to process after grouping. The final file was not created.")

if __name__ == '__main__':
    process_blast_results() 
import os
import pandas as pd
import json
import re
import time
import mygene
from tqdm import tqdm
from pathlib import Path
import argparse
# 🔹 Function to clean ENSEMBL IDs (remove version number only)
def clean_ensembl_id(ensembl_id):
    """
    Cleans ENSEMBL gene IDs by removing the version number.
    
    Example:
        ENSG00000223972.5 → ENSG00000223972
    """
    return ensembl_id.split('.')[0] if '.' in ensembl_id else ensembl_id

# 🔹 Main function to map ENSEMBL IDs to HGNC symbols
def map_ensembl_to_hgnc(input_file, output_file, batch_size=50, log_progress=True):
    """
    Maps ENSEMBL IDs to HGNC gene symbols and filters out transcripts without valid symbols.
    
    Parameters:
    -----------
    input_file: str
        Path to input counts file with ENSEMBL IDs as row indices.
    output_file: str
        Path to output file with mapped HGNC symbols.
    batch_size: int
        Number of IDs to query at once.
    log_progress: bool
        Whether to print progress logs.
    """
    if log_progress:
        print(f"Processing {input_file}...")

    # Read the input file
    try:
        df = pd.read_csv(input_file, index_col=0, dtype={'index': str})
        df.index = df.index.map(str)  # Ensure indices are strings
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return False

    ensembl_ids = df.index.tolist()
    total_transcripts = len(ensembl_ids)
    if log_progress:
        print(f"Found {total_transcripts} ENSEMBL transcripts to map")

    # Clean the ENSEMBL IDs (remove version numbers)
    clean_ids = [clean_ensembl_id(eid) for eid in ensembl_ids if eid is not None]

    # Dictionary to store ENSEMBL-to-HGNC mappings
    ensembl_to_hgnc = {}

    # Use local mapping file if available
    local_mapping_file = "ensembl_to_hgnc_mapping.json"
    if os.path.exists(local_mapping_file):
        with open(local_mapping_file, 'r') as f:
            ensembl_to_hgnc = json.load(f)
    else:
        mg = mygene.MyGeneInfo()
        if log_progress:
            print("Querying gene symbols from MyGene.info API...")
        batches = [clean_ids[i:i+batch_size] for i in range(0, len(clean_ids), batch_size)]
        for batch in tqdm(batches, desc="Processing batches"):
            batch = [b for b in batch if b is not None]
            if not batch:
                continue
            try:
                results = mg.getgenes(batch, fields='symbol', species='human')
                for result in results:
                    if 'symbol' in result and 'query' in result:
                        ensembl_to_hgnc[result['query']] = result['symbol']
                time.sleep(0.2)
            except Exception as e:
                print(f"Error querying batch: {e}")
                time.sleep(1)
                continue
        # Save the mapping for future use
        with open(local_mapping_file, 'w') as f:
            json.dump(ensembl_to_hgnc, f)

    # Create a new index list using the mapped HGNC symbols
    new_index = []
    for eid in df.index:
        clean_id = clean_ensembl_id(eid)
        symbol = ensembl_to_hgnc.get(clean_id, None)
        # If no valid symbol is returned or it's empty, assign an empty string
        if symbol is None or symbol.strip() == "":
            new_index.append("")
        else:
            new_index.append(symbol)

    # Assign the new index
    df.index = new_index
    # Drop rows where the index is missing or empty (after stripping whitespace)
    df = df.dropna()
    df = df[df.index.str.strip() != ""]

    # Rename the index column header to "Geneid"
    df.index.name = "Geneid"

    # Save the mapped data
    df.to_csv(output_file)

    mapped_count = df.shape[0]
    unmapped_count = total_transcripts - mapped_count
    print(f"Mapped {mapped_count} genes to HGNC symbols")
    print(f"{unmapped_count} transcripts could not be mapped (or were removed)")
    print(f"Saved mapped data to {output_file}")

    return True

# Process all files in the filtered_counts directory
def hngc_mapping(input_dir="filtered_counts", output_dir="hgnc_mapped"):
    """
    Process all filtered count files and map to HGNC symbols.
    
    Parameters:
    -----------
    input_dir: str
        Directory containing filtered count files.
    output_dir: str
        Directory for saving HGNC mapped count files.
        
    Returns:
        bool: True if processing was successful, False otherwise
    """
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        filtered_files = [f for f in os.listdir(input_dir) if f.endswith("_filtered_counts.csv")]
        if not filtered_files:
            print(f"No filtered count files found in {input_dir}")
            return False
            
        print(f"Found {len(filtered_files)} filtered count files to process")
        skipped_files = []
        success = True
        
        for filename in filtered_files:
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename.replace("_filtered_counts.csv", "_hgnc_mapped.csv"))
            
            # Skip if output file already exists
            if os.path.exists(output_path):
                file_id = filename.replace("_filtered_counts.csv", "")
                print(f"Skipping {file_id}: Output file already exists")
                skipped_files.append(file_id)
                continue
                
            if not map_ensembl_to_hgnc(input_path, output_path):
                success = False
                break
        
        # Print summary
        processed_count = len(filtered_files) - len(skipped_files)
        if processed_count == 0:
            print("All files have already been processed.")
        else:
            print(f"HGNC mapping complete! Processed {processed_count} files, skipped {len(skipped_files)} files.")
            
        return success
        
    except Exception as e:
        print(f"Error in HGNC mapping: {str(e)}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process filtered count files.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing filtered count files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for saving HGNC mapped count files')
    args = parser.parse_args()
    
    hngc_mapping(args.input_dir, args.output_dir)

import os
import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def filter_raw_counts(input_file, output_file, gene_exp_cutoff=0.6):
    """
    Filter raw count data similar to TCGAanalyze_Preprocessing with a gene expression cutoff of 0.6.

    Parameters:
    -----------
    input_file: str
        Path to input raw counts CSV file.
    output_file: str
        Path to output filtered counts CSV file.
    gene_exp_cutoff: float
        Expression threshold for filtering genes (default = 0.6).
    filter_fold: int
        Minimum number of samples a gene must be expressed in.
    """
    print(f"Processing {input_file}...")

    # Read raw counts file
    try:
        counts_df = pd.read_csv(input_file, index_col=0)
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return False

    if counts_df.shape[0] == 0 or counts_df.shape[1] == 0:
        print(f"File {input_file} has no valid data")
        return False

    print(f"Original dimensions: {counts_df.shape}")

    # 1️⃣ Apply the Gene Expression Cutoff (Main Step from Paper)
    genes_to_keep = counts_df[counts_df.mean(axis=1) > gene_exp_cutoff].index
    counts_df_filtered = counts_df.loc[genes_to_keep]
    print(f"Genes retained after expression cutoff ({gene_exp_cutoff}): {counts_df_filtered.shape[0]}")

    # 2️⃣ Filter genes based on expression in multiple samples
    # non_zero_counts = (counts_df_filtered > 0).sum(axis=1)
    # min_samples = max(filter_fold, int(counts_df_filtered.shape[1] / 3))
    # genes_to_keep_final = non_zero_counts[non_zero_counts >= min_samples].index
    # counts_df_filtered = counts_df_filtered.loc[genes_to_keep_final]
    
    print(f"Filtered dimensions: {counts_df_filtered.shape}")
    print(f"Removed {counts_df.shape[0] - counts_df_filtered.shape[0]} genes")

    # Save filtered data
    counts_df_filtered.to_csv(output_file)
    print(f"Saved filtered data to {output_file}")

    return True

def preprocess_all_files(input_dir="counts", output_dir="filtered_counts"):
    """
    Process all count files in the specified directory.

    Parameters:
    -----------
    input_dir: str
        Directory containing raw count files.
    output_dir: str
        Directory for saving filtered count files.
        
    Returns:
        bool: True if processing was successful, False otherwise
    """
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        raw_files = [f for f in os.listdir(input_dir) if f.endswith("_counts.csv")]

        if not raw_files:
            print(f"No raw count files found in {input_dir}")
            return False

        print(f"Found {len(raw_files)} raw count files to process")
        skipped_files = []
        success = True

        for filename in raw_files:
            file_id = filename.replace("_counts.csv", "")
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, f"{file_id}_filtered_counts.csv")

            # Skip if output file already exists
            if os.path.exists(output_path):
                print(f"Skipping {file_id}: Output file already exists")
                skipped_files.append(file_id)
                continue

            if not filter_raw_counts(input_path, output_path):
                success = False
                break

        # Print summary
        processed_count = len(raw_files) - len(skipped_files)
        if processed_count == 0:
            print("All files have already been processed.")
        else:
            print(f"Processing complete! Processed {processed_count} files, skipped {len(skipped_files)} files.")
            
        return success
        
    except Exception as e:
        print(f"Error in preprocessing: {str(e)}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process raw count files.')    
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing raw count files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for saving filtered count files')
    args = parser.parse_args()
    
    process_all_files(args.input_dir, args.output_dir)
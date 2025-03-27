import os
import pandas as pd
import glob
import argparse

def merge_feature_counts(input_dir, output_dir):
    """
    Merge feature count files from a specified directory into a single DataFrame.
    
    Args:
        input_dir (str): Parent directory path containing project folders with _counts.tsv files
        output_dir (str): Directory where output CSV files are stored
        
    Returns:
        pd.DataFrame: Merged count matrix with genes as rows and samples as columns
    """
    # Find all project directories (PRJNA*) in the input directory
    project_dirs = glob.glob(os.path.join(input_dir, "PRJNA*"))
    
    if not project_dirs:
        raise ValueError(f"No project directories (PRJNA*) found in: {input_dir}")
    
    # Initialize an empty list to store dataframes
    dfs = []
    skipped_projects = []
    
    # Process each project directory
    for project_dir in project_dirs:
        # Get project ID from directory name
        project_id = os.path.basename(project_dir)
        output_file = os.path.join(output_dir, f"{project_id}_counts.csv")
        
        # Skip if output file already exists
        if os.path.exists(output_file):
            print(f"Skipping {project_id}: Output file already exists")
            skipped_projects.append(project_id)
            continue
            
        # Find all _counts.tsv files in the project directory
        count_files = glob.glob(os.path.join(project_dir, "SRR*_counts.tsv"))
        
        if not count_files:
            print(f"Warning: No count files found in project directory: {project_dir}")
            continue
            
        # Process each count file in the project directory
        for file in count_files:
            # Extract SRR number from filename (remove _counts.tsv and .fastq_trimmed)
            srr_id = os.path.basename(file).split("_counts.tsv")[0].split(".fastq_trimmed")[0]
            
            # Read the count file, skipping comment lines
            df = pd.read_csv(file, sep="\t", comment="#")
            
            # Set Geneid as the index
            df.set_index("Geneid", inplace=True)
            
            # Keep only the relevant count column (assumed to be the last column)
            df = df.iloc[:, [-1]]
            
            # Rename the count column to simplified SRR ID
            df.columns = [srr_id]
            
            # Store the cleaned dataframe
            dfs.append(df)
    
    # If all projects were skipped, return None
    if len(skipped_projects) == len(project_dirs):
        print("All projects have already been processed.")
        return None
    
    # If no valid count files were found in any unprocessed project
    if not dfs:
        raise ValueError(f"No valid count files found in any unprocessed project directory under: {input_dir}")
    
    # Merge all dataframes on Geneid using an outer join
    merged_counts = pd.concat(dfs, axis=1, join="outer")
    
    # Fill missing values with 0
    merged_counts.fillna(0, inplace=True)
    
    # Convert counts to integers
    merged_counts = merged_counts.astype(int)
    
    # Save the merged counts for each project
    for project_dir in project_dirs:
        project_id = os.path.basename(project_dir)
        output_file = os.path.join(output_dir, f"{project_id}_counts.csv")
        if project_id not in skipped_projects:
            merged_counts.to_csv(output_file)
            print(f"Saved merged counts to {output_file}")
    
    return merged_counts

def save_merged_counts(input_dir, output_dir):
    """
    Process feature count files and save the merged counts to a CSV file.
    
    Args:
        input_dir (str): Directory path containing _counts.tsv files
        output_dir (str): Directory where the merged counts CSV file should be saved
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get merged counts
    merged_counts = merge_feature_counts(input_dir, output_dir)
    
    # Return True if either:
    # 1. All projects were already processed (merged_counts is None)
    # 2. New counts were successfully processed (merged_counts is not None)
    return True

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Merge feature count files into a single CSV file.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing feature count files (*_counts.tsv)')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory where the merged counts CSV file should be saved')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the counts
    merged_counts = save_merged_counts(args.input_dir, args.output_dir)
    #print(merged_counts.head())  # Display first few rows

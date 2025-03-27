import pandas as pd
import os
from pathlib import Path
import argparse

def filter_genes_by_quantile(counts_df, quantile_cutoff=0.25):
    """
    Filters genes by removing those with average expression below the specified quantile cutoff.
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Normalized gene expression counts (genes as rows, samples as columns).
    quantile_cutoff : float
        The quantile threshold (default 0.25).
    
    Returns:
    --------
    pd.DataFrame
        A DataFrame with only genes that meet or exceed the threshold.
    """
    # Calculate the average expression for each gene across all samples
    gene_avg_expr = counts_df.mean(axis=1)
    
    # Determine the threshold value corresponding to the given quantile
    threshold = gene_avg_expr.quantile(quantile_cutoff)
    print(f"Expression threshold at {quantile_cutoff} quantile: {threshold}")
    
    # Filter: keep genes whose average expression is >= threshold
    filtered_counts = counts_df.loc[gene_avg_expr >= threshold]
    
    return filtered_counts

def filter_counts(input_file, output_file, quantile_cutoff=0.25):
    """
    Filter normalized count data using quantile-based filtering.
    
    Parameters:
    -----------
    input_file : str
        Path to input normalized count file.
    output_file : str
        Path to output filtered count file.
    quantile_cutoff : float
        The quantile threshold for filtering (default 0.25).
        
    Returns:
    --------
    bool
        True if filtering was successful, False otherwise.
    """
    try:
        # Load normalized counts
        counts_df = pd.read_csv(input_file, index_col=0)
        
        # Apply quantile filtering
        filtered_counts = filter_genes_by_quantile(counts_df, quantile_cutoff)
        
        # Report the number of genes remaining
        num_genes_remaining = filtered_counts.shape[0]
        print(f"Filtering complete. {num_genes_remaining} genes remain after filtering.")
        
        # Save the filtered counts
        filtered_counts.to_csv(output_file)
        print(f"Saved filtered counts to {output_file}")
        
        return True
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return False

def filter_all_files(input_dir, output_dir, quantile_cutoff=0.25):
    """
    Process normalized count files for a specific project.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing normalized count files.
    output_dir : str
        Directory for saving filtered count files.
    quantile_cutoff : float
        The quantile threshold for filtering (default 0.25).
        
    Returns:
    --------
    bool
        True if filtering was successful, False otherwise.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Find all normalized count files
    input_files = [f for f in os.listdir(input_dir) if f.endswith("_normalized.csv")]
    
    if not input_files:
        print(f"No normalized count files found in {input_dir}")
        return False
        
    print(f"Found {len(input_files)} normalized count files to process")
    skipped_files = []
    all_successful = True
    
    for filename in input_files:
        # Extract project ID from filename
        project_id = filename.replace("_normalized.csv", "")
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, f"{project_id}_quantile.csv")
        
        # Skip if output file already exists
        if os.path.exists(output_path):
            print(f"Skipping {project_id}: Output file already exists")
            skipped_files.append(project_id)
            continue
            
        # Process the file
        success = filter_counts(input_path, output_path, quantile_cutoff)
        
        if success:
            print(f"Filtering complete for project {project_id}")
        else:
            print(f"Failed to filter counts for project {project_id}")
            all_successful = False
    
    # Print summary
    processed_count = len(input_files) - len(skipped_files)
    if processed_count == 0:
        print("All files have already been processed.")
    else:
        print(f"Filtering complete! Processed {processed_count} files, skipped {len(skipped_files)} files.")
        
    return all_successful

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter normalized count files using quantile-based filtering.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing normalized count files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for saving filtered count files')
    parser.add_argument('-q', '--quantile_cutoff', type=float, default=0.25, help='Quantile threshold for filtering (default: 0.25)')
    args = parser.parse_args()
    
    filter_all_files(args.input_dir, args.output_dir, args.quantile_cutoff)

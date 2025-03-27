import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import os
from pathlib import Path
import argparse

def within_lane_normalization(counts_df, gc_series, frac=0.3):
    """
    Adjust counts for gene-level effects (GC-content) within each sample using LOESS.
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Raw gene counts (genes as rows, samples as columns).
    gc_series : pd.Series
        GC content for each gene (index must match counts_df.index).
    frac : float (default=0.3)
        Fraction of data used for each LOESS fit point.
        
    Returns:
    --------
    pd.DataFrame
        Counts after within-lane normalization.
    """
    norm_counts = counts_df.copy()
    for sample in counts_df.columns:
        # Log-transform counts (adding 1 to avoid log(0))
        log_counts = np.log(counts_df[sample] + 1)
        # Create a DataFrame for the LOESS fit
        data = pd.DataFrame({'log_counts': log_counts, 'gc': gc_series})
        data = data.dropna()  # Remove any missing values
        
        # Fit LOESS: estimate the effect of GC-content on log(count)
        fitted = lowess(endog=data['log_counts'], exog=data['gc'], frac=frac, return_sorted=False)
        
        # Remove the GC-content trend: compute residuals and add back the overall mean
        adjusted_log_counts = data['log_counts'] - fitted + data['log_counts'].mean()
        
        # Transform back to the original count scale and convert to integers
        norm_counts.loc[data.index, sample] = np.round(np.exp(adjusted_log_counts) - 1).astype(int)
    return norm_counts

def between_lane_normalization(counts_df):
    """
    Adjust counts for differences in sequencing depth (library size) between samples.
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Gene counts (raw or after within-lane normalization).
        
    Returns:
    --------
    pd.DataFrame
        Counts scaled to account for sequencing depth differences.
    """
    # Compute library sizes (total counts per sample)
    lib_sizes = counts_df.sum(axis=0)
    # Determine the median library size
    median_lib_size = np.median(lib_sizes)
    # Compute scaling factors for each sample
    scaling_factors = lib_sizes / median_lib_size
    # Divide each sample's counts by its scaling factor
    norm_counts = counts_df.div(scaling_factors, axis=1)
    return norm_counts

def filter_low_expression_genes(counts_df, min_avg_count=1):
    """
    Filter out genes with very low average expression across samples.
    
    Parameters:
    -----------
    counts_df : pd.DataFrame
        Gene counts (after normalization).
    min_avg_count : float (default=1)
        Minimum average count required to keep a gene.
        
    Returns:
    --------
    pd.DataFrame
        Filtered counts DataFrame.
    """
    filtered = counts_df[counts_df.mean(axis=1) >= min_avg_count]
    return filtered

def normalize_counts(input_file, output_file, gc_content_file):
    """
    Normalize count data using within-lane and between-lane normalization.
    
    Parameters:
    -----------
    input_file : str
        Path to input count file.
    output_file : str
        Path to output normalized count file.
    gc_content_file : str
        Path to GC content file.
        
    Returns:
    --------
    bool
        True if normalization was successful, False otherwise.
    """
    try:
        # Load raw counts
        counts_df = pd.read_csv(input_file, index_col=0)
        
        # Load gene GC content with Gene_name as index
        gc_content_df = pd.read_csv(gc_content_file, index_col=1)  # Use Gene_name as index
        
        # Check for duplicates
        dup_count_index = counts_df.index.duplicated().sum()
        dup_gc_index = gc_content_df.index.duplicated().sum()
        
        print(f"Found {dup_count_index} duplicate genes in count data")
        print(f"Found {dup_gc_index} duplicate genes in GC content data")
        
        # Handle duplicates in count data if needed
        if dup_count_index > 0:
            print("Removing duplicate genes from count data (keeping first occurrence)...")
            counts_df = counts_df.loc[~counts_df.index.duplicated(keep='first')]
        
        # Handle duplicates in GC content data
        if dup_gc_index > 0:
            print("Removing duplicate genes from GC content data (keeping first occurrence)...")
            gc_content_df = gc_content_df.loc[~gc_content_df.index.duplicated(keep='first')]
        
        # Extract GC content as a series (using the GC_content column)
        gc_series = gc_content_df['GC_content']
        
        # Now safely align GC content with count data
        common_genes = set(counts_df.index).intersection(set(gc_series.index))
        print(f"Found {len(common_genes)} genes in common between count data and GC content")
        
        # Filter both datasets to include only common genes
        counts_df = counts_df.loc[counts_df.index.isin(common_genes)]
        gc_series = gc_series.loc[gc_series.index.isin(common_genes)]
        
        # Step 1: Within-lane normalization (adjust for GC-content)
        norm_within = within_lane_normalization(counts_df, gc_series, frac=0.3)
        
        # Step 2: Between-lane normalization (adjust for sequencing depth)
        norm_between = between_lane_normalization(norm_within)
        
        # Step 3: Filter lowly expressed genes
        normalized_counts = filter_low_expression_genes(norm_between, min_avg_count=1)
        
        # Report the number of genes remaining
        num_genes_remaining = normalized_counts.shape[0]
        print(f"Normalization complete. {num_genes_remaining} genes remain after normalization.")
        
        # Save the normalized counts
        normalized_counts.to_csv(output_file)
        print(f"Saved normalized counts to {output_file}")
        
        return True
        
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return False

def normalize_all_files(input_dir, output_dir, gc_content_file):
    """
    Process all HGNC mapped count files in the input directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing HGNC mapped count files.
    output_dir : str
        Directory for saving normalized count files.
    gc_content_file : str
        Path to GC content file.
        
    Returns:
    --------
    bool
        True if normalization was successful, False otherwise
    """
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Find all HGNC mapped count files
        input_files = [f for f in os.listdir(input_dir) if f.endswith("_hgnc_mapped.csv")]
        
        if not input_files:
            print(f"No HGNC mapped count files found in {input_dir}")
            return False
            
        print(f"Found {len(input_files)} HGNC mapped count files to process")
        skipped_files = []
        success = True
        
        # First check which files need processing
        files_to_process = []
        for filename in input_files:
            output_path = os.path.join(output_dir, filename.replace("_hgnc_mapped.csv", "_normalized.csv"))
            if os.path.exists(output_path):
                file_id = filename.replace("_hgnc_mapped.csv", "")
                print(f"Skipping {file_id}: Output file already exists")
                skipped_files.append(file_id)
            else:
                files_to_process.append(filename)
        
        # If all files are already processed, return early
        if not files_to_process:
            print("All files have already been processed. No new files to normalize.")
            return True
            
        print(f"Processing {len(files_to_process)} new files...")
        
        # Process only the files that need normalization
        for filename in files_to_process:
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename.replace("_hgnc_mapped.csv", "_normalized.csv"))
            
            if not normalize_counts(input_path, output_path, gc_content_file):
                success = False
                break
        
        # Print summary
        processed_count = len(files_to_process)
        if processed_count == 0:
            print("All files have already been processed.")
        else:
            print(f"Normalization complete! Processed {processed_count} files, skipped {len(skipped_files)} files.")
            
        return success
        
    except Exception as e:
        print(f"Error in normalization: {str(e)}")
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Normalize HGNC mapped count files.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing HGNC mapped count files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for saving normalized count files')
    parser.add_argument('-g', '--gc_content_file', required=True, help='Path to GC content file')
    args = parser.parse_args()
    
    normalize_all_files(args.input_dir, args.output_dir, args.gc_content_file)

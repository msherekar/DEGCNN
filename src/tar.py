import os
import pandas as pd
import scipy.io
import scipy.sparse
import gzip
import numpy as np

# Path where the extracted GSM files are stored
data_dir = "GSE161045_RAW/"

# Identify all sample files
mtx_files = [f for f in os.listdir(data_dir) if f.endswith("matrix.mtx.gz")]

# Initialize an empty dictionary to store sample counts
expression_data = {}
gene_names = None

# Process each sample
for mtx_file in mtx_files:
    # Extract the base sample name (e.g., "GSM4888887_CTRL1")
    sample_name = mtx_file.replace("_matrix.mtx.gz", "")

    # Define full file paths
    mtx_path = os.path.join(data_dir, sample_name + "_matrix.mtx.gz")
    genes_path = os.path.join(data_dir, sample_name + "_features.tsv.gz")
    barcodes_path = os.path.join(data_dir, sample_name + "_barcodes.tsv.gz")

    # Load gene names (features)
    with gzip.open(genes_path, "rt") as f:
        genes = [line.strip().split("\t")[0] for line in f]

    # Load sparse count matrix (MTX format)
    counts = scipy.io.mmread(gzip.open(mtx_path, "rb")).todense()

    # Check number of barcodes
    with gzip.open(barcodes_path, "rt") as f:
        barcodes = [line.strip() for line in f]

    # If multiple barcodes, sum across all columns (bulk RNA-seq samples should have one column)
    if len(barcodes) > 1:
        counts = np.sum(counts, axis=1)  # Sum across all barcodes

    # Ensure counts is 1D
    counts = np.array(counts).flatten()

    # Store only one sample per file
    expression_data[sample_name] = counts

    # Store gene names once (ensuring consistency across all samples)
    if gene_names is None:
        gene_names = genes

# Convert to a Pandas DataFrame
counts_df = pd.DataFrame(expression_data, index=gene_names)

# Save as CSV
counts_df.to_csv("raw_counts/GSE161045_raw_counts.csv")

print("Final merged count matrix saved as 'geo_raw_counts.csv'.")

import os
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import re
import numpy as np

# Activate pandas conversion
pandas2ri.activate()

# Load TCGAanalyze_DEA and limma DEA function in R
ro.r('''
library(TCGAbiolinks)
library(limma)

run_TCGA_DEA <- function(count_data, metadata) {
    library(edgeR)

    # Ensure sample names match by setting rownames from the Sample column
    rownames(metadata) <- metadata$Sample  
    metadata$Condition <- factor(metadata$Condition, levels = c("Control", "Disease"))

    # Convert count data to matrix
    count_matrix <- as.matrix(count_data)

    # Ensure metadata and count matrix have matching samples
    common_samples <- intersect(colnames(count_matrix), rownames(metadata))
    if (length(common_samples) == 0) {
        stop("Error: No matching samples found between metadata and count matrix. Check input data.")
    }

    # Subset both dataframes to include only matching samples
    count_matrix <- count_matrix[, common_samples]
    metadata <- metadata[common_samples, , drop=FALSE]

    # Remove genes with zero variance
    count_matrix <- count_matrix[apply(count_matrix, 1, var) > 0, ]

    # Ensure no zero library sizes (samples with all zero counts)
    sample_sums <- colSums(count_matrix)
    if (any(sample_sums == 0)) {
        stop("Error: At least one sample has zero total counts. Check input data.")
    }

    # Run DEA using TCGAanalyze_DEA
    res <- TCGAanalyze_DEA(
        mat1 = count_matrix[, metadata$Condition == "Control"], 
        mat2 = count_matrix[, metadata$Condition == "Disease"], 
        Cond1type = "Control", 
        Cond2type = "Disease", 
        method = "glmLRT",  
        fdr.cut = 0.01
    )

    # Convert results to DataFrame
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)

    return(res_df)
}


run_limma_DEA <- function(expr_data, metadata) {
    # Convert input data
    expr_matrix <- as.matrix(expr_data)
    colData <- data.frame(condition = metadata$Condition)
    rownames(colData) <- metadata$Sample

    # Ensure Condition is a factor with exactly two levels
    colData$condition <- factor(colData$condition)

    # Error check: Ensure there are exactly two conditions
    if (nlevels(colData$condition) != 2) {
        stop("Error: Metadata must contain exactly two conditions (Control and Disease). Check your input data.")
    }

    # Create design matrix
    design <- model.matrix(~condition, data=colData)

    # Remove intercept column if necessary
    if (colnames(design)[1] == "(Intercept)") {
        design <- design[, -1, drop=FALSE]  # Remove intercept
    }

    # Fit linear model
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)

    # Get results
    results <- topTable(fit, number=Inf, sort.by="p")
    results$Gene <- rownames(results)

    return(results)
}
''')

def preprocess_counts(count_data):
    """
    Check and fix negative values and zero-count issues in the count data.

    Args:
        count_data (pd.DataFrame): The count matrix.

    Returns:
        pd.DataFrame: Processed count matrix.
    """
    # Check for negative values
    if (count_data < 0).sum().sum() > 0:
        print("⚠️ Negative values detected! Checking if data is log-transformed...")

        # Check if the data is log2(TPM+1) transformed (raw counts are usually higher)
        if count_data.max().max() < 20:
            print("🔄 Converting log2-transformed values back to raw counts...")
            count_data = (2 ** count_data) - 1
            count_data[count_data < 0] = 0  # Ensure no negatives
            count_data = count_data.round().astype(int)  # Convert to integer counts
        else:
            print("⚠️ Data contains negative values but is NOT log-transformed. Setting negatives to 0.")
            count_data[count_data < 0] = 0  # Replace negatives with 0

    # Remove genes with zero counts across all samples
    print("📉 Removing zero-count genes...")
    count_data = count_data.loc[(count_data.sum(axis=1) > 0)]

    # Remove samples with zero total counts
    sample_sums = count_data.sum(axis=0)
    if (sample_sums == 0).sum() > 0:
        print(f"⚠️ {sum(sample_sums == 0)} samples have zero counts and will be removed.")
        count_data = count_data.loc[:, sample_sums > 0]

    # Transpose the count data if needed (depends on the expected format in R)
    # count_data = count_data.T  # Uncomment if R expects samples as rows

    print(f"✅ Final processed count matrix: {count_data.shape[0]} genes x {count_data.shape[1]} samples")
    return count_data


def perform_tcga_dea(count_file, metadata_file, output_folder):
    """
    Perform differential expression analysis using TCGAanalyze_DEA.
    """
    # Load count data
    count_data = pd.read_csv(count_file, index_col=0)
    print(f"Original count data shape: {count_data.shape}")
    
    # Preprocessing: Undo log2 transformation if necessary
    if (count_data < 0).values.any():
        count_data = (2 ** count_data) - 1
        count_data[count_data < 0] = 0
        count_data = count_data.round().astype(int)
    
    # Remove genes with zero counts
    count_data = count_data.loc[count_data.sum(axis=1) > 0]
    
    # Load metadata without setting the index (retain 'Sample' column)
    metadata = pd.read_csv(metadata_file)
    
    print(f"Count data columns (samples): {count_data.columns[:3]}...")
    print(f"Metadata Sample column: {metadata['Sample'][:3]}...")
    
    # Ensure sample names match
    common_samples = set(count_data.columns).intersection(set(metadata['Sample']))
    print(f"Matching samples found: {len(common_samples)}")
    if len(common_samples) == 0:
        raise ValueError("No matching sample names between count data and metadata")
    
    # Filter both datasets to include only common samples
    count_data = count_data[list(common_samples)]
    metadata = metadata[metadata['Sample'].isin(common_samples)]
    
    # Do not set metadata index here – let R use the 'Sample' column
    print(f"Final count data shape: {count_data.shape[0]} genes x {count_data.shape[1]} samples")
    print(f"Final metadata shape: {metadata.shape[0]} samples")
    
    # Convert to R objects
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        count_r = ro.conversion.py2rpy(count_data)
        meta_r = ro.conversion.py2rpy(metadata)
    
    # Run DEA using TCGAanalyze_DEA
    dea_results = ro.globalenv['run_TCGA_DEA'](count_r, meta_r)
    
    # Convert R output to Pandas DataFrame
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        results_df = ro.conversion.rpy2py(dea_results)
    
    # Save results
    dataset_name = os.path.basename(count_file).replace(".csv", "")
    output_path = os.path.join(output_folder, f"{dataset_name}_TCGA_DEA_results.csv")
    results_df.to_csv(output_path, index=False)
    
    print(f"✅ DEA results saved for {dataset_name} at {output_path}")


def perform_dea(count_file, metadata_file, output_folder):
    """
    Perform differential expression analysis using limma (for normalized data).

    Args:
        count_file (str): Path to gene expression matrix (FPKM/TPM normalized).
        metadata_file (str): Path to metadata file (Sample, Condition).
        output_folder (str): Folder to save DEA results.
    """
    # Load count data
    count_data = pd.read_csv(count_file, index_col=0)
    
    # Undo log2 transformation if needed
    count_data = (2 ** count_data) - 1
    count_data[count_data < 0] = 0
    count_data = count_data.round().astype(int)
    
    # Load metadata without setting the index
    metadata = pd.read_csv(metadata_file)
    metadata["Condition"] = metadata["Condition"].astype(str)

    # Convert to R objects
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        count_r = ro.conversion.py2rpy(count_data)
        meta_r = ro.conversion.py2rpy(metadata)

    # Run limma DEA in R
    dea_results = ro.globalenv['run_limma_DEA'](count_r, meta_r)

    # Convert R output to Pandas DataFrame
    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        results_df = ro.conversion.rpy2py(dea_results)

    # Save results
    dataset_name = os.path.basename(count_file).replace(".csv", "")
    output_path = os.path.join(output_folder, f"{dataset_name}_limma_DEA_results.csv")
    results_df.to_csv(output_path, index=False)

    print(f"✅ DEA results saved for {dataset_name} at {output_path}")


# Example Usage
counts_folder = "ex_filtered/"  
metadata_folder = "ex_meta_data/"  
output_folder = "deg/"  

os.makedirs(output_folder, exist_ok=True)

# Helper function to extract GSE ID from filename
def extract_gse_id(filename):
    match = re.search(r'(GSE\d+)', filename)
    return match.group(1) if match else None

# Get list of count and metadata files
count_files = {extract_gse_id(f): os.path.join(counts_folder, f)
               for f in os.listdir(counts_folder) if f.endswith(".csv")}
metadata_files = {extract_gse_id(f): os.path.join(metadata_folder, f)
                  for f in os.listdir(metadata_folder) if f.endswith(".csv")}

# Process each dataset
processed = 0
for gse_id in count_files.keys():
    if gse_id in metadata_files:
        print(f"🚀 Processing DEA for {gse_id}...")
        perform_tcga_dea(count_files[gse_id], metadata_files[gse_id], output_folder)
        # Alternatively, call perform_dea for limma-based analysis:
        # perform_dea(count_files[gse_id], metadata_files[gse_id], output_folder)
        processed += 1
    else:
        print(f"⚠️ Warning: No matching metadata file found for {gse_id}, skipping.")

print(f"✅ Successfully processed {processed}/{len(count_files)} datasets")

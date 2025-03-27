import os
import pandas as pd
import numpy as np
import re
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from pathlib import Path
import argparse

# Activate pandas conversion for rpy2
pandas2ri.activate()

# Embedded R code: functions for TCGA DEA using TCGAanalyze_DEA and limma
ro.r('''
library(TCGAbiolinks)
library(limma)
library(edgeR)

run_TCGA_DEA <- function(count_data, metadata) {
    # Print metadata to inspect condition values
    print("Checking metadata structure:")
    print(head(metadata))
    print("Condition values:")
    print(unique(metadata$Condition))
    
    # Print sample names from both datasets
    print("Count matrix column names:")
    print(colnames(count_data))
    print("Metadata sample names:")
    print(metadata$Sample)
    
    # Check for NAs in count data
    if(any(is.na(count_data))) {
        print("WARNING: NA values found in count_data. Replacing with zeros.")
        count_data[is.na(count_data)] <- 0
    }
    
    # Ensure counts are integers (edgeR requirement)
    count_data <- round(count_data)
    
    # Clean condition values to ensure exact match
    metadata$Condition <- trimws(metadata$Condition)
    
    # Verify we have the expected condition values
    if(!all(c("Control", "Disease") %in% unique(metadata$Condition))) {
        stop("Error: Metadata must contain exactly 'Control' and 'Disease' conditions.")
    }
    
    # Convert count data to matrix
    count_matrix <- as.matrix(count_data)
    
    # Additional NA check after conversion
    if(any(is.na(count_matrix))) {
        print("WARNING: NA values found after matrix conversion. Replacing with zeros.")
        count_matrix[is.na(count_matrix)] <- 0
    }
    
    # Create a data frame for metadata with rownames
    metadata_df <- data.frame(
        Condition = metadata$Condition,
        row.names = metadata$Sample
    )
    
    # Ensure metadata has the same order as count matrix
    metadata_df <- metadata_df[colnames(count_matrix), , drop=FALSE]
    
    # Convert condition to factor with specific levels
    metadata_df$Condition <- factor(metadata_df$Condition, levels = c("Control", "Disease"))
    
    # Debug messages to track condition assignments
    print("Control samples:")
    print(rownames(metadata_df)[metadata_df$Condition == "Control"])
    print("Disease samples:")
    print(rownames(metadata_df)[metadata_df$Condition == "Disease"])
    
    # Final check to ensure no NAs before DGEList creation
    if(any(is.na(count_matrix))) {
        print("CRITICAL: NAs still present. Replacing with zeros as last resort.")
        count_matrix[is.na(count_matrix)] <- 0
    }
    
    # Remove genes with zero variance
    var_by_gene <- apply(count_matrix, 1, var)
    if(any(is.na(var_by_gene))) {
        print("WARNING: NAs found in variance calculation. Removing problematic genes.")
        var_by_gene[is.na(var_by_gene)] <- 0
    }
    count_matrix <- count_matrix[var_by_gene > 0, , drop=FALSE]
    
    # Check that no sample has zero total counts
    sample_sums <- colSums(count_matrix)
    if (any(sample_sums == 0)) {
        stop("Error: At least one sample has zero total counts. Check input data.")
    }
    
    # Print some diagnostics
    print(paste("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples"))
    print(paste("Number of Control samples:", sum(metadata_df$Condition == "Control")))
    print(paste("Number of Disease samples:", sum(metadata_df$Condition == "Disease")))
    
    # Safely create DGEList with error handling
    dge <- tryCatch({
        DGEList(counts = count_matrix, 
                group = metadata_df$Condition)
    }, error = function(e) {
        message("Error in DGEList creation: ", e$message)
        message("Attempting to identify problematic values...")
        
        # Try to identify problematic values
        problem_rows <- apply(count_matrix, 1, function(x) any(is.na(x) | is.infinite(x) | x < 0))
        if(any(problem_rows)) {
            message("Found ", sum(problem_rows), " problematic genes. Removing them.")
            count_matrix <- count_matrix[!problem_rows, , drop=FALSE]
        }
        
        # Try again with cleaned matrix
        DGEList(counts = count_matrix, 
                group = metadata_df$Condition)
    })
    
    # Step 1: Perform DEA with FDR cutoff of 0.01 to get differentially expressed genes
    res <- tryCatch({
        TCGAanalyze_DEA(
            mat1 = count_matrix[, metadata_df$Condition == "Control", drop=FALSE], 
            mat2 = count_matrix[, metadata_df$Condition == "Disease", drop=FALSE], 
            Cond1type = "Control", 
            Cond2type = "Disease", 
            method = "glmLRT",  
            fdr.cut = 0.01  # This will give us only differentially expressed genes
        )
    }, error = function(e) {
        message("Error in TCGAanalyze_DEA: ", e$message)
        message("Falling back to simplified DEA using edgeR directly...")
        
        # Fallback to direct edgeR analysis if TCGAanalyze_DEA fails
        y <- DGEList(counts = count_matrix, group = metadata_df$Condition)
        y <- calcNormFactors(y)
        design <- model.matrix(~metadata_df$Condition)
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design)
        qlf <- glmQLFTest(fit, coef = 2)
        
        # Convert edgeR results to similar format as TCGAanalyze_DEA
        tt <- topTags(qlf, n=Inf)$table
        tt$Gene <- rownames(tt)
        return(tt)
    })
    
    # Convert results to data frame and add gene names as a column
    if (is.null(res) || nrow(res) == 0) {
        message("No differentially expressed genes found. Creating empty result with all genes.")
        res_df <- data.frame(
            Gene = rownames(count_matrix),
            logFC = rep(0, nrow(count_matrix)),
            PValue = rep(1, nrow(count_matrix)),
            FDR = rep(1, nrow(count_matrix))
        )
    } else {
        res_df <- as.data.frame(res)
        if (!"Gene" %in% colnames(res_df)) {
            res_df$Gene <- rownames(res_df)
        }
    }
    
    # Step 2: Create a data frame with all genes and their counts
    all_genes_df <- data.frame(
        Gene = rownames(count_matrix),
        stringsAsFactors = FALSE
    )
    
    # Add mean counts for each condition
    control_samples <- metadata_df$Condition == "Control"
    disease_samples <- metadata_df$Condition == "Disease"
    
    all_genes_df$Control_mean <- rowMeans(count_matrix[, control_samples, drop=FALSE])
    all_genes_df$Disease_mean <- rowMeans(count_matrix[, disease_samples, drop=FALSE])
    
    # Merge with DEA results
    final_df <- merge(all_genes_df, res_df, by = "Gene", all = TRUE)
    
    # Add labels:
    # - For differentially expressed genes (FDR < 0.01):
    #   - If logFC < 0, label "0" (down-regulated)
    #   - If logFC > 0, label "1" (up-regulated)
    # - For non-differentially expressed genes (FDR >= 0.01 or NA), label "2" (neutral)
    final_df$Label <- ifelse(is.na(final_df$FDR) | final_df$FDR >= 0.01, "2",
                            ifelse(final_df$logFC < 0, "0", "1"))
    
    # Fill NA values in DEA results with appropriate values for non-DE genes
    final_df$logFC[is.na(final_df$logFC)] <- 0
    final_df$FDR[is.na(final_df$FDR)] <- 1
    final_df$PValue[is.na(final_df$PValue)] <- 1
    
    return(final_df)
}

run_limma_DEA <- function(expr_data, metadata) {
    expr_matrix <- as.matrix(expr_data)
    colData <- data.frame(condition = metadata$Condition)
    rownames(colData) <- metadata$Sample
    colData$condition <- factor(colData$condition)
    
    if (nlevels(colData$condition) != 2) {
        stop("Error: Metadata must contain exactly two conditions (Control and Disease).")
    }
    
    design <- model.matrix(~condition, data = colData)
    if (colnames(design)[1] == "(Intercept)") {
        design <- design[, -1, drop=FALSE]
    }
    
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)
    results <- topTable(fit, number=Inf, sort.by="p")
    results$Gene <- rownames(results)
    
    results$Label <- ifelse(results$`adj.P.Val` < 0.01,
                            ifelse(results$logFC < 0, "0", ifelse(results$logFC > 0, "1", NA)),
                            "2")
    
    return(results)
}
''')

def preprocess_counts(count_data):
    """
    Preprocess the count data: fixes negative values, removes genes with zero counts, and removes zero-sum samples.
    """
    print("🔍 Starting count data preprocessing...")
    
    # Check data types and ensure numeric
    try:
        count_data = count_data.astype(float)
        print("✅ Converted data to float type")
    except Exception as e:
        print(f"⚠️ Error converting data to float: {e}")
        # Try to fix non-numeric values
        for col in count_data.columns:
            count_data[col] = pd.to_numeric(count_data[col], errors='coerce')
    
    # Clean up column names and index
    count_data.columns = count_data.columns.str.strip()
    count_data.index = count_data.index.str.strip()
    
    # Strip the '_1_trimmed.fastq.gz' suffix from column names
    original_columns = count_data.columns.tolist()
    count_data.columns = [col.replace('_1_trimmed.fastq.gz', '') for col in count_data.columns]
    print(f"✅ Cleaned up column names and index")
    if original_columns != count_data.columns.tolist():
        print("✅ Stripped '_1_trimmed.fastq.gz' suffix from column names")
    
    # Check for NA values
    na_count = count_data.isna().sum().sum()
    if na_count > 0:
        print(f"⚠️ {na_count} NA values detected! Replacing with 0...")
        count_data = count_data.fillna(0)
    
    # Check for infinite values
    inf_count = np.isinf(count_data).sum().sum()
    if inf_count > 0:
        print(f"⚠️ {inf_count} infinite values detected! Replacing with 0...")
        count_data = count_data.replace([np.inf, -np.inf], 0)
    
    # Convert floating point numbers to integers
    print("🔄 Converting floating point numbers to integers...")
    count_data = count_data.round().astype(int)
    
    # Check for negative values and if counts appear log-transformed
    neg_count = (count_data < 0).sum().sum()
    if neg_count > 0:
        print(f"⚠️ {neg_count} negative values detected! Replacing with 0...")
        count_data[count_data < 0] = 0
    
    # Final NA check
    if count_data.isna().sum().sum() > 0:
        print("⚠️ NA values still present after preprocessing! Applying final cleanup...")
        count_data = count_data.fillna(0)
    
    print("📉 Removing genes with zero counts...")
    count_data = count_data.loc[count_data.sum(axis=1) > 0]
    
    sample_sums = count_data.sum(axis=0)
    if (sample_sums == 0).sum() > 0:
        print(f"⚠️ {sum(sample_sums == 0)} samples have zero counts and will be removed.")
        count_data = count_data.loc[:, sample_sums > 0]
    
    print(f"✅ Final processed count matrix: {count_data.shape[0]} genes x {count_data.shape[1]} samples")
    print("Sample names after preprocessing:")
    print(count_data.columns.tolist())
    return count_data

def perform_tcga_dea(count_file, metadata_file, output_file):
    """
    Load the count data and metadata, preprocess the count data, then perform DE analysis using run_TCGA_DEA.
    """
    try:
        print(f"\n🔍 Processing file: {os.path.basename(count_file)}")
        
        # Load count data with error handling
        try:
            count_data = pd.read_csv(count_file, index_col=0)
            print(f"Count data shape: {count_data.shape}")
            print("Count data samples:", sorted(count_data.columns.tolist()))
        except Exception as e:
            print(f"⚠️ Error reading count file: {e}")
            print("Trying to read file with different settings...")
            count_data = pd.read_csv(count_file, index_col=0, on_bad_lines='skip')
            print(f"Count data shape after skipping bad lines: {count_data.shape}")
            print("Count data samples:", sorted(count_data.columns.tolist()))
        
        # Preprocess the count data
        count_data = preprocess_counts(count_data)
        
        # Load metadata with error handling
        try:
            metadata = pd.read_csv(metadata_file)
            print(f"Metadata shape: {metadata.shape}")
            print("Metadata samples:", sorted(metadata['Sample'].tolist()))
        except Exception as e:
            print(f"⚠️ Error reading metadata file: {e}")
            return False
        
        # Check if required columns exist in metadata
        if 'Sample' not in metadata.columns or 'Condition' not in metadata.columns:
            print("⚠️ Error: Metadata file must contain 'Sample' and 'Condition' columns.")
            return False
            
        # Clean metadata
        for col in metadata.columns:
            if metadata[col].dtype == 'object':  # Only clean string columns
                metadata[col] = metadata[col].str.strip()
        
        # Check condition values - must be exactly "Control" and "Disease"
        valid_conditions = set(metadata['Condition'].unique())
        if "Control" not in valid_conditions or "Disease" not in valid_conditions:
            print(f"⚠️ Warning: Expected conditions 'Control' and 'Disease', found {valid_conditions}")
            
            # Map similar conditions to expected values
            condition_map = {}
            for cond in valid_conditions:
                if cond.lower().strip() in ["control", "ctrl", "normal", "healthy"]:
                    condition_map[cond] = "Control"
                elif cond.lower().strip() in ["disease", "tumor", "cancer", "sick"]:
                    condition_map[cond] = "Disease"
            
            if condition_map:
                print(f"🔄 Mapping conditions: {condition_map}")
                metadata['Condition'] = metadata['Condition'].replace(condition_map)
        
        # Find common samples between count data and metadata
        count_samples = set(count_data.columns)
        meta_samples = set(metadata['Sample'])
        common_samples = count_samples.intersection(meta_samples)
        
        if len(common_samples) == 0:
            print("⚠️ Error: No matching sample names between count data and metadata.")
            print("Count data samples:", sorted(list(count_samples)))
            print("Metadata samples:", sorted(list(meta_samples)))
            return False
        
        # Subset data to common samples
        count_data = count_data[list(common_samples)]
        metadata = metadata[metadata['Sample'].isin(common_samples)]
        
        # Print condition counts
        condition_counts = metadata['Condition'].value_counts()
        print(f"Condition counts: {condition_counts.to_dict()}")
        
        # Make sure we have both conditions
        if len(condition_counts) < 2:
            print("⚠️ Error: Need at least two conditions (Control and Disease) in metadata.")
            return False
            
        # Final check to ensure conditions are exactly "Control" and "Disease"
        if "Control" not in condition_counts.index or "Disease" not in condition_counts.index:
            print("⚠️ Error: Conditions must be exactly 'Control' and 'Disease'.")
            return False
        
        # Check for NAs before conversion to R
        if count_data.isna().any().any():
            print("⚠️ Warning: NA values present before R conversion. Fixing...")
            count_data = count_data.fillna(0)
        
        # Convert to R objects with verbose errors
        try:
            with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
                count_r = ro.conversion.py2rpy(count_data)
                meta_r = ro.conversion.py2rpy(metadata)
        except Exception as e:
            print(f"⚠️ Error during Python to R conversion: {e}")
            return False
        
        # Run the DEA analysis
        try:
            print("🧬 Running differential expression analysis...")
            dea_results = ro.globalenv['run_TCGA_DEA'](count_r, meta_r)
            print("✅ DEA analysis completed")
        except Exception as e:
            print(f"⚠️ Error during DEA analysis: {e}")
            return False
        
        # Convert results back to Python
        try:
            with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
                results_df = ro.conversion.rpy2py(dea_results)
        except Exception as e:
            print(f"⚠️ Error converting results back to Python: {e}")
            return False
        
        # Check results before saving
        if results_df.empty:
            print("⚠️ Warning: Empty results dataframe.")
            return False
        
        # Count labels for summary
        label_counts = results_df['Label'].value_counts().to_dict()
        up_regulated = label_counts.get('1', 0)
        down_regulated = label_counts.get('0', 0)
        neutral = label_counts.get('2', 0)
        
        print(f"\n📊 DEA Results Summary:")
        print(f"  - Up-regulated genes (Label 1): {up_regulated}")
        print(f"  - Down-regulated genes (Label 0): {down_regulated}")
        print(f"  - Neutral genes (Label 2): {neutral}")
        print(f"  - Total genes: {len(results_df)}")
        
        # Save results to CSV
        results_df.to_csv(output_file, index=False)
        print(f"✅ DEA results saved to {output_file}")
        return True
        
    except Exception as e:
        print(f"❌ Error processing {count_file}: {e}")
        import traceback
        traceback.print_exc()
        return False

def extract_identifier(filename):
    """
    Extract a general identifier from the filename by taking the first token before an underscore.
    """
    basename = os.path.splitext(filename)[0]
    parts = basename.split("_")
    return parts[0] if parts else basename

def differential_expression_analysis(input_dir, metadata_dir, output_dir):
    """
    Process filtered count files with metadata files from a directory.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing quantile count files.
    metadata_dir : str
        Directory containing metadata files.
    output_dir : str
        Directory for saving DEA results.
        
    Returns:
    --------
    bool
        True if DEA was successful, False otherwise
    """
    try:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Find all filtered count files
        input_files = {extract_identifier(f): os.path.join(input_dir, f)
                      for f in os.listdir(input_dir) if f.endswith("_quantile.csv")}
        
        if not input_files:
            print(f"No filtered count files found in {input_dir}")
            return False
            
        print(f"Found {len(input_files)} filtered count files to process")
        skipped_files = []
        failed_files = []
        success = True
        
        # Find all metadata files
        metadata_files = {extract_identifier(f): os.path.join(metadata_dir, f)
                         for f in os.listdir(metadata_dir) if f.endswith("_metadata.csv")}
        
        if not metadata_files:
            print(f"No metadata files found in {metadata_dir}")
            return False
            
        print(f"Found {len(metadata_files)} metadata files")
        
        # First check which files need processing
        files_to_process = []
        for identifier in input_files:
            if identifier not in metadata_files:
                print(f"⚠️ No matching metadata file found for {identifier}")
                failed_files.append(identifier)
                continue
                
            output_file = os.path.join(output_dir, f"{identifier}_sdeg.csv")
            if os.path.exists(output_file):
                print(f"Skipping {identifier}: Output file already exists")
                skipped_files.append(identifier)
            else:
                files_to_process.append(identifier)
        
        # If all files are already processed, return early
        if not files_to_process:
            print("All files have already been processed. No new files to analyze.")
            return True
            
        print(f"Processing {len(files_to_process)} new files...")
        
        # Process only the files that need DEA
        for identifier in files_to_process:
            output_file = os.path.join(output_dir, f"{identifier}_sdeg.csv")
            print(f"🚀 Processing DEA for {identifier}...")
            try:
                if not perform_tcga_dea(input_files[identifier], metadata_files[identifier], output_file):
                    print(f"❌ Failed to process {identifier}")
                    failed_files.append(identifier)
                    continue
            except Exception as e:
                print(f"❌ Error processing {identifier}: {str(e)}")
                failed_files.append(identifier)
                continue
        
        # Print summary
        processed_count = len(files_to_process) - len(failed_files)
        print("\n📊 Processing Summary:")
        print(f"  - Total files found: {len(input_files)}")
        print(f"  - Files skipped (already processed): {len(skipped_files)}")
        print(f"  - Files successfully processed: {processed_count}")
        print(f"  - Files failed: {len(failed_files)}")
        
        if failed_files:
            print("\n❌ Failed files:")
            for f in failed_files:
                print(f"  - {f}")
        
        # Return True if at least one file was processed successfully
        return processed_count > 0
        
    except Exception as e:
        print(f"Error in differential expression analysis: {str(e)}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform differential expression analysis on quantile filtered count files.')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing quantile filtered count files')
    parser.add_argument('-m', '--metadata_dir', required=True, help='Directory containing metadata files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for saving DEA results')
    args = parser.parse_args()
    
    differential_expression_analysis(args.input_dir, args.metadata_dir, args.output_dir)

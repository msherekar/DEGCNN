import os
import pickle as pkl
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from statsmodels.stats.multitest import multipletests

# 📌 **Step 1: Load Count Matrix**
counts_df = pd.read_csv(
    "/Users/mukulsherekar/pythonProject/DEGCNN/hgnc_mapped/PRJNA767074_hgnc_mapped.csv", 
    index_col=0
)  # Genes as rows, samples as columns
counts_df = counts_df.T  # ✅ Transpose so samples are rows, genes are columns

# 📌 **Step 2: Load Metadata**
metadata = pd.read_csv("/Users/mukulsherekar/pythonProject/DEGCNN/meta_data/GSE184942_metadata.csv")  
metadata.set_index("Sample", inplace=True)  # ✅ Ensure "Sample" column is set as index

# ✅ **Ensure all samples in metadata exist in count matrix**
common_samples = metadata.index.intersection(counts_df.index)
counts_df = counts_df.loc[common_samples]  # Keep only matched samples
metadata = metadata.loc[common_samples]  # Keep only matched metadata

print(f"\n✅ Matched Samples: {len(common_samples)}")

# 📌 **Step 3: Remove Low-Count Genes**
genes_before_filtering = counts_df.shape[1]
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]  # Remove genes with <10 total counts
counts_df = counts_df[genes_to_keep]
print(f"🔹 Genes before filtering: {genes_before_filtering}, after filtering: {counts_df.shape[1]}")

# 📌 **Step 4: Use Control Samples for Normalization**
# ✅ Ensure "Condition" column is categorical and set Control as the reference
metadata["Condition"] = metadata["Condition"].astype("category")
metadata["Condition"] = metadata["Condition"].cat.reorder_categories(["Control", "Disease"], ordered=True)

# 📌 **Step 5: Create PyDESeq2 Object (No control_genes, just reference Control samples)**
inference = DefaultInference(n_cpus=8)  # Parallel processing

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~Condition",  # ✅ Ensures Control is the reference
    refit_cooks=True,
    inference=inference,
)

# 📌 **Step 6: Run Differential Expression Analysis (DEA)**
dds.deseq2()
print(dds)
print(dds.varm["dispersions"])
print(dds.varm["LFC"])

# 📌 **Step 7: Compute Statistical Significance**
ds = DeseqStats(dds, contrast=["Condition", "Control", "Disease"], inference=inference)
ds.summary()

# 📌 **Step 8: Extract DEGs**
results_df = ds.results_df  # Get results (log2FC, p-values, padj)
#print("🔹 First 5 DESeq2 results:\n", results_df.head())
#print(f"✅ Total genes in results: {len(results_df)}")

# 📌 **Step 9: Apply Quantile Filtering (Keep Top 75% Expressed Genes)**
quantile_cutoff = results_df["baseMean"].quantile(0.25)  # 25% quantile cutoff
filtered_results = results_df[results_df["baseMean"] > quantile_cutoff]

print(f"🔹 Genes after quantile filtering: {len(filtered_results)}")
print(filtered_results.head())

# ================================
# Step 10: Apply False Discovery Rate (FDR) Correction
# ================================
filtered_results = filtered_results.copy()  # Fix SettingWithCopyWarning
filtered_results["padj"] = multipletests(filtered_results["pvalue"], method="fdr_bh")[1]

# Remove NaN values
filtered_results = filtered_results.dropna(subset=["pvalue", "padj"])

# ================================
# Step 11: Label Genes Based on FDR and log2FoldChange
# ================================
# Labeling scheme:
#   - 0: Significant and downregulated (padj < 0.05 and log2FoldChange < 0)
#   - 1: Significant and upregulated (padj < 0.05 and log2FoldChange > 0)
#   - 2: Not significant (padj >= 0.05)
def label_gene(row, fdr_cutoff=0.1):
    if row["padj"] < fdr_cutoff:
        return 0 if row["log2FoldChange"] < 0 else 1
    else:
        return 2

filtered_results["DEG_label"] = filtered_results.apply(label_gene, axis=1)

# ================================
# Step 12: Save Final Results Based on the New Cutoff
# ================================
# Now using an FDR cutoff of 0.05 (change to 0.1 if preferred)
final_results = filtered_results[filtered_results["padj"] < 0.1]

print(f"🔹 Genes passing FDR < 0.1: {len(final_results)}")
print(final_results.head())

if len(final_results) == 0:
    print("⚠️ No DEGs found with padj < 0.1. Consider adjusting the cutoff (e.g., to 0.1) or reviewing the data.")
else:
    final_results.to_csv("PyDESeq2_filtered_DEGs.csv", index=True)
    print(f"✅ Done! {len(final_results)} DEGs saved in PyDESeq2_filtered_DEGs.csv")


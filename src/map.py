# map hgnc genes to ensembl ids
import pandas as pd
import json

# Load the Ensembl-to-HGNC mapping
with open("ensembl_to_hgnc_mapping.json", "r") as f:
    ensembl_to_hgnc = json.load(f)  # Mapping of {ENSEMBL_ID: HGNC_NAME}

# Load GC-content file
gc_content = pd.read_csv("gc_content.csv")  # Your existing GC-content file
print(gc_content.head(5))
# Map Ensembl IDs to HGNC names
gc_content["HGNC_gene_name"] = gc_content["Gene_stable_ID"].map(ensembl_to_hgnc)
print(gc_content.head(5))
# Save the new file
gc_content.to_csv("gc_content_with_hgnc.csv", index=False)

print("New GC-content file saved as gc_content_with_hgnc.csv")

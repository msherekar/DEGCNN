import os
import pandas as pd
import glob

dataset_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/training_datasets"

# Get all CSV files
csv_files = glob.glob(os.path.join(dataset_dir, "*.csv"))

for file_path in csv_files:
    filename = os.path.basename(file_path)
    
    # Skip *_genes.csv files
    if filename.endswith("_genes.csv"):
        print(f"Skipping gene list file: {filename}")
        continue
    
    print(f"Processing: {filename}")
    
    try:
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Remove the first column (column at index 0)
        df = df.iloc[:, 1:]
        
        # Write back to the same file without index
        df.to_csv(file_path, index=False)
        
        print(f"Successfully removed first column from {filename}")
    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")

print("First column removal completed for all applicable CSV files") 
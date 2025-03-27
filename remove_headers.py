import os
import glob

dataset_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/training_datasets"

# Get all CSV files
csv_files = glob.glob(os.path.join(dataset_dir, "*.csv"))
label_files = glob.glob(os.path.join(dataset_dir, "*_label.txt"))

all_files = csv_files + label_files

for file_path in all_files:
    print(f"Processing: {os.path.basename(file_path)}")
    
    # Read all lines
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Skip the first line (header)
    if len(lines) > 1:
        with open(file_path, 'w') as f:
            f.writelines(lines[1:])
    
    print(f"Removed header from {os.path.basename(file_path)}")

print("Header removal completed for all files") 
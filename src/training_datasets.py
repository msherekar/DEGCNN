import pandas as pd
import os
import argparse
from datetime import datetime
from colorama import init, Fore, Style
from sklearn.model_selection import train_test_split

# Initialize colorama for colored output
init()

def save_data_and_labels_and_genes(df, output_base_path):
    """Save data, label, and genes files"""
    # Save the main data file (all columns)
    df.to_csv(f"{output_base_path}.csv", index=False)
    
    # Save the label file (just the Label column)
    df['Label'].to_csv(f"{output_base_path}_label.txt", index=False, header=False)
    
    # Save the genes file (just the Gene_id column)
    df[['Geneid']].to_csv(f"{output_base_path}_genes.csv", index=False)
    
    return len(df)

def split_and_save_data(df, base_name, output_dir, is_bio_data=False):
    """Split data into train/test sets and save data, label, and genes files"""
    # Split the data
    if is_bio_data:
        # For bio data: 80% test (T3), 20% fine-tune (F1)
        test_df, finetune_df = train_test_split(df, train_size=0.8, random_state=42, stratify=df['Label'])
        
        # Get label counts
        test_label_counts = test_df['Label'].value_counts().to_dict()
        finetune_label_counts = finetune_df['Label'].value_counts().to_dict()
        
        # Save test data (T3)
        test_path = os.path.join(output_dir, f"{base_name}_bio_test")
        test_size = save_data_and_labels_and_genes(test_df, test_path)
        
        # Save fine-tune data (F1)
        finetune_path = os.path.join(output_dir, f"{base_name}_bio_fine_tune")
        finetune_size = save_data_and_labels_and_genes(finetune_df, finetune_path)
        
        return test_size, finetune_size, test_label_counts, finetune_label_counts
    else:
        # For non-bio data: 80% train (T1), 20% test (T2)
        train_df, test_df = train_test_split(df, train_size=0.8, random_state=42, stratify=df['Label'])
        
        # Get label counts
        train_label_counts = train_df['Label'].value_counts().to_dict()
        test_label_counts = test_df['Label'].value_counts().to_dict()
        
        # Save train data (T1)
        train_path = os.path.join(output_dir, f"{base_name}_non_bio_train")
        train_size = save_data_and_labels_and_genes(train_df, train_path)
        
        # Save test data (T2)
        test_path = os.path.join(output_dir, f"{base_name}_non_bio_test")
        test_size = save_data_and_labels_and_genes(test_df, test_path)
        
        return train_size, test_size, train_label_counts, test_label_counts

def process_dataset(input_dir, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create report file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = open(os.path.join(output_dir, f"training_split_report_{timestamp}.txt"), "w")
    report_file.write(f"Training Dataset Split Report\n")
    report_file.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    report_file.write(f"{'='*50}\n\n")
    
    # Process each P and Q file
    for filename in os.listdir(input_dir):
        if filename.endswith('_P.csv'):  # Process non-bio data
            print(f"\n{Fore.CYAN}Processing non-bio data {filename}...{Style.RESET_ALL}")
            input_path = os.path.join(input_dir, filename)
            df = pd.read_csv(input_path)
            
            # Check for NaN values
            nan_count = df.isna().sum().sum()
            rows_with_nan = df[df.isna().any(axis=1)]
            nan_rows_count = len(rows_with_nan)
            
            print(f"Found {Fore.YELLOW}{nan_count}{Style.RESET_ALL} NaN values in {Fore.YELLOW}{nan_rows_count}{Style.RESET_ALL} rows")
            report_file.write(f"\nNaN values in {filename}:\n")
            report_file.write(f"Total NaN values: {nan_count}\n")
            report_file.write(f"Rows with NaN values: {nan_rows_count}\n")
            
            if nan_rows_count > 0:
                # Remove rows with NaN values
                original_count = len(df)
                df = df.dropna()
                removed_count = original_count - len(df)
                print(f"{Fore.RED}Removed {removed_count} rows with NaN values{Style.RESET_ALL}")
                report_file.write(f"Removed {removed_count} rows with NaN values\n")
            
            base_name = os.path.splitext(filename)[0].replace('_P', '')
            
            # Split into T1 (train) and T2 (test)
            train_size, test_size, train_label_counts, test_label_counts = split_and_save_data(df, base_name, output_dir, is_bio_data=False)
            
            # Write report for non-bio split
            report_file.write(f"\nNon-bio data split ({filename}):\n")
            report_file.write(f"T1 (train) size: {train_size}\n")
            report_file.write(f"T2 (test) size: {test_size}\n")
            
            # Write label distribution
            report_file.write(f"\nLabel distribution in train set:\n")
            for label, count in train_label_counts.items():
                report_file.write(f"  Label {label}: {count} genes ({count/train_size*100:.1f}%)\n")
                
            report_file.write(f"\nLabel distribution in test set:\n")
            for label, count in test_label_counts.items():
                report_file.write(f"  Label {label}: {count} genes ({count/test_size*100:.1f}%)\n")
                
            report_file.write(f"\nOutput files:\n")
            report_file.write(f"- {base_name}_non_bio_train.csv, _label.txt, and _genes.csv\n")
            report_file.write(f"- {base_name}_non_bio_test.csv, _label.txt, and _genes.csv\n")
            
            print(f"Non-bio split complete:")
            print(f"T1 (train) size: {Fore.GREEN}{train_size}{Style.RESET_ALL}")
            print(f"T2 (test) size: {Fore.YELLOW}{test_size}{Style.RESET_ALL}")
            
            # Print label distribution
            print(f"\nLabel distribution in train set:")
            for label, count in train_label_counts.items():
                print(f"  Label {label}: {Fore.BLUE}{count}{Style.RESET_ALL} genes ({Fore.BLUE}{count/train_size*100:.1f}%{Style.RESET_ALL})")
                
            print(f"\nLabel distribution in test set:")
            for label, count in test_label_counts.items():
                print(f"  Label {label}: {Fore.BLUE}{count}{Style.RESET_ALL} genes ({Fore.BLUE}{count/test_size*100:.1f}%{Style.RESET_ALL})")
            
            print(f"Output files:")
            print(f"- {base_name}_non_bio_train.csv, _label.txt, and _genes.csv")
            print(f"- {base_name}_non_bio_test.csv, _label.txt, and _genes.csv")
            
        elif filename.endswith('_Q.csv'):  # Process bio data
            print(f"\n{Fore.CYAN}Processing bio data {filename}...{Style.RESET_ALL}")
            input_path = os.path.join(input_dir, filename)
            df = pd.read_csv(input_path)
            
            # Check for NaN values
            nan_count = df.isna().sum().sum()
            rows_with_nan = df[df.isna().any(axis=1)]
            nan_rows_count = len(rows_with_nan)
            
            print(f"Found {Fore.YELLOW}{nan_count}{Style.RESET_ALL} NaN values in {Fore.YELLOW}{nan_rows_count}{Style.RESET_ALL} rows")
            report_file.write(f"\nNaN values in {filename}:\n")
            report_file.write(f"Total NaN values: {nan_count}\n")
            report_file.write(f"Rows with NaN values: {nan_rows_count}\n")
            
            if nan_rows_count > 0:
                # Remove rows with NaN values
                original_count = len(df)
                df = df.dropna()
                removed_count = original_count - len(df)
                print(f"{Fore.RED}Removed {removed_count} rows with NaN values{Style.RESET_ALL}")
                report_file.write(f"Removed {removed_count} rows with NaN values\n")
            
            base_name = os.path.splitext(filename)[0].replace('_Q', '')
            
            # Split into T3 (test) and F1 (fine-tune)
            test_size, finetune_size, test_label_counts, finetune_label_counts = split_and_save_data(df, base_name, output_dir, is_bio_data=True)
            
            # Write report for bio split
            report_file.write(f"\nBio data split ({filename}):\n")
            report_file.write(f"T3 (test) size: {test_size}\n")
            report_file.write(f"F1 (fine-tune) size: {finetune_size}\n")
            
            # Write label distribution
            report_file.write(f"\nLabel distribution in test set:\n")
            for label, count in test_label_counts.items():
                report_file.write(f"  Label {label}: {count} genes ({count/test_size*100:.1f}%)\n")
                
            report_file.write(f"\nLabel distribution in fine-tune set:\n")
            for label, count in finetune_label_counts.items():
                report_file.write(f"  Label {label}: {count} genes ({count/finetune_size*100:.1f}%)\n")
                
            report_file.write(f"\nOutput files:\n")
            report_file.write(f"- {base_name}_bio_test.csv, _label.txt, and _genes.csv\n")
            report_file.write(f"- {base_name}_bio_fine_tune.csv, _label.txt, and _genes.csv\n")
            
            print(f"Bio split complete:")
            print(f"T3 (test) size: {Fore.GREEN}{test_size}{Style.RESET_ALL}")
            print(f"F1 (fine-tune) size: {Fore.YELLOW}{finetune_size}{Style.RESET_ALL}")
            
            # Print label distribution
            print(f"\nLabel distribution in test set:")
            for label, count in test_label_counts.items():
                print(f"  Label {label}: {Fore.BLUE}{count}{Style.RESET_ALL} genes ({Fore.BLUE}{count/test_size*100:.1f}%{Style.RESET_ALL})")
                
            print(f"\nLabel distribution in fine-tune set:")
            for label, count in finetune_label_counts.items():
                print(f"  Label {label}: {Fore.BLUE}{count}{Style.RESET_ALL} genes ({Fore.BLUE}{count/finetune_size*100:.1f}%{Style.RESET_ALL})")
            
            print(f"Output files:")
            print(f"- {base_name}_bio_test.csv, _label.txt, and _genes.csv")
            print(f"- {base_name}_bio_fine_tune.csv, _label.txt, and _genes.csv")
    
    report_file.close()
    print(f"\n{Fore.GREEN}Report saved to: {os.path.join(output_dir, f'training_split_report_{timestamp}.txt')}{Style.RESET_ALL}")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Split P and Q datasets into training, test, and fine-tune sets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add arguments
    parser.add_argument(
        "-i", "--input-dir", 
        default="/Users/mukulsherekar/pythonProject/DEGCNN/split_datasets",
        help="Directory containing the P and Q files"
    )
    
    parser.add_argument(
        "-o", "--output-dir", 
        default="/Users/mukulsherekar/pythonProject/DEGCNN/training_datasets",
        help="Directory where the training files will be saved"
    )
    
    args = parser.parse_args()
    
    # Process datasets
    process_dataset(args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
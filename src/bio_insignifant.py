import os
import pandas as pd
import argparse
from datetime import datetime
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

def process_files(bio_validated_dir, filtered_quantile_dir, split_datasets_dir, report_dir=None):
    """
    Process files to find biologically insignificant genes and add them to destination files.
    
    Args:
        bio_validated_dir (str): Directory containing source files with biologically_validated column
        filtered_quantile_dir (str): Directory containing filtered quantile files with gene data
        split_datasets_dir (str): Directory containing destination P.csv files
        report_dir (str, optional): Directory for saving reports.
        
    Returns:
        str: Path to the generated report file
    """
    # Set up reporting
    if report_dir is None:
        report_dir = bio_validated_dir
    
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(report_dir, f"bio_insignificant_report_{timestamp}.txt")
    
    report_data = []
    genes_summary = {}
    
    # Get all CSV files in the bio_validated directory
    csv_files = [f for f in os.listdir(bio_validated_dir) if f.endswith('.csv')]
    
    if not csv_files:
        print(f"{Fore.RED}No CSV files found in {bio_validated_dir}")
        return None
    
    # Process each file
    for csv_file in csv_files:
        project_id = csv_file.split('_')[0]  # Extract PRJNA ID
        genes_summary[project_id] = {"yes": 0, "no": 0, "total": 0}
        
        # Define paths for all related files
        bioval_path = os.path.join(bio_validated_dir, csv_file)
        quantile_file = f"{project_id}_quantile.csv"
        quantile_path = os.path.join(filtered_quantile_dir, quantile_file)
        q_file = f"{project_id}_Q.csv"
        q_path = os.path.join(split_datasets_dir, q_file)
        p_file = f"{project_id}_P.csv"
        p_path = os.path.join(split_datasets_dir, p_file)
        
        # Check if all required files exist
        if not os.path.exists(quantile_path):
            message = f"Quantile file {quantile_file} not found in {filtered_quantile_dir}"
            report_data.append({"file": csv_file, "status": "ERROR", "message": message})
            print(f"{Fore.RED}{message}")
            continue
            
        if not os.path.exists(p_path):
            message = f"Destination P file {p_file} not found in {split_datasets_dir}"
            report_data.append({"file": csv_file, "status": "ERROR", "message": message})
            print(f"{Fore.RED}{message}")
            continue
            
        if not os.path.exists(q_path):
            message = f"Q file {q_file} not found in {split_datasets_dir}"
            report_data.append({"file": csv_file, "status": "ERROR", "message": message})
            print(f"{Fore.RED}{message}")
            continue
        
        try:
            # Step 1: Read biologically validated file and find genes marked as "no"
            bioval_df = pd.read_csv(bioval_path)
            
            # Check if biologically_validated column exists
            if 'biologically_validated' not in bioval_df.columns:
                message = f"Column 'biologically_validated' not found in {csv_file}"
                report_data.append({"file": csv_file, "status": "ERROR", "message": message})
                print(f"{Fore.RED}{message}")
                continue
            
            # Count genes with "yes" and "no" values
            genes_summary[project_id]["total"] = len(bioval_df)
            genes_summary[project_id]["yes"] = len(bioval_df[bioval_df['biologically_validated'] == 'yes'])
            genes_summary[project_id]["no"] = len(bioval_df[bioval_df['biologically_validated'] == 'no'])
            
            # Filter rows where biologically_validated is "no"
            filtered_bioval_df = bioval_df[bioval_df['biologically_validated'] == 'no']
            
            # Get the list of genes to process
            genes_to_find = filtered_bioval_df['Geneid'].tolist()
            
            if len(genes_to_find) == 0:
                message = f"No genes with biologically_validated='no' found in {csv_file}"
                report_data.append({"file": csv_file, "status": "INFO", "message": message})
                print(f"{Fore.YELLOW}{message}")
                continue
            
            print(f"{Fore.CYAN}Found {len(genes_to_find)} genes with biologically_validated='no' in {csv_file}")
            
            # Step 2: Find these genes in the quantile file
            # Read the quantile file
            quantile_df = pd.read_csv(quantile_path)
            
            # Filter quantile data to only include the genes we want
            filtered_quantile_df = quantile_df[quantile_df['Geneid'].isin(genes_to_find)]
            
            genes_found = filtered_quantile_df['Geneid'].nunique()
            
            if genes_found == 0:
                message = f"None of the {len(genes_to_find)} genes were found in {quantile_file}"
                report_data.append({"file": csv_file, "status": "WARNING", "message": message})
                print(f"{Fore.YELLOW}{message}")
                continue
            
            if genes_found < len(genes_to_find):
                missing_count = len(genes_to_find) - genes_found
                print(f"{Fore.YELLOW}Warning: {missing_count} genes were not found in {quantile_file}")
            
            # Step 3: Read the Q file to get labels for the genes
            q_df = pd.read_csv(q_path)
            
            # Create a mapping from gene ID to label
            gene_to_label = dict(zip(q_df['Geneid'], q_df['Label']))
            
            # Initialize label column if it doesn't exist in filtered_quantile_df
            if 'Label' not in filtered_quantile_df.columns:
                filtered_quantile_df['Label'] = 0  # Default value
            
            # Update labels based on Q file
            labels_updated = 0
            genes_with_labels = []
            
            # Loop through rows and update labels
            for index, row in filtered_quantile_df.iterrows():
                gene_id = row['Geneid']
                if gene_id in gene_to_label:
                    filtered_quantile_df.at[index, 'Label'] = gene_to_label[gene_id]
                    labels_updated += 1
                    genes_with_labels.append(gene_id)
                else:
                    # If gene not found in Q file, use label from bioval file if available
                    bioval_gene_label = filtered_bioval_df.loc[filtered_bioval_df['Geneid'] == gene_id, 'Label'].values
                    if len(bioval_gene_label) > 0:
                        filtered_quantile_df.at[index, 'Label'] = bioval_gene_label[0]
                        labels_updated += 1
                        genes_with_labels.append(gene_id)
            
            # Report on label updates
            label_match_rate = (labels_updated / genes_found) * 100 if genes_found > 0 else 0
            print(f"{Fore.CYAN}Updated labels for {labels_updated} out of {genes_found} genes ({label_match_rate:.1f}%)")
            
            # Ensure columns match destination file's expected format
            p_header_df = pd.read_csv(p_path, nrows=0)
            p_columns = p_header_df.columns.tolist()
            
            # Make sure we have all required columns
            for col in p_columns:
                if col not in filtered_quantile_df.columns:
                    message = f"Missing column '{col}' required for {p_file}"
                    report_data.append({"file": csv_file, "status": "ERROR", "message": message})
                    print(f"{Fore.RED}{message}")
                    continue
            
            # Reorder columns to match P file
            filtered_quantile_df = filtered_quantile_df[p_columns]
            
            # Step 4: Add these rows to the destination P.csv file
            # Append to destination file without headers
            filtered_quantile_df.to_csv(p_path, mode='a', header=False, index=False)
            
            message = f"Added {genes_found} genes from {quantile_file} to {p_file} with label data from {q_file}"
            report_data.append({"file": csv_file, "status": "SUCCESS", "message": message})
            print(f"{Fore.GREEN}{message}")
            
        except Exception as e:
            message = f"Error processing {csv_file}: {str(e)}"
            report_data.append({"file": csv_file, "status": "ERROR", "message": message})
            print(f"{Fore.RED}{message}")
    
    if report_data:
        generate_report(report_data, report_path, timestamp, genes_summary)
    
    return report_path

def generate_report(report_data, report_path, timestamp, genes_summary):
    """Generate a color-coded report of the processing results.
    
    Args:
        report_data (list): List of dictionaries with processing results
        report_path (str): Path to save the report
        timestamp (str): Timestamp for the report header
        genes_summary (dict): Summary of genes with biologically_validated="yes" and "no"
    """
    with open(report_path, 'w') as f:
        f.write(f"Bio Insignificant Data Processing Report - {timestamp}\n")
        f.write("="*50 + "\n\n")
        
        # Write gene count summary
        f.write("GENE VALIDATION SUMMARY\n")
        f.write("-"*50 + "\n")
        f.write(f"{'Project ID':<15} {'Total Genes':<15} {'Yes':<10} {'No':<10}\n")
        f.write("-"*50 + "\n")
        
        for project_id, counts in genes_summary.items():
            f.write(f"{project_id:<15} {counts['total']:<15} {counts['yes']:<10} {counts['no']:<10}\n")
        
        f.write("\n\n")
        f.write("PROCESSING RESULTS\n")
        f.write("-"*50 + "\n\n")
        
        for item in report_data:
            status_color = ""
            if item["status"] == "SUCCESS":
                status_color = "\033[92m"  # Green in ANSI
            elif item["status"] == "ERROR":
                status_color = "\033[91m"  # Red in ANSI
            elif item["status"] == "WARNING":
                status_color = "\033[93m"  # Yellow in ANSI
            elif item["status"] == "INFO":
                status_color = "\033[96m"  # Cyan in ANSI
            
            # Write to file (with color codes for terminal viewing)
            f.write(f"{status_color}[{item['status']}]\033[0m {item['file']}: {item['message']}\n")
            
    print(f"\n{Fore.CYAN}Report saved to: {report_path}")

def main():
    """Main function to parse arguments and run the processing"""
    parser = argparse.ArgumentParser(
        description='Process biologically insignificant data from CSV files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--bioval', '-b', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/bio_validated',
                        help='Directory containing biologically validated CSV files')
    
    parser.add_argument('--quantile', '-q', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/filtered_quantile',
                        help='Directory containing filtered quantile CSV files')
    
    parser.add_argument('--output', '-o', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/split_datasets',
                        help='Directory containing destination P.csv files')
    
    parser.add_argument('--report', '-r', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/bio_validated',
                        help='Directory to save reports')
    
    args = parser.parse_args()
    
    print(f"{Fore.CYAN}Starting processing of biologically insignificant data...")
    print(f"{Fore.CYAN}Bio Validated directory: {args.bioval}")
    print(f"{Fore.CYAN}Filtered Quantile directory: {args.quantile}")
    print(f"{Fore.CYAN}Output directory: {args.output}")
    print(f"{Fore.CYAN}Report directory: {args.report}")
    
    report_path = process_files(args.bioval, args.quantile, args.output, args.report)
    
    if report_path:
        print(f"{Fore.CYAN}Processing complete! Report saved at: {report_path}")
    else:
        print(f"{Fore.RED}Processing complete with errors. No report generated.")

if __name__ == "__main__":
    main() 
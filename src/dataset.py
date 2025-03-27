import pandas as pd
import os
from datetime import datetime
from colorama import init, Fore, Style

# Initialize colorama for colored output
init()

def merge_datasets(sdeg_file, quantile_file, output_file):
    """Merge sdeg and quantile files into a single dataset"""
    # Read the input files
    sdeg_df = pd.read_csv(sdeg_file)
    quantile_df = pd.read_csv(quantile_file)
    
    # Merge based on Gene_id
    merged_df = pd.merge(quantile_df, sdeg_df[['Geneid', 'Label']], on='Geneid', how='inner')
    
    # Save the merged dataset
    merged_df.to_csv(output_file, index=False)
    
    return len(merged_df), len(sdeg_df), len(quantile_df)

def main():
    # Directory containing the files
    sdeg_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/sdeg"
    quantile_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/filtered_quantile"
    
    # Create output directory if it doesn't exist
    output_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/datasets_alzheimers"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create report file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = open(os.path.join(output_dir, f"merge_report_{timestamp}.txt"), "w")
    report_file.write(f"Dataset Merge Report\n")
    report_file.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    report_file.write(f"{'='*50}\n\n")
    
    # Get all PRJNA identifiers from sdeg files
    prjna_ids = set()
    for filename in os.listdir(sdeg_dir):
        if filename.endswith('_sdeg.csv'):
            prjna_id = filename.split('_sdeg.csv')[0]
            prjna_ids.add(prjna_id)
    
    # Process each PRJNA dataset
    for prjna_id in prjna_ids:
        print(f"\n{Fore.CYAN}Processing {prjna_id}...{Style.RESET_ALL}")
        
        # Define file paths
        sdeg_file = os.path.join(sdeg_dir, f"{prjna_id}_sdeg.csv")
        quantile_file = os.path.join(quantile_dir, f"{prjna_id}_quantile.csv")
        output_file = os.path.join(output_dir, f"{prjna_id}_dataset.csv")
        
        # Check if both files exist
        if not os.path.exists(sdeg_file):
            print(f"{Fore.RED}Error: {sdeg_file} not found{Style.RESET_ALL}")
            continue
        if not os.path.exists(quantile_file):
            print(f"{Fore.RED}Error: {quantile_file} not found{Style.RESET_ALL}")
            continue
        
        # Merge the datasets
        merged_size, sdeg_size, quantile_size = merge_datasets(sdeg_file, quantile_file, output_file)
        
        # Write report
        report_file.write(f"\n{'='*50}\n")
        report_file.write(f"Report for {prjna_id}\n")
        report_file.write(f"{'='*50}\n\n")
        report_file.write(f"SDEG file: {sdeg_file} ({sdeg_size} rows)\n")
        report_file.write(f"Quantile file: {quantile_file} ({quantile_size} rows)\n")
        report_file.write(f"Merged file: {output_file} ({merged_size} rows)\n")
        report_file.write(f"Genes retained: {merged_size}/{sdeg_size} ({(merged_size/sdeg_size)*100:.2f}%)\n")
        report_file.write(f"\n{'='*50}\n")
        
        # Print report
        print(f"\n{Fore.GREEN}Merge complete:{Style.RESET_ALL}")
        print(f"SDEG file: {sdeg_file} ({Fore.YELLOW}{sdeg_size}{Style.RESET_ALL} rows)")
        print(f"Quantile file: {quantile_file} ({Fore.YELLOW}{quantile_size}{Style.RESET_ALL} rows)")
        print(f"Merged file: {output_file} ({Fore.GREEN}{merged_size}{Style.RESET_ALL} rows)")
        print(f"Genes retained: {Fore.GREEN}{merged_size}/{sdeg_size} ({(merged_size/sdeg_size)*100:.2f}%){Style.RESET_ALL}")
    
    report_file.close()
    print(f"\n{Fore.GREEN}Report saved to: {os.path.join(output_dir, f'merge_report_{timestamp}.txt')}{Style.RESET_ALL}")

if __name__ == "__main__":
    main()
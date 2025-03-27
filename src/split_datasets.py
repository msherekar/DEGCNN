import pandas as pd
import os
import sys
from datetime import datetime
from colorama import init, Fore, Style
import argparse

# Initialize colorama for colored output
init()

def get_statistics(df):
    total_rows = len(df)
    non_bio_rows = len(df[df['Label'] == 2])
    bio_rows = len(df[df['Label'].isin([0, 1])])
    up_regulated = len(df[df['Label'] == 0])
    down_regulated = len(df[df['Label'] == 1])
    
    return {
        'total_rows': total_rows,
        'non_bio_rows': non_bio_rows,
        'bio_rows': bio_rows,
        'up_regulated': up_regulated,
        'down_regulated': down_regulated
    }

def write_report(report_file, stats, filename):
    report_file.write(f"\n{'='*50}\n")
    report_file.write(f"Report for {filename}\n")
    report_file.write(f"{'='*50}\n\n")
    report_file.write(f"Total number of rows: {stats['total_rows']}\n")
    report_file.write(f"Biologically relevant rows: {stats['bio_rows']}\n")
    report_file.write(f"  - Up-regulated: {stats['up_regulated']}\n")
    report_file.write(f"  - Down-regulated: {stats['down_regulated']}\n")
    report_file.write(f"Non-biologically relevant rows: {stats['non_bio_rows']}\n")
    report_file.write(f"\nPercentage distribution:\n")
    report_file.write(f"  - Biologically relevant: {(stats['bio_rows']/stats['total_rows'])*100:.2f}%\n")
    report_file.write(f"  - Non-biologically relevant: {(stats['non_bio_rows']/stats['total_rows'])*100:.2f}%\n")
    report_file.write(f"  - Up-regulated: {(stats['up_regulated']/stats['total_rows'])*100:.2f}%\n")
    report_file.write(f"  - Down-regulated: {(stats['down_regulated']/stats['total_rows'])*100:.2f}%\n")
    report_file.write(f"\n{'='*50}\n")

def print_colored_report(stats, filename):
    print(f"\n{Fore.CYAN}{'='*50}")
    print(f"Report for {filename}")
    print(f"{'='*50}{Style.RESET_ALL}\n")
    print(f"Total number of rows: {Fore.GREEN}{stats['total_rows']}{Style.RESET_ALL}")
    print(f"Biologically relevant rows: {Fore.GREEN}{stats['bio_rows']}{Style.RESET_ALL}")
    print(f"  - Up-regulated: {Fore.RED}{stats['up_regulated']}{Style.RESET_ALL}")
    print(f"  - Down-regulated: {Fore.BLUE}{stats['down_regulated']}{Style.RESET_ALL}")
    print(f"Non-biologically relevant rows: {Fore.YELLOW}{stats['non_bio_rows']}{Style.RESET_ALL}")
    print(f"\nPercentage distribution:")
    print(f"  - Biologically relevant: {Fore.GREEN}{(stats['bio_rows']/stats['total_rows'])*100:.2f}%{Style.RESET_ALL}")
    print(f"  - Non-biologically relevant: {Fore.YELLOW}{(stats['non_bio_rows']/stats['total_rows'])*100:.2f}%{Style.RESET_ALL}")
    print(f"  - Up-regulated: {Fore.RED}{(stats['up_regulated']/stats['total_rows'])*100:.2f}%{Style.RESET_ALL}")
    print(f"  - Down-regulated: {Fore.BLUE}{(stats['down_regulated']/stats['total_rows'])*100:.2f}%{Style.RESET_ALL}")
    print(f"\n{Fore.CYAN}{'='*50}{Style.RESET_ALL}")

def split_dataset(input_file, output_dir, report_file):
    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Get statistics before splitting
    stats = get_statistics(df)
    
    # Split based on Label values
    non_bio_df = df[df['Label'] == 2]  # P (non-bio) data
    bio_df = df[df['Label'].isin([0, 1])]  # Q (bio) data
    
    # Create output filenames
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    # Remove '_dataset' from the base name if it exists
    base_name = base_name.replace('_dataset', '')
    non_bio_output = os.path.join(output_dir, f"{base_name}_P.csv")
    bio_output = os.path.join(output_dir, f"{base_name}_Q.csv")
    
    # Save the split datasets
    non_bio_df.to_csv(non_bio_output, index=False)
    bio_df.to_csv(bio_output, index=False)
    
    # Write and print reports
    write_report(report_file, stats, base_name)
    print_colored_report(stats, base_name)
    
    print(f"\nSplit files saved as:")
    print(f"- {Fore.GREEN}{non_bio_output}{Style.RESET_ALL} (non-bio data, Label=2)")
    print(f"- {Fore.GREEN}{bio_output}{Style.RESET_ALL} (bio data, Label=0,1)")

def main(input_dir=None, output_dir=None):
    # Directory containing the CSV files - use arguments if provided, otherwise use defaults
    if input_dir is None:
        input_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/datasets_alzheimers"
    
    # Create output directory if it doesn't exist
    if output_dir is None:
        output_dir = "/Users/mukulsherekar/pythonProject/DEGCNN/split_datasets"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Create report file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = open(os.path.join(output_dir, f"split_report_{timestamp}.txt"), "w")
    report_file.write(f"Dataset Split Report\n")
    report_file.write(f"Input directory: {input_dir}\n")
    report_file.write(f"Output directory: {output_dir}\n")
    report_file.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    report_file.write(f"{'='*50}\n\n")
    
    # Process each CSV file in the directory
    file_count = 0
    for filename in os.listdir(input_dir):
        if filename.endswith(".csv"):
            input_file = os.path.join(input_dir, filename)
            print(f"\n{Fore.CYAN}Processing {filename}...{Style.RESET_ALL}")
            split_dataset(input_file, output_dir, report_file)
            file_count += 1
    
    report_file.write(f"\nTotal files processed: {file_count}\n")
    report_file.close()
    print(f"\n{Fore.GREEN}Total files processed: {file_count}{Style.RESET_ALL}")
    print(f"{Fore.GREEN}Report saved to: {os.path.join(output_dir, f'split_report_{timestamp}.txt')}{Style.RESET_ALL}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Split CSV datasets into biologically relevant and non-biologically relevant files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add arguments
    parser.add_argument(
        "-i", "--input-dir", 
        default="/Users/mukulsherekar/pythonProject/DEGCNN/datasets_alzheimers",
        help="Directory containing the input CSV files"
    )
    
    parser.add_argument(
        "-o", "--output-dir", 
        default="/Users/mukulsherekar/pythonProject/DEGCNN/split_datasets",
        help="Directory where the split files will be saved"
    )
    
    # Parse arguments and call main function
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
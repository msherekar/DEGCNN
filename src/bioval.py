import os
import pandas as pd
import numpy as np
from gprofiler import GProfiler
import glob
import sys
import argparse
import logging

def setup_logger(log_level):
    """Setup logger with appropriate log level"""
    logging_level = getattr(logging, log_level.upper())
    logging.basicConfig(
        level=logging_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    return logging.getLogger('bioval')

def biological_validation(input_dir, output_dir, organism, logger):
    """
    Process each file in the input directory, extract genes with labels 0 or 1,
    use gprofiler to validate biological significance, and create output CSV files.
    
    Args:
        input_dir (str): Directory containing CSV files in SDEG folder
        output_dir (str): Directory to save output files
        organism (str): Organism to use for GProfiler (e.g., 'hsapiens')
        logger: Logger object for output messages
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    
    # Initialize GProfiler
    logger.info("Initializing GProfiler...")
    gp = GProfiler(return_dataframe=True)
    
    # Process each file in the input folder
    sdeg_files = glob.glob(os.path.join(input_dir, '*.csv'))
    
    logger.info(f"Looking for CSV files in: {input_dir}")
    logger.info(f"Found {len(sdeg_files)} files: {[os.path.basename(f) for f in sdeg_files]}")
    
    if not sdeg_files:
        logger.error(f"Error: No CSV files found in the input directory: {input_dir}")
        return False
    
    success_count = 0
    for file_path in sdeg_files:
        logger.info(f"\nProcessing {file_path}...")
        
        # Extract filename without extension
        base_filename = os.path.basename(file_path).split('.')[0]
        
        try:
            # Read the sdeg file
            logger.info(f"Reading file {file_path}...")
            df = pd.read_csv(file_path)
            logger.info(f"File loaded with {len(df)} rows.")
            
            # Check if Label column exists
            if 'Label' not in df.columns:
                logger.error(f"Error: 'Label' column not found in {file_path}. Available columns: {df.columns.tolist()}")
                continue
            
            # Filter genes with labels 0 or 1
            filtered_df = df[df['Label'].isin([0, 1])]
            logger.info(f"Filtered to {len(filtered_df)} genes with labels 0 or 1.")
            
            if filtered_df.empty:
                logger.warning(f"No genes with labels 0 or 1 found in {file_path}")
                continue
            
            # Get list of genes
            genes = filtered_df['Geneid'].tolist()
            logger.info(f"Running GProfiler on {len(genes)} genes...")
            
            # Run g:Profiler on the gene list using specified sources
            try:
                result = gp.profile(
                    organism=organism,
                    query=genes,
                    sources=["KEGG", "REAC", "WP", "GO:BP"]
                )
                logger.info(f"GProfiler analysis complete. Found {len(result) if not result.empty else 0} results.")
                
                # Filter for Alzheimer-related pathways using the term name
                alzheimers_results = result[result['name'].str.contains("alzheimer", case=False, na=False)]
                logger.info(f"Found {len(alzheimers_results)} Alzheimer-related pathways.")
                
                if alzheimers_results.empty:
                    filtered_df['biologically_validated'] = 'no'
                    logger.info("No Alzheimer-specific pathways found. Marking all genes as not validated.")
                else:
                    # Extract validated genes from the 'intersections' column
                    validated_genes = set()
                    for _, row in alzheimers_results.iterrows():
                        if 'intersections' in row and isinstance(row['intersections'], str):
                            validated_genes.update(row['intersections'].split(','))
                    
                    logger.info(f"Found {len(validated_genes)} biologically validated genes based on Alzheimer-specific pathways.")
                    
                    # Mark genes as validated or not
                    filtered_df['biologically_validated'] = filtered_df['Geneid'].apply(
                        lambda x: 'yes' if x in validated_genes else 'no'
                    )
                    
            except Exception as e:
                logger.error(f"Error in gprofiler analysis for {file_path}: {str(e)}")
                # If error occurs, mark all as not validated
                filtered_df['biologically_validated'] = 'no'
            
            # Save the result to output file
            output_path = os.path.join(output_dir, f"{base_filename}_bioval.csv")
            filtered_df.to_csv(output_path, index=False)
            logger.info(f"Results saved to {output_path}")
            success_count += 1
            
        except Exception as e:
            logger.error(f"Error processing file {file_path}: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
    
    logger.info(f"Successfully processed {success_count} out of {len(sdeg_files)} files.")
    return success_count > 0

def main():
    """Main function with command-line argument parsing"""
    parser = argparse.ArgumentParser(
        description='Perform biological validation on SDEG files using GProfiler',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/sdeg',
                        help='Directory containing SDEG CSV files')
    
    parser.add_argument('-o', '--output', 
                        default='/Users/mukulsherekar/pythonProject/DEGCNN/bio_validated',
                        help='Directory to save output files')
    
    parser.add_argument('-g', '--organism', 
                        default='hsapiens',
                        help='Organism to use for GProfiler analysis')
    
    parser.add_argument('-l', '--log-level',
                        choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='info',
                        help='Set the logging level')
    
    args = parser.parse_args()
    
    # Setup logger
    logger = setup_logger(args.log_level)
    
    # Convert relative paths to absolute paths
    current_dir = os.getcwd()
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    
    logger.info("Starting biological validation script...")
    logger.info(f"Current working directory: {current_dir}")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Organism: {args.organism}")
    
    success = biological_validation(input_dir, output_dir, args.organism, logger)
    
    if success:
        logger.info("Script execution completed successfully.")
    else:
        logger.error("Script execution failed.")

if __name__ == "__main__":
    main()

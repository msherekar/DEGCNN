# Pipeline for data engineering for DEGCNN
# qc->trimming->qc->allignment->featurecounts->counts.py(to combine)->preprocessing.py->hngc.py->normalize.py->filter.py(by quantile)->sdeg->datasplit->degcnn.py

# Combine feature counts from all SRR files
from counts import save_merged_counts

# Preprocess the combined feature counts
from preprocessing import preprocess_all_files

# Map HGNC symbols to gene IDs
from hngc import hngc_mapping   

# Normalize the counts
from normalize import normalize_all_files

# Filter the counts by quantile
from filter import filter_all_files

# Perform DEA analysis and label the data -0,1,2
from sdeg import differential_expression_analysis

import os
import logging
from pathlib import Path
import argparse
from typing import Dict, Optional, Tuple

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('data_engineering.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def setup_directories(base_dir: str) -> Dict[str, str]:
    """
    Create and return a dictionary of all required directories for the pipeline.
    This structure is designed to be compatible with Snakemake workflows.
    
    Args:
        base_dir (str): Path to feature_counts directory or a project-specific subdirectory
        
    Returns:
        Dict[str, str]: Dictionary containing paths to all required directories
    """
    # Check if base_dir is the feature_counts directory or a project-specific directory
    base_dir_name = os.path.basename(base_dir)
    
    if base_dir_name == "feature_counts":
        # If the base_dir is the feature_counts directory itself
        project_root = os.path.dirname(base_dir)
        feature_counts_dir = base_dir
    else:
        # If the base_dir is a project-specific directory within feature_counts
        project_root = os.path.dirname(os.path.dirname(base_dir))
        feature_counts_dir = os.path.dirname(base_dir)
    
    # Define common output directories at the project root level
    directories = {
        'feature_counts': feature_counts_dir,
        'counts': os.path.join(project_root, 'counts'),
        'preprocessed': os.path.join(project_root, 'preprocessed'),
        'hgnc_mapped': os.path.join(project_root, 'hgnc_mapped'),
        'normalized': os.path.join(project_root, 'normalized'),
        'filtered_counts': os.path.join(project_root, 'filtered_counts'),
        'filtered_quantile': os.path.join(project_root, 'filtered_quantile'),
        'sdeg': os.path.join(project_root, 'sdeg'),
        'logs': os.path.join(project_root, 'logs')
    }
    
    # Create directories if they don't exist
    for dir_path in directories.values():
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        logger.info(f"Created/verified directory: {dir_path}")
    
    return directories

def data_engineering(base_dir, metadata_dir, gc_content_file, skip_existing=False):
    """
    Run the complete data engineering pipeline.
    
    Args:
        base_dir (str): Base directory containing project folders
        metadata_dir (str): Directory containing metadata files
        gc_content_file (str): Path to GC content file
        skip_existing (bool): Whether to skip steps that have already been completed
    """
    # Set up directories
    dirs = setup_directories(base_dir)
    
    # Set up logging
    log_file = os.path.join(dirs['logs'], 'pipeline.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    # Step 1: Merge feature counts
    logging.info("Step 1: Merging feature counts...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['counts'], f"{os.path.basename(base_dir)}_counts.csv")):
        if not save_merged_counts(base_dir, dirs['counts']):
            logging.error("Failed to merge feature counts")
            return False
    else:
        logging.info("Skipping feature count merging - output already exists")
    
    # Step 2: Preprocess counts
    logging.info("Step 2: Preprocessing counts...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['filtered_counts'], f"{os.path.basename(base_dir)}_filtered_counts.csv")):
        if not preprocess_all_files(dirs['counts'], dirs['filtered_counts']):
            logging.error("Failed to preprocess counts")
            return False
    else:
        logging.info("Skipping preprocessing - output already exists")
    
    # Step 3: Map to HGNC
    logging.info("Step 3: Mapping to HGNC...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['hgnc_mapped'], f"{os.path.basename(base_dir)}_hgnc_mapped.csv")):
        if not hngc_mapping(dirs['filtered_counts'], dirs['hgnc_mapped']):
            logging.error("Failed to map to HGNC")
            return False
    else:
        logging.info("Skipping HGNC mapping - output already exists")
    
    # Step 4: Normalize counts
    logging.info("Step 4: Normalizing counts...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['normalized'], f"{os.path.basename(base_dir)}_normalized.csv")):
        if not normalize_all_files(dirs['hgnc_mapped'], dirs['normalized'], gc_content_file):
            logging.error("Failed to normalize counts")
            return False
    else:
        logging.info("Skipping normalization - output already exists")
    
    # Step 5: Filter counts
    logging.info("Step 5: Filtering counts...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['filtered_quantile'], f"{os.path.basename(base_dir)}_filtered_counts_quantile.csv")):
        if not filter_all_files(dirs['normalized'], dirs['filtered_quantile']):
            logging.error("Failed to filter counts")
            return False
    else:
        logging.info("Skipping filtering - output already exists")
    
    # Step 6: Differential expression analysis
    logging.info("Step 6: Performing differential expression analysis...")
    if not skip_existing or not os.path.exists(os.path.join(dirs['sdeg'], f"{os.path.basename(base_dir)}_sdeg.csv")):
        if not differential_expression_analysis(dirs['filtered_quantile'], metadata_dir, dirs['sdeg']):
            logging.error("Failed to perform differential expression analysis")
            return False
    else:
        logging.info("Skipping differential expression analysis - output already exists")
    
    logging.info("Pipeline completed successfully!")
    return True

def main():
    """Command-line interface for the data engineering pipeline."""
    parser = argparse.ArgumentParser(description='Data Engineering Pipeline for DEGCNN')
    parser.add_argument('-b', '--base_dir', required=True, help='Base directory for all data processing')
    parser.add_argument('-m', '--metadata_dir', required=True, help='Directory containing metadata files')
    parser.add_argument('-g', '--gc_content_file', required=True, help='Path to GC content file')
    parser.add_argument('--skip_existing', action='store_true', help='Skip processing if output files exist')
    
    args = parser.parse_args()
    
    success = data_engineering(
        base_dir=args.base_dir,
        metadata_dir=args.metadata_dir,
        gc_content_file=args.gc_content_file,
        skip_existing=args.skip_existing
    )
    
    if success:
        logger.info("Data engineering completed successfully!")
    else:
        logger.error("Data engineering failed.")
        exit(1)

if __name__ == "__main__":
    main()


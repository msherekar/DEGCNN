# Reproduction of DEGNext
# RRNA-seq Data Engineering and Analysis Pipeline using CNN

This directory contains a comprehensive pipeline for processing RNA-seq data, performing differential expression analysis, and preparing datasets for machine learning applications.

## 📁 File Structure

### 🔧 Core Data Processing

#### `counts.py`
**Purpose**: Merges feature count files from multiple samples into a single count matrix.
- **Input**: Individual `*_counts.tsv` files from featureCounts
- **Output**: Combined `*_counts.csv` files per project
- **Usage**: `python counts.py -i /path/to/feature_counts -o /path/to/output`

#### `preprocessing.py`
**Purpose**: Filters raw count data to remove low-expression genes.
- **Input**: Raw count matrices
- **Output**: Filtered count matrices
- **Usage**: `python preprocessing.py -i counts/ -o filtered_counts/`

#### `hngc.py`
**Purpose**: Maps Ensembl gene IDs to HGNC gene symbols for standardization.
- **Input**: Count files with Ensembl IDs
- **Output**: Count files with HGNC symbols
- **Usage**: `python hngc.py -i filtered_counts/ -o hgnc_mapped/`

#### `normalize.py`
**Purpose**: Performs within-lane and between-lane normalization of count data.
- **Input**: Filtered count matrices
- **Output**: Normalized expression matrices
- **Usage**: `python normalize.py -i hgnc_mapped/ -o normalized/ -g gc_content.csv`

#### `filter.py`
**Purpose**: Filters genes based on expression quantiles to retain highly expressed genes.
- **Input**: Normalized expression matrices
- **Output**: Quantile-filtered expression matrices
- **Usage**: `python filter.py -i normalized/ -o filtered_quantile/ -q 0.25`

### 🔬 Differential Expression Analysis

#### `sdeg.py`
**Purpose**: Performs differential expression analysis using TCGAanalyze_DEA (edgeR-based).
- **Input**: Filtered count matrices + metadata
- **Output**: DEA results with gene labels (0, 1, 2)
- **Usage**: `python sdeg.py -i filtered_quantile/ -m meta_data/ -o sdeg/`

#### `sdeg_GEO.py`
**Purpose**: GEO-specific version of differential expression analysis.
- **Input**: GEO count matrices + metadata
- **Output**: DEA results for GEO datasets
- **Usage**: `python sdeg_GEO.py -i filtered_quantile/ -m meta_data/ -o sdeg_GEO/`

#### `deg.py`
**Purpose**: Alternative DEA implementation using limma and TCGAanalyze_DEA.
- **Input**: Count matrices + metadata
- **Output**: DEA results
- **Usage**: `python deg.py` (configured for specific directories)

#### `pydea.py`
**Purpose**: Python-based differential expression analysis implementation.
- **Input**: Expression matrices + metadata
- **Output**: DEA results with statistical significance
- **Usage**: `python pydea.py` (configured for specific directories)

### 🧬 Biological Validation

#### `bioval.py`
**Purpose**: Performs biological validation using GProfiler to identify Alzheimer's-related genes.
- **Input**: SDEG files with labeled genes
- **Output**: Biologically validated gene lists
- **Usage**: `python bioval.py -i datasets_alzheimers/ -o bio_validated/ -g hsapiens`

#### `bio_insignifant.py`
**Purpose**: Processes biologically insignificant genes and adds them to training datasets.
- **Input**: Bio-validated files + quantile files
- **Output**: Enhanced training datasets
- **Usage**: `python bio_insignifant.py -b bio_validated/ -q filtered_quantile/ -o split_datasets/`

#### `enrichr.py`
**Purpose**: Performs gene set enrichment analysis using Enrichr API.
- **Input**: Gene lists from SDEG analysis
- **Output**: Enrichment results and validated genes
- **Usage**: `python enrichr.py` (configured for specific files)

#### `significance.py`
**Purpose**: Comprehensive biological significance analysis using GProfiler.
- **Input**: SDEG files
- **Output**: Significance analysis results and visualizations
- **Usage**: `python significance.py` (configured for specific directories)

### 📊 Dataset Management

#### `dataset.py`
**Purpose**: Merges SDEG results with quantile-filtered data into complete datasets.
- **Input**: SDEG files + quantile files
- **Output**: Complete datasets with labels
- **Usage**: `python dataset.py` (configured for specific directories)

#### `split_datasets.py`
**Purpose**: Splits datasets into P and Q subsets based on biological validation.
- **Input**: Complete datasets
- **Output**: P and Q dataset files
- **Usage**: `python split_datasets.py -i datasets_alzheimers_GEO/ -o split_datasets_GEO/`

#### `split.py`
**Purpose**: Splits datasets and generates label distribution visualizations.
- **Input**: Dataset files with Dataset column
- **Output**: P/Q splits + distribution plots
- **Usage**: `python split.py` (configured for specific directories)

#### `combined_datasetprocessing.py`
**Purpose**: Comprehensive dataset processing pipeline for training data preparation.
- **Input**: P and Q dataset files
- **Output**: Training, test, and fine-tune datasets
- **Usage**: `python combined_datasetprocessing.py -i split_datasets/ -o training_datasets_GEO/`

#### `training_datasets.py`
**Purpose**: Creates training datasets with proper formatting for ML models.
- **Input**: Split datasets
- **Output**: Formatted training datasets
- **Usage**: `python training_datasets.py -i split_datasets/ -o training_datasets/`

### 🔄 Data Engineering Pipeline

#### `data_engineering.py`
**Purpose**: Orchestrates the complete data engineering pipeline from raw counts to training datasets.
- **Input**: Raw count files + metadata + GC content
- **Output**: Complete processed datasets
- **Usage**: `python data_engineering.py -b feature_counts/ -m meta_data/ -g gc_content.csv`

### 📥 Data Acquisition

#### `download.py`
**Purpose**: Downloads RNA-seq data from SRA using prefetch and fasterq-dump.
- **Input**: Project IDs (PRJNA*)
- **Output**: FASTQ files
- **Usage**: `python download.py` (configured for specific project IDs)

#### `tar.py`
**Purpose**: Processes compressed MTX files from GEO datasets.
- **Input**: MTX format files
- **Output**: Count matrices
- **Usage**: `python tar.py` (configured for specific directories)

### 🔧 Utility Scripts

#### `map.py`
**Purpose**: Maps Ensembl IDs to HGNC gene symbols using pre-built mapping.
- **Input**: GC content file with Ensembl IDs
- **Output**: GC content file with HGNC symbols
- **Usage**: `python map.py` (configured for specific files)

#### `mygenescript.py`
**Purpose**: Converts NCBI gene IDs to Ensembl IDs using MyGene.info API.
- **Input**: Count files with NCBI IDs
- **Output**: Count files with Ensembl IDs
- **Usage**: `python mygenescript.py` (configured for specific files)

#### `validation.py`
**Purpose**: Extracts biologically validated data for transfer learning.
- **Input**: Q dataset files
- **Output**: Biologically validated training files
- **Usage**: `python validation.py` (configured for specific directories)

## 🚀 Quick Start

### 1. Data Download
```bash
python download.py  # Configure project IDs in the script
```

### 2. Run Complete Pipeline
```bash
python data_engineering.py -b feature_counts/ -m meta_data/ -g gc_content.csv
```

### 3. Biological Validation
```bash
python bioval.py -i datasets_alzheimers/ -o bio_validated/ -g hsapiens
```

### 4. Create Training Datasets
```bash
python combined_datasetprocessing.py -i split_datasets/ -o training_datasets_GEO/
```

## 📋 Dependencies

### Core Dependencies
- `pandas` - Data manipulation
- `numpy` - Numerical operations
- `matplotlib` - Plotting
- `scikit-learn` - Machine learning utilities

### Bioinformatics Dependencies
- `rpy2` - R integration for DEA
- `gprofiler` - Biological pathway analysis
- `mygene` - Gene ID conversion

### R Dependencies (via rpy2)
- `TCGAbiolinks`
- `limma`
- `edgeR`

## 📊 Output Structure

```
project_root/
├── counts/                    # Merged count matrices
├── filtered_counts/           # Preprocessed counts
├── hgnc_mapped/              # HGNC-mapped counts
├── normalized/                # Normalized expression
├── filtered_quantile/         # Quantile-filtered data
├── sdeg/                      # Differential expression results
├── datasets_alzheimers/       # Complete datasets
├── split_datasets/            # P/Q dataset splits
├── training_datasets/         # ML-ready training data
└── bio_validated/             # Biologically validated genes
```

## 🔧 Configuration

Most scripts use relative paths and can be configured by modifying the default paths in each script. Key configuration files:

- `gc_content.csv` - GC content data for normalization
- `ensembl_to_hgnc_mapping.json` - Gene ID mapping
- Metadata files in `meta_data/` directory

## 📈 Pipeline Flow

1. **Data Acquisition**: Download FASTQ files from SRA
2. **Quality Control**: FastQC and trimming
3. **Alignment**: Map reads to reference genome
4. **Feature Counting**: Generate count matrices
5. **Data Engineering**: Preprocess and normalize
6. **Differential Expression**: Identify DEGs
7. **Biological Validation**: Pathway enrichment
8. **Dataset Preparation**: Create ML-ready datasets



## Citation
Kakati, T., Bhattacharyya, D.K., Kalita, J.K., and Norden-Krichmar, T.M. (2022) “DEGNext: Classification of Differentially Expressed Genes from RNA-seq data using a Convolutional Neural Network with Transfer Learning”, BMC Bioinformatics, 23:17, doi:10.1186/s12859-021-04527-4.
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04527-4

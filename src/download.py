import os
import subprocess
import requests
import time
import defusedxml.ElementTree as ET  # Secure XML parsing

# Define base directories
BASE_DIR = "/Users/mukulsherekar/pythonProject/DEGCNN"
SRA_DOWNLOAD_DIR = os.path.join(BASE_DIR, "downloads")    # Where .sra files will be stored
FASTQ_OUTPUT_DIR = os.path.join(BASE_DIR, "fastq_files")  # Where FASTQ files will be saved

# Ensure directories exist
os.makedirs(SRA_DOWNLOAD_DIR, exist_ok=True)
os.makedirs(FASTQ_OUTPUT_DIR, exist_ok=True)

def get_sra_run_ids(project_id):
    """
    Fetch numeric SRA IDs given a GEO Project ID (PRJNAxxxxxx).
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={project_id}&retmode=xml"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching SRA Run IDs for {project_id}: {e}")
        return []

    try:
        root = ET.fromstring(response.content)
        numeric_ids = [id_elem.text for id_elem in root.findall(".//Id")]
        if not numeric_ids:
            print(f"No SRA Run IDs found for {project_id}")
            return []
        return numeric_ids
    except ET.ParseError:
        print(f"Error parsing XML for {project_id}")
        return []

def convert_numeric_to_srr(numeric_id):
    """
    Convert an internal numeric SRA ID to a public SRR accession.
    Includes delay to prevent HTTP 429 errors (Too Many Requests).
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={numeric_id}&retmode=xml"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching efetch for {numeric_id}: {e}")
        return None

    try:
        root = ET.fromstring(response.content)
        run_elem = root.find(".//RUN")
        if run_elem is not None:
            srr_accession = run_elem.get("accession")
            print(f"Converted numeric ID {numeric_id} to accession {srr_accession}")
            return srr_accession
        else:
            print(f"No RUN element found for numeric ID {numeric_id}")
            return None
    except ET.ParseError as e:
        print(f"Error parsing efetch XML for numeric ID {numeric_id}: {e}")
        return None
    finally:
        time.sleep(1)  # Prevent hitting API rate limits

def download_fastq(srr_id, output_dir=FASTQ_OUTPUT_DIR):
    """
    Download and convert SRA files directly to compressed FASTQ, handling paired-end reads correctly.
    """
    print(f"Downloading {srr_id} using prefetch...")

    # ✅ Prefetch command (unchanged)
    prefetch_command = f"prefetch {srr_id} --output-directory {SRA_DOWNLOAD_DIR} --max-size 100G"
    try:
        subprocess.run(prefetch_command, shell=True, check=True)
        print(f"Prefetch completed for {srr_id}")

        # ✅ Locate the downloaded .sra file
        sra_file_path = os.path.join(SRA_DOWNLOAD_DIR, srr_id, f"{srr_id}.sra")
        if not os.path.exists(sra_file_path):
            print(f"Error: Expected SRA file not found at {sra_file_path}")
            return
        print(f"Using SRA file: {sra_file_path}")

        # ✅ Convert SRA directly to compressed FASTQ files
        fastq_command = f"""
        fasterq-dump --split-files --threads 4 -O {output_dir} {sra_file_path} && \
        pigz -p 4 {output_dir}/{srr_id}_*.fastq
        """

        subprocess.run(fastq_command, shell=True, check=True)
        print(f"Finished downloading {srr_id} directly as FASTQ.GZ")

    except subprocess.CalledProcessError as e:
        print(f"Error processing {srr_id}: {e}")


def process_project(project_id):
    """
    Process a single GEO Project ID.
    Fetch numeric SRA IDs, convert them to SRR accessions, and download FASTQ files.
    """
    print(f"Fetching SRA Run IDs for Project: {project_id}")
    numeric_ids = get_sra_run_ids(project_id)
    if not numeric_ids:
        print(f"No SRA data found for project {project_id}")
        return

    for num_id in numeric_ids:
        srr_id = convert_numeric_to_srr(num_id)
        if srr_id:
            download_fastq(srr_id)
        else:
            print(f"Skipping numeric ID {num_id} due to conversion error.")

def process_multiple_projects(project_ids):
    """
    Process multiple GEO Project IDs sequentially.
    """
    for project_id in project_ids:
        process_project(project_id)
        time.sleep(2)  # Prevent overloading NCBI

# Example Usage
if __name__ == "__main__":
    project_ids = ["PRJNA683625"]
    process_multiple_projects(project_ids)

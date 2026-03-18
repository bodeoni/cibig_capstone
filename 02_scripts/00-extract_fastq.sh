#!/bin/bash
#SBATCH --job-name=fasterq_dump      
#SBATCH --output=99_logs/00_fasterq_dump_%j.log                        
#SBATCH --cpus-per-task=16            
#SBATCH --mem=64G                    
#SBATCH --partition=normal
#SBATCH --nodelist=node06

# Paths
DATA_DIR="/scratch/onilee/capstone"
TMP_DIR="/lscratch/onilee/tmp_space"
INPUT="$DATA_DIR/01_data/00_sra"
OUTPUT="$DATA_DIR/01_data/01_raw_fastq"

# move into data directory
cd "$DATA_DIR"
# Create temp directory for fasterq-dump
mkdir -p "$TMP_DIR"

echo "Job started at: $(date)"

for dir in $INPUT/SRR*/; do
    srr_id=$(basename "$dir") # Extract the SRR ID from the directory path
    outdir="$OUTPUT/$srr_id" # Output directory for FASTQ files
    mkdir -p "$outdir" # Create output directory if it doesn't exist
    
    echo "--------------------------------------------"
    echo "Processing $srr_id"
    echo "--------------------------------------------"
    
    # 1. Extraction: Uses 16 threads to generate raw FASTQ
    fasterq-dump --split-3 "$srr_id" \
                 --threads 16 \
                 --outdir "$outdir" \
                 --progress \
                 --temp "$TMP_DIR"

    # 2. Background Compression: Gzip the FASTQ files in the background
    # We use 'nice' to ensure compression doesn't starve the next extraction for I/O
    #echo "Starting compression for $srr_id..."
    #gzip -f "$outdir/"*.fastq
done

echo "Job finished at: $(date)"
rm -rf "$TMP_DIR" # Clean up temp directory
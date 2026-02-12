#!/bin/bash
#SBATCH --job-name=QC      
#SBATCH --output=99_logs/01_qc_%j.log                        
#SBATCH --cpus-per-task=16                              
#SBATCH --partition=normal
#SBATCH --nodelist=node06

# Paths
## Input data directory
DATA_DIR="/data/onilee/capstone"

## Output directories
fastqc_dir="$DATA_DIR/03_analysis/01_qc/a_fastqc"
multiqc_dir="$DATA_DIR/03_analysis/01_qc/b_multiqc"

echo "Job started at: $(date)"

echo "----------------------------------------------"
echo "-----------------FASTQC-----------------------"
echo "----------------------------------------------"

# Load FastQC module
module load FastQC/0.12.1

# Record FastQC version
fastqc -v >> "$DATA_DIR/00_meta/software_versions.txt"

for dir in "$DATA_DIR"/01_data/01_raw_fastq/SRR*/; do
    srr_id=$(basename "$dir") # Extract SRR ID from directory name
    
    # make output directory for FastQC reports
    outdir="$fastqc_dir/$srr_id"
    mkdir -p "$outdir"

    echo "--------------------------------------------"
    echo "Running FastQC for $srr_id"
    echo "--------------------------------------------"
    
    # Run FastQC on all FASTQ files for the current SRR ID
    fastqc -t $SLURM_CPUS_PER_TASK \
           -o "$outdir" \
           "$dir"/"$srr_id"_1.fastq "$dir"/"$srr_id"_2.fastq
done

echo "----------------------------------------------"
echo "-----------------MULTIQC----------------------"
echo "----------------------------------------------"

# Load MultiQC module
module load MultiQC/1.9

# Record MultiQC version
multiqc --version >> "$DATA_DIR/00_meta/software_versions.txt"

mkdir -p "$multiqc_dir"

multiqc -o "$multiqc_dir" \
 -n "multiqc_raw_fastqc.html" "$fastqc_dir"

echo "Job ended at: $(date)"

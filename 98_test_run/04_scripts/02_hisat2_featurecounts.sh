#!/bin/bash
#SBATCH --job-name=hisat2_fc
#SBATCH --output=/projects/Whitefly_RNASeq/98_test_run/99_logs/approach2_hisat2_%j.log
#SBATCH --error=/projects/Whitefly_RNASeq/98_test_run/99_logs/approach2_hisat2_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

# ================================================================================
# APPROACH 2: HISAT2 Alignment + featureCounts
# ================================================================================
# Input: Trimmed clean reads
# Output: BAM alignments + gene counts from featureCounts
# ================================================================================

set -euo pipefail

# --- INITIALIZE WAVE MODULE ENVIRONMENT ---
source /usr/local/bioinfo/modules-5.4.0/init/bash
module use /usr/local/bioinfo/modulefiles/environment/

# Module loads
module load bioinfo-shared apptainer

# ================================================================================
# CONFIGURATION
# ================================================================================

# Directories
PROJECT_DIR="/projects/Whitefly_RNASeq/98_test_run"
INPUT_DIR="${PROJECT_DIR}/01_data/01_trimmed"
OUTPUT_DIR="${PROJECT_DIR}/02_analysis/approach2_hisat2"
REF_DIR="${PROJECT_DIR}/98_refs"
LOG_DIR="${PROJECT_DIR}/99_logs"

# Reference files
GENOME_FASTA="${REF_DIR}/MEAM1_scaffold_v1.2.fa.gz"
GTF_FILE="${REF_DIR}/MEAM1_v1.2.gff3.gz"
HISAT2_INDEX="${REF_DIR}/hisat2_index/MEAM1"  # Prefix for HISAT2 index

# Containers
HISAT2_CONTAINER="/projects/onilee/software/containers/hisat2_2.2.2.sif"  
SUBREAD_CONTAINER="/projects/onilee/software/containers/subread_2.1.1.sif" 
SAMTOOLS_CONTAINER="/projects/onilee/software/containers/samtools_1.23.sif"

# Resources
THREADS=32

# Sample IDs
SAMPLES=("SRR28578498" "SRR28578499" "SRR28578500")

# ================================================================================
# CREATE OUTPUT DIRECTORIES
# ================================================================================

mkdir -p "${OUTPUT_DIR}"/{bam,counts,logs}
mkdir -p "${REF_DIR}/hisat2_index"

echo "=========================================="
echo "APPROACH 2: HISAT2 + featureCounts"
echo "=========================================="
echo "Started: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# ================================================================================
# STEP 1: BUILD HISAT2 INDEX (if not exists)
# ================================================================================

if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "[$(date)] Building HISAT2 index..."
    
    # Uncompress reference if needed
    gunzip -c "${GENOME_FASTA}" > "${REF_DIR}/MEAM1_scaffold_v1.2.fa"
    
    apptainer exec "${HISAT2_CONTAINER}" hisat2-build \
        -p "${THREADS}" \
        "${REF_DIR}/MEAM1_scaffold_v1.2.fa" \
        "${HISAT2_INDEX}"
    
    echo "[$(date)] HISAT2 index built successfully"
else
    echo "[$(date)] HISAT2 index already exists, skipping..."
fi

# ================================================================================
# STEP 2: ALIGN WITH HISAT2
# ================================================================================

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================="
    
    # Input files
    READ1="${INPUT_DIR}/${SAMPLE}_1.clean.fastq.gz"
    READ2="${INPUT_DIR}/${SAMPLE}_2.clean.fastq.gz"
    
    # Output files
    SAM_FILE="${OUTPUT_DIR}/bam/${SAMPLE}.sam"
    BAM_FILE="${OUTPUT_DIR}/bam/${SAMPLE}.bam"
    SORTED_BAM="${OUTPUT_DIR}/bam/${SAMPLE}.sorted.bam"
    
    # Check input files exist
    if [[ ! -f "${READ1}" ]] || [[ ! -f "${READ2}" ]]; then
        echo "ERROR: Input files not found for ${SAMPLE}"
        continue
    fi
    
    echo "[$(date)] Aligning ${SAMPLE} with HISAT2..."
    
    apptainer exec "${HISAT2_CONTAINER}" hisat2 \
        -x "${HISAT2_INDEX}" \
        -1 "${READ1}" \
        -2 "${READ2}" \
        -S "${SAM_FILE}" \
        -p "${THREADS}" \
        --dta \
        --summary-file "${OUTPUT_DIR}/logs/${SAMPLE}_hisat2_summary.txt" \
        2> "${OUTPUT_DIR}/logs/${SAMPLE}_hisat2.log"
    
    echo "[$(date)] Converting to BAM and sorting..."
    
    # Convert SAM to BAM (use system samtools or container)
    apptainer exec "${SAMTOOLS_CONTAINER}" samtools view -@ "${THREADS}" -bS "${SAM_FILE}" > "${BAM_FILE}"
    
    # Sort BAM
    apptainer exec "${SAMTOOLS_CONTAINER}" samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" "${BAM_FILE}"
    
    # Index BAM
    apptainer exec "${SAMTOOLS_CONTAINER}" samtools index "${SORTED_BAM}"
    
    # Remove intermediate files
    rm "${SAM_FILE}" "${BAM_FILE}"
    
    echo "[$(date)] ${SAMPLE} alignment completed"
    echo "  BAM: ${SORTED_BAM}"
done

# ================================================================================
# STEP 3: QUANTIFY WITH FEATURECOUNTS
# ================================================================================

echo ""
echo "=========================================="
echo "Running featureCounts"
echo "=========================================="

# Uncompress GTF/GFF if needed
gunzip -c "${GTF_FILE}" > "${REF_DIR}/MEAM1_v1.2.gff3"

# Prepare list of BAM files
BAM_FILES=""
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES="${BAM_FILES} ${OUTPUT_DIR}/bam/${SAMPLE}.sorted.bam"
done

echo "[$(date)] Counting features..."

apptainer exec "${SUBREAD_CONTAINER}" featureCounts \
    -T "${THREADS}" \
    -p \
    -t mRNA \
    -g Parent \
    -a "${REF_DIR}/MEAM1_v1.2.gff3" \
    -o "${OUTPUT_DIR}/counts/all_samples_counts.txt" \
    ${BAM_FILES}

# Create individual count files for consistency
echo "[$(date)] Extracting individual sample counts..."

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    COL=$((i + 7))  # featureCounts output: cols 1-6 are annotation, 7+ are samples
    
    awk -v col="$COL" 'NR==1 {print} NR>1 {print $1"\t"$col}' \
        "${OUTPUT_DIR}/counts/all_samples_counts.txt" \
        > "${OUTPUT_DIR}/counts/${SAMPLE}_counts.txt"
done

# ================================================================================
# STEP 4: SUMMARY
# ================================================================================

echo ""
echo "=========================================="
echo "APPROACH 2 COMPLETED"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output locations:"
echo "  BAM files: ${OUTPUT_DIR}/bam/"
echo "  Count files: ${OUTPUT_DIR}/counts/"
echo "  Logs: ${OUTPUT_DIR}/logs/"
echo ""
echo "Alignment summary:"
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ -f "${OUTPUT_DIR}/logs/${SAMPLE}_hisat2_summary.txt" ]]; then
        echo "  ${SAMPLE}:"
        grep "overall alignment rate" "${OUTPUT_DIR}/logs/${SAMPLE}_hisat2_summary.txt"
    fi
done
echo ""
echo "featureCounts summary:"
cat "${OUTPUT_DIR}/counts/all_samples_counts.txt.summary"
echo ""

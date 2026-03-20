#!/bin/bash
#SBATCH --job-name=star_only
#SBATCH --output=/projects/Whitefly_RNASeq/98_test_run/99_logs/approach3_star_%j.log
#SBATCH --error=/projects/Whitefly_RNASeq/98_test_run/99_logs/approach3_star_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G


# ================================================================================
# APPROACH 3: STAR Alignment Only (from trimmed reads, no decontamination)
# ================================================================================
# Input: Trimmed clean reads
# Output: BAM alignments + gene counts
# ================================================================================

set -euo pipefail

# Module loads (if needed)
module load bioinfo-shared apptainer

# ================================================================================
# CONFIGURATION
# ================================================================================

# Directories
PROJECT_DIR="/projects/Whitefly_RNASeq/98_test_run"
INPUT_DIR="${PROJECT_DIR}/01_data/01_trimmed"
OUTPUT_DIR="${PROJECT_DIR}/02_analysis/approach3_star_only"
REF_DIR="${PROJECT_DIR}/98_refs"
LOG_DIR="${PROJECT_DIR}/99_logs"

# Reference files
GENOME_FASTA="${REF_DIR}/MEAM1_scaffold_v1.2.fa.gz"
GTF_FILE="${REF_DIR}/MEAM1_v1.2.gff3.gz"
STAR_INDEX="${REF_DIR}/star_index"  # Directory for STAR index (shared with Approach 1)

# Containers
STAR_CONTAINER="/projects/onilee/software/containers/star_2.7.11.sif" 
SUBREAD_CONTAINER="/projects/onilee/software/containers/subread_2.1.1.sif"

# Resources
THREADS=16

# Sample IDs
SAMPLES=("SRR28578498" "SRR28578499" "SRR28578500")

# ================================================================================
# CREATE OUTPUT DIRECTORIES
# ================================================================================

mkdir -p "${OUTPUT_DIR}"/{bam,counts,logs}
mkdir -p "${STAR_INDEX}"

echo "=========================================="
echo "APPROACH 3: STAR Only (No Decontamination)"
echo "=========================================="
echo "Started: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# ================================================================================
# STEP 1: BUILD STAR INDEX (if not exists)
# ================================================================================

if [ ! -f "${STAR_INDEX}/SAindex" ]; then
    echo "[$(date)] Building STAR index..."
    
    # Uncompress reference if needed
    gunzip -c "${GENOME_FASTA}" > "${REF_DIR}/MEAM1_scaffold_v1.2.fa"
    gunzip -c "${GTF_FILE}" > "${REF_DIR}/MEAM1_v1.2.gff3"
    
    apptainer exec "${STAR_CONTAINER}" STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${REF_DIR}/MEAM1_scaffold_v1.2.fa" \
        --sjdbGTFfile "${REF_DIR}/MEAM1_v1.2.gff3" \
        --sjdbGTFtagExonParentTranscript Parent \
        --runThreadN "${THREADS}" \
        --genomeSAindexNbases 12 \
        --sjdbOverhang 99
    
    echo "[$(date)] STAR index built successfully"
else
    echo "[$(date)] STAR index already exists, skipping..."
fi

# ================================================================================
# STEP 2: ALIGN WITH STAR AND QUANTIFY
# ================================================================================

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================="
    
    # Input files
    READ1="${INPUT_DIR}/${SAMPLE}_1.clean.fastq.gz"
    READ2="${INPUT_DIR}/${SAMPLE}_2.clean.fastq.gz"
    
    # Output prefix
    OUT_PREFIX="${OUTPUT_DIR}/bam/${SAMPLE}_"
    
    # Check input files exist
    if [[ ! -f "${READ1}" ]] || [[ ! -f "${READ2}" ]]; then
        echo "ERROR: Input files not found for ${SAMPLE}"
        continue
    fi
    
    echo "[$(date)] Aligning ${SAMPLE} with STAR..."
    
    apptainer exec "${STAR_CONTAINER}" STAR \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${OUT_PREFIX}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --limitBAMsortRAM 30000000000
    
    # Move log files
    mv "${OUT_PREFIX}"Log.* "${OUTPUT_DIR}/logs/"
    
    echo "[$(date)] ${SAMPLE} alignment completed"
    echo "  BAM: ${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
done

# ================================================================================
# STEP 3: QUANTIFY WITH FEATURECOUNTS
# ================================================================================

echo ""
echo "=========================================="
echo "Running featureCounts"
echo "=========================================="

# Prepare list of BAM files
BAM_FILES=""
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES="${BAM_FILES} ${OUTPUT_DIR}/bam/${SAMPLE}_Aligned.sortedByCoordinate.out.bam"
done

echo "[$(date)] Counting features..."

# Use featureCounts with proper hierarchy handling
apptainer exec "${SUBREAD_CONTAINER}" featureCounts \
    -T "${THREADS}" \
    -p \
    -t mRNA \
    -g Parent \
    -a "${REF_DIR}/MEAM1_v1.2.gff3" \
    -o "${OUTPUT_DIR}/counts/all_samples_genelevel.txt" \
    ${BAM_FILES}

# Create individual count files for each sample
echo "[$(date)] Extracting individual sample counts..."

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    COL=$((i + 7))  # featureCounts output: cols 1-6 are annotation, 7+ are samples
    
    awk -v col="$COL" 'NR==1 {print "GeneID\t"$col} NR>1 {print $1"\t"$col}' \
        "${OUTPUT_DIR}/counts/all_samples_genelevel.txt" \
        > "${OUTPUT_DIR}/counts/${SAMPLE}_genelevel.txt"
done

# ================================================================================
# STEP 4: SUMMARY
# ================================================================================

echo ""
echo "=========================================="
echo "APPROACH 3 COMPLETED"
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
    if [[ -f "${OUTPUT_DIR}/logs/${SAMPLE}_Log.final.out" ]]; then
        echo "  ${SAMPLE}:"
        grep "Uniquely mapped reads %" "${OUTPUT_DIR}/logs/${SAMPLE}_Log.final.out"
        grep "% of reads mapped to multiple loci" "${OUTPUT_DIR}/logs/${SAMPLE}_Log.final.out"
    fi
done
echo ""
echo "featureCounts summary:"
cat "${OUTPUT_DIR}/counts/all_samples_genelevel.txt.summary"
echo ""
echo "Count files summary:"
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ -f "${OUTPUT_DIR}/counts/${SAMPLE}_genelevel.txt" ]]; then
        TOTAL=$(tail -n +2 "${OUTPUT_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '{sum+=$2} END {print sum}')
        GENES=$(tail -n +2 "${OUTPUT_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '$2 > 0' | wc -l)
        echo "  ${SAMPLE}: ${TOTAL} reads assigned to ${GENES} genes"
    fi
done
echo ""

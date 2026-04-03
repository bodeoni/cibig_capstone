#!/bin/bash
#SBATCH --job-name=star_align_counts
#SBATCH --output=/scratch/onilee/capstone/99_logs/04_star_align_%j.log
#SBATCH --error=/scratch/onilee/capstone/99_logs/04_star_align_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on: $(hostname)"

module load bioinfo-shared apptainer

# ================================================================================
# CONFIGURATION
# ================================================================================

BASE_DIR="/scratch/onilee/capstone"

# Input: HISAT2-decontaminated unmapped reads
# HISAT2 --un-conc-gz with prefix "<sample>_unmapped.fastq.gz" produces
# paired files: <sample>_unmapped.fastq.gz.1 and <sample>_unmapped.fastq.gz.2
DECONTAM_DIR="${BASE_DIR}/03_analysis/03_hisat2_decontam/unmapped_fastq"

# Reference files (downloaded by 03a_download_genomes.sh)
REF_DIR="${BASE_DIR}/01_data/03_references"
GENOME_FASTA="${REF_DIR}/MEAM1_scaffold_v1.2.fa"
GFF_FILE="${REF_DIR}/MEAM1_v1.2.gff3"
STAR_INDEX="${REF_DIR}/star_index"

# Output directories
OUT_BASE="${BASE_DIR}/03_analysis/04_star_align"
BAM_DIR="${OUT_BASE}/bam"
COUNTS_DIR="${OUT_BASE}/counts"
LOG_DIR="${OUT_BASE}/logs"

# Apptainer containers
STAR_CONTAINER="/projects/onilee/software/containers/star_2.7.11.sif"
SUBREAD_CONTAINER="/projects/onilee/software/containers/subread_2.1.1.sif"

THREADS=${SLURM_CPUS_PER_TASK:-10}

# ================================================================================
# SETUP
# ================================================================================

mkdir -p "${BAM_DIR}" "${COUNTS_DIR}" "${LOG_DIR}"

mkdir -p "${BASE_DIR}/00_meta"
{
  echo "=== $(date) ==="
  echo -n "STAR: "; apptainer exec "${STAR_CONTAINER}" STAR --version 2>&1
  echo -n "featureCounts: "; apptainer exec "${SUBREAD_CONTAINER}" featureCounts -v 2>&1 | head -n 1
} >> "${BASE_DIR}/00_meta/software_versions.txt"

# ================================================================================
# STEP 1: DECOMPRESS REFERENCE FILES
# ================================================================================

if [ ! -f "${GENOME_FASTA}" ]; then
  if [ -f "${GENOME_FASTA}.gz" ]; then
    echo "[$(date)] Decompressing genome FASTA..."
    gunzip -c "${GENOME_FASTA}.gz" > "${GENOME_FASTA}"
  else
    echo "ERROR: Genome FASTA not found: ${GENOME_FASTA}(.gz)"
    exit 1
  fi
fi

if [ ! -f "${GFF_FILE}" ]; then
  if [ -f "${GFF_FILE}.gz" ]; then
    echo "[$(date)] Decompressing GFF3 annotation..."
    gunzip -c "${GFF_FILE}.gz" > "${GFF_FILE}"
  else
    echo "ERROR: GFF3 annotation not found: ${GFF_FILE}(.gz)"
    exit 1
  fi
fi

# ================================================================================
# STEP 2: BUILD STAR INDEX (if not already built)
# ================================================================================

if [ ! -d "${STAR_INDEX}" ]; then
  echo "[$(date)] Building STAR index..."
  mkdir -p "${STAR_INDEX}"
  apptainer exec "${STAR_CONTAINER}" STAR \
    --runMode genomeGenerate \
    --genomeDir "${STAR_INDEX}" \
    --genomeFastaFiles "${GENOME_FASTA}" \
    --sjdbGTFfile "${GFF_FILE}" \
    --sjdbGTFtagExonParentTranscript Parent \
    --runThreadN "${THREADS}" \
    --genomeSAindexNbases 12 \
    --sjdbOverhang 99
  echo "[$(date)] STAR index built successfully"
else
  echo "[$(date)] STAR index already exists, skipping build"
fi

# ================================================================================
# STEP 3: DISCOVER SAMPLES
# ================================================================================

# Detect samples from decontaminated FASTQ files.
# HISAT2 --un-conc-gz appends .1 and .2 to the supplied prefix, so files are
# named: <sample>_unmapped.fastq.gz.1  and  <sample>_unmapped.fastq.gz.2

shopt -s nullglob
r1_files=("${DECONTAM_DIR}"/*_unmapped.fastq.1.gz)

if [ ${#r1_files[@]} -eq 0 ]; then
  echo "ERROR: No decontaminated FASTQ files found in ${DECONTAM_DIR}"
  echo "  Expected pattern: *_unmapped.fastq.1.gz"
  exit 1
fi

echo "[$(date)] Found ${#r1_files[@]} sample(s)"

# ================================================================================
# STEP 4: ALIGN WITH STAR
# ================================================================================

BAM_FILES=()

for r1 in "${r1_files[@]}"; do
  SAMPLE=$(basename "${r1}" _unmapped.fastq.1.gz)
  r2="${DECONTAM_DIR}/${SAMPLE}_unmapped.fastq.2.gz"

  if [ ! -f "${r2}" ]; then
    echo "WARNING: R2 not found for ${SAMPLE}, skipping"
    echo "  Expected: ${r2}"
    continue
  fi

  echo ""
  echo "=========================================="
  echo "Aligning sample: ${SAMPLE}"
  echo "=========================================="
  echo "[$(date)] R1: ${r1}"
  echo "[$(date)] R2: ${r2}"

  OUT_PREFIX="${BAM_DIR}/${SAMPLE}_"

  apptainer exec "${STAR_CONTAINER}" STAR \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --runThreadN "${THREADS}" \
    --outFileNamePrefix "${OUT_PREFIX}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --limitBAMsortRAM 60000000000

  # Move STAR logs to log directory
  mv "${OUT_PREFIX}"Log.* "${LOG_DIR}/"

  BAM="${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
  BAM_FILES+=("${BAM}")

  echo "[$(date)] ${SAMPLE} alignment complete: ${BAM}"
done

if [ ${#BAM_FILES[@]} -eq 0 ]; then
  echo "ERROR: No BAM files produced — check alignment logs in ${LOG_DIR}"
  exit 1
fi

# ================================================================================
# STEP 5: QUANTIFY WITH FEATURECOUNTS
# ================================================================================

echo ""
echo "=========================================="
echo "Running featureCounts on all samples"
echo "=========================================="
echo "[$(date)] BAM files: ${BAM_FILES[*]}"

apptainer exec "${SUBREAD_CONTAINER}" featureCounts \
  -T "${THREADS}" \
  -p \
  -t mRNA \
  -g Parent \
  -a "${GFF_FILE}" \
  -o "${COUNTS_DIR}/all_samples_genelevel.txt" \
  "${BAM_FILES[@]}"

echo "[$(date)] featureCounts complete"

# ================================================================================
# STEP 6: EXTRACT INDIVIDUAL SAMPLE COUNT FILES
# ================================================================================

echo "[$(date)] Extracting per-sample count files..."

# Re-read the ordered sample list from the featureCounts header (line 2)
# Columns 1-6 are annotation fields; samples start at column 7
mapfile -t HEADER < <(head -n 2 "${COUNTS_DIR}/all_samples_genelevel.txt" | tail -n 1 | cut -f7-)

for i in "${!HEADER[@]}"; do
  BAM_PATH="${HEADER[$i]}"
  SAMPLE=$(basename "${BAM_PATH}" _Aligned.sortedByCoord.out.bam)
  COL=$((i + 7))

  awk -v col="${COL}" 'NR==1 {next} {print $1"\t"$col}' \
    "${COUNTS_DIR}/all_samples_genelevel.txt" \
    > "${COUNTS_DIR}/${SAMPLE}_genelevel.txt"

  echo "  Extracted: ${COUNTS_DIR}/${SAMPLE}_genelevel.txt"
done

# ================================================================================
# SUMMARY
# ================================================================================

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output locations:"
echo "  BAM files:   ${BAM_DIR}/"
echo "  Count files: ${COUNTS_DIR}/"
echo "  STAR logs:   ${LOG_DIR}/"
echo ""
echo "featureCounts assignment summary:"
cat "${COUNTS_DIR}/all_samples_genelevel.txt.summary"
echo ""
echo "Per-sample read assignment:"
for bam in "${BAM_FILES[@]}"; do
  SAMPLE=$(basename "${bam}" _Aligned.sortedByCoord.out.bam)
  COUNT_FILE="${COUNTS_DIR}/${SAMPLE}_genelevel.txt"
  if [ -f "${COUNT_FILE}" ]; then
    TOTAL=$(awk 'NR>1 {sum+=$2} END {print sum}' "${COUNT_FILE}")
    GENES=$(awk 'NR>1 && $2>0' "${COUNT_FILE}" | wc -l)
    echo "  ${SAMPLE}: ${TOTAL} reads assigned across ${GENES} genes"
  fi
done
echo ""
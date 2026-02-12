#!/bin/bash
#SBATCH --job-name=hisat2_decontam
#SBATCH --output=99_logs/03_decontam_hisat2_%j.log
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on: $(hostname)"

BASE_DIR="/data/onilee/capstone"

# ==============================
# Activate conda environment
# ==============================
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate bioinfo

mkdir -p "$BASE_DIR/00_meta"
{
  echo "=== $(date) ==="
  echo -n "hisat2: "; hisat2 --version | head -n 1
  echo -n "samtools: "; samtools --version | head -n 1
} >> "$BASE_DIR/00_meta/software_versions.txt"

echo "Using HISAT2 at: $(which hisat2)"
hisat2 --version | head -n 1

# ==============================
# Directories
# ==============================
CLEAN_DIR="$BASE_DIR/01_data/02_trimmed_fastq"
REF="$BASE_DIR/01_data/03_references/MEAM1_contiminants.fa"

# Recommended: processing output folder (not analysis), but keep yours if you want
OUT_BASE="$BASE_DIR/03_analysis/03_hisat2_decontam"
INDEX_DIR="$OUT_BASE/index"
BAM_DIR="$OUT_BASE/bam"
UNMAP_DIR="$OUT_BASE/unmapped_fastq"
SUMMARY_DIR="$OUT_BASE/summaries"

mkdir -p "$OUT_BASE" "$INDEX_DIR" "$BAM_DIR" "$UNMAP_DIR" "$SUMMARY_DIR"

# ==============================
# Build HISAT2 index (only once)
# ==============================
if [ ! -f "$INDEX_DIR/contam.1.ht2" ]; then
  echo "Building HISAT2 index..."
  hisat2-build -p "$SLURM_CPUS_PER_TASK" "$REF" "$INDEX_DIR/contam"
fi

# ==============================
# Run HISAT2 filtering
# ==============================
shopt -s nullglob
r1_files=("$CLEAN_DIR"/*_1.clean.fastq.gz)

if [ ${#r1_files[@]} -eq 0 ]; then
  echo "ERROR: No files matching *_1.clean.fastq.gz found in $CLEAN_DIR"
  exit 1
fi

for r1 in "${r1_files[@]}"; do
  sample=$(basename "$r1" _1.clean.fastq.gz)
  r2="$CLEAN_DIR/${sample}_2.clean.fastq.gz"

  if [ ! -f "$r2" ]; then
    echo "WARNING: Missing R2 for $sample"
    echo "  Expected: $r2"
    echo "  Skipping..."
    continue
  fi

  echo "Processing sample: $sample"
  
  #  Path for unmapped reads output (paired-end)
  unpref="$UNMAP_DIR/${sample}_unmapped.fastq.gz"

  hisat2 \
    -p "$SLURM_CPUS_PER_TASK" \
    -x "$INDEX_DIR/contam" \
    -1 "$r1" -2 "$r2" \
    --very-sensitive \
    --no-unal \
    --summary-file "$SUMMARY_DIR/${sample}.hisat2.summary.txt" \
    --un-conc-gz "$unpref" \
    -S /dev/stdout \
  | samtools view -@ "$SLURM_CPUS_PER_TASK" -b \
  | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM_DIR/${sample}.sorted.bam"

  samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM_DIR/${sample}.sorted.bam"
done

echo "Job finished at: $(date)"
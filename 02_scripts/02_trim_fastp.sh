#!/bin/bash
#SBATCH --job-name=fastp_clean
#SBATCH --output=99_logs/02_trim_fastp_%j.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

DATA_DIR="/data/onilee/capstone" # main project directory

RAW_DIR="$DATA_DIR/01_data/01_raw_fastq" # Directory containing raw FASTQ files (SRR*/ subdirectories)
OUT_REPORTS="$DATA_DIR/03_analysis/02_fastp" # Directory to save fastp reports and MultiQC results
CLEAN_FASTQ_DIR="$DATA_DIR/01_data/02_trimmed_fastq" # Directory to save trimmed FASTQ files

REPORT_DIR="$OUT_REPORTS/4-fastp_reports"
MULTIQC_DIR="$OUT_REPORTS/5-multiqc_post"

mkdir -p "$CLEAN_FASTQ_DIR" "$REPORT_DIR" "$MULTIQC_DIR"

echo "Job started at: $(date)"
echo "RAW_DIR: $RAW_DIR"
echo "CLEAN_FASTQ_DIR: $CLEAN_FASTQ_DIR"
echo "OUT_REPORTS: $OUT_REPORTS"

# Load modules
module load fastp/0.20.1
module load MultiQC/1.9

# Record software versions"
echo "Software versions for fastp step:" >> "$DATA_DIR/00_meta/software_versions.txt"
fastp --version >> "$DATA_DIR/00_meta/software_versions.txt"
multiqc --version >> "$DATA_DIR/00_meta/software_versions.txt"

echo "----------------------------------------------"
echo "----------------- RUNNING FASTP --------------"
echo "----------------------------------------------"

for dir in "$RAW_DIR"/SRR*/; do
  srr_id=$(basename "$dir")
  echo ""
  echo "Processing sample: $srr_id"
  echo "Directory: $dir"

  # Try to find paired fastqs (common fasterq-dump naming)
  r1=$(ls "$dir"/*_1*.fastq 2>/dev/null | head -n 1 || true)
  r2=$(ls "$dir"/*_2*.fastq 2>/dev/null | head -n 1 || true)

  if [[ -z "${r1:-}" || -z "${r2:-}" ]]; then
    echo "  WARNING: Could not find paired *_1/*.fastq and *_2/*.fastq for $srr_id"
    echo "  Files present:"
    ls -lh "$dir" || true
    continue
  fi

  out_r1="$CLEAN_FASTQ_DIR/${srr_id}_1.clean.fastq.gz"
  out_r2="$CLEAN_FASTQ_DIR/${srr_id}_2.clean.fastq.gz"
  html="$REPORT_DIR/${srr_id}.fastp.html"
  json="$REPORT_DIR/${srr_id}.fastp.json"

  echo "  R1: $r1"
  echo "  R2: $r2"
  echo "  OUT R1: $out_r1"
  echo "  OUT R2: $out_r2"

  fastp \
    -i "$r1" -I "$r2" \
    -o "$out_r1" -O "$out_r2" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 75 \
    --thread "$SLURM_CPUS_PER_TASK" \
    --html "$html" \
    --json "$json"

done

echo "----------------------------------------------"
echo "----------------- RUNNING MULTIQC ------------"
echo "----------------------------------------------"

multiqc -o "$MULTIQC_DIR" "$REPORT_DIR"

echo "Job finished at: $(date)"

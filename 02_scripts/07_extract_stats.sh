#!/bin/bash
#SBATCH --job-name=extract_stats
#SBATCH --output=/scratch/onilee/capstone/99_logs/07_extract_stats_%j.log
#SBATCH --error=/scratch/onilee/capstone/99_logs/07_extract_stats_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on:     $(hostname)"

module load bioinfo-wave miniconda3
eval "$(conda shell.bash hook)"
conda activate rnaseq

BASE_DIR="/scratch/onilee/capstone"
R_SCRIPT="${BASE_DIR}/02_scripts/07_extract_stats.R"

# ================================================================================
# PREFLIGHT CHECKS
# ================================================================================

echo ""
echo "=========================================="
echo "Preflight checks"
echo "=========================================="

# Check R script exists
if [ ! -f "${R_SCRIPT}" ]; then
  echo "ERROR: R script not found: ${R_SCRIPT}"
  exit 1
fi
echo "  OK: ${R_SCRIPT}"

# Check jsonlite is available (required for fastp JSON parsing)
Rscript -e 'if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite not installed")' || {
  echo "ERROR: R package 'jsonlite' is not installed in the rnaseq environment."
  echo "  Install with: install.packages('jsonlite')"
  exit 1
}
echo "  OK: jsonlite available"

# Warn (not fail) for each expected input directory
for dir in \
  "${BASE_DIR}/03_analysis/01_qc/c_fastp_reports" \
  "${BASE_DIR}/03_analysis/03_hisat2_decontam/summaries" \
  "${BASE_DIR}/03_analysis/04_star_align/logs"; do
  if [ ! -d "${dir}" ]; then
    echo "  WARNING: Directory not found (will be skipped): ${dir}"
  else
    echo "  OK: ${dir}"
  fi
done

fc_summary="${BASE_DIR}/03_analysis/04_star_align/counts/all_samples_genelevel.txt.summary"
if [ ! -f "${fc_summary}" ]; then
  echo "  WARNING: featureCounts summary not found (will be skipped): ${fc_summary}"
else
  echo "  OK: ${fc_summary}"
fi

# ================================================================================
# RUN
# ================================================================================

echo ""
echo "=========================================="
echo "Extracting pipeline statistics"
echo "=========================================="
echo "  Base directory: ${BASE_DIR}"
echo "  Rscript:        $(which Rscript)"
echo ""

Rscript "${R_SCRIPT}" "${BASE_DIR}"

# ================================================================================
# SUMMARY
# ================================================================================

OUT_DIR="${BASE_DIR}/04_results/tables/qc_stats"

echo ""
echo "=========================================="
echo "EXTRACTION COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files in ${OUT_DIR}/:"
for f in "${OUT_DIR}"/*.csv; do
  rows=$(tail -n +2 "${f}" | wc -l)
  echo "  $(basename "${f}")  (${rows} samples)"
done
echo ""

#!/bin/bash
#SBATCH --job-name=plot_qc_stats
#SBATCH --output=/scratch/onilee/capstone/99_logs/08_plot_qc_stats_%j.log
#SBATCH --error=/scratch/onilee/capstone/99_logs/08_plot_qc_stats_%j.err
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
R_SCRIPT="${BASE_DIR}/02_scripts/08_plot_qc_stats.R"

# ================================================================================
# PREFLIGHT CHECKS
# ================================================================================

echo ""
echo "=========================================="
echo "Preflight checks"
echo "=========================================="

if [ ! -f "${R_SCRIPT}" ]; then
  echo "ERROR: R script not found: ${R_SCRIPT}"
  exit 1
fi
echo "  OK: ${R_SCRIPT}"

# Check required R packages
REQUIRED_PKGS=(ggplot2 dplyr tidyr scales ggrepel pheatmap RColorBrewer grid)
for pkg in "${REQUIRED_PKGS[@]}"; do
  Rscript -e "if (!requireNamespace('${pkg}', quietly = TRUE)) quit(status = 1)" >/dev/null 2>&1 || {
    echo "ERROR: R package '${pkg}' is not installed in the rnaseq environment."
    echo "  Install with: install.packages('${pkg}')"
    exit 1
  }
  echo "  OK: R package ${pkg}"
done

# Check metadata
if [ ! -f "${BASE_DIR}/00_meta/sample_metadata.csv" ]; then
  echo "ERROR: sample_metadata.csv not found: ${BASE_DIR}/00_meta/sample_metadata.csv"
  exit 1
fi
echo "  OK: sample_metadata.csv"

# Check all four QC CSV inputs produced by 07_extract_stats
for csv in \
  "01_trimming_stats.csv" \
  "02_decontamination_stats.csv" \
  "03_star_alignment_stats.csv" \
  "04_featurecounts_stats.csv"; do
  path="${BASE_DIR}/04_results/tables/qc_stats/${csv}"
  if [ ! -f "${path}" ]; then
    echo "ERROR: Required input not found: ${path}"
    echo "  Run 07_extract_stats.sh first."
    exit 1
  fi
  echo "  OK: ${csv}"
done

# ================================================================================
# RUN
# ================================================================================

echo ""
echo "=========================================="
echo "Generating QC plots"
echo "=========================================="
echo "  Base directory: ${BASE_DIR}"
echo "  Rscript:        $(which Rscript)"
echo ""

Rscript "${R_SCRIPT}" "${BASE_DIR}"

# ================================================================================
# SUMMARY
# ================================================================================

OUT_DIR="${BASE_DIR}/04_results/figures/qc_plots"

echo ""
echo "=========================================="
echo "PLOTTING COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files in ${OUT_DIR}/:"
for f in "${OUT_DIR}"/*.{png,pdf}; do
  if [ -f "${f}" ]; then
    echo "  $(basename "${f}")"
  fi
done
echo ""

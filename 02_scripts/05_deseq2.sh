#!/bin/bash
#SBATCH --job-name=deseq2_analysis
#SBATCH --output=/scratch/onilee/capstone/99_logs/05_deseq2_%j.log
#SBATCH --error=/scratch/onilee/capstone/99_logs/05_deseq2_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on:     $(hostname)"

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate rnaseq

# ================================================================================
# CONFIGURATION
# ================================================================================

BASE_DIR="/scratch/onilee/capstone"

# ── Inputs ────────────────────────────────────────────────────────────────────
COUNT_FILE="${BASE_DIR}/03_analysis/04_star_align/counts/all_samples_genelevel.txt"
META_FILE="${BASE_DIR}/00_meta/sample_metadata.csv"
R_SCRIPT="${BASE_DIR}/02_scripts/05_deseq2.R"

# ── Outputs ───────────────────────────────────────────────────────────────────
OUT_DIR="${BASE_DIR}/03_analysis/05_deseq2"

# ================================================================================
# PREFLIGHT CHECKS
# ================================================================================

echo ""
echo "=========================================="
echo "Preflight checks"
echo "=========================================="

for f in "${COUNT_FILE}" "${META_FILE}" "${R_SCRIPT}"; do
  if [ ! -f "${f}" ]; then
    echo "ERROR: Required file not found: ${f}"
    exit 1
  fi
  echo "  OK: ${f}"
done

# ================================================================================
# SETUP
# ================================================================================

mkdir -p "${OUT_DIR}/plots" "${OUT_DIR}/tables"

# Copy R script alongside outputs for reproducibility
cp "${R_SCRIPT}" "${OUT_DIR}/05_deseq2.R"

echo ""
echo "=========================================="
echo "Running DESeq2 analysis"
echo "=========================================="
echo "  Count matrix:  ${COUNT_FILE}"
echo "  Metadata:      ${META_FILE}"
echo "  Output dir:    ${OUT_DIR}"
echo "  R environment: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "  Rscript:       $(which Rscript)"
echo ""

# ================================================================================
# RUN DESEQ2
# ================================================================================

Rscript "${R_SCRIPT}" \
  "${COUNT_FILE}" \
  "${META_FILE}" \
  "${OUT_DIR}"

# ================================================================================
# COPY KEY RESULTS TO 04_results/
# ================================================================================

echo ""
echo "=========================================="
echo "Copying results to 04_results/"
echo "=========================================="

RESULTS_DIR="${BASE_DIR}/04_results"
mkdir -p "${RESULTS_DIR}/figures" "${RESULTS_DIR}/tables"

cp "${OUT_DIR}"/plots/*.png "${RESULTS_DIR}/figures/" 2>/dev/null && \
  echo "  Copied plots to ${RESULTS_DIR}/figures/" || \
  echo "  WARNING: No PNG plots found to copy"

cp "${OUT_DIR}"/tables/*_sig.csv "${RESULTS_DIR}/tables/" 2>/dev/null && \
  echo "  Copied sig DEG tables to ${RESULTS_DIR}/tables/" || \
  echo "  WARNING: No sig tables found to copy"

# ================================================================================
# SUMMARY
# ================================================================================

echo ""
echo "=========================================="
echo "DESEQ2 COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Full analysis outputs: ${OUT_DIR}/"
echo "  tables/  - CSV files (all genes + sig DEGs, overall + per time point)"
echo "  plots/   - PCA, sample distance heatmap, volcano, MA plot, DEG heatmap"
echo ""
echo "Key results copied to: ${RESULTS_DIR}/"
echo ""

echo "Significant DEG summary:"
for f in "${OUT_DIR}"/tables/*_sig.csv; do
  label=$(basename "${f}" .csv)
  count=$(tail -n +2 "${f}" | wc -l)
  echo "  ${label}: ${count} DEGs"
done

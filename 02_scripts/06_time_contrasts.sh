#!/bin/bash
#SBATCH --job-name=time_contrasts
#SBATCH --output=/scratch/onilee/capstone/99_logs/06_time_contrasts_%j.log
#SBATCH --error=/scratch/onilee/capstone/99_logs/06_time_contrasts_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on:     $(hostname)"

module load bioinfo-wave miniconda3
eval "$(conda shell.bash hook)"
conda activate rnaseq

# ================================================================================
# CONFIGURATION
# ================================================================================

BASE_DIR="/scratch/onilee/capstone"

COUNT_FILE="${BASE_DIR}/03_analysis/04_star_align/counts/all_samples_genelevel.txt"
META_FILE="${BASE_DIR}/00_meta/sample_metadata.csv"
R_SCRIPT="${BASE_DIR}/02_scripts/06_time_contrasts.R"

OUT_DIR="${BASE_DIR}/03_analysis/06_time_contrasts"

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

cp "${R_SCRIPT}" "${OUT_DIR}/06_time_contrasts.R"

echo ""
echo "=========================================="
echo "Running time-point contrast analysis"
echo "=========================================="
echo "  Count matrix: ${COUNT_FILE}"
echo "  Metadata:     ${META_FILE}"
echo "  Output dir:   ${OUT_DIR}"
echo "  Rscript:      $(which Rscript)"
echo ""

# ================================================================================
# RUN
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
mkdir -p "${RESULTS_DIR}/figures/time_contrasts" "${RESULTS_DIR}/tables/time_contrasts"

cp "${OUT_DIR}"/plots/*.png "${RESULTS_DIR}/figures/time_contrasts/" 2>/dev/null && \
  echo "  Copied plots" || echo "  WARNING: No PNG plots found"

cp "${OUT_DIR}"/tables/*_sig_lfc1.csv "${RESULTS_DIR}/tables/time_contrasts/" 2>/dev/null && \
  echo "  Copied lfc1 sig tables" || echo "  WARNING: No sig tables found"

# ================================================================================
# SUMMARY
# ================================================================================

echo ""
echo "=========================================="
echo "TIME CONTRASTS COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Full outputs: ${OUT_DIR}/"
echo ""

echo "Significant DEG counts (|LFC|>1):"
for f in "${OUT_DIR}"/tables/*_sig_lfc1.csv; do
  label=$(basename "${f}" _sig_lfc1.csv)
  count=$(tail -n +2 "${f}" | wc -l)
  echo "  ${label}: ${count} DEGs"
done

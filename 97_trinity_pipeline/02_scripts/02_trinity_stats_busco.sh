#!/bin/bash
#SBATCH --job-name=trinity_stats_busco
#SBATCH --output=/scratch/onilee/capstone/97_trinity_pipeline/99_logs/02_trinity_stats_busco_%j.log
#SBATCH --error=/scratch/onilee/capstone/97_trinity_pipeline/99_logs/02_trinity_stats_busco_%j.err
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

# =============================================================================
# 02_trinity_stats_busco.sh
# Assembly quality assessment for the Trinity de novo transcriptome
#
# Two complementary assessments:
#   1. TrinityStats.pl — basic assembly statistics (N50, contig count, total
#      bases, ExN50 profile). Ships inside the Trinity container so no extra
#      install needed.
#   2. BUSCO (transcriptome mode, insecta_odb10) — benchmarks against single-
#      copy orthologs shared across insects to estimate assembly completeness.
#      Bemisia tabaci is an insect (Hemiptera: Aleyrodidae); insecta_odb10 is
#      the standard lineage used in published B. tabaci transcriptome papers.
#      If you want a more specific lineage, swap in hemiptera_odb10.
#
# Input:  Trinity.fasta from 01_trinity_assembly.sh
# Output: assembly_stats.txt  — TrinityStats.pl output
#         busco_insecta/       — full BUSCO run directory
#         busco_insecta/short_summary*.txt — the key result file
#
# Run after:  01_trinity_assembly.sh (Trinity.fasta must exist)
# Run before: 02_trinotate.sh (only proceed if BUSCO looks reasonable,
#             i.e. >70% complete for a de novo whitefly transcriptome)
# =============================================================================

set -euo pipefail

echo "============================================================"
echo "  Trinity Assembly Statistics + BUSCO"
echo "  Started:  $(date)"
echo "  Node:     $(hostname)"
echo "  Job ID:   ${SLURM_JOB_ID}"
echo "============================================================"

# =============================================================================
# MODULES / CONTAINERS
# =============================================================================

module load bioinfo-wave apptainer busco

# Trinity container — provides TrinityStats.pl
SIF="/projects/onilee/software/containers/trinity_2.15.2.sif"
APPTAINER_BINDS="/scratch,/projects"


# =============================================================================
# CONFIGURATION
# =============================================================================

CAPSTONE_DIR="/scratch/onilee/capstone"
TRINITY_DIR="${CAPSTONE_DIR}/97_trinity_pipeline"
ASSEMBLY_DIR="${TRINITY_DIR}/03_analysis/01_trinity_assembly"
STATS_DIR="${TRINITY_DIR}/03_analysis/01_trinity_assembly/qc"
LOG_DIR="${TRINITY_DIR}/99_logs"

FASTA="${ASSEMBLY_DIR}/Trinity.fasta"

# BUSCO settings
BUSCO_LINEAGE="insecta_odb12"    # swap to hemiptera_odb10 for a finer-grained check
BUSCO_OUT_NAME="busco_${BUSCO_LINEAGE%_odb12}"   # → busco_insecta
BUSCO_OUT_DIR="${STATS_DIR}/${BUSCO_OUT_NAME}"
BUSCO_LINEAGE_PATH="/projects/onilee/databases/lineages/${BUSCO_LINEAGE}"

THREADS=12

mkdir -p "${STATS_DIR}" "${LOG_DIR}"

# =============================================================================
# PREFLIGHT
# =============================================================================

echo ""
echo "=========================================="
echo "Preflight checks"
echo "=========================================="

if [ ! -f "${FASTA}" ]; then
    echo "ERROR: Trinity.fasta not found at:"
    echo "  ${FASTA}"
    echo ""
    echo "  If you moved the assembly to a different location, update"
    echo "  ASSEMBLY_DIR in this script and re-submit."
    exit 1
fi

FASTA_SIZE=$(du -sh "${FASTA}" | cut -f1)
N_TRANSCRIPTS=$(grep -c "^>" "${FASTA}")
echo "  OK: Trinity.fasta found"
echo "    Path:        ${FASTA}"
echo "    Size:        ${FASTA_SIZE}"
echo "    Transcripts: ${N_TRANSCRIPTS}"
echo ""

if [ ! -d "${BUSCO_LINEAGE_PATH}" ]; then
    echo "ERROR: BUSCO lineage not found at:"
    echo "  ${BUSCO_LINEAGE_PATH}"
    echo "  Check that the lineage directory exists under /projects/onilee/databases/lineages/"
    exit 1
fi
echo "  OK: BUSCO lineage found: ${BUSCO_LINEAGE_PATH}"
echo ""

# =============================================================================
# SECTION 1 — TRINITYSTATISTICS.PL
# =============================================================================

echo "=========================================="
echo "Section 1: TrinityStats.pl"
echo "=========================================="
echo ""

STATS_FILE="${STATS_DIR}/assembly_stats_${SLURM_JOB_ID}.txt"

echo "--- TrinityStats.pl output ---"
apptainer exec --bind "${APPTAINER_BINDS}" "${SIF}" \
    TrinityStats.pl "${FASTA}" | tee "${STATS_FILE}"

echo ""
echo "Stats saved to: ${STATS_FILE}"

# =============================================================================
# SECTION 2 — ExN50 profile
# =============================================================================
# ExN50 is a more informative metric than N50 for transcriptomes because it
# considers transcript expression level. It requires the abundance estimates
# from Salmon/RSEM, so we skip it here (Salmon runs in step 04). A note is
# left so this can be revisited after quantification.

echo ""
echo "=========================================="
echo "Section 2: ExN50 (deferred — needs Salmon quant)"
echo "=========================================="
echo "  ExN50 requires per-transcript abundance estimates."
echo "  Re-run after 04_salmon_quant.sh with:"
echo ""
echo "    apptainer exec --bind \"${APPTAINER_BINDS}\" \"${SIF}\" \\"
echo "      \$(which contig_ExN50_statistic.pl) \\"
echo "      <salmon_quant_dir>/quant.sf \\"
echo "      ${FASTA} | tee ${STATS_DIR}/ExN50_profile.txt"
echo ""

# =============================================================================
# SECTION 3 — BUSCO
# =============================================================================

echo "=========================================="
echo "Section 3: BUSCO (${BUSCO_LINEAGE})"
echo "=========================================="
echo ""

# BUSCO writes to the current directory by default; cd to STATS_DIR so all
# output lands in one place.  The --out_path flag (BUSCO ≥5) makes this explicit.

if command -v busco &>/dev/null; then
    echo "  BUSCO version: $(busco --version 2>&1 | head -1)"
else
    echo "  WARNING: busco not found in PATH. Check module load above."
    echo "  Skipping BUSCO section."
    # Don't exit — TrinityStats output is already saved above
fi

if command -v busco &>/dev/null; then
    # Remove a previous partial BUSCO run for this lineage to avoid confusion
    if [ -d "${BUSCO_OUT_DIR}" ]; then
        echo "  Removing previous BUSCO output directory: ${BUSCO_OUT_DIR}"
        rm -rf "${BUSCO_OUT_DIR}"
    fi

    echo ""
    echo "  Running BUSCO..."
    echo "  Lineage:    ${BUSCO_LINEAGE}"
    echo "  Mode:       transcriptome"
    echo "  Output dir: ${BUSCO_OUT_DIR}"
    echo "  CPUs:       ${THREADS}"
    echo ""

    busco \
        --in              "${FASTA}" \
        --out             "${BUSCO_OUT_NAME}" \
        --out_path        "${STATS_DIR}" \
        --mode            transcriptome \
        --lineage_dataset "${BUSCO_LINEAGE_PATH}" \
        --offline \
        --cpu             "${THREADS}" \
        --force

    echo ""
    echo "BUSCO finished: $(date)"

    # Print the short summary
    SUMMARY=$(find "${BUSCO_OUT_DIR}" -name "short_summary*.txt" | head -1)
    if [ -n "${SUMMARY}" ]; then
        echo ""
        echo "--- BUSCO short summary ---"
        cat "${SUMMARY}"
        # Copy to logs for easy reference
        cp "${SUMMARY}" "${LOG_DIR}/busco_short_summary_${SLURM_JOB_ID}.txt"
        echo ""
        echo "Summary copied to: ${LOG_DIR}/busco_short_summary_${SLURM_JOB_ID}.txt"
    else
        echo "WARNING: short_summary file not found in ${BUSCO_OUT_DIR}"
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "============================================================"
echo "  ASSEMBLY QC COMPLETE"
echo "  Finished: $(date)"
echo "============================================================"
echo ""
echo "  Outputs:"
echo "    TrinityStats:    ${STATS_FILE}"
echo "    BUSCO directory: ${BUSCO_OUT_DIR}"
echo "    BUSCO summary:   ${LOG_DIR}/busco_short_summary_${SLURM_JOB_ID}.txt"
echo ""
echo "  Interpretation guidance:"
echo "    N50 ≥ 1 kb        — reasonable contig length for a de novo assembly"
echo "    BUSCO complete ≥ 70%  — acceptable for de novo whitefly transcriptome"
echo "    BUSCO duplicated  — high duplication can reflect alternative splicing"
echo "                        or assembly artifacts; review alongside N50"
echo ""
echo "  Next steps (only proceed if BUSCO results are satisfactory):"
echo "    2. Functional annotation:  sbatch 02_trinotate.sh"
echo "    3. Salmon index:           sbatch 03_salmon_index.sh"
echo "    4. Salmon quantification:  sbatch 04_salmon_quant.sh  (one job per sample)"
echo "    5. Differential expression: sbatch 05_deseq2.sh"
echo ""

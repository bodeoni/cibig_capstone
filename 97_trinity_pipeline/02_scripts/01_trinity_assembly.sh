#!/bin/bash
#SBATCH --job-name=trinity_assembly
#SBATCH --output=/scratch/onilee/capstone/97_trinity_pipeline/99_logs/01_trinity_assembly_%j.log
#SBATCH --error=/scratch/onilee/capstone/97_trinity_pipeline/99_logs/01_trinity_assembly_%j.err
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --partition=normal
#SBATCH --nodelist=node06

# =============================================================================
# 01_trinity_assembly.sh
# De novo transcriptome assembly using Trinity
#
# Strategy: all 18 decontaminated samples are assembled together into a single
# reference transcriptome. Using all samples maximises the chances of capturing
# lowly expressed transcripts and transcripts that are condition-specific
# (e.g. only expressed in viruliferous whiteflies at 72h).
#
# Input:  Decontaminated paired-end FASTQs from HISAT2 step
#         (03_analysis/03_hisat2_decontam/unmapped_fastq/)
# Output: Trinity.fasta  — the assembled transcriptome
#         Trinity.fasta.gene_trans_map — gene-to-transcript mapping (for tximport)
#
# After this script:
#   → 02_trinotate.sh   — functional annotation of Trinity.fasta
#   → 03_salmon_index.sh + 04_salmon_quant.sh — quantification
#
# Resume behaviour: Trinity automatically resumes an interrupted job if the
# output directory already exists and contains partial progress. To start
# fresh, delete the output directory first:
#   rm -rf ${ASSEMBLY_DIR}
# =============================================================================

set -euo pipefail

echo "=============================================="
echo "  Trinity De Novo Transcriptome Assembly"
echo "  Started: $(date)"
echo "  Node:    $(hostname)"
echo "  Job ID:  ${SLURM_JOB_ID}"
echo "=============================================="

# =============================================================================
# MODULES
# =============================================================================

module load bioinfo-shared trinity

echo ""
echo "Trinity version: $(Trinity --version 2>&1 | head -1)"
echo ""

# =============================================================================
# CONFIGURATION
# =============================================================================

CAPSTONE_DIR="/scratch/onilee/capstone"
TRINITY_DIR="${CAPSTONE_DIR}/97_trinity_pipeline"

# Input: decontaminated reads from the main pipeline HISAT2 step

DECONTAM_DIR="${CAPSTONE_DIR}/03_analysis/03_hisat2_decontam/unmapped_fastq"

# Output directory — Trinity writes Trinity.fasta here
ASSEMBLY_DIR="${TRINITY_DIR}/03_analysis/01_trinity_assembly"
LOG_DIR="${TRINITY_DIR}/99_logs"

# Resources (keep in sync with #SBATCH headers above)
THREADS=12
MAX_MEM="100G"

# Strand specificity
# If the library was prepared with a strand-preserving protocol (e.g. Illumina
# TruSeq Stranded, dUTP method), set SS_LIB_TYPE to RF. If unstranded, leave
# it empty. Check the BioProject metadata or SRA record for library layout.
#
#   RF  →  reverse-stranded (most common for dUTP/TruSeq Stranded)
#   FR  →  forward-stranded
#   ""  →  unstranded (comment out the flag below)
#
#SS_LIB_TYPE="RF"

# In-silico normalisation cap (max coverage per k-mer)
# Trinity's built-in normalisation reduces the read set to ≤50× k-mer coverage
# before assembly. This dramatically lowers memory and runtime for large
# datasets with no meaningful loss in assembly quality. Strongly recommended
# when total input exceeds ~300 M read pairs.
NORMALIZE_MAX_COV=50

mkdir -p "${LOG_DIR}"

# =============================================================================
# PREFLIGHT CHECKS
# =============================================================================

echo "=========================================="
echo "Preflight checks"
echo "=========================================="

# Check decontamination directory
if [ ! -d "${DECONTAM_DIR}" ]; then
    echo "ERROR: Decontamination output directory not found:"
    echo "  ${DECONTAM_DIR}"
    echo ""
    echo "  Check whether the main pipeline used 02_ or 03_ as the prefix:"
    echo "    ls ${CAPSTONE_DIR}/03_analysis/"
    exit 1
fi
echo "  OK: decontamination directory found"

# Discover R1 (*.1) and R2 (*.2) files
# HISAT2 --un-conc-gz with a .gz prefix produces files ending in .gz.1 / .gz.2
r1_files=("${DECONTAM_DIR}"/*_unmapped.fastq.1.gz)
r2_files=("${DECONTAM_DIR}"/*_unmapped.fastq.2.gz)

if [ "${#r1_files[@]}" -eq 0 ] || [ ! -f "${r1_files[0]}" ]; then
    echo "ERROR: No R1 files (*_unmapped.fastq.1.gz) found in:"
    echo "  ${DECONTAM_DIR}"
    exit 1
fi

if [ "${#r1_files[@]}" -ne "${#r2_files[@]}" ]; then
    echo "ERROR: R1 and R2 file counts do not match."
    echo "  R1 files: ${#r1_files[@]}"
    echo "  R2 files: ${#r2_files[@]}"
    exit 1
fi

echo "  OK: found ${#r1_files[@]} R1/R2 file pairs"
echo ""
echo "  Samples to assemble:"
for r1 in "${r1_files[@]}"; do
    sample=$(basename "${r1}" _unmapped.fastq.1.gz)
    echo "    ${sample}"
done

# Check for existing output (Trinity resume behaviour)
if [ -d "${ASSEMBLY_DIR}" ] && [ -f "${ASSEMBLY_DIR}/Trinity.fasta" ]; then
    echo ""
    echo "  NOTE: Trinity.fasta already exists in output directory."
    echo "    The assembly appears complete. To re-assemble from scratch:"
    echo "      rm -rf ${ASSEMBLY_DIR}"
    echo "    Exiting — delete the directory if a fresh run is intended."
    exit 0
elif [ -d "${ASSEMBLY_DIR}" ]; then
    echo ""
    echo "  NOTE: Output directory exists but Trinity.fasta not found."
    echo "    Trinity will attempt to resume from partial progress."
fi

echo ""

# =============================================================================
# BUILD INPUT FILE LISTS
# =============================================================================
# Trinity --left and --right take comma-separated lists of file paths.

LEFT_FILES=$(printf "%s," "${r1_files[@]}")
LEFT_FILES="${LEFT_FILES%,}"   # strip trailing comma

RIGHT_FILES=$(printf "%s," "${r2_files[@]}")
RIGHT_FILES="${RIGHT_FILES%,}"

echo "=========================================="
echo "Input summary"
echo "=========================================="
echo "  Samples:        ${#r1_files[@]}"
echo "  Decontam dir:   ${DECONTAM_DIR}"
echo "  Output dir:     ${ASSEMBLY_DIR}"
echo "  CPUs:           ${THREADS}"
echo "  Max memory:     ${MAX_MEM}"
#echo "  Strand type:    ${SS_LIB_TYPE:-unstranded}"
echo "  Norm. max cov:  ${NORMALIZE_MAX_COV}x"
echo ""

# =============================================================================
# RUN TRINITY
# =============================================================================

echo "=========================================="
echo "Running Trinity assembly"
echo "Started: $(date)"
echo "=========================================="
echo ""

Trinity \
    --seqType      fq \
    --left         "${LEFT_FILES}" \
    --right        "${RIGHT_FILES}" \
    #--SS_lib_type  "${SS_LIB_TYPE}" \
    --CPU          "${THREADS}" \
    --max_memory   "${MAX_MEM}" \
    --output       "${ASSEMBLY_DIR}" \
    --normalize_reads \
    --normalize_max_cov "${NORMALIZE_MAX_COV}" \
    --verbose

echo ""
echo "Trinity assembly finished: $(date)"

# =============================================================================
# POST-ASSEMBLY STATISTICS
# =============================================================================

echo ""
echo "=========================================="
echo "Assembly statistics"
echo "=========================================="

FASTA="${ASSEMBLY_DIR}/Trinity.fasta"

if [ ! -f "${FASTA}" ]; then
    echo "ERROR: Trinity.fasta not found — assembly may have failed."
    echo "  Check the log: ${LOG_DIR}/01_trinity_assembly_${SLURM_JOB_ID}.log"
    exit 1
fi

# TrinityStats.pl ships with Trinity and reports N50, total bases,
# transcript and gene counts
echo ""
echo "--- TrinityStats.pl ---"
TrinityStats.pl "${FASTA}"

# Save stats to file for records
TrinityStats.pl "${FASTA}" > "${TRINITY_DIR}/99_logs/trinity_assembly_stats_${SLURM_JOB_ID}.txt" 2>&1
echo ""
echo "Stats saved to: ${TRINITY_DIR}/99_logs/trinity_assembly_stats_${SLURM_JOB_ID}.txt"

# Basic file info
echo ""
echo "--- Assembly file sizes ---"
ls -lh "${ASSEMBLY_DIR}/Trinity.fasta"
ls -lh "${ASSEMBLY_DIR}/Trinity.fasta.gene_trans_map" 2>/dev/null || \
    echo "  (gene_trans_map not found — check Trinity version)"

# Count transcripts and genes
N_TRANSCRIPTS=$(grep -c "^>" "${FASTA}")
echo ""
echo "  Total transcripts assembled: ${N_TRANSCRIPTS}"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=============================================="
echo "  TRINITY ASSEMBLY COMPLETE"
echo "  Finished: $(date)"
echo "=============================================="
echo ""
echo "  Output:"
echo "    Transcriptome FASTA:   ${ASSEMBLY_DIR}/Trinity.fasta"
echo "    Gene-transcript map:   ${ASSEMBLY_DIR}/Trinity.fasta.gene_trans_map"
echo "    Assembly stats:        ${TRINITY_DIR}/99_logs/trinity_assembly_stats_${SLURM_JOB_ID}.txt"
echo ""
echo "  Next steps:"
echo "    1. Review assembly stats (N50, number of transcripts, % BUSCO completeness)"
echo "    2. Run Trinotate annotation:  sbatch 02_trinotate.sh"
echo "    3. Build Salmon index:        sbatch 03_salmon_index.sh"
echo ""
echo "  BUSCO completeness assessment (recommended before proceeding):"
echo "    module load bioinfo-shared busco"
echo "    busco -i ${ASSEMBLY_DIR}/Trinity.fasta \\"
echo "          -o busco_trinity -m transcriptome \\"
echo "          -l insecta_odb10 --cpu ${THREADS}"
echo ""

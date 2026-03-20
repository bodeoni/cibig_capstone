#!/bin/bash

# ================================================================================
# COMPARE RESULTS FROM ALL THREE APPROACHES - UPDATED
# ================================================================================
# This script compares mapping rates and count statistics across all approaches
# ================================================================================

set -euo pipefail

PROJECT_DIR="/projects/Whitefly_RNASeq/98_test_run"
RESULTS_DIR="${PROJECT_DIR}/03_results"

mkdir -p "${RESULTS_DIR}"

echo "=========================================="
echo "COMPARISON OF THREE APPROACHES"
echo "=========================================="
echo "Generated: $(date)"
echo ""

# ================================================================================
# APPROACH 1: STAR with decontaminated reads
# ================================================================================

echo "=========================================="
echo "APPROACH 1: STAR (decontaminated input)"
echo "=========================================="

APP1_DIR="${PROJECT_DIR}/02_analysis/approach1_star"

if [[ -d "${APP1_DIR}" ]]; then
    echo ""
    echo "Mapping Statistics:"
    for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
        LOG_FILE="${APP1_DIR}/logs/${SAMPLE}_Log.final.out"
        if [[ -f "${LOG_FILE}" ]]; then
            echo ""
            echo "  ${SAMPLE}:"
            echo "    Input reads: $(grep "Number of input reads" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Uniquely mapped: $(grep "Uniquely mapped reads %" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Multi-mapped: $(grep "% of reads mapped to multiple loci" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Unmapped (too short): $(grep "% of reads unmapped: too short" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Unmapped (other): $(grep "% of reads unmapped: other" "${LOG_FILE}" | awk '{print $NF}')"
        fi
    done
    
    echo ""
    echo "Gene Assignment (featureCounts):"
    if [[ -f "${APP1_DIR}/counts/all_samples_genelevel.txt.summary" ]]; then
        cat "${APP1_DIR}/counts/all_samples_genelevel.txt.summary"
    fi
    
    echo ""
    for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
        COUNT_FILE="${APP1_DIR}/counts/${SAMPLE}_genelevel.txt"
        if [[ -f "${COUNT_FILE}" ]]; then
            TOTAL=$(tail -n +2 "${COUNT_FILE}" | awk '{sum+=$2} END {print sum}')
            GENES=$(tail -n +2 "${COUNT_FILE}" | awk '$2 > 0' | wc -l)
            echo "  ${SAMPLE}: ${TOTAL} reads assigned to ${GENES} genes"
        fi
    done
else
    echo "  Not run yet"
fi

# ================================================================================
# APPROACH 2: HISAT2 + featureCounts
# ================================================================================

echo ""
echo "=========================================="
echo "APPROACH 2: HISAT2 + featureCounts"
echo "=========================================="

APP2_DIR="${PROJECT_DIR}/02_analysis/approach2_hisat2"

if [[ -d "${APP2_DIR}" ]]; then
    echo ""
    echo "Mapping Statistics:"
    for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
        LOG_FILE="${APP2_DIR}/logs/${SAMPLE}_hisat2_summary.txt"
        if [[ -f "${LOG_FILE}" ]]; then
            echo ""
            echo "  ${SAMPLE}:"
            cat "${LOG_FILE}"
        fi
    done
    
    echo ""
    echo "Gene Assignment (featureCounts):"
    if [[ -f "${APP2_DIR}/counts/all_samples_counts.txt.summary" ]]; then
        cat "${APP2_DIR}/counts/all_samples_counts.txt.summary"
    fi
else
    echo "  Not run yet"
fi

# ================================================================================
# APPROACH 3: STAR only (no decontamination)
# ================================================================================

echo ""
echo "=========================================="
echo "APPROACH 3: STAR Only (no decontam)"
echo "=========================================="

APP3_DIR="${PROJECT_DIR}/02_analysis/approach3_star_only"

if [[ -d "${APP3_DIR}" ]]; then
    echo ""
    echo "Mapping Statistics:"
    for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
        LOG_FILE="${APP3_DIR}/logs/${SAMPLE}_Log.final.out"
        if [[ -f "${LOG_FILE}" ]]; then
            echo ""
            echo "  ${SAMPLE}:"
            echo "    Input reads: $(grep "Number of input reads" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Uniquely mapped: $(grep "Uniquely mapped reads %" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Multi-mapped: $(grep "% of reads mapped to multiple loci" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Unmapped (too short): $(grep "% of reads unmapped: too short" "${LOG_FILE}" | awk '{print $NF}')"
            echo "    Unmapped (other): $(grep "% of reads unmapped: other" "${LOG_FILE}" | awk '{print $NF}')"
        fi
    done
    
    echo ""
    echo "Gene Assignment (featureCounts):"
    if [[ -f "${APP3_DIR}/counts/all_samples_genelevel.txt.summary" ]]; then
        cat "${APP3_DIR}/counts/all_samples_genelevel.txt.summary"
    fi
    
    echo ""
    for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
        COUNT_FILE="${APP3_DIR}/counts/${SAMPLE}_genelevel.txt"
        if [[ -f "${COUNT_FILE}" ]]; then
            TOTAL=$(tail -n +2 "${COUNT_FILE}" | awk '{sum+=$2} END {print sum}')
            GENES=$(tail -n +2 "${COUNT_FILE}" | awk '$2 > 0' | wc -l)
            echo "  ${SAMPLE}: ${TOTAL} reads assigned to ${GENES} genes"
        fi
    done
else
    echo "  Not run yet"
fi

# ================================================================================
# SIDE-BY-SIDE COMPARISON TABLE
# ================================================================================

echo ""
echo "=========================================="
echo "SIDE-BY-SIDE COMPARISON"
echo "=========================================="

# Create comparison CSV
COMP_FILE="${RESULTS_DIR}/approach_comparison.csv"

echo "Sample,Approach,Input_Reads,Uniquely_Mapped_%,Multi_Mapped_%,Reads_Assigned,Genes_Detected" > "${COMP_FILE}"

for SAMPLE in SRR28578498 SRR28578499 SRR28578500; do
    
    # Approach 1
    if [[ -f "${APP1_DIR}/logs/${SAMPLE}_Log.final.out" ]]; then
        INPUT=$(grep "Number of input reads" "${APP1_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        UNIQUE=$(grep "Uniquely mapped reads %" "${APP1_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        MULTI=$(grep "% of reads mapped to multiple loci" "${APP1_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        
        if [[ -f "${APP1_DIR}/counts/${SAMPLE}_genelevel.txt" ]]; then
            ASSIGNED=$(tail -n +2 "${APP1_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '{sum+=$2} END {print sum}')
            GENES=$(tail -n +2 "${APP1_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '$2 > 0' | wc -l)
        else
            ASSIGNED="NA"
            GENES="NA"
        fi
        
        echo "${SAMPLE},Approach1_STAR_decontam,${INPUT},${UNIQUE},${MULTI},${ASSIGNED},${GENES}" >> "${COMP_FILE}"
    fi
    
    # Approach 3
    if [[ -f "${APP3_DIR}/logs/${SAMPLE}_Log.final.out" ]]; then
        INPUT=$(grep "Number of input reads" "${APP3_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        UNIQUE=$(grep "Uniquely mapped reads %" "${APP3_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        MULTI=$(grep "% of reads mapped to multiple loci" "${APP3_DIR}/logs/${SAMPLE}_Log.final.out" | awk '{print $NF}')
        
        if [[ -f "${APP3_DIR}/counts/${SAMPLE}_genelevel.txt" ]]; then
            ASSIGNED=$(tail -n +2 "${APP3_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '{sum+=$2} END {print sum}')
            GENES=$(tail -n +2 "${APP3_DIR}/counts/${SAMPLE}_genelevel.txt" | awk '$2 > 0' | wc -l)
        else
            ASSIGNED="NA"
            GENES="NA"
        fi
        
        echo "${SAMPLE},Approach3_STAR_only,${INPUT},${UNIQUE},${MULTI},${ASSIGNED},${GENES}" >> "${COMP_FILE}"
    fi
    
done

if [[ -f "${COMP_FILE}" ]]; then
    echo ""
    echo "Comparison table saved to: ${COMP_FILE}"
    echo ""
    column -t -s',' "${COMP_FILE}"
fi

# ================================================================================
# GENE-LEVEL COMPARISON
# ================================================================================

echo ""
echo "=========================================="
echo "GENE-LEVEL DIFFERENCES"
echo "=========================================="

# Compare gene counts between Approach 1 and 3
SAMPLE="SRR28578498"  # Use first sample as example

if [[ -f "${APP1_DIR}/counts/${SAMPLE}_genelevel.txt" ]] && [[ -f "${APP3_DIR}/counts/${SAMPLE}_genelevel.txt" ]]; then
    
    echo ""
    echo "Analyzing ${SAMPLE}..."
    
    # Join the count files
    join -t $'\t' -1 1 -2 1 \
        <(tail -n +2 "${APP1_DIR}/counts/${SAMPLE}_genelevel.txt" | sort -k1,1) \
        <(tail -n +2 "${APP3_DIR}/counts/${SAMPLE}_genelevel.txt" | sort -k1,1) | \
        awk 'BEGIN{OFS="\t"} {
            diff = $3 - $2
            fold = ($2 > 0) ? $3/$2 : "inf"
            print $1, $2, $3, diff, fold
        }' | \
        sort -k4,4rn > "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt"
    
    # Add header
    echo -e "GeneID\tApp1_Count\tApp3_Count\tDifference\tFold_Change" | \
        cat - "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt" > "${RESULTS_DIR}/${SAMPLE}_gene_comparison_header.txt"
    
    mv "${RESULTS_DIR}/${SAMPLE}_gene_comparison_header.txt" "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt"
    
    echo ""
    echo "Top 20 genes with MOST EXTRA counts in Approach 3 (no decontam):"
    head -21 "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt" | column -t
    
    # Summary stats
    echo ""
    TOTAL_GENES=$(tail -n +2 "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt" | wc -l)
    INFLATED=$(awk '$5 > 2' "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt" | tail -n +2 | wc -l)
    SIMILAR=$(awk '$5 >= 0.8 && $5 <= 1.2' "${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt" | tail -n +2 | wc -l)
    
    echo "Gene-level summary:"
    echo "  Total genes compared: ${TOTAL_GENES}"
    echo "  Similar counts (±20%): ${SIMILAR} ($(echo "scale=1; ${SIMILAR}*100/${TOTAL_GENES}" | bc)%)"
    echo "  Inflated in App3 (>2x): ${INFLATED} ($(echo "scale=1; ${INFLATED}*100/${TOTAL_GENES}" | bc)%)"
fi

echo ""
echo "=========================================="
echo "COMPARISON COMPLETE"
echo "=========================================="
echo ""
echo "Files created:"
echo "  ${RESULTS_DIR}/approach_comparison.csv"
echo "  ${RESULTS_DIR}/${SAMPLE}_gene_comparison.txt"
echo ""

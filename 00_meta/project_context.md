# Project Context: CiBIG Capstone — Viruliferous vs Non-viruliferous *Bemisia tabaci* RNA-seq

> **For the CLI agent reading this file:**
> This document is your complete mental map of the project. Most raw data, intermediate files, and large outputs live on a remote HPC cluster (`/scratch/onilee/capstone` on the WAVE cluster). You do **not** need cluster access to write a report. Everything required — results tables, figures, and analysis outputs — has been downloaded to the **local working directory** you are currently in. All file paths in this document are relative to that local root unless stated otherwise.

---

## 1. Project Overview

### Biological Question
Does begomovirus acquisition alter the transcriptome of the whitefly vector *Bemisia tabaci*, and if so, does the effect depend on how long the insect has been feeding on an infected plant?

### Organism
*Bemisia tabaci* MEAM1 biotype (silverleaf whitefly) — a major agricultural pest and the primary vector of begomoviruses (family *Geminiviridae*). The MEAM1 genome (scaffold v1.2) is the reference used throughout.

### Experimental Design
- **18 samples** total: 3 time points × 2 groups × 3 biological replicates
- **Groups:** `virus` (viruliferous — fed on begomovirus-infected plants) vs `control` (non-viruliferous)
- **Time points:** 24h, 48h, 72h post-acquisition access period
- **Library type:** Bulk RNA-seq, paired-end, 150 bp reads
- **Data source:** BioProject PRJNA1096732 (public SRA data)
- Full sample table: `00_meta/sample_metadata.csv`

### Analytical Goals
1. Identify differentially expressed genes (DEGs) between viruliferous and control whiteflies at each time point and overall
2. Characterise the temporal transcriptional response (how gene expression changes over 24h → 48h → 72h) independently in the virus and control groups
3. Determine whether the dominant signal in the data is driven by virus acquisition or by time

---

## 2. Technical Stack

### Execution Environment
- **Cluster:** WAVE HPC, SLURM scheduler, node `node06`
- **Local OS:** Windows with WSL2 (Ubuntu)
- **Working directory (cluster):** `/scratch/onilee/capstone`
- **Working directory (local):** current directory (this repo root)

### Tools by Pipeline Stage

| Stage | Tool | Version | Execution method |
|---|---|---|---|
| SRA download | fasterq-dump (SRA Toolkit) | — | module |
| Raw QC | FastQC | v0.12.1 | module |
| QC aggregation | MultiQC | v1.9 | module |
| Adapter trimming | fastp | v0.20.1 | module |
| Decontamination | HISAT2 | v2.2.2 | conda: `bioinfo` |
| BAM handling | samtools | v1.23 | conda: `bioinfo` |
| Alignment | STAR | v2.7.11b | Apptainer container |
| Quantification | featureCounts (Subread) | v2.1.1 | Apptainer container |
| Differential expression | DESeq2 | v1.50+ | conda: `rnaseq` (R) |
| Visualisation | ggplot2, pheatmap, ggrepel | — | conda: `rnaseq` (R) |
| Stats extraction | Custom R script | — | conda: `rnaseq` (R) |

---

## 3. Data Processing Pipeline

The pipeline runs in 9 sequential steps. Steps 1–5 produce intermediate files that remain on the cluster. Steps 6–9 produce the analysis outputs that have been downloaded locally.

---

### Step 1 — SRA to FASTQ
**Script:** `02_scripts/00-extract_fastq.sh`
**Tool:** fasterq-dump `--split-3`
**What it does:** Downloads SRA archives and splits them into paired R1/R2 FASTQ files.
**Output (cluster only):** `01_data/01_raw_fastq/SRR*/{SRR*_1.fastq, SRR*_2.fastq}`

---

### Step 2 — Raw Read QC
**Script:** `02_scripts/01_qc_fastq.sh`
**Tools:** FastQC + MultiQC
**What it does:** Assesses base quality, GC content, adapter content, and sequence duplication on raw reads before any processing.
**Output (cluster only):** `03_analysis/01_qc/`

---

### Step 3 — Adapter Trimming
**Script:** `02_scripts/02_trim_fastp.sh`
**Tool:** fastp
**What it does:** Removes sequencing adapters, performs sliding-window quality trimming from both ends, and discards reads shorter than 75 bp post-trimming. Also generates per-sample JSON reports used later by the stats extraction step.
**Key parameters:** `--detect_adapter_for_pe`, `--cut_front`, `--cut_tail`, `--cut_window_size 4`, `--cut_mean_quality 20`, `--length_required 75`
**Output (cluster only):** `01_data/02_trimmed_fastq/` (clean FASTQs), `03_analysis/02_fastp/` (JSON/HTML reports)

---

### Step 4 — Endosymbiont Decontamination
**Script:** `02_scripts/03b_decontam_hisat2.sh`
**Tool:** HISAT2 (used here as a filter, not for transcript alignment)
**What it does:** Aligns trimmed reads to a combined contaminant reference containing three endosymbiont genomes (Portiera, Hamiltonella, Rickettsia) and the mitochondrial genome. Reads that **do not map** to this contaminant reference are retained as putative host-derived reads. This step is critical — endosymbionts are obligate intracellular bacteria in whitefly cells and contribute very heavily to the total RNA pool.
**Output (cluster only):** `03_analysis/03_hisat2_decontam/unmapped_fastq/` (decontaminated FASTQs ending in `_unmapped.fastq.gz.1` / `.gz.2`)

---

### Step 5 — STAR Alignment + featureCounts Quantification
**Script:** `02_scripts/04_star_align_counts.sh`
**Tools:** STAR (alignment), featureCounts (quantification)
**What it does:** Aligns decontaminated reads to the MEAM1 host genome. STAR is used instead of HISAT2 here because it is better suited for splice-aware alignment with an annotated reference and produces the alignment statistics needed for QC. featureCounts then counts reads over exonic regions and assigns them to genes.
**Key parameters:**
- STAR: `--outSAMtype BAM SortedByCoordinate`, `--sjdbGTFtagExonParentTranscript Parent`, `--genomeSAindexNbases 12`
- featureCounts: `-p` (paired-end), `-t exon`, `-g Parent` (groups by transcript via the GFF3 `Parent` attribute; genes are predominantly single-isoform so this is effectively gene-level)
**Output (cluster only):** `03_analysis/04_star_align/` (BAMs + count matrix)
**Key output file:** `03_analysis/04_star_align/counts/all_samples_genelevel.txt` — the count matrix fed into DESeq2

---

### Step 6 — Differential Expression: Virus vs Control
**Script:** `02_scripts/05_deseq2.sh` → calls `02_scripts/05_deseq2.R`
**Tool:** DESeq2
**What it does:** Tests for differential expression between virus and control groups. Runs four parallel models:

| Analysis | DESeq2 design | Contrast | Purpose |
|---|---|---|---|
| Overall | `~ Time_Point + Group` | virus vs control | Estimates virus effect across all time points, controlling for time |
| 24h subset | `~ Group` | virus vs control | Virus effect at 24h only (6 samples) |
| 48h subset | `~ Group` | virus vs control | Virus effect at 48h only (6 samples) |
| 72h subset | `~ Group` | virus vs control | Virus effect at 72h only (6 samples) |

**Reference levels:** `control` (Group), `24h` (Time_Point)
**Low-count pre-filter:** genes with fewer than 10 total counts removed before fitting
**P-value adjustment:** Benjamini-Hochberg (BH) FDR

**Output (available locally):** `03_analysis/05_deseq2/tables/` and `03_analysis/05_deseq2/plots/`

---

### Step 7 — Differential Expression: Time-point Contrasts
**Script:** `02_scripts/06_time_contrasts.sh` → calls `02_scripts/06_time_contrasts.R`
**Tool:** DESeq2
**What it does:** Tests how gene expression changes over time (48h vs 24h and 72h vs 24h), run in three parallel models:

| Analysis label | DESeq2 design | Samples used | Purpose |
|---|---|---|---|
| `overall_*` | `~ Group + Time_Point` | All 18 | Time effect with group as covariate |
| `virus_*` | `~ Time_Point` | 9 virus only | Time effect within viruliferous flies |
| `control_*` | `~ Time_Point` | 9 control only | Time effect within control flies |

**Output (available locally):** `03_analysis/06_time_contrasts/tables/` and `03_analysis/06_time_contrasts/plots/`

---

### Step 8 — QC Statistics Extraction
**Script:** `02_scripts/07_extract_stats.sh` → calls `02_scripts/07_extract_stats.R`
**What it does:** Parses fastp JSON reports, HISAT2 summaries, STAR logs, and featureCounts summary files into structured CSV tables for reporting.
**Output (available locally):** `04_results/tables/qc_stats/`

---

### Step 9 — QC Visualisation
**Script:** `02_scripts/08_plot_qc_stats.sh` → calls `02_scripts/08_plot_qc_stats.R`
**What it does:** Generates 12 publication-quality figures from the QC tables.
**Output (available locally):** `04_results/figures/qc_plots/`

---

## 4. Significance Cutoffs

Three cutoffs are used throughout and appear consistently in all file names:

| Label suffix | Criteria | Interpretation |
|---|---|---|
| `_lfc1` | padj < 0.05 AND \|log2FoldChange\| > 1 | Strict: ≥ 2-fold change |
| `_lfc058` | padj < 0.05 AND \|log2FoldChange\| > 0.58 | Moderate: ≥ 1.5-fold change |
| `_pval` | padj < 0.05 only | Permissive: any significant change regardless of magnitude |

When discussing results, use `lfc1` as the primary cutoff and reference others for context.

---

## 5. File System Guide — What is Local vs Cluster

### What lives on the cluster only (not needed for the report)
- Raw SRA files: `01_data/00_sra/`
- Raw and trimmed FASTQs: `01_data/01_raw_fastq/`, `01_data/02_trimmed_fastq/`
- Reference genomes and STAR index: `01_data/03_references/`
- QC reports (FastQC/MultiQC HTML): `03_analysis/01_qc/`, `03_analysis/02_fastp/`
- Decontamination FASTQs and BAMs: `03_analysis/03_hisat2_decontam/`
- STAR BAM files: `03_analysis/04_star_align/bam/`
- Raw featureCounts count matrix: `03_analysis/04_star_align/counts/`

### What is available locally (everything needed for the report)

#### QC Statistics Tables → `04_results/tables/qc_stats/`
| File | Contents |
|---|---|
| `01_trimming_stats.csv` | Per-sample fastp metrics: raw reads, clean reads, Q30%, GC%, duplication rate |
| `02_decontamination_stats.csv` | Per-sample HISAT2 alignment rates to the contaminant reference; retained unmapped pairs |
| `03_star_alignment_stats.csv` | Per-sample STAR alignment rates: uniquely mapped, multi-mapped, unmapped |
| `04_featurecounts_stats.csv` | Per-sample assignment rates: assigned, no feature, ambiguous, multi-mapping |
| `05_read_funnel.csv` | Read pairs at each stage: raw → post-trim → post-decontam → uniquely aligned; includes % retention at each step |
| `06_deg_summary.csv` | DEG counts (total/up/down) for every comparison × every cutoff; 30 rows total |
| `00_pipeline_summary.csv` | Wide-format merge of all per-sample metrics across all stages |

#### DESeq2 Results — Virus vs Control → `03_analysis/05_deseq2/tables/`
Files follow the naming pattern `deseq2_<comparison>_<type>.csv`:
- `<comparison>` is one of: `virus_vs_control`, `24h`, `48h`, `72h`
- `<type>` is one of: `all` (all tested genes), `sig_lfc1`, `sig_lfc058`, `sig_pval`
- Each `_all.csv` file contains columns: `GeneID`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`
- Additional stat files in this directory: `genes_tested.csv` (filter counts), `pca_variance.csv` (PC1/PC2 variance explained)

#### DESeq2 Results — Time Contrasts → `03_analysis/06_time_contrasts/tables/`
Files follow the naming pattern `<contrast>_<type>.csv`:
- `<contrast>` is one of: `overall_48h_vs_24h`, `overall_72h_vs_24h`, `virus_48h_vs_24h`, `virus_72h_vs_24h`, `control_48h_vs_24h`, `control_72h_vs_24h`
- `<type>` is one of: `all`, `sig_lfc1`, `sig_lfc058`, `sig_pval`
- Additional stat file: `genes_tested_A.csv`

#### Figures — QC Plots → `04_results/figures/qc_plots/`
12 figures, each available as PNG and PDF:
| File | What it shows |
|---|---|
| `01_pipeline_read_funnel` | Read retention across all 4 pipeline stages |
| `02_raw_read_counts` | Per-sample raw read counts |
| `03_trimming_filter_reasons` | Reads removed by fastp, broken down by reason |
| `04_q30_before_after` | Q30 quality score before vs after trimming |
| `05_gc_content_before_after` | GC content before vs after trimming |
| `06_decontam_rate_per_sample` | Endosymbiont contamination rate per sample |
| `07_decontam_stacked_bar` | Stacked: host-retained vs contaminant-mapped reads |
| `08_star_alignment_stacked` | Stacked: STAR alignment categories per sample |
| `09_star_uniquely_mapped_scatter` | Scatter plot of uniquely mapped read counts |
| `10_featurecounts_stacked` | Stacked: featureCounts assignment categories |
| `11_featurecounts_assigned_dot` | Fraction of reads assigned to annotated genes |
| `12_qc_summary_heatmap` | Heatmap of all scaled QC metrics across samples |

#### Figures — DESeq2 Main Analysis → `03_analysis/05_deseq2/plots/`
| File pattern | What it shows |
|---|---|
| `pca_plot.png/.pdf` | PCA of all 18 samples (VST-normalised) |
| `sample_distance_heatmap.png/.pdf` | Euclidean sample-to-sample distance heatmap |
| `volcano_virus_vs_control_<cutoff>.png/.pdf` | Volcano: overall virus vs control |
| `volcano_<tp>_<cutoff>.png/.pdf` | Volcano: virus vs control at 24h / 48h / 72h |
| `ma_virus_vs_control_<cutoff>.png/.pdf` | MA plot: overall virus vs control |

#### Figures — Time Contrast Analysis → `03_analysis/06_time_contrasts/plots/` and `04_results/figures/time_contrasts/`
| File pattern | What it shows |
|---|---|
| `<contrast>_volcano_<cutoff>.png/.pdf` | Volcano for each of the 6 time contrasts |
| `<contrast>_ma_<cutoff>.png/.pdf` | MA plot for each of the 6 time contrasts |
| `<contrast>_lfc1_heatmap.png/.pdf` | Top DEG heatmap (lfc1 cutoff) for each contrast |

---

## 6. Analytical Framework and Key Interpretive Points

### Two independent analyses, both important
The project runs two complementary DE analyses that answer different questions:
- **Analysis 1 (virus vs control):** Is there a transcriptional signature of virus acquisition? Tests each time point and overall.
- **Analysis 2 (time contrasts):** How does gene expression change over time? Run separately in the virus group, control group, and both combined.

Both analyses are of equal importance to the report. Do not treat one as secondary.

### What the DESeq2 design formulas mean
- `~ Time_Point + Group`: estimates the effect of Group (virus vs control) after accounting for the time point. This is the correct overall model because time point is a known source of variance.
- `~ Group + Time_Point`: estimates the effect of Time_Point after accounting for Group. Used for the overall time contrast analysis.
- `~ Group` (per-timepoint): estimates virus effect within a single time point. Simpler model because time is constant.
- `~ Time_Point` (per-group): estimates time effect within one group only.

### What the PCA means for this dataset
PCA was run on VST-normalised counts of all 18 samples. The proportion of variance explained by PC1 and PC2 is in `03_analysis/05_deseq2/tables/pca_variance.csv`. The separation pattern (whether samples cluster by time point or by group) directly indicates which factor drives the dominant transcriptional signal in the dataset.

### Endosymbiont contamination
A high proportion of reads in each sample map to the combined contaminant reference (endosymbionts + mitochondria). This is biologically expected — *Bemisia tabaci* harbours obligate endosymbionts at high abundance in bacteriocytes. The decontamination step is not a sign of poor library quality; it is a necessary biological filter. The post-decontamination read counts are the appropriate baseline for all downstream statistics.

### Gene ID format
All gene identifiers follow the pattern `BtaXXXXX` (e.g., `Bta01566`). These are MEAM1 gene model IDs from the scaffold v1.2 annotation. Functional annotation is outside the scope of this capstone; do not attempt to assign gene names or GO terms.

---

## 7. Future Work (Not Executed)

`97_trinity_pipeline/` contains scripts for an alternative reference-free approach:
1. De novo transcriptome assembly with Trinity (all 18 samples combined)
2. Functional annotation with Trinotate
3. Quantification with Salmon
4. DE analysis with tximport + DESeq2

Scripts are written and ready (`97_trinity_pipeline/02_scripts/`) but the pipeline was not run as part of the capstone. It is a planned next step. Do not reference it as completed work.

---

## 8. Report Writing Guide

### Suggested results chapter structure

| Section | Key figures | Key tables |
|---|---|---|
| Read-level QC | `01_pipeline_read_funnel`, `02_raw_read_counts`, `04_q30_before_after` | `01_trimming_stats.csv` |
| Decontamination | `06_decontam_rate_per_sample`, `07_decontam_stacked_bar` | `02_decontamination_stats.csv` |
| Alignment and quantification | `08_star_alignment_stacked`, `10_featurecounts_stacked` | `03_star_alignment_stats.csv`, `04_featurecounts_stats.csv` |
| Overall read funnel summary | `01_pipeline_read_funnel` | `05_read_funnel.csv` |
| Sample relationships (PCA + clustering) | `pca_plot`, `sample_distance_heatmap` | `pca_variance.csv` |
| Virus vs control DE | `volcano_virus_vs_control_*`, `ma_virus_vs_control_*`, per-timepoint volcanos | `deseq2_virus_vs_control_all.csv`, `06_deg_summary.csv` (rows where Source = virus_vs_control) |
| Time-point DE | `overall_*_volcano_*`, `virus_*_volcano_*`, `control_*_volcano_*`, heatmaps | `06_deg_summary.csv` (rows where Source = time_contrasts), individual `*_sig_lfc1.csv` tables |
| DEG summary overview | — | `06_deg_summary.csv` (all 30 rows) |

### Tone and framing notes
- The audience is a CiBIG supervisor who wants to see the student's thought process, not just the results
- Explain *why* each analytical decision was made, not just what was done
- The fact that the dominant signal is temporal (not viral) is itself a meaningful biological finding — frame it as such, not as a negative result
- When reporting DEG numbers, always state the cutoff used (`lfc1`, `lfc058`, or `pval`)
- All p-values are BH-adjusted FDR. Do not use the word "significant" without specifying the cutoff and that it is FDR-corrected

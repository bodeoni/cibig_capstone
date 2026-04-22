# CiBIG Capstone: Viruliferous vs Non-viruliferous Whiteflies RNA-seq

**Author:** Olabode Onile-ere

**Program:** Certificate in Bioinformatics and Genomics (CiBIG)

**Cluster:** WAVE HPC (SLURM scheduler)

**Working directory (cluster):** `/scratch/onilee/capstone` inside node06

---

## Biological Background

*Bemisia tabaci* (silverleaf whitefly, MEAM1 biotype) is a major agricultural pest and vector of begomoviruses (family *Geminiviridae*). This project investigates transcriptional changes in *B. tabaci* following acquisition of an Old World begomovirus (TYLCV-related), comparing viruliferous (virus-carrying) and non-viruliferous (control) flies across three time points post-acquisition access period (24h, 48h, 72h).

---

## Source Data

| Field | Value |
|---|---|
| Publication | "Gene expression differences in *Bemisia tabaci* following acquisition of an Old World begomovirus" |
| DOI | https://doi.org/10.1038/s41597-025-06417-3 |
| BioProject | [PRJNA1096732](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1096732) |
| Organism | *Bemisia tabaci* MEAM1 |
| Library type | Bulk RNA-seq, paired-end |
| Genome reference | MEAM1 scaffold v1.2 (whiteflygenomics.org) |

---

## Sample Metadata

18 samples total: 3 time points × 2 groups × 3 biological replicates.

| Run | Time Point | Group |
|---|---|---|
| SRR28578498 | 24h | virus |
| SRR28578499 | 72h | control |
| SRR28578500 | 48h | control |
| SRR28578501 | 24h | control |
| SRR28578502 | 72h | virus |
| SRR28578503 | 48h | virus |
| SRR28578504 | 24h | virus |
| SRR28578505 | 72h | control |
| SRR28578506 | 72h | virus |
| SRR28578507 | 48h | virus |
| SRR28578508 | 24h | virus |
| SRR28578509 | 72h | control |
| SRR28578510 | 48h | control |
| SRR28578511 | 24h | control |
| SRR28578512 | 72h | virus |
| SRR28578513 | 48h | virus |
| SRR28578514 | 48h | control |
| SRR28578515 | 24h | control |

Full metadata: `00_meta/sample_metadata.csv`

**Sample breakdown:**

| | 24h | 48h | 72h |
|---|---|---|---|
| virus | 3 | 3 | 3 |
| control | 3 | 3 | 3 |

---

## Pipeline Overview

The pipeline runs in 9 sequential steps. All SLURM scripts are in `02_scripts/` and should be submitted from the cluster working directory `/scratch/onilee/capstone`.

```
Raw SRA files
    │
    ▼ 00-extract_fastq.sh
Raw FASTQ (paired-end)
    │
    ▼ 01_qc_fastq.sh
FastQC / MultiQC reports
    │
    ▼ 02_trim_fastp.sh
Trimmed FASTQ
    │
    ▼ 03b_decontam_hisat2.sh
Decontaminated FASTQ (endosymbionts + mitochondria removed)
    │
    ▼ 04_star_align_counts.sh
STAR BAM files + featureCounts gene count matrix
    │
    ▼ 05_deseq2.sh  (→ 05_deseq2.R)
Differential expression: virus vs control (overall + per time point)
    │
    ▼ 06_time_contrasts.sh  (→ 06_time_contrasts.R)
Differential expression: time-point effects (overall, virus-only, control-only)
    │
    ▼ 07_extract_stats.sh  (→ 07_extract_stats.R)
QC statistics tables (trimming, decontamination, alignment, DEG summary)
    │
    ▼ 08_plot_qc_stats.sh  (→ 08_plot_qc_stats.R)
QC visualisation plots (12 figures)
```

---

## Step-by-Step Instructions

### Step 0: Prerequisites

Ensure the following are available on the cluster before running any step:

**Conda environments:**
- `bioinfo` — contains HISAT2 v2.2.2, samtools v1.23
- `rnaseq` — contains R with DESeq2 v1.50+, ggplot2, pheatmap, RColorBrewer

**Apptainer containers** (in `/projects/onilee/software/containers/`):
- `star_2.7.11b.sif` — STAR aligner
- `subread_2.1.1.sif` — featureCounts

**Modules required:** `bioinfo-wave`, `miniconda3`, `apptainer`

---

### Step 1: Extract FASTQ from SRA

```bash
sbatch 02_scripts/00-extract_fastq.sh
```

- Tool: `fasterq-dump` (SRA Toolkit)
- Input: `01_data/00_sra/SRR*/`
- Output: `01_data/01_raw_fastq/SRR*/{SRR*_1.fastq, SRR*_2.fastq}`
- Resources: 16 CPUs, 64G RAM

---

### Step 2: Raw Read QC

```bash
sbatch 02_scripts/01_qc_fastq.sh
```

- Tools: FastQC v0.12.1, MultiQC v1.9
- Input: `01_data/01_raw_fastq/SRR*/`
- Output:
  - `03_analysis/01_qc/a_fastqc/` — per-sample FastQC reports
  - `03_analysis/01_qc/b_multiqc/multiqc_raw_fastqc.html` — summary report
- Resources: 16 CPUs

---

### Step 3: Adapter Trimming and Post-trim QC

```bash
sbatch 02_scripts/02_trim_fastp.sh
```

- Tool: fastp v0.20.1
- Input: `01_data/01_raw_fastq/SRR*/`
- Output:
  - `01_data/02_trimmed_fastq/SRR*_1.clean.fastq.gz` / `SRR*_2.clean.fastq.gz`
  - `03_analysis/02_fastp/` — HTML/JSON reports and MultiQC summary
- Key parameters: `--detect_adapter_for_pe`, `--cut_front`, `--cut_tail`, `--cut_window_size 4`, `--cut_mean_quality 20`, `--length_required 75`
- Resources: 24 CPUs, 100G RAM

---

### Step 4: Decontamination (HISAT2)

```bash
sbatch 02_scripts/03b_decontam_hisat2.sh
```

Removes reads mapping to endosymbionts (Portiera, Hamiltonella, Rickettsia) and mitochondria. Only unmapped (host-derived) reads are retained.

- Tool: HISAT2 v2.2.2 (conda: `bioinfo`)
- Reference: `01_data/03_references/MEAM1_contiminants.fa` (combined contaminant reference)
- Input: `01_data/02_trimmed_fastq/`
- Output:
  - `03_analysis/03_hisat2_decontam/unmapped_fastq/<sample>_unmapped.fastq.gz.1` / `.2`
  - `03_analysis/03_hisat2_decontam/bam/` — alignments to contaminants
  - `03_analysis/03_hisat2_decontam/summaries/` — per-sample alignment summaries
- Key parameters: `--very-sensitive`, `--un-conc-gz`
- Resources: 24 CPUs, 100G RAM

> **Reference downloads** (`03a_download_genomes.sh`): MEAM1 genome and annotation from whiteflygenomics.org; endosymbiont genomes from the same source; mitochondrial genome (NC_006279.1) from NCBI.

---

### Step 5: STAR Alignment and featureCounts Quantification

```bash
sbatch 02_scripts/04_star_align_counts.sh
```

- Tools: STAR v2.7.11b (apptainer), featureCounts / Subread v2.1.1 (apptainer)
- Input: decontaminated FASTQ from Step 4
- Reference genome: `01_data/03_references/MEAM1_scaffold_v1.2.fa`
- Reference annotation: `01_data/03_references/MEAM1_v1.2.gff3`
- STAR index: `01_data/03_references/star_index/` (built on first run; skipped if directory exists)
- Output:
  - `03_analysis/04_star_align/bam/` — coordinate-sorted BAM files
  - `03_analysis/04_star_align/counts/all_samples_genelevel.txt` — combined count matrix
  - `03_analysis/04_star_align/counts/<sample>_genelevel.txt` — per-sample counts
  - `03_analysis/04_star_align/logs/` — STAR alignment logs
- Key parameters:
  - STAR: `--outSAMtype BAM SortedByCoordinate`, `--sjdbGTFtagExonParentTranscript Parent`, `--genomeSAindexNbases 12`, `--sjdbOverhang 99`
  - featureCounts: `-p` (paired-end), `-t exon`, `-g Parent`
- Resources: 24 CPUs, 64G RAM

> **Note on featureCounts:** `-t exon -g Parent` counts reads over exonic regions only. The `Parent` attribute in the MEAM1 GFF3 links each exon to its transcript ID (e.g., `Bta00001-mRNA`). Genes in this annotation are predominantly single-isoform, so transcript-level counts are effectively gene-level.

---

### Step 6: Differential Expression — Virus vs Control

```bash
sbatch 02_scripts/05_deseq2.sh
```

Calls `05_deseq2.R` via the `rnaseq` conda environment.

- Tool: DESeq2 (R)
- Input: `03_analysis/04_star_align/counts/all_samples_genelevel.txt`, `00_meta/sample_metadata.csv`
- Design: `~ Time_Point + Group` (overall); `~ Group` per time-point subset
- Low-count filter: genes with < 10 total counts removed
- Output directory: `03_analysis/05_deseq2/`

**Tables** (in `tables/`):

| File | Contents |
|---|---|
| `deseq2_virus_vs_control_all.csv` | All tested genes (overall) |
| `deseq2_virus_vs_control_sig_lfc1.csv` | padj < 0.05 & \|LFC\| > 1 |
| `deseq2_virus_vs_control_sig_lfc058.csv` | padj < 0.05 & \|LFC\| > 0.58 |
| `deseq2_virus_vs_control_sig_pval.csv` | padj < 0.05 only |
| `deseq2_<tp>_all.csv` | All genes, per time point (24h/48h/72h) |
| `deseq2_<tp>_sig_lfc1/lfc058/pval.csv` | Significant DEGs per time point |

**Plots** (in `plots/`): PCA, sample distance heatmap, volcano plots × 3 cutoffs, MA plots × 3 cutoffs, top DEG heatmap. Also produced per time point: 3 volcano plots.

---

### Step 7: Differential Expression — Time-point Contrasts

```bash
sbatch 02_scripts/06_time_contrasts.sh
```

Calls `06_time_contrasts.R` via the `rnaseq` conda environment.

Three analyses, each with contrasts 48h vs 24h and 72h vs 24h:

| Analysis | Design | Samples |
|---|---|---|
| A. Overall time effect | `~ Group + Time_Point` | All 18 |
| B. Time effect in virus group | `~ Time_Point` | 9 virus only |
| C. Time effect in control group | `~ Time_Point` | 9 control only |

- Output directory: `03_analysis/06_time_contrasts/`
- Per contrast: 4 tables (`_all`, `_sig_lfc1`, `_sig_lfc058`, `_sig_pval`), 3 volcano plots, 3 MA plots, 1 heatmap (top DEGs at lfc1 cutoff)

---

### Step 8: QC Statistics Extraction

```bash
sbatch 02_scripts/07_extract_stats.sh
```

Calls `07_extract_stats.R` via the `rnaseq` conda environment. Parses outputs from all upstream steps and consolidates them into summary tables.

- Input: fastp JSON reports, HISAT2 summaries, STAR logs, featureCounts summary, DESeq2 tables
- Output: `04_results/tables/qc_stats/`

| File | Contents |
|---|---|
| `01_trimming_stats.csv` | Per-sample fastp metrics (reads, Q30, GC, duplication rate) |
| `02_decontamination_stats.csv` | Per-sample HISAT2 alignment to contaminant reference |
| `03_star_alignment_stats.csv` | Per-sample STAR alignment statistics |
| `04_featurecounts_stats.csv` | Per-sample featureCounts assignment rates |
| `05_read_funnel.csv` | Read pairs at each pipeline stage (raw → trimmed → decontaminated → aligned) |
| `06_deg_summary.csv` | DEG counts (total / up / down) per comparison and significance cutoff |
| `00_pipeline_summary.csv` | Merged per-sample summary across all steps |

---

### Step 9: QC Visualisation

```bash
sbatch 02_scripts/08_plot_qc_stats.sh
```

Calls `08_plot_qc_stats.R` via the `rnaseq` conda environment. Generates 12 publication-quality figures from the QC statistics tables.

- Input: tables from Step 8 (`04_results/tables/qc_stats/`)
- Output: `04_results/figures/qc_plots/` (PNG + PDF for each figure)

| Figure | Description |
|---|---|
| `01_pipeline_read_funnel` | Read pairs retained at each pipeline stage |
| `02_raw_read_counts` | Per-sample raw read counts |
| `03_trimming_filter_reasons` | Breakdown of reads removed by fastp |
| `04_q30_before_after` | Q30 percentage before vs after trimming |
| `05_gc_content_before_after` | GC content before vs after trimming |
| `06_decontam_rate_per_sample` | Endosymbiont contamination rate per sample |
| `07_decontam_stacked_bar` | Stacked bar: host vs contaminant reads |
| `08_star_alignment_stacked` | Stacked bar: STAR alignment categories |
| `09_star_uniquely_mapped_scatter` | Uniquely mapped reads per sample (scatter) |
| `10_featurecounts_stacked` | Stacked bar: featureCounts assignment categories |
| `11_featurecounts_assigned_dot` | Fraction of reads assigned to genes (dot plot) |
| `12_qc_summary_heatmap` | Heatmap of all QC metrics across samples |

---

## Significance Cutoffs Used Throughout

| Label | Criteria | Interpretation |
|---|---|---|
| `lfc1` | padj < 0.05 & \|LFC\| > 1 | 2-fold change |
| `lfc058` | padj < 0.05 & \|LFC\| > 0.58 | 1.5-fold change |
| `pval` | padj < 0.05 only | No fold-change filter |

All p-values are Benjamini-Hochberg adjusted. Reference level: `control` for Group; `24h` for Time_Point.

---

## Project Structure

```
cibig_capstone/
├── 00_meta/
│   ├── methods.md              # Detailed methods for each pipeline step
│   ├── sample_metadata.csv     # Sample-level metadata (Run, Time_Point, Group)
│   └── software_versions.txt   # Tool versions logged during each run
├── 01_data/
│   ├── 00_sra/                 # Downloaded SRA files (SRR*/)
│   ├── 01_raw_fastq/           # Raw paired FASTQ (SRR*/)
│   ├── 02_trimmed_fastq/       # fastp-trimmed FASTQ
│   └── 03_references/          # Genome FASTA, GFF3, STAR index, contaminant reference
├── 02_scripts/
│   ├── 00-extract_fastq.sh     # SRA → FASTQ
│   ├── 01_qc_fastq.sh          # FastQC + MultiQC
│   ├── 02_trim_fastp.sh        # fastp trimming
│   ├── 03a_download_genomes.sh # Download reference genomes
│   ├── 03b_decontam_hisat2.sh  # HISAT2 decontamination
│   ├── 04_star_align_counts.sh # STAR alignment + featureCounts
│   ├── 05_deseq2.sh / .R       # DESeq2: virus vs control
│   ├── 06_time_contrasts.sh / .R  # DESeq2: time-point contrasts
│   ├── 07_extract_stats.sh / .R   # QC statistics extraction
│   └── 08_plot_qc_stats.sh / .R   # QC visualisation (12 figures)
├── 03_analysis/
│   ├── 01_qc/                  # FastQC and MultiQC outputs
│   ├── 02_fastp/               # fastp reports
│   ├── 03_hisat2_decontam/     # Decontamination outputs
│   ├── 04_star_align/          # BAM files and count matrix
│   ├── 05_deseq2/              # DESeq2 results — virus vs control (tables + plots)
│   └── 06_time_contrasts/      # DESeq2 results — time-point contrasts (tables + plots)
├── 04_results/
│   ├── figures/
│   │   ├── qc_plots/           # 12 QC figures (PNG + PDF)
│   │   └── time_contrasts/     # Time-point contrast plots
│   ├── tables/
│   │   ├── qc_stats/           # QC statistics CSVs (steps 8–9 output)
│   │   └── time_contrasts/     # Significant DEG tables for time contrasts
│   └── reports/                # Alignment summary reports
├── 97_trinity_pipeline/        # Future work: Trinity de novo assembly pipeline (scripts prepared)
├── 98_test_run/                # Preliminary work: approach comparison on 3-sample subset
└── 99_logs/                    # SLURM log files
```

---

## Software Versions

| Tool | Version | Used in | Execution |
|---|---|---|---|
| fasterq-dump (SRA Toolkit) | — | Step 1 | module |
| FastQC | v0.12.1 | Step 2 | module |
| MultiQC | v1.9 | Steps 2–3 | module |
| fastp | v0.20.1 | Step 3 | module |
| HISAT2 | v2.2.2 | Step 4 | conda: bioinfo |
| samtools | v1.23 | Step 4 | conda: bioinfo |
| STAR | v2.7.11b | Step 5 | apptainer |
| featureCounts (Subread) | v2.1.1 | Step 5 | apptainer |
| DESeq2 | v1.50+ | Steps 6–7 | conda: rnaseq |
| ggplot2 | — | Steps 6–9 | conda: rnaseq |
| pheatmap | — | Steps 6–9 | conda: rnaseq |
| RColorBrewer | — | Steps 6–9 | conda: rnaseq |
| ggrepel | — | Steps 6–9 | conda: rnaseq |
| scales | — | Steps 8–9 | conda: rnaseq |

Full version log: `00_meta/software_versions.txt`

---

## Reproducibility Notes

- All SLURM scripts use `set -euo pipefail` — the job will abort on any error.
- Samples are discovered dynamically from directory contents (no hardcoded lists) in Steps 5–7.
- The STAR index is built once and reused. If you change the reference genome, delete `01_data/03_references/star_index/` to trigger a rebuild.
- R scripts are copied into their output directory at runtime for a snapshot of the exact code used.
- `session_info.txt` is written at the end of each R script recording exact R and package versions.
- All significance testing uses Benjamini-Hochberg FDR correction (`pAdjustMethod = "BH"`).
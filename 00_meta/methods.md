# Methods

## Overview

Bulk RNA-seq data from viruliferous and non-viruliferous *Bemisia tabaci* (MEAM1 biotype) were processed to identify differentially expressed genes following begomovirus acquisition. The pipeline covers raw data extraction, quality control, trimming, endosymbiont decontamination, alignment, quantification, and differential expression analysis. All steps are run on the Negishi HPC cluster using SLURM. Containerised tools are executed via Apptainer; R-based steps use a dedicated conda environment.

---

## Step 1: SRA to FASTQ extraction

Script: [02_scripts/00-extract_fastq.sh](../02_scripts/00-extract_fastq.sh)

- Tool: `fasterq-dump` (NCBI SRA Toolkit)
- Input: SRA directories under `01_data/00_sra/SRR*/`
- Output: Paired FASTQ files under `01_data/01_raw_fastq/SRR*/`
- Parameters:
  - `--split-3` — split paired-end reads into R1 and R2
  - `--threads 16`
  - `--outdir` set per sample
  - `--progress` enabled
  - `--temp /lscratch/onilee/tmp_space`
- SLURM resources: 16 CPUs, 64G RAM, partition `normal`, node `node06`

---

## Step 2: Raw read QC

Script: [02_scripts/01_qc_fastq.sh](../02_scripts/01_qc_fastq.sh)

- Tools: FastQC v0.12.1, MultiQC v1.9
- Input: Raw FASTQ files in `01_data/01_raw_fastq/SRR*/`
- Output:
  - Per-sample FastQC reports: `03_analysis/01_qc/a_fastqc/SRR*/`
  - Combined MultiQC report: `03_analysis/01_qc/b_multiqc/multiqc_raw_fastqc.html`
- Parameters:
  - `fastqc -t $SLURM_CPUS_PER_TASK -o <outdir> <R1> <R2>`
  - `multiqc -o <multiqc_dir> -n multiqc_raw_fastqc.html <fastqc_dir>`
- SLURM resources: 16 CPUs, partition `normal`, node `node06`

---

## Step 3: Adapter trimming and post-trim QC

Script: [02_scripts/02_trim_fastp.sh](../02_scripts/02_trim_fastp.sh)

- Tools: fastp v0.20.1, MultiQC v1.9
- Input: Raw FASTQ files in `01_data/01_raw_fastq/SRR*/`
- Output:
  - Trimmed FASTQs: `01_data/02_trimmed_fastq/SRR*_1.clean.fastq.gz`, `SRR*_2.clean.fastq.gz`
  - fastp HTML/JSON reports: `03_analysis/02_fastp/4-fastp_reports/`
  - Post-trim MultiQC summary: `03_analysis/02_fastp/5-multiqc_post/`
- Parameters:
  - `--detect_adapter_for_pe` — automatic adapter detection for paired-end data
  - `--cut_front --cut_tail` — sliding window quality trimming from both ends
  - `--cut_window_size 4 --cut_mean_quality 20` — window size and quality threshold
  - `--length_required 75` — discard reads shorter than 75 bp post-trimming
  - `--thread $SLURM_CPUS_PER_TASK`
- SLURM resources: 24 CPUs, 100G RAM, partition `normal`, node `node06`

---

## Step 4a: Reference genome download

Script: [02_scripts/03a_download_genomes.sh](../02_scripts/03a_download_genomes.sh)

Downloads all reference sequences to `01_data/03_references/`:

| File | Source | Description |
|---|---|---|
| `MEAM1_scaffold_v1.2.fa.gz` | whiteflygenomics.org | MEAM1 whitefly genome |
| `MEAM1_v1.2.gff3.gz` | whiteflygenomics.org | Gene annotation (GFF3) |
| `NC_006279.1.fa` | NCBI | Mitochondrial genome |
| `Portiera_meam1_genome_v2.0.fa` | whiteflygenomics.org | Portiera endosymbiont |
| `Hamiltonella_meam1_genome_v2_0.fa` | whiteflygenomics.org | Hamiltonella endosymbiont |
| `Rickettsia_meam1_genome_v2.0.fa` | whiteflygenomics.org | Rickettsia endosymbiont |
| `MEAM1_contiminants.fa` | merged locally | All contaminant genomes combined |

---

## Step 4b: Endosymbiont decontamination (HISAT2)

Script: [02_scripts/03b_decontam_hisat2.sh](../02_scripts/03b_decontam_hisat2.sh)

Trimmed reads are aligned to a combined contaminant reference (mitochondria + three endosymbionts). Reads that do **not** map (i.e., putative host-derived reads) are retained for downstream analysis.

- Tool: HISAT2 v2.2.2, samtools v1.23 (conda environment: `bioinfo`)
- Input: `01_data/02_trimmed_fastq/SRR*_1.clean.fastq.gz`, `SRR*_2.clean.fastq.gz`
- Reference index: built from `01_data/03_references/MEAM1_contiminants.fa`
- Output:
  - Decontaminated FASTQ: `03_analysis/03_hisat2_decontam/unmapped_fastq/<sample>_unmapped.fastq.gz.1` / `.2`
  - Contaminant BAMs: `03_analysis/03_hisat2_decontam/bam/<sample>.sorted.bam`
  - Per-sample alignment summaries: `03_analysis/03_hisat2_decontam/summaries/<sample>.hisat2.summary.txt`
- Key parameters:
  - `--very-sensitive` — highest sensitivity preset
  - `--no-unal` — suppress unaligned reads from SAM output
  - `--un-conc-gz <prefix>` — write unmapped read pairs to gzipped FASTQ
- SLURM resources: 24 CPUs, 100G RAM, partition `normal`, node `node06`

---

## Step 5: STAR alignment and featureCounts quantification

Script: [02_scripts/04_star_align_counts.sh](../02_scripts/04_star_align_counts.sh)

Decontaminated reads are aligned to the MEAM1 host genome using STAR and quantified with featureCounts.

- Tools: STAR v2.7.11b (apptainer), featureCounts/Subread v2.1.1 (apptainer)
- Input: `03_analysis/03_hisat2_decontam/unmapped_fastq/<sample>_unmapped.fastq.gz.1` / `.2`
- Reference genome: `01_data/03_references/MEAM1_scaffold_v1.2.fa`
- Reference annotation: `01_data/03_references/MEAM1_v1.2.gff3`
- STAR index: `01_data/03_references/star_index/` (built on first run if directory absent)
- Output:
  - BAM files: `03_analysis/04_star_align/bam/<sample>_Aligned.sortedByCoord.out.bam`
  - Combined count matrix: `03_analysis/04_star_align/counts/all_samples_genelevel.txt`
  - Per-sample counts: `03_analysis/04_star_align/counts/<sample>_genelevel.txt`
  - STAR logs: `03_analysis/04_star_align/logs/`
- STAR parameters:
  - `--outSAMtype BAM SortedByCoordinate`
  - `--sjdbGTFtagExonParentTranscript Parent` — GFF3 uses `Parent` to link exons to transcripts
  - `--genomeSAindexNbases 12` — reduced suffix array index for the relatively small MEAM1 genome
  - `--sjdbOverhang 99` — read length minus 1 for splice junction detection
  - `--readFilesCommand zcat` — decompress gzipped input on the fly
- featureCounts parameters:
  - `-p` — paired-end mode
  - `-t exon` — count reads overlapping exon features only (avoids inflating counts with intronic reads)
  - `-g Parent` — group exons by transcript ID (the `Parent` attribute in MEAM1 GFF3 resolves to transcript IDs such as `Bta00001-mRNA`; genes are predominantly single-isoform so transcript counts are effectively gene-level)
- SLURM resources: 24 CPUs, 64G RAM, partition `normal`, node `node06`

---

## Step 6: Differential expression — virus vs control

Scripts: [02_scripts/05_deseq2.sh](../02_scripts/05_deseq2.sh), [02_scripts/05_deseq2.R](../02_scripts/05_deseq2.R)

- Tool: DESeq2 (R), executed via conda environment `rnaseq`
- Input: `03_analysis/04_star_align/counts/all_samples_genelevel.txt`, `00_meta/sample_metadata.csv`
- Output: `03_analysis/05_deseq2/`

### Statistical design

| Analysis | DESeq2 design | Contrast |
|---|---|---|
| Overall virus effect | `~ Time_Point + Group` | virus vs control |
| 24h subset | `~ Group` | virus vs control |
| 48h subset | `~ Group` | virus vs control |
| 72h subset | `~ Group` | virus vs control |

- Reference levels: `control` (Group), `24h` (Time_Point)
- Low-count pre-filter: genes with < 10 total counts across all samples removed before fitting
- P-value adjustment: Benjamini-Hochberg (BH) FDR

### Significance cutoffs applied

| Label | Criteria |
|---|---|
| `lfc1` | padj < 0.05 & \|log2FoldChange\| > 1 (≥ 2-fold) |
| `lfc058` | padj < 0.05 & \|log2FoldChange\| > 0.58 (≥ 1.5-fold) |
| `pval` | padj < 0.05 (no fold-change filter) |

### Outputs per analysis

- Tables: `_all.csv` (all genes), `_sig_lfc1.csv`, `_sig_lfc058.csv`, `_sig_pval.csv`
- Plots (overall only): PCA, sample distance heatmap, top DEG heatmap
- Plots (all analyses): volcano plots × 3 cutoffs, MA plots × 3 cutoffs

- SLURM resources: 4 CPUs, 32G RAM, partition `normal`, node `node06`

---

## Step 7: Differential expression — time-point contrasts

Scripts: [02_scripts/06_time_contrasts.sh](../02_scripts/06_time_contrasts.sh), [02_scripts/06_time_contrasts.R](../02_scripts/06_time_contrasts.R)

- Tool: DESeq2 (R), executed via conda environment `rnaseq`
- Input: same count matrix and metadata as Step 6
- Output: `03_analysis/06_time_contrasts/`

### Statistical design

Three parallel analyses, each testing 48h vs 24h and 72h vs 24h (24h as reference):

| Analysis | DESeq2 design | Samples used |
|---|---|---|
| A. Overall time effect | `~ Group + Time_Point` | All 18 (Group as covariate) |
| B. Virus group time effect | `~ Time_Point` | 9 virus samples only |
| C. Control group time effect | `~ Time_Point` | 9 control samples only |

- For analyses B and C, `fitType = "local"` is used as a safeguard for the smaller per-group sample size
- Row-level count filter applied per subset: genes with < 10 total counts within that subset removed

### Contrasts and output labels

| Label | Comparison |
|---|---|
| `overall_48h_vs_24h` | 48h vs 24h, all samples |
| `overall_72h_vs_24h` | 72h vs 24h, all samples |
| `virus_48h_vs_24h` | 48h vs 24h, virus group |
| `virus_72h_vs_24h` | 72h vs 24h, virus group |
| `control_48h_vs_24h` | 48h vs 24h, control group |
| `control_72h_vs_24h` | 72h vs 24h, control group |

### Outputs per contrast

- Tables: `<label>_all.csv`, `<label>_sig_lfc1.csv`, `<label>_sig_lfc058.csv`, `<label>_sig_pval.csv`
- Plots: volcano plots × 3 cutoffs, MA plots × 3 cutoffs, top DEG heatmap (lfc1 cutoff)

- SLURM resources: 4 CPUs, 32G RAM, partition `normal`, node `node06`

---

## Step 8: QC statistics extraction

Scripts: [02_scripts/07_extract_stats.sh](../02_scripts/07_extract_stats.sh), [02_scripts/07_extract_stats.R](../02_scripts/07_extract_stats.R)

- Tool: R, executed via conda environment `rnaseq`
- Input: fastp JSON reports, HISAT2 alignment summaries, STAR logs, featureCounts summary file, DESeq2 sig tables
- Output: `04_results/tables/qc_stats/`

Parses all upstream pipeline outputs and consolidates them into structured CSV tables for reporting and visualisation:

| Output file | Contents |
|---|---|
| `01_trimming_stats.csv` | Per-sample fastp metrics including read counts, Q30 %, GC %, duplication rate |
| `02_decontamination_stats.csv` | Per-sample HISAT2 alignment rates to the combined contaminant reference |
| `03_star_alignment_stats.csv` | Per-sample STAR alignment statistics (uniquely mapped, multi-mapped, unmapped) |
| `04_featurecounts_stats.csv` | Per-sample featureCounts assignment rates |
| `05_read_funnel.csv` | Read pairs remaining after each pipeline stage (raw → trimmed → decontaminated → uniquely aligned) |
| `06_deg_summary.csv` | DEG counts (total, up-regulated, down-regulated) per comparison and cutoff |
| `00_pipeline_summary.csv` | Wide-format merge of all per-sample metrics |

---

## Step 9: QC visualisation

Scripts: [02_scripts/08_plot_qc_stats.sh](../02_scripts/08_plot_qc_stats.sh), [02_scripts/08_plot_qc_stats.R](../02_scripts/08_plot_qc_stats.R)

- Tool: R (ggplot2, pheatmap, ggrepel, scales, grid), executed via conda environment `rnaseq`
- Input: tables from Step 8
- Output: `04_results/figures/qc_plots/` — PNG and PDF versions of each figure

Generates 12 publication-quality figures covering every stage of the pipeline:

| Figure | Type | Description |
|---|---|---|
| `01_pipeline_read_funnel` | Bar | Read pairs retained at each stage |
| `02_raw_read_counts` | Bar | Per-sample raw read counts |
| `03_trimming_filter_reasons` | Stacked bar | Reads removed by fastp by reason |
| `04_q30_before_after` | Paired dot | Q30 % before vs after trimming |
| `05_gc_content_before_after` | Paired dot | GC % before vs after trimming |
| `06_decontam_rate_per_sample` | Bar | Endosymbiont mapping rate per sample |
| `07_decontam_stacked_bar` | Stacked bar | Host vs contaminant read proportions |
| `08_star_alignment_stacked` | Stacked bar | STAR alignment category proportions |
| `09_star_uniquely_mapped_scatter` | Scatter | Uniquely mapped read count per sample |
| `10_featurecounts_stacked` | Stacked bar | featureCounts assignment categories |
| `11_featurecounts_assigned_dot` | Dot | Fraction of reads assigned to genes |
| `12_qc_summary_heatmap` | Heatmap | All QC metrics scaled across samples |

Samples are ordered by time point then group (control before virus) throughout. A manifest (`00_plot_manifest.csv`) is written on completion.

- SLURM resources: 4 CPUs, 16G RAM, partition `normal`, node `node06`

---

## Software versions

All module versions are appended to [00_meta/software_versions.txt](software_versions.txt) at runtime during each step.

| Tool | Version | Steps | Environment |
|---|---|---|---|
| fasterq-dump | — | 1 | module |
| FastQC | v0.12.1 | 2 | module |
| MultiQC | v1.9 | 2–3 | module |
| fastp | v0.20.1 | 3 | module |
| HISAT2 | v2.2.2 | 4 | conda: bioinfo |
| samtools | v1.23 | 4 | conda: bioinfo |
| STAR | v2.7.11b | 5 | apptainer |
| featureCounts (Subread) | v2.1.1 | 5 | apptainer |
| DESeq2 | v1.50+ | 6–7 | conda: rnaseq |
| ggplot2 | — | 6–9 | conda: rnaseq |
| pheatmap | — | 6–9 | conda: rnaseq |
| RColorBrewer | — | 6–9 | conda: rnaseq |
| ggrepel | — | 6–9 | conda: rnaseq |
| scales | — | 8–9 | conda: rnaseq |

---

## Future work: Trinity de novo assembly pipeline

Scripts prepared in `97_trinity_pipeline/02_scripts/` but not executed as part of the capstone. The pipeline was designed as an alternative approach that avoids dependence on the MEAM1 reference genome by assembling a transcriptome directly from the experimental reads.

Planned pipeline:
1. `01_trinity_assembly.sh` — de novo assembly of all 18 decontaminated samples with Trinity; `--normalize_reads --normalize_max_cov 50`
2. `02_trinotate.sh` — functional annotation of the assembled transcriptome with Trinotate
3. `03_salmon_index.sh` — build a Salmon index from the Trinity assembly
4. `04_salmon_quant.sh` — quantify each sample against the assembly with Salmon
5. `05_deseq2.R` — differential expression analysis using tximport + DESeq2

This approach would allow comparison of reference-based and reference-free DE results, and would enable functional annotation of differentially expressed transcripts. It remains as a planned next step beyond the capstone submission.

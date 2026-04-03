# Methods

## Overview
This project processes bulk RNA-seq data from viruliferous and non-viruliferous flies. The workflow includes SRA extraction, read-level QC, and adapter/quality trimming. All steps are run on the Negishi cluster using SLURM.

## Step 1: SRA to FASTQ extraction
Script: [02_scripts/00-extract_fastq.sh](02_scripts/00-extract_fastq.sh)

- Tool: `fasterq-dump`
- Input: SRA directories under `01_data/00_sra/SRR*/`
- Output: Paired FASTQ files under `01_data/01_raw_fastq/SRR*/`
- Parameters:
	- `--split-3` to split paired-end reads
	- `--threads 16` for parallel extraction
	- `--outdir` set per sample directory
	- `--progress` enabled
	- `--temp` set to `/lscratch/onilee/tmp_space`
- SLURM resources: 16 CPUs, 64G RAM, partition `normal`, node `node06`

## Step 2: Raw read QC
Script: [02_scripts/01_qc_fastq.sh](02_scripts/01_qc_fastq.sh)

- Tools: `FastQC` v0.12.1, `MultiQC` v1.9
- Input: Raw FASTQ files in `01_data/01_raw_fastq/SRR*/`
- Output:
	- Per-sample FastQC reports in `03_analysis/01_qc/a_fastqc/SRR*/`
	- Combined MultiQC report in `03_analysis/01_qc/b_multiqc/multiqc_raw_fastqc.html`
- Parameters:
	- `fastqc -t $SLURM_CPUS_PER_TASK -o <outdir> <R1> <R2>`
	- `multiqc -o <multiqc_dir> -n multiqc_raw_fastqc.html <fastqc_dir>`
- SLURM resources: 16 CPUs, partition `normal`, node `node06`

## Step 3: Trimming and post-trim QC
Script: [02_scripts/02_trim_fastp.sh](02_scripts/02_trim_fastp.sh)

- Tools: `fastp` v0.20.1, `MultiQC` v1.9
- Input: Raw FASTQ files in `01_data/01_raw_fastq/SRR*/`
- Output:
	- Trimmed FASTQs in `01_data/02_trimmed_fastq/`
	- fastp HTML/JSON reports in `03_analysis/02_fastp/4-fastp_reports/`
	- MultiQC summary in `03_analysis/02_fastp/5-multiqc_post/`
- Parameters:
	- Adapter detection: `--detect_adapter_for_pe`
	- Quality trimming: `--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20`
	- Minimum length: `--length_required 75`
	- Threads: `--thread $SLURM_CPUS_PER_TASK`
	- Reports: `--html <sample>.fastp.html --json <sample>.fastp.json`
- SLURM resources: 24 CPUs, partition `normal`, node `node06`

## Step 4 (cont.): STAR alignment and featureCounts quantification
Script: [02_scripts/04_star_align_counts.sh](02_scripts/04_star_align_counts.sh)

- Tools: STAR v2.7.11 (apptainer), featureCounts / Subread v2.1.1 (apptainer)
- Input: HISAT2-decontaminated unmapped reads in `03_analysis/03_hisat2_decontam/unmapped_fastq/`
  - File naming: `<sample>_unmapped.fastq.gz.1` / `<sample>_unmapped.fastq.gz.2` (from `--un-conc-gz`)
- Reference:
  - Genome: `01_data/03_references/MEAM1_scaffold_v1.2.fa`
  - Annotation: `01_data/03_references/MEAM1_v1.2.gff3`
  - STAR index: `01_data/03_references/star_index/`
- STAR parameters:
  - `--outSAMtype BAM SortedByCoordinate`
  - `--sjdbGTFtagExonParentTranscript Parent`
  - `--genomeSAindexNbases 12` (reduced for small genome)
  - `--sjdbOverhang 99`
- featureCounts parameters:
  - `-p` (paired-end)
  - `-t mRNA` (feature type)
  - `-g Parent` (attribute for gene ID)
- Output:
  - BAM files: `03_analysis/04_star_align/bam/`
  - Combined counts: `03_analysis/04_star_align/counts/all_samples_genelevel.txt`
  - Per-sample counts: `03_analysis/04_star_align/counts/<sample>_genelevel.txt`
  - STAR logs: `03_analysis/04_star_align/logs/`
- Samples: discovered dynamically from decontam output directory (no hardcoding)
- SLURM resources: 24 CPUs, 64G RAM, partition `normal`, node `node06`

## Software versions
All module versions are appended to [00_meta/software_versions.txt](00_meta/software_versions.txt) during each step.
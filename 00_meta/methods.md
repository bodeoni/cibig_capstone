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

## Software versions
All module versions are appended to [00_meta/software_versions.txt](00_meta/software_versions.txt) during each step.
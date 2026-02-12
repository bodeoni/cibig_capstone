# CiBIG Capstone: Viruliferous vs Non-viruliferous Flies RNA-seq

**Author:** Olabode Onile-ere

**Program:** Certificate in Bioinformatics and Genomics (CiBIG)

## Goal
Identify differentially expressed genes between viruliferous and non-viruliferous flies using bulk RNA-seq.

## Samples
Raw sample identifiers are stored under `01_data/00_sra/` and corresponding FASTQs are under `01_data/01_raw_fastq/`.
Each sample folder is named with the SRA run ID (e.g., `SRR28578498`).

## Pipeline
1. Extract FASTQ from SRA: `02_scripts/00-extract_fastq.sh`
2. QC: FastQC + MultiQC: `02_scripts/01_qc_fastq.sh`
3. Trimming: fastp: `02_scripts/02_trim_fastp.sh`

## Software versions
See `00_meta/software_versions.txt` for complete module listings and tool versions.

## Project structure
- `00_meta/`: Project notes and software versions.
- `01_data/`: Raw and processed data.
- `01_data/00_sra/`: Original SRA files.
- `01_data/01_raw_fastq/`: Raw paired-end FASTQ files per sample.
- `01_data/02_trimmed_fastq/`: FASTQ files after trimming.
- `02_scripts/`: Pipeline scripts for extraction, QC, and trimming.
- `03_analysis/`: Analysis outputs by stage (QC, fastp reports).
- `04_results/`: Final figures, reports, and tables.
- `99_logs/`: Log files from pipeline runs.

## Notes
- Add a sample manifest (e.g., `sample_manifest.tsv`) to document condition assignments and replicate structure.
- Update this README with alignment, counting, and DE analysis details once those steps are completed.

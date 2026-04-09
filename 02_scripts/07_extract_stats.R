#!/usr/bin/env Rscript
# =============================================================================
# Pipeline QC Statistics Extraction
# Project: CIBIG Capstone - Bemisia tabaci begomovirus response
#
# Extracts and collates per-sample statistics from:
#   1. fastp JSON reports   → pre- and post-trimming stats
#   2. HISAT2 summary files → decontamination stats
#   3. STAR Log.final.out   → alignment stats
#   4. featureCounts summary → quantification stats
#
# Outputs clean CSV tables to 04_results/tables/qc_stats/
# =============================================================================

suppressPackageStartupMessages(library(jsonlite))

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]   # /scratch/onilee/capstone

# Input directories
fastp_dir    <- file.path(base_dir, "03_analysis/01_qc/c_fastp_reports")
hisat2_dir   <- file.path(base_dir, "03_analysis/03_hisat2_decontam/summaries")
star_dir     <- file.path(base_dir, "03_analysis/04_star_align/logs")
fc_summary   <- file.path(base_dir, "03_analysis/04_star_align/counts/all_samples_genelevel.txt.summary")

# Output directory
out_dir <- file.path(base_dir, "04_results/tables/qc_stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== QC Statistics Extraction Started:", format(Sys.time()), "===\n\n")

# =============================================================================
# HELPER: safe numeric extraction
# =============================================================================

pct <- function(x, total) round(100 * x / total, 2)

# =============================================================================
# 1. FASTP — pre- and post-trimming stats
# =============================================================================

cat("[1/4] Parsing fastp JSON reports...\n")

json_files <- list.files(fastp_dir, pattern = "\\.fastp\\.json$", full.names = TRUE)

if (length(json_files) == 0) {
  cat("  WARNING: No fastp JSON files found in", fastp_dir, "\n")
} else {
  rows <- lapply(json_files, function(f) {
    sample <- sub("\\.fastp\\.json$", "", basename(f))
    j      <- fromJSON(f)

    bf <- j$summary$before_filtering
    af <- j$summary$after_filtering
    fr <- j$filtering_result

    total_in  <- bf$total_reads
    total_out <- af$total_reads

    data.frame(
      Sample                   = sample,
      # Pre-trimming
      Raw_Reads                = total_in,
      Raw_Bases                = bf$total_bases,
      Raw_Q30_pct              = round(bf$q30_rate * 100, 2),
      Raw_GC_pct               = round(bf$gc_content * 100, 2),
      Raw_R1_mean_len          = bf$read1_mean_length,
      Raw_R2_mean_len          = bf$read2_mean_length,
      # Post-trimming
      Clean_Reads              = total_out,
      Clean_Bases              = af$total_bases,
      Clean_Q30_pct            = round(af$q30_rate * 100, 2),
      Clean_GC_pct             = round(af$gc_content * 100, 2),
      Clean_R1_mean_len        = af$read1_mean_length,
      Clean_R2_mean_len        = af$read2_mean_length,
      # Filtering summary
      Reads_passed_filter      = fr$passed_filter_reads,
      Reads_passed_filter_pct  = pct(fr$passed_filter_reads, total_in),
      Reads_low_quality        = fr$low_quality_reads,
      Reads_too_short          = fr$too_short_reads,
      Reads_too_many_N         = fr$too_many_N_reads,
      stringsAsFactors         = FALSE
    )
  })

  fastp_df <- do.call(rbind, rows)
  fastp_df <- fastp_df[order(fastp_df$Sample), ]

  write.csv(fastp_df,
    file.path(out_dir, "01_trimming_stats.csv"),
    row.names = FALSE)

  cat(sprintf("  Parsed %d samples → 01_trimming_stats.csv\n", nrow(fastp_df)))
}

# =============================================================================
# 2. HISAT2 — decontamination stats
# =============================================================================

cat("[2/4] Parsing HISAT2 decontamination summaries...\n")

# HISAT2 summary format (key lines):
#   N reads; of these:
#   N (X%) were paired; of these:
#     N (X%) aligned concordantly 0 times
#     N (X%) aligned concordantly exactly 1 time
#     N (X%) aligned concordantly >1 times
#   X% overall alignment rate

hisat2_files <- list.files(hisat2_dir, pattern = "\\.hisat2\\.summary\\.txt$", full.names = TRUE)

if (length(hisat2_files) == 0) {
  cat("  WARNING: No HISAT2 summary files found in", hisat2_dir, "\n")
} else {
  rows <- lapply(hisat2_files, function(f) {
    sample <- sub("\\.hisat2\\.summary\\.txt$", "", basename(f))
    lines  <- readLines(f)

    grab_n   <- function(pattern) {
      ln <- grep(pattern, lines, value = TRUE)
      if (length(ln) == 0) return(NA_real_)
      as.numeric(gsub(",", "", regmatches(ln[1], regexpr("^\\s*([0-9,]+)", ln[1], perl = TRUE))[[1]]))
    }
    grab_pct <- function(pattern) {
      ln <- grep(pattern, lines, value = TRUE)
      if (length(ln) == 0) return(NA_real_)
      as.numeric(regmatches(ln[1], regexpr("[0-9.]+(?=%)", ln[1], perl = TRUE))[[1]])
    }

    total        <- grab_n("reads; of these")
    conc_0       <- grab_n("aligned concordantly 0 times")
    conc_1       <- grab_n("aligned concordantly exactly 1 time")
    conc_gt1     <- grab_n("aligned concordantly >1 times")
    overall_rate <- grab_pct("overall alignment rate")

    mapped_pairs   <- total - conc_0
    unmapped_pairs <- conc_0

    data.frame(
      Sample                        = sample,
      Input_Read_Pairs              = total,
      Aligned_Concordantly_1        = conc_1,
      Aligned_Concordantly_1_pct    = pct(conc_1, total),
      Aligned_Concordantly_gt1      = conc_gt1,
      Aligned_Concordantly_gt1_pct  = pct(conc_gt1, total),
      Mapped_to_contaminants_pairs  = mapped_pairs,
      Mapped_to_contaminants_pct    = round(overall_rate, 2),
      Retained_unmapped_pairs       = unmapped_pairs,
      Retained_unmapped_pct         = round(100 - overall_rate, 2),
      stringsAsFactors              = FALSE
    )
  })

  hisat2_df <- do.call(rbind, rows)
  hisat2_df <- hisat2_df[order(hisat2_df$Sample), ]

  write.csv(hisat2_df,
    file.path(out_dir, "02_decontamination_stats.csv"),
    row.names = FALSE)

  cat(sprintf("  Parsed %d samples → 02_decontamination_stats.csv\n", nrow(hisat2_df)))
}

# =============================================================================
# 3. STAR — alignment stats
# =============================================================================

cat("[3/4] Parsing STAR Log.final.out files...\n")

# STAR Log.final.out uses "field label |<tab>value" format

star_files <- list.files(star_dir, pattern = "Log\\.final\\.out$", full.names = TRUE)

if (length(star_files) == 0) {
  cat("  WARNING: No STAR Log.final.out files found in", star_dir, "\n")
} else {
  rows <- lapply(star_files, function(f) {
    # Sample name from prefix before _Log.final.out
    sample <- sub("_Log\\.final\\.out$", "", basename(f))
    lines  <- readLines(f)

    # Extract value after the | delimiter
    get_val <- function(pattern) {
      ln <- grep(pattern, lines, value = TRUE, fixed = TRUE)
      if (length(ln) == 0) return(NA_character_)
      trimws(sub(".*\\|\\s*", "", ln[1]))
    }

    to_num <- function(x) suppressWarnings(as.numeric(gsub("%", "", x)))
    to_int <- function(x) suppressWarnings(as.integer(gsub(",", "", x)))

    data.frame(
      Sample                          = sample,
      Input_Reads                     = to_int(get_val("Number of input reads")),
      Avg_Input_Read_Length           = to_num(get_val("Average input read length")),
      Uniquely_Mapped_Reads           = to_int(get_val("Uniquely mapped reads number")),
      Uniquely_Mapped_pct             = to_num(get_val("Uniquely mapped reads %")),
      Avg_Mapped_Length               = to_num(get_val("Average mapped length")),
      Splices_Total                   = to_int(get_val("Number of splices: Total")),
      Multi_Mapped_Reads              = to_int(get_val("Number of reads mapped to multiple loci")),
      Multi_Mapped_pct                = to_num(get_val("% of reads mapped to multiple loci")),
      Unmapped_TooShort_pct           = to_num(get_val("% of reads unmapped: too short")),
      Unmapped_Mismatch_pct           = to_num(get_val("% of reads unmapped: too many mismatches")),
      Unmapped_Other_pct              = to_num(get_val("% of reads unmapped: other")),
      stringsAsFactors                = FALSE
    )
  })

  star_df <- do.call(rbind, rows)
  star_df  <- star_df[order(star_df$Sample), ]

  write.csv(star_df,
    file.path(out_dir, "03_star_alignment_stats.csv"),
    row.names = FALSE)

  cat(sprintf("  Parsed %d samples → 03_star_alignment_stats.csv\n", nrow(star_df)))
}

# =============================================================================
# 4. FEATURECOUNTS — quantification stats
# =============================================================================

cat("[4/4] Parsing featureCounts summary...\n")

if (!file.exists(fc_summary)) {
  cat("  WARNING: featureCounts summary not found:", fc_summary, "\n")
} else {
  fc_raw <- read.table(fc_summary, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)

  # Rows are status categories; columns are Status + BAM file paths
  # Transpose so each sample is a row
  status_col <- fc_raw[, 1]
  bam_cols   <- fc_raw[, -1]

  # Clean sample names from BAM paths
  sample_names <- sub(".*/(SRR[0-9]+)_.*", "\\1", colnames(bam_cols))

  fc_t           <- as.data.frame(t(bam_cols))
  colnames(fc_t) <- status_col
  fc_t$Sample    <- sample_names

  # Calculate totals and percentages for key categories
  fc_out <- data.frame(
    Sample                      = fc_t$Sample,
    Assigned                    = as.integer(fc_t$Assigned),
    Unassigned_NoFeatures       = as.integer(fc_t$Unassigned_NoFeatures),
    Unassigned_Ambiguity        = as.integer(fc_t$Unassigned_Ambiguity),
    Unassigned_MultiMapping     = as.integer(fc_t$Unassigned_MultiMapping),
    stringsAsFactors            = FALSE
  )

  fc_out$Total_Reads  <- rowSums(fc_t[, status_col, drop = FALSE] |>
                           lapply(as.integer) |> as.data.frame())
  fc_out$Assigned_pct            <- pct(fc_out$Assigned,                  fc_out$Total_Reads)
  fc_out$Unassigned_NoFeatures_pct   <- pct(fc_out$Unassigned_NoFeatures,  fc_out$Total_Reads)
  fc_out$Unassigned_Ambiguity_pct    <- pct(fc_out$Unassigned_Ambiguity,   fc_out$Total_Reads)
  fc_out$Unassigned_MultiMapping_pct <- pct(fc_out$Unassigned_MultiMapping, fc_out$Total_Reads)

  fc_out <- fc_out[order(fc_out$Sample), ]

  write.csv(fc_out,
    file.path(out_dir, "04_featurecounts_stats.csv"),
    row.names = FALSE)

  cat(sprintf("  Parsed %d samples → 04_featurecounts_stats.csv\n", nrow(fc_out)))
}

# =============================================================================
# 5. COMBINED SUMMARY TABLE — one row per sample, key stats only
# =============================================================================

cat("\nBuilding combined pipeline summary...\n")

all_dfs <- list()

if (exists("fastp_df"))  all_dfs[["fastp"]]  <- fastp_df[,  c("Sample", "Raw_Reads", "Reads_passed_filter_pct", "Clean_Q30_pct")]
if (exists("hisat2_df")) all_dfs[["hisat2"]] <- hisat2_df[, c("Sample", "Mapped_to_contaminants_pct", "Retained_unmapped_pairs")]
if (exists("star_df"))   all_dfs[["star"]]   <- star_df[,   c("Sample", "Input_Reads", "Uniquely_Mapped_pct", "Multi_Mapped_pct")]
if (exists("fc_out"))    all_dfs[["fc"]]     <- fc_out[,    c("Sample", "Assigned_pct", "Unassigned_NoFeatures_pct", "Unassigned_Ambiguity_pct")]

if (length(all_dfs) > 0) {
  summary_df <- Reduce(function(a, b) merge(a, b, by = "Sample", all = TRUE), all_dfs)
  summary_df <- summary_df[order(summary_df$Sample), ]

  write.csv(summary_df,
    file.path(out_dir, "00_pipeline_summary.csv"),
    row.names = FALSE)

  cat(sprintf("  Combined summary → 00_pipeline_summary.csv (%d samples)\n", nrow(summary_df)))

  cat("\n--- Pipeline summary (key metrics) ---\n")
  print(summary_df, row.names = FALSE)
}

cat("\n=== Extraction Complete:", format(Sys.time()), "===\n")
cat("Output directory:", out_dir, "\n\n")
cat("Files written:\n")
for (f in list.files(out_dir, pattern = "\\.csv$", full.names = FALSE)) {
  cat("  ", f, "\n")
}

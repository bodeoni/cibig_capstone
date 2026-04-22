#!/usr/bin/env Rscript
# =============================================================================
# Time-point contrast analysis
# Project: CIBIG Capstone - Bemisia tabaci begomovirus response
#
# Analyses:
#   A. Overall time effects    (~ Group + Time_Point; controls for group)
#      - 48h vs 24h, 72h vs 24h
#   B. Time effects in virus group   (subset; ~ Time_Point)
#      - 48h vs 24h, 72h vs 24h
#   C. Time effects in control group (subset; ~ Time_Point)
#      - 48h vs 24h, 72h vs 24h
#
# Three LFC cutoffs per contrast:
#   lfc1   : padj < 0.05 & |LFC| > 1    (2-fold change)
#   lfc058 : padj < 0.05 & |LFC| > 0.58 (1.5-fold change)
#   pval   : padj < 0.05 only
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

args       <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]   # all_samples_genelevel.txt  (from featureCounts)
meta_file  <- args[2]   # sample_metadata.csv
out_dir    <- args[3]   # output directory

dir.create(file.path(out_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

cat("=== Time-point Contrast Analysis Started:", format(Sys.time()), "===\n\n")

# =============================================================================
# HELPERS
# =============================================================================

# Volcano plot for a given LFC threshold (lfc_thresh = 0 → padj only)
make_volcano <- function(res_df, lfc_thresh, title, out_path_base) {
  vol <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
  vol$sig <- "Not significant"
  if (lfc_thresh > 0) {
    vol$sig[vol$padj < 0.05 & vol$log2FoldChange >  lfc_thresh] <- "Up"
    vol$sig[vol$padj < 0.05 & vol$log2FoldChange < -lfc_thresh] <- "Down"
  } else {
    vol$sig[vol$padj < 0.05 & vol$log2FoldChange > 0] <- "Up"
    vol$sig[vol$padj < 0.05 & vol$log2FoldChange < 0] <- "Down"
  }
  vol$sig <- factor(vol$sig, levels = c("Up", "Down", "Not significant"))
  top10   <- head(vol[vol$sig != "Not significant", ], 10)
  n_up    <- sum(vol$sig == "Up")
  n_dn    <- sum(vol$sig == "Down")
  cut_lbl <- if (lfc_thresh > 0) paste0("padj<0.05, |LFC|>", lfc_thresh) else "padj<0.05"

  p <- ggplot(vol, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_point(data = top10, size = 2, alpha = 0.9) +
    geom_text(data = top10, aes(label = GeneID),
      vjust = -0.5, hjust = 0.5, size = 2.8, colour = "black", check_overlap = TRUE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("Up" = "#D6604D", "Down" = "#2166AC", "Not significant" = "grey70")) +
    labs(
      title    = title,
      subtitle = paste0(n_up + n_dn, " DEGs (", cut_lbl, ")  |  ", n_up, " up, ", n_dn, " down"),
      x        = expression(log[2]~"fold change"),
      y        = expression(-log[10]~"adjusted p-value"),
      colour   = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  if (lfc_thresh > 0)
    p <- p + geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", colour = "grey40")

  ggsave(paste0(out_path_base, ".pdf"), p, width = 7, height = 6, dpi = 300)
  ggsave(paste0(out_path_base, ".png"), p, width = 7, height = 6, dpi = 300)
  invisible(p)
}

# Heatmap of top DEGs (row z-score of VST counts)
make_heatmap <- function(sig_df, vsd, meta, ann_colors, label, title, out_dir) {
  if (nrow(sig_df) < 2) {
    cat("    Skipping heatmap: fewer than 2 DEGs\n")
    return(invisible(NULL))
  }
  top_genes <- head(sig_df$GeneID, min(50, nrow(sig_df)))
  top_genes <- top_genes[top_genes %in% rownames(assay(vsd))]
  if (length(top_genes) < 2) {
    cat("    Skipping heatmap: genes not found in vsd\n")
    return(invisible(NULL))
  }
  mat        <- assay(vsd)[top_genes, , drop = FALSE]
  mat_scaled <- t(scale(t(mat)))

  # Build annotation: always include Time_Point; include Group only if >1 level
  ann_col <- data.frame(Time_Point = meta$Time_Point, row.names = rownames(meta))
  if ("Group" %in% colnames(meta) && nlevels(factor(meta$Group)) > 1)
    ann_col <- cbind(data.frame(Group = meta$Group, row.names = rownames(meta)), ann_col)
  ann_colors_use <- ann_colors[names(ann_colors) %in% colnames(ann_col)]

  hm_title <- paste0(title, "\n(top ", length(top_genes), " DEGs, row z-score of VST counts)")

  for (ext in c("pdf", "png")) {
    out_file <- file.path(out_dir, "plots", paste0(label, "_heatmap.", ext))
    if (ext == "pdf") pdf(out_file, width = 10, height = 12) else
      png(out_file, width = 1000, height = 1200, res = 100)
    pheatmap(
      mat_scaled,
      annotation_col    = ann_col,
      annotation_colors = ann_colors_use,
      color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_rows      = TRUE,
      cluster_cols      = TRUE,
      show_rownames     = length(top_genes) <= 60,
      show_colnames     = TRUE,
      fontsize_row      = 7,
      fontsize_col      = 9,
      main              = hm_title
    )
    dev.off()
  }
}

# MA plot (lfc_thresh = 0 → padj only colouring)
make_ma <- function(res_df, lfc_thresh, title, out_path_base) {
  ma <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
  ma$sig <- "Not significant"
  if (lfc_thresh > 0) {
    ma$sig[ma$padj < 0.05 & ma$log2FoldChange >  lfc_thresh] <- "Up"
    ma$sig[ma$padj < 0.05 & ma$log2FoldChange < -lfc_thresh] <- "Down"
  } else {
    ma$sig[ma$padj < 0.05 & ma$log2FoldChange > 0] <- "Up"
    ma$sig[ma$padj < 0.05 & ma$log2FoldChange < 0] <- "Down"
  }
  ma$sig       <- factor(ma$sig, levels = c("Up", "Down", "Not significant"))
  ma$mean_expr <- log10(ma$baseMean + 1)
  cut_lbl      <- if (lfc_thresh > 0) paste0("padj<0.05, |LFC|>", lfc_thresh) else "padj<0.05"

  p <- ggplot(ma, aes(mean_expr, log2FoldChange, colour = sig)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "black") +
    scale_colour_manual(values = c("Up" = "#D6604D", "Down" = "#2166AC", "Not significant" = "grey70")) +
    labs(
      title    = paste0("MA plot: ", title),
      subtitle = cut_lbl,
      x        = expression(log[10]~"mean expression"),
      y        = expression(log[2]~"fold change"),
      colour   = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

  if (lfc_thresh > 0)
    p <- p + geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", colour = "grey40")

  ggsave(paste0(out_path_base, ".pdf"), p, width = 7, height = 5, dpi = 300)
  ggsave(paste0(out_path_base, ".png"), p, width = 7, height = 5, dpi = 300)
  invisible(p)
}

# Filter at 3 cutoffs, write tables, produce 3 volcanos + MA plot + 1 heatmap (lfc1)
run_contrast <- function(res_df, label, title, vsd, meta, ann_colors, out_dir) {
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

  sig_lfc1   <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,    ]
  sig_lfc058 <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.58, ]
  sig_pval   <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

  write.csv(res_df,    file.path(out_dir, "tables", paste0(label, "_all.csv")),        row.names = FALSE)
  write.csv(sig_lfc1,  file.path(out_dir, "tables", paste0(label, "_sig_lfc1.csv")),   row.names = FALSE)
  write.csv(sig_lfc058,file.path(out_dir, "tables", paste0(label, "_sig_lfc058.csv")), row.names = FALSE)
  write.csv(sig_pval,  file.path(out_dir, "tables", paste0(label, "_sig_pval.csv")),   row.names = FALSE)

  cat(sprintf("    |LFC|>1:    %d DEGs [%d up, %d down]\n",
    nrow(sig_lfc1),   sum(sig_lfc1$log2FoldChange > 0),   sum(sig_lfc1$log2FoldChange < 0)))
  cat(sprintf("    |LFC|>0.58: %d DEGs [%d up, %d down]\n",
    nrow(sig_lfc058), sum(sig_lfc058$log2FoldChange > 0), sum(sig_lfc058$log2FoldChange < 0)))
  cat(sprintf("    padj only:  %d DEGs [%d up, %d down]\n",
    nrow(sig_pval),   sum(sig_pval$log2FoldChange > 0),   sum(sig_pval$log2FoldChange < 0)))

  for (cut in list(list(thresh = 1, tag = "lfc1"), list(thresh = 0.58, tag = "lfc058"), list(thresh = 0, tag = "pval"))) {
    make_volcano(res_df, cut$thresh, title,
      out_path_base = file.path(out_dir, "plots", paste0(label, "_volcano_", cut$tag)))
  }

  for (cut in list(list(thresh = 1, tag = "lfc1"), list(thresh = 0.58, tag = "lfc058"), list(thresh = 0, tag = "pval"))) {
    make_ma(res_df, cut$thresh, title,
      out_path_base = file.path(out_dir, "plots", paste0(label, "_ma_", cut$tag)))
  }

  # Heatmap uses lfc1 DEGs; falls back to lfc058 if lfc1 is empty
  hm_sig <- if (nrow(sig_lfc1) >= 2) sig_lfc1 else sig_lfc058
  make_heatmap(hm_sig, vsd, meta, ann_colors,
    label   = paste0(label, "_lfc1"),
    title   = title,
    out_dir = out_dir)

  invisible(list(lfc1 = sig_lfc1, lfc058 = sig_lfc058, pval = sig_pval))
}

# =============================================================================
# 1. LOAD COUNT MATRIX
# =============================================================================

cat("[1/5] Loading count matrix...\n")

raw_header    <- readLines(count_file, n = 2)[2]
header_fields <- strsplit(raw_header, "\t")[[1]]
bam_paths     <- header_fields[7:length(header_fields)]
sample_names  <- sub(".*/(SRR[0-9]+)_.*", "\\1", bam_paths)

cat("  Detected", length(sample_names), "samples:", paste(sample_names, collapse = ", "), "\n")

counts_raw <- read.table(count_file, header = FALSE, sep = "\t", skip = 2, stringsAsFactors = FALSE)
gene_ids   <- counts_raw[, 1]
counts     <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
colnames(counts) <- sample_names
rownames(counts) <- gene_ids

cat("  Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# =============================================================================
# 2. LOAD METADATA
# =============================================================================

cat("[2/5] Loading metadata...\n")

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[meta$Run != "", ]
rownames(meta) <- meta$Run
meta <- meta[sample_names, ]
meta$Group      <- factor(meta$Group,      levels = c("control", "virus"))
meta$Time_Point <- factor(meta$Time_Point, levels = c("24h", "48h", "72h"))

cat("  Sample breakdown:\n")
print(table(Group = meta$Group, Time_Point = meta$Time_Point))
stopifnot(all(rownames(meta) == colnames(counts)))

ann_colors <- list(
  Group      = c(control = "#2166AC", virus = "#D6604D"),
  Time_Point = c(`24h` = "#FEE08B", `48h` = "#FDAE61", `72h` = "#D73027")
)

# =============================================================================
# 3. ANALYSIS A: OVERALL TIME EFFECTS (~ Group + Time_Point)
# 24h is the reference level. Contrasts: 48h vs 24h, 72h vs 24h
# Group is included as a covariate to increase power.
# =============================================================================

cat("\n[3/5] Analysis A: Overall time effects (~ Group + Time_Point)...\n")

keep_A <- rowSums(counts) >= 10
dds_A  <- DESeqDataSetFromMatrix(counts[keep_A, ], meta, design = ~ Group + Time_Point)
dds_A  <- DESeq(dds_A)
vsd_A  <- vst(dds_A, blind = TRUE)

write.csv(
  data.frame(
    Analysis               = "A_overall_time",
    Total_Genes_Quantified = nrow(counts),
    Genes_Passing_Filter   = sum(keep_A),
    Genes_Tested           = nrow(dds_A)
  ),
  file.path(out_dir, "tables", "genes_tested_A.csv"),
  row.names = FALSE
)
cat("  Genes passing filter (Analysis A):", sum(keep_A), "\n")

for (ctp in list(c("48h", "24h"), c("72h", "24h"))) {
  tp_num <- ctp[1]; tp_ref <- ctp[2]
  label  <- paste0("overall_", tp_num, "_vs_", tp_ref)
  title  <- paste0("Time effect: ", tp_num, " vs ", tp_ref, " (all samples)")

  cat(sprintf("\n  %s vs %s\n", tp_num, tp_ref))

  res    <- results(dds_A, contrast = c("Time_Point", tp_num, tp_ref),
                    alpha = 0.05, pAdjustMethod = "BH")
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)

  run_contrast(res_df, label, title, vsd_A, meta, ann_colors, out_dir)
}

# =============================================================================
# 4. ANALYSES B & C: WITHIN-GROUP TIME EFFECTS
# For each group (virus, control): subset samples, fit ~ Time_Point,
# contrast 48h vs 24h and 72h vs 24h.
# fitType = "local" is used as a safeguard with smaller sample sizes.
# =============================================================================

cat("\n[4/5] Analyses B & C: Within-group time effects...\n")

for (grp in c("virus", "control")) {
  cat(sprintf("\n  Group: %s\n", grp))

  grp_idx    <- meta$Group == grp
  grp_counts <- counts[rowSums(counts[, grp_idx]) >= 10, grp_idx]
  grp_meta   <- meta[grp_idx, , drop = FALSE]
  grp_meta$Time_Point <- droplevels(grp_meta$Time_Point)

  dds_grp <- DESeqDataSetFromMatrix(grp_counts, grp_meta, design = ~ Time_Point)
  dds_grp <- DESeq(dds_grp, fitType = "local")
  vsd_grp <- vst(dds_grp, blind = TRUE)

  for (ctp in list(c("48h", "24h"), c("72h", "24h"))) {
    tp_num <- ctp[1]; tp_ref <- ctp[2]
    label  <- paste0(grp, "_", tp_num, "_vs_", tp_ref)
    title  <- paste0("Time effect: ", tp_num, " vs ", tp_ref, " (", grp, " group)")

    cat(sprintf("    %s vs %s\n", tp_num, tp_ref))

    res    <- results(dds_grp, contrast = c("Time_Point", tp_num, tp_ref),
                      alpha = 0.05, pAdjustMethod = "BH")
    res_df <- as.data.frame(res)
    res_df$GeneID <- rownames(res_df)

    run_contrast(res_df, label, title, vsd_grp, grp_meta, ann_colors, out_dir)
  }
}

# =============================================================================
# 5. SESSION INFO
# =============================================================================

cat("\n[5/5] Writing session info...\n")

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Analysis Complete:", format(Sys.time()), "===\n")
cat("Output directory:", out_dir, "\n\n")

cat("Contrasts run:\n")
cat("  Analysis A (overall):  overall_48h_vs_24h, overall_72h_vs_24h\n")
cat("  Analysis B (virus):    virus_48h_vs_24h,   virus_72h_vs_24h\n")
cat("  Analysis C (control):  control_48h_vs_24h, control_72h_vs_24h\n\n")

cat("Per-contrast outputs (tables/ and plots/):\n")
cat("  *_all.csv          all tested genes\n")
cat("  *_sig_lfc1.csv     padj<0.05 & |LFC|>1\n")
cat("  *_sig_lfc058.csv   padj<0.05 & |LFC|>0.58\n")
cat("  *_sig_pval.csv     padj<0.05 only\n")
cat("  *_volcano_lfc1/lfc058/pval.png  volcano plots\n")
cat("  *_lfc1_heatmap.png              top DEG heatmap\n\n")

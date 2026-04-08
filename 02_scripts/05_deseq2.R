#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis
# Project: CIBIG Capstone - Bemisia tabaci begomovirus response
# Design:  ~ Time_Point + Group  (virus vs control, accounting for time)
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Paths (passed from SLURM script) ─────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
count_file  <- args[1]   # all_samples_genelevel.txt
meta_file   <- args[2]   # sample_metadata.csv
out_dir     <- args[3]   # output directory

dir.create(file.path(out_dir, "plots"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

cat("=== DESeq2 Analysis Started:", format(Sys.time()), "===\n\n")

# =============================================================================
# 1. LOAD COUNT MATRIX
# =============================================================================

cat("[1/7] Loading count matrix...\n")

# Read the raw header line (line 2, skipping featureCounts comment on line 1)
# to extract SRR IDs from the full BAM paths before R mangling
raw_header <- readLines(count_file, n = 2)[2]
header_fields <- strsplit(raw_header, "\t")[[1]]
# Columns 1-6: Geneid Chr Start End Strand Length; 7+ are BAM paths
bam_paths    <- header_fields[7:length(header_fields)]
sample_names <- sub(".*/(SRR[0-9]+)_.*", "\\1", bam_paths)

cat("  Detected", length(sample_names), "samples:", paste(sample_names, collapse = ", "), "\n")

# Read the count matrix (skip comment line 1; suppress header since we set it)
counts_raw <- read.table(
  count_file,
  header = FALSE,
  sep    = "\t",
  skip   = 2,           # skip comment + header rows
  stringsAsFactors = FALSE
)

# Column 1 = GeneID; columns 7+ = counts; drop annotation columns 2-6
gene_ids <- counts_raw[, 1]
counts   <- counts_raw[, 7:ncol(counts_raw)]
colnames(counts) <- sample_names
rownames(counts) <- gene_ids
counts <- as.matrix(counts)

cat("  Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# =============================================================================
# 2. LOAD METADATA & MATCH TO COUNT MATRIX
# =============================================================================

cat("[2/7] Loading metadata...\n")

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[meta$Run != "", ]   # drop any blank trailing rows
rownames(meta) <- meta$Run

# Reorder metadata to match count matrix column order
meta <- meta[sample_names, ]

# Set factor levels: control is reference for Group; 24h reference for Time_Point
meta$Group      <- factor(meta$Group,      levels = c("control", "virus"))
meta$Time_Point <- factor(meta$Time_Point, levels = c("24h", "48h", "72h"))

cat("  Sample breakdown:\n")
print(table(Group = meta$Group, Time_Point = meta$Time_Point))

# Confirm row/col alignment
stopifnot(all(rownames(meta) == colnames(counts)))

# =============================================================================
# 3. BUILD DESEQ2 DATASET & RUN
# =============================================================================

cat("[3/7] Running DESeq2 (design: ~ Time_Point + Group)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ Time_Point + Group
)

# Pre-filter: remove genes with < 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("  Genes retained after low-count filter:", nrow(dds), "\n")

dds <- DESeq(dds)

cat("  DESeq2 complete\n")

# =============================================================================
# 4. EXTRACT & SAVE RESULTS  (virus vs control, overall effect)
# =============================================================================

cat("[4/7] Extracting results (virus vs control)...\n")

res <- results(dds,
  contrast    = c("Group", "virus", "control"),
  alpha       = 0.05,
  pAdjustMethod = "BH"
)

res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# All genes
write.csv(res_df,
  file      = file.path(out_dir, "tables", "deseq2_virus_vs_control_all.csv"),
  row.names = FALSE)

# Significant DEGs (padj < 0.05 & |LFC| > 1)
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

write.csv(sig,
  file      = file.path(out_dir, "tables", "deseq2_virus_vs_control_sig.csv"),
  row.names = FALSE)

cat("  Total tested genes:", nrow(res_df), "\n")
cat("  Significant DEGs (padj<0.05, |LFC|>1):", nrow(sig), "\n")
cat("  Upregulated in virus:", sum(sig$log2FoldChange > 0), "\n")
cat("  Downregulated in virus:", sum(sig$log2FoldChange < 0), "\n")

# =============================================================================
# 5. PLOT 1: PCA
# =============================================================================

cat("[5/7] Generating plots...\n")
cat("  - PCA\n")

vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("Group", "Time_Point"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, colour = Group, shape = Time_Point)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_colour_manual(values = c("control" = "#2166AC", "virus" = "#D6604D")) +
  labs(
    title    = expression(paste("PCA of VST-normalised counts (", italic("B. tabaci"), ")")),
    subtitle = "Virus acquisition vs control across three time points",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)"),
    colour   = "Group",
    shape    = "Time Point"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "plots", "pca_plot.pdf"),
  pca_plot, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "plots", "pca_plot.png"),
  pca_plot, width = 7, height = 5, dpi = 300)

# ── Plot 2: Sample Distance Heatmap ──────────────────────────────────────────

cat("  - Sample distance heatmap\n")

samp_dists <- dist(t(assay(vsd)))
dist_mat   <- as.matrix(samp_dists)

annotation_col <- data.frame(
  Group      = meta$Group,
  Time_Point = meta$Time_Point,
  row.names  = rownames(meta)
)

ann_colors <- list(
  Group      = c(control = "#2166AC", virus = "#D6604D"),
  Time_Point = c(`24h` = "#FEE08B", `48h` = "#FDAE61", `72h` = "#D73027")
)

pdf(file.path(out_dir, "plots", "sample_distance_heatmap.pdf"), width = 9, height = 7)
pheatmap(
  dist_mat,
  annotation_col  = annotation_col,
  annotation_row  = annotation_col,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  main            = "Sample-to-sample Euclidean distances (VST counts)",
  fontsize        = 10
)
dev.off()

# also save PNG
png(file.path(out_dir, "plots", "sample_distance_heatmap.png"), width = 900, height = 700, res = 100)
pheatmap(
  dist_mat,
  annotation_col  = annotation_col,
  annotation_row  = annotation_col,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  main            = "Sample-to-sample Euclidean distances (VST counts)",
  fontsize        = 10
)
dev.off()

# ── Plot 3: Volcano Plot ──────────────────────────────────────────────────────

cat("  - Volcano plot\n")

volcano_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
volcano_df$significance <- "Not significant"
volcano_df$significance[volcano_df$padj < 0.05 & volcano_df$log2FoldChange >  1] <- "Up in virus"
volcano_df$significance[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < -1] <- "Down in virus"
volcano_df$significance <- factor(volcano_df$significance,
  levels = c("Up in virus", "Down in virus", "Not significant"))

# Top 10 genes to label
top_label <- head(volcano_df[volcano_df$significance != "Not significant", ], 10)

volcano_plot <- ggplot(volcano_df, aes(log2FoldChange, -log10(padj), colour = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_point(data = top_label, size = 2, alpha = 0.9) +
  geom_text(data = top_label, aes(label = GeneID),
    vjust = -0.5, hjust = 0.5, size = 2.8, colour = "black", check_overlap = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Up in virus"      = "#D6604D",
    "Down in virus"    = "#2166AC",
    "Not significant"  = "grey70"
  )) +
  labs(
    title    = "Volcano plot: virus vs control",
    subtitle = paste0(
      nrow(sig), " DEGs (padj<0.05, |LFC|>1)  |  ",
      sum(sig$log2FoldChange > 0), " up, ",
      sum(sig$log2FoldChange < 0), " down in virus"
    ),
    x      = expression(log[2]~"fold change (virus / control)"),
    y      = expression(-log[10]~"adjusted p-value"),
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "plots", "volcano_plot.pdf"),
  volcano_plot, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "plots", "volcano_plot.png"),
  volcano_plot, width = 7, height = 6, dpi = 300)

# ── Plot 4: MA Plot ───────────────────────────────────────────────────────────

cat("  - MA plot\n")

ma_df <- volcano_df
ma_df$mean_expr <- log10(ma_df$baseMean + 1)

ma_plot <- ggplot(ma_df, aes(mean_expr, log2FoldChange, colour = significance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"),
             colour = c("grey40", "black", "grey40")) +
  scale_colour_manual(values = c(
    "Up in virus"      = "#D6604D",
    "Down in virus"    = "#2166AC",
    "Not significant"  = "grey70"
  )) +
  labs(
    title  = "MA plot: virus vs control",
    x      = expression(log[10]~"mean expression"),
    y      = expression(log[2]~"fold change (virus / control)"),
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "plots", "ma_plot.pdf"),
  ma_plot, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "plots", "ma_plot.png"),
  ma_plot, width = 7, height = 5, dpi = 300)

# ── Plot 5: Heatmap of Top 50 DEGs ───────────────────────────────────────────

cat("  - Heatmap of top 50 DEGs\n")

if (nrow(sig) >= 2) {
  # Use top 50 by adjusted p-value (or all if fewer)
  top_genes <- head(sig$GeneID, min(50, nrow(sig)))
  mat       <- assay(vsd)[top_genes, ]
  # Z-score per gene (row-scale)
  mat_scaled <- t(scale(t(mat)))

  annotation_col_hm <- data.frame(
    Group      = meta$Group,
    Time_Point = meta$Time_Point,
    row.names  = colnames(mat)
  )

  pdf(file.path(out_dir, "plots", "heatmap_top_degs.pdf"), width = 10, height = 12)
  pheatmap(
    mat_scaled,
    annotation_col    = annotation_col_hm,
    annotation_colors = ann_colors,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 9,
    main              = paste0("Top ", length(top_genes), " DEGs (row z-score of VST counts)")
  )
  dev.off()

  png(file.path(out_dir, "plots", "heatmap_top_degs.png"), width = 1000, height = 1200, res = 100)
  pheatmap(
    mat_scaled,
    annotation_col    = annotation_col_hm,
    annotation_colors = ann_colors,
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 9,
    main              = paste0("Top ", length(top_genes), " DEGs (row z-score of VST counts)")
  )
  dev.off()
} else {
  cat("  WARNING: fewer than 2 significant DEGs, skipping heatmap\n")
}

# =============================================================================
# 6. TIME-POINT-SPECIFIC RESULTS
# =============================================================================

cat("[6/7] Time-point-specific results...\n")

dds_tp <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ Time_Point + Group + Time_Point:Group
)
dds_tp <- dds_tp[rowSums(counts(dds_tp)) >= 10, ]
dds_tp <- DESeq(dds_tp)

for (tp in c("24h", "48h", "72h")) {
  tp_label <- gsub("h", "h", tp)

  # The interaction term tests whether the virus effect differs at this time point
  # For direct virus vs control within a time point, subset and rerun
  sub_idx <- meta$Time_Point == tp
  sub_counts <- counts[rowSums(counts) >= 10, sub_idx]
  sub_meta   <- meta[sub_idx, ]
  sub_meta$Group <- droplevels(sub_meta$Group)

  dds_sub <- DESeqDataSetFromMatrix(sub_counts, sub_meta, design = ~ Group)
  dds_sub <- DESeq(dds_sub, fitType = "local")

  res_sub <- results(dds_sub,
    contrast = c("Group", "virus", "control"),
    alpha    = 0.05, pAdjustMethod = "BH"
  )
  res_sub_df <- as.data.frame(res_sub)
  res_sub_df$GeneID <- rownames(res_sub_df)
  res_sub_df <- res_sub_df[order(res_sub_df$padj, na.last = TRUE), ]

  sig_sub <- res_sub_df[!is.na(res_sub_df$padj) & res_sub_df$padj < 0.05 &
                        abs(res_sub_df$log2FoldChange) > 1, ]

  write.csv(res_sub_df,
    file      = file.path(out_dir, "tables", paste0("deseq2_", tp, "_all.csv")),
    row.names = FALSE)
  write.csv(sig_sub,
    file      = file.path(out_dir, "tables", paste0("deseq2_", tp, "_sig.csv")),
    row.names = FALSE)

  cat(sprintf("  %s: %d DEGs (padj<0.05, |LFC|>1)  [%d up, %d down in virus]\n",
    tp, nrow(sig_sub),
    sum(sig_sub$log2FoldChange > 0),
    sum(sig_sub$log2FoldChange < 0)
  ))
}

# =============================================================================
# 7. SESSION INFO
# =============================================================================

cat("[7/7] Writing session info...\n")

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Analysis Complete:", format(Sys.time()), "===\n")
cat("Output directory:", out_dir, "\n")
cat("  tables/ - CSV results files\n")
cat("  plots/  - PCA, volcano, MA, heatmap (PDF + PNG)\n\n")

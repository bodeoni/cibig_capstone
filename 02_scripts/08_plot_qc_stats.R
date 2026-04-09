#!/usr/bin/env Rscript
# =============================================================================
# 08_plot_qc_stats.R
# Pipeline QC Visualisation
# Project: CIBIG Capstone – Bemisia tabaci begomovirus response
#
# Plots generated (all saved as PNG + PDF in 04_results/figures/qc_plots/):
#
#   01_pipeline_read_funnel         – read-pair retention through pipeline
#   02_raw_read_counts              – sequencing depth per sample
#   03_trimming_filter_reasons      – reads removed by discard category
#   04_q30_before_after             – Q30 rate dumbbell (pre vs post trimming)
#   05_gc_content_before_after      – GC content dumbbell (pre vs post trimming)
#   06_decontam_rate_per_sample     – % reads mapping to contaminants
#   07_decontam_stacked_bar         – % retained vs contaminant (stacked)
#   08_star_alignment_stacked       – STAR mapping fate (stacked %)
#   09_star_uniquely_mapped_scatter – uniquely mapped % dot plot
#   10_featurecounts_stacked        – featureCounts assignment categories (stacked %)
#   11_featurecounts_assigned_dot   – assigned % dot plot per sample
#   12_qc_summary_heatmap           – pheatmap of all key QC metrics
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else "/scratch/onilee/capstone"

stats_dir <- file.path(base_dir, "04_results/tables/qc_stats")
meta_file <- file.path(base_dir, "00_meta/sample_metadata.csv")
out_dir   <- file.path(base_dir, "04_results/figures/qc_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== QC Plots Started:", format(Sys.time()), "===\n\n")
cat("  Stats directory:", stats_dir, "\n")
cat("  Output directory:", out_dir, "\n\n")

# --------------------------------------------------------------------------- #
# Load data
# --------------------------------------------------------------------------- #

meta  <- read.csv(meta_file, stringsAsFactors = FALSE)
trim  <- read.csv(file.path(stats_dir, "01_trimming_stats.csv"),       stringsAsFactors = FALSE)
decon <- read.csv(file.path(stats_dir, "02_decontamination_stats.csv"), stringsAsFactors = FALSE)
star  <- read.csv(file.path(stats_dir, "03_star_alignment_stats.csv"),  stringsAsFactors = FALSE)
fc    <- read.csv(file.path(stats_dir, "04_featurecounts_stats.csv"),   stringsAsFactors = FALSE)

# Rename Run → Sample for join
names(meta)[names(meta) == "Run"] <- "Sample"

add_meta <- function(df) left_join(df, meta, by = "Sample")
trim  <- add_meta(trim)
decon <- add_meta(decon)
star  <- add_meta(star)
fc    <- add_meta(fc)

# --------------------------------------------------------------------------- #
# Sample ordering: time point (24h → 48h → 72h) then group (control, virus)
# Short label strips the "SRR28578" common prefix
# --------------------------------------------------------------------------- #

sample_order <- meta %>%
  mutate(
    Time_Point = factor(Time_Point, levels = c("24h", "48h", "72h")),
    Group      = factor(Group,      levels = c("control", "virus"))
  ) %>%
  arrange(Time_Point, Group) %>%
  pull(Sample)

short_id <- function(x) sub("SRR28578", "", as.character(x))

reorder_df <- function(df) {
  df$Sample <- factor(df$Sample, levels = sample_order)
  df$Short  <- factor(short_id(df$Sample), levels = short_id(sample_order))
  df$Time_Point <- factor(df$Time_Point, levels = c("24h", "48h", "72h"))
  df$Group      <- factor(df$Group,      levels = c("control", "virus"))
  df
}

trim  <- reorder_df(trim)
decon <- reorder_df(decon)
star  <- reorder_df(star)
fc    <- reorder_df(fc)

# --------------------------------------------------------------------------- #
# Colour palettes
# --------------------------------------------------------------------------- #

col_group <- c(virus = "#C0392B", control = "#2980B9")
col_time  <- c("24h" = "#E67E22", "48h" = "#27AE60", "72h" = "#8E44AD")
shape_tp  <- c("24h" = 16, "48h" = 17, "72h" = 15)

# --------------------------------------------------------------------------- #
# Publication theme
# --------------------------------------------------------------------------- #

theme_qc <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92"),
      strip.background  = element_rect(fill = "grey88", colour = "grey70"),
      strip.text        = element_text(face = "bold", size = base_size - 1),
      axis.title        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", size = base_size + 2,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(colour = "grey45", size = base_size - 1,
                                       margin = margin(b = 8)),
      plot.caption      = element_text(colour = "grey55", size = 9,
                                       hjust = 0, margin = margin(t = 8)),
      legend.position   = "right",
      legend.title      = element_text(face = "bold", size = base_size - 1),
      legend.key.size   = unit(0.9, "lines")
    )
}

# --------------------------------------------------------------------------- #
# Helper: save PNG + PDF and log
# --------------------------------------------------------------------------- #

plot_log <- character(0)   # accumulate a manifest of plots

save_plot <- function(p, name, width = 14, height = 7, desc = "") {
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, p, width = width, height = height, bg = "white")
  cat("  Saved:", name, "\n")
  plot_log <<- c(plot_log, sprintf("%-45s  %s", paste0(name, ".png"), desc))
  invisible(p)
}

# pheatmap version (returns a grob)
save_pheatmap <- function(ph, name, width = 12, height = 8, desc = "") {
  png(file.path(out_dir, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300, bg = "white")
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
  pdf(file.path(out_dir, paste0(name, ".pdf")),
      width = width, height = height, bg = "white")
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
  cat("  Saved:", name, "\n")
  plot_log <<- c(plot_log, sprintf("%-45s  %s", paste0(name, ".png"), desc))
}

# =============================================================================
# PLOT 1 — Pipeline read-retention funnel
# =============================================================================

cat("[1/12] Pipeline read-retention funnel...\n")

funnel <- trim %>%
  select(Sample, Short, Group, Time_Point, Raw_Reads, Clean_Reads) %>%
  left_join(decon %>% select(Sample, Retained_unmapped_pairs), by = "Sample") %>%
  left_join(star  %>% select(Sample, Input_Reads, Uniquely_Mapped_pct,
                              Uniquely_Mapped_Reads),           by = "Sample") %>%
  mutate(
    Step1_Raw          = Raw_Reads / 2 / 1e6,          # pairs
    Step2_Trimmed      = Clean_Reads / 2 / 1e6,
    Step3_Decontam     = Retained_unmapped_pairs / 1e6,
    Step4_UniqueAligned = Uniquely_Mapped_Reads / 1e6
  ) %>%
  select(Short, Group, Time_Point,
         `1. Raw` = Step1_Raw,
         `2. Trimmed` = Step2_Trimmed,
         `3. Decontaminated` = Step3_Decontam,
         `4. Uniquely aligned` = Step4_UniqueAligned) %>%
  pivot_longer(cols = c(`1. Raw`, `2. Trimmed`,
                         `3. Decontaminated`, `4. Uniquely aligned`),
               names_to = "Step", values_to = "Pairs_M") %>%
  mutate(Step = factor(Step, levels = c("1. Raw", "2. Trimmed",
                                         "3. Decontaminated", "4. Uniquely aligned")))

p1 <- ggplot(funnel, aes(x = Step, y = Pairs_M, group = Short,
                          colour = Group, linetype = Time_Point)) +
  geom_line(linewidth = 0.85, alpha = 0.85) +
  geom_point(aes(shape = Time_Point), size = 3.2) +
  scale_colour_manual("Group", values = col_group) +
  scale_shape_manual("Time point", values = shape_tp) +
  scale_linetype_manual("Time point",
                        values = c("24h" = "solid", "48h" = "dashed", "72h" = "dotted")) +
  scale_y_continuous(labels = label_number(suffix = " M"),
                     expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    title    = "Read-Pair Retention Through the Pipeline",
    subtitle = "Each line = one sample; note the large drop at decontamination (endosymbiont removal)",
    x        = "Pipeline Step",
    y        = "Read Pairs (millions)",
    caption  = "Raw and Trimmed are divided by 2 (paired-end); Decontaminated and Aligned are already pair-level"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(size = 11))

save_plot(p1, "01_pipeline_read_funnel", width = 11, height = 6,
          desc = "Read-pair retention at each pipeline step")

# =============================================================================
# PLOT 2 — Sequencing depth per sample
# =============================================================================

cat("[2/12] Sequencing depth per sample...\n")

mean_depth <- mean(trim$Raw_Reads / 2 / 1e6)

p2 <- trim %>%
  mutate(Raw_M = Raw_Reads / 2 / 1e6) %>%
  ggplot(aes(y = reorder(Short, as.integer(Short)), x = Raw_M, fill = Group)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = mean_depth,
             linetype = "dashed", colour = "grey25", linewidth = 0.7) +
  annotate("text", x = mean_depth + 0.3, y = 0.4,
           label = sprintf("Mean = %.1f M", mean_depth),
           hjust = 0, size = 3.2, colour = "grey25") +
  geom_text(aes(label = sprintf("%.1f M", Raw_M)),
            hjust = -0.1, size = 3, fontface = "plain") +
  facet_grid(Time_Point ~ Group, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = col_group) +
  scale_x_continuous(labels = label_number(suffix = " M"),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "Sequencing Depth per Sample",
    subtitle = "Raw read pairs in millions; grouped by time point and experimental group",
    x        = "Raw Read Pairs (millions)",
    y        = "Sample (SRR suffix)",
    fill     = "Group"
  ) +
  theme_qc()

save_plot(p2, "02_raw_read_counts", width = 11, height = 9,
          desc = "Raw sequencing depth per sample")

# =============================================================================
# PLOT 3 — Reads removed by trimming category
# =============================================================================

cat("[3/12] Trimming filter reasons...\n")

filter_long <- trim %>%
  mutate(
    `Low quality`  = Reads_low_quality / 1e3,
    `Too short`    = Reads_too_short   / 1e3,
    `Too many N`   = Reads_too_many_N  / 1e3
  ) %>%
  select(Short, Group, Time_Point, `Low quality`, `Too short`, `Too many N`) %>%
  pivot_longer(c(`Low quality`, `Too short`, `Too many N`),
               names_to = "Reason", values_to = "Reads_K") %>%
  mutate(Reason = factor(Reason,
                          levels = c("Low quality", "Too short", "Too many N")))

p3 <- ggplot(filter_long, aes(x = Short, y = Reads_K, fill = Reason)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  facet_grid(. ~ Time_Point + Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(
    "Low quality" = "#E74C3C",
    "Too short"   = "#F39C12",
    "Too many N"  = "#8E44AD"
  )) +
  scale_y_continuous(labels = label_number(suffix = "K"),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Reads Removed During Adapter Trimming",
    subtitle = "Counts in thousands; classified by discard reason (fastp)",
    x        = "Sample (SRR suffix)",
    y        = "Reads Discarded (thousands)",
    fill     = "Discard Reason",
    caption  = "Parameters: --cut_mean_quality 20, --cut_window_size 4, --length_required 75"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_plot(p3, "03_trimming_filter_reasons", width = 15, height = 6,
          desc = "Reads discarded by trimming category")

# =============================================================================
# PLOT 4 — Q30 rate before and after trimming (dumbbell)
# =============================================================================

cat("[4/12] Q30 dumbbell plot...\n")

q30_long <- trim %>%
  select(Short, Group, Time_Point, `Before trimming` = Raw_Q30_pct,
         `After trimming` = Clean_Q30_pct) %>%
  pivot_longer(c(`Before trimming`, `After trimming`),
               names_to = "Stage", values_to = "Q30") %>%
  mutate(Stage = factor(Stage, levels = c("Before trimming", "After trimming")))

p4 <- ggplot() +
  geom_segment(data = trim,
               aes(y = Short, yend = Short,
                   x = Raw_Q30_pct, xend = Clean_Q30_pct,
                   colour = Group),
               linewidth = 1.2, alpha = 0.45) +
  geom_point(data = q30_long,
             aes(y = Short, x = Q30, colour = Group, shape = Stage),
             size = 3.8) +
  geom_text(data = trim,
            aes(y = Short, x = Clean_Q30_pct,
                label = sprintf("+%.2f%%", Clean_Q30_pct - Raw_Q30_pct)),
            hjust = -0.25, size = 2.9, colour = "grey30") +
  facet_grid(Time_Point ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual("Group", values = col_group) +
  scale_shape_manual("Stage",
                     values = c("Before trimming" = 1, "After trimming" = 16)) +
  scale_x_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title    = "Q30 Base Quality Before and After Trimming",
    subtitle = "Open circle = before trimming; filled circle = after; label = improvement (Δ%)",
    x        = "Q30 Rate (%)",
    y        = "Sample (SRR suffix)"
  ) +
  theme_qc()

save_plot(p4, "04_q30_before_after", width = 10, height = 9,
          desc = "Q30 quality improvement from trimming (dumbbell)")

# =============================================================================
# PLOT 5 — GC content before and after trimming (dumbbell)
# =============================================================================

cat("[5/12] GC content dumbbell plot...\n")

gc_long <- trim %>%
  select(Short, Group, Time_Point, `Before trimming` = Raw_GC_pct,
         `After trimming` = Clean_GC_pct) %>%
  pivot_longer(c(`Before trimming`, `After trimming`),
               names_to = "Stage", values_to = "GC") %>%
  mutate(Stage = factor(Stage, levels = c("Before trimming", "After trimming")))

p5 <- ggplot() +
  geom_segment(data = trim,
               aes(y = Short, yend = Short,
                   x = Raw_GC_pct, xend = Clean_GC_pct,
                   colour = Group),
               linewidth = 1.2, alpha = 0.45) +
  geom_point(data = gc_long,
             aes(y = Short, x = GC, colour = Group, shape = Stage),
             size = 3.8) +
  facet_grid(Time_Point ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual("Group", values = col_group) +
  scale_shape_manual("Stage",
                     values = c("Before trimming" = 1, "After trimming" = 16)) +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  labs(
    title    = "GC Content Before and After Trimming",
    subtitle = "Open circle = before trimming; filled circle = after; minimal change expected",
    x        = "GC Content (%)",
    y        = "Sample (SRR suffix)"
  ) +
  theme_qc()

save_plot(p5, "05_gc_content_before_after", width = 10, height = 9,
          desc = "GC content before and after trimming (dumbbell)")

# =============================================================================
# PLOT 6 — Endosymbiont/mitochondrial contamination rate
# =============================================================================

cat("[6/12] Decontamination rate per sample...\n")

mean_contam <- mean(decon$Mapped_to_contaminants_pct)

p6 <- decon %>%
  ggplot(aes(x = Short, y = Mapped_to_contaminants_pct, fill = Group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_hline(yintercept = mean_contam,
             linetype = "dashed", colour = "grey25", linewidth = 0.7) +
  annotate("text",
           x    = 0.4,
           y    = mean_contam + 1.2,
           label = sprintf("Mean = %.1f%%", mean_contam),
           hjust = 0, size = 3.2, colour = "grey25") +
  geom_text(aes(label = sprintf("%.1f%%", Mapped_to_contaminants_pct)),
            vjust = -0.4, size = 2.9, fontface = "plain") +
  facet_grid(. ~ Time_Point, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_group) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     limits = c(0, 65),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Endosymbiont and Mitochondrial Contamination Rate",
    subtitle = paste0("% read pairs mapping to Portiera, Hamiltonella, Rickettsia or mitochondria (HISAT2); ",
                      "mean = ", sprintf("%.1f%%", mean_contam)),
    x        = "Sample (SRR suffix)",
    y        = "Contaminant Read Pairs (%)",
    fill     = "Group",
    caption  = "Reads mapping to contaminant genomes were removed; retained reads used for host genome alignment"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p6, "06_decontam_rate_per_sample", width = 13, height = 6,
          desc = "Endosymbiont contamination rate per sample")

# =============================================================================
# PLOT 7 — Decontamination stacked bar (% composition)
# =============================================================================

cat("[7/12] Decontamination stacked bar...\n")

decon_stack <- decon %>%
  select(Short, Group, Time_Point,
         `Contaminants (removed)` = Mapped_to_contaminants_pct,
         `Retained (host-derived)` = Retained_unmapped_pct) %>%
  pivot_longer(c(`Contaminants (removed)`, `Retained (host-derived)`),
               names_to = "Fate", values_to = "Pct") %>%
  mutate(Fate = factor(Fate,
                        levels = c("Contaminants (removed)", "Retained (host-derived)")))

p7 <- ggplot(decon_stack, aes(x = Short, y = Pct, fill = Fate)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(data = decon,
            aes(x = Short, y = 50, label = sprintf("%.0f%%", Retained_unmapped_pct)),
            inherit.aes = FALSE, size = 3, colour = "white", fontface = "bold") +
  facet_grid(. ~ Time_Point + Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(
    "Contaminants (removed)"  = "#E74C3C",
    "Retained (host-derived)" = "#2980B9"
  )) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Read Composition After Decontamination",
    subtitle = "Stacked bars show % retained vs % removed; white labels = % retained per sample",
    x        = "Sample (SRR suffix)",
    y        = "Read Pairs (%)",
    fill     = "Fate"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_plot(p7, "07_decontam_stacked_bar", width = 15, height = 6,
          desc = "% reads retained vs removed during decontamination")

# =============================================================================
# PLOT 8 — STAR alignment fate stacked bar
# =============================================================================

cat("[8/12] STAR alignment stacked bar...\n")

star_stack <- star %>%
  mutate(
    Unmapped_pct = Unmapped_TooShort_pct + Unmapped_Mismatch_pct + Unmapped_Other_pct
  ) %>%
  select(Short, Group, Time_Point,
         `Uniquely mapped` = Uniquely_Mapped_pct,
         `Multi-mapped`    = Multi_Mapped_pct,
         `Unmapped`        = Unmapped_pct) %>%
  pivot_longer(c(`Uniquely mapped`, `Multi-mapped`, `Unmapped`),
               names_to = "Fate", values_to = "Pct") %>%
  mutate(Fate = factor(Fate,
                        levels = c("Unmapped", "Multi-mapped", "Uniquely mapped")))

p8 <- ggplot(star_stack, aes(x = Short, y = Pct, fill = Fate)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(data = star,
            aes(x = Short, y = Uniquely_Mapped_pct / 2,
                label = sprintf("%.1f%%", Uniquely_Mapped_pct)),
            inherit.aes = FALSE, size = 2.8, colour = "white", fontface = "bold") +
  facet_grid(. ~ Time_Point + Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(
    "Uniquely mapped" = "#27AE60",
    "Multi-mapped"    = "#F39C12",
    "Unmapped"        = "#E74C3C"
  )) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "STAR Alignment Outcome per Sample",
    subtitle = "White labels = unique alignment rate; stack = proportional fate of all input reads",
    x        = "Sample (SRR suffix)",
    y        = "Read Pairs (%)",
    fill     = "Mapping Fate",
    caption  = "Unmapped = too short + too many mismatches + other"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_plot(p8, "08_star_alignment_stacked", width = 15, height = 6,
          desc = "STAR alignment fate stacked bar")

# =============================================================================
# PLOT 9 — Uniquely mapped % dot plot
# =============================================================================

cat("[9/12] STAR uniquely mapped scatter...\n")

mean_uniq <- mean(star$Uniquely_Mapped_pct)

p9 <- star %>%
  ggplot(aes(x = Short, y = Uniquely_Mapped_pct,
             colour = Group, shape = Time_Point)) +
  geom_hline(yintercept = mean_uniq,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_point(size = 4.5, stroke = 0.5) +
  geom_text_repel(aes(label = sprintf("%.1f%%", Uniquely_Mapped_pct)),
                  size = 3, max.overlaps = 20, box.padding = 0.4,
                  show.legend = FALSE) +
  annotate("text", x = 0.4, y = mean_uniq + 0.4,
           label = sprintf("Mean = %.1f%%", mean_uniq),
           hjust = 0, size = 3.2, colour = "grey30") +
  facet_grid(. ~ Time_Point, scales = "free_x", space = "free_x") +
  scale_colour_manual("Group", values = col_group) +
  scale_shape_manual("Time point", values = shape_tp) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    title    = "STAR Unique Alignment Rate per Sample",
    subtitle = paste0("Dashed line = overall mean (", sprintf("%.1f%%", mean_uniq), "); ",
                      "coloured by group; high rates indicate good reference coverage"),
    x        = "Sample (SRR suffix)",
    y        = "Uniquely Mapped (%)"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p9, "09_star_uniquely_mapped_scatter", width = 13, height = 6,
          desc = "Unique alignment rate per sample (dot plot)")

# =============================================================================
# PLOT 10 — featureCounts assignment stacked bar
# =============================================================================

cat("[10/12] featureCounts stacked bar...\n")

fc_total_cols <- c("Assigned", "Unassigned_NoFeatures",
                   "Unassigned_Ambiguity", "Unassigned_MultiMapping")

fc_stack <- fc %>%
  mutate(
    Other = pmax(0, Total_Reads - Assigned - Unassigned_NoFeatures -
                    Unassigned_Ambiguity - Unassigned_MultiMapping),
    `Assigned`           = Assigned_pct,
    `No feature`         = Unassigned_NoFeatures_pct,
    `Ambiguous`          = Unassigned_Ambiguity_pct,
    `Multi-mapping`      = Unassigned_MultiMapping_pct,
    `Other / unmapped`   = pmax(0, 100 - Assigned_pct - Unassigned_NoFeatures_pct -
                                         Unassigned_Ambiguity_pct - Unassigned_MultiMapping_pct)
  ) %>%
  select(Short, Group, Time_Point,
         `Assigned`, `No feature`, `Ambiguous`, `Multi-mapping`, `Other / unmapped`) %>%
  pivot_longer(c(`Assigned`, `No feature`, `Ambiguous`, `Multi-mapping`, `Other / unmapped`),
               names_to = "Category", values_to = "Pct") %>%
  mutate(Category = factor(Category,
                            levels = c("Other / unmapped", "Multi-mapping",
                                       "Ambiguous", "No feature", "Assigned")))

p10 <- ggplot(fc_stack, aes(x = Short, y = Pct, fill = Category)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(data = fc,
            aes(x = Short, y = Assigned_pct / 2,
                label = sprintf("%.1f%%", Assigned_pct)),
            inherit.aes = FALSE, size = 2.8, colour = "white", fontface = "bold") +
  facet_grid(. ~ Time_Point + Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(
    "Assigned"          = "#27AE60",
    "No feature"        = "#E74C3C",
    "Ambiguous"         = "#F39C12",
    "Multi-mapping"     = "#8E44AD",
    "Other / unmapped"  = "#BDC3C7"
  )) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "featureCounts Read Assignment by Category",
    subtitle = "White labels = % assigned per sample; 'No feature' = reads falling outside annotated exons",
    x        = "Sample (SRR suffix)",
    y        = "Aligned Reads (%)",
    fill     = "Assignment",
    caption  = "Parameters: -t exon -g Parent (GFF3); read pairs assigned to MEAM1 transcript features"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_plot(p10, "10_featurecounts_stacked", width = 15, height = 6,
          desc = "featureCounts read assignment stacked bar")

# =============================================================================
# PLOT 11 — featureCounts assignment % dot plot
# =============================================================================

cat("[11/12] featureCounts assignment rate scatter...\n")

mean_assigned <- mean(fc$Assigned_pct)

p11 <- fc %>%
  ggplot(aes(x = Short, y = Assigned_pct,
             colour = Group, shape = Time_Point)) +
  geom_hline(yintercept = mean_assigned,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_point(size = 4.5, stroke = 0.5) +
  geom_text_repel(aes(label = sprintf("%.1f%%", Assigned_pct)),
                  size = 3, max.overlaps = 20, box.padding = 0.4,
                  show.legend = FALSE) +
  annotate("text", x = 0.4, y = mean_assigned + 0.7,
           label = sprintf("Mean = %.1f%%", mean_assigned),
           hjust = 0, size = 3.2, colour = "grey30") +
  facet_grid(. ~ Time_Point, scales = "free_x", space = "free_x") +
  scale_colour_manual("Group", values = col_group) +
  scale_shape_manual("Time point", values = shape_tp) +
  scale_y_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title    = "featureCounts Assignment Rate per Sample",
    subtitle = paste0("Proportion of aligned reads assigned to an annotated feature; ",
                      "mean = ", sprintf("%.1f%%", mean_assigned)),
    x        = "Sample (SRR suffix)",
    y        = "Assigned Reads (%)"
  ) +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p11, "11_featurecounts_assigned_dot", width = 13, height = 6,
          desc = "featureCounts assignment rate dot plot")

# =============================================================================
# PLOT 12 — QC summary heatmap (pheatmap, all key metrics)
# =============================================================================

cat("[12/12] QC summary heatmap...\n")

# Join all key % metrics into one matrix
hm_df <- meta %>%
  left_join(trim  %>% select(Sample, Reads_passed_filter_pct,
                              Raw_Q30_pct, Clean_Q30_pct, Raw_GC_pct),
            by = "Sample") %>%
  left_join(decon %>% select(Sample, Mapped_to_contaminants_pct,
                              Retained_unmapped_pct),
            by = "Sample") %>%
  left_join(star  %>% select(Sample, Uniquely_Mapped_pct, Multi_Mapped_pct,
                              Unmapped_TooShort_pct),
            by = "Sample") %>%
  left_join(fc    %>% select(Sample, Assigned_pct, Unassigned_NoFeatures_pct,
                              Unassigned_Ambiguity_pct),
            by = "Sample")

# Sample order
hm_df <- hm_df %>%
  mutate(
    Time_Point = factor(Time_Point, levels = c("24h", "48h", "72h")),
    Group      = factor(Group,      levels = c("control", "virus"))
  ) %>%
  arrange(Time_Point, Group)

hm_mat <- hm_df %>%
  select(
    `Pass filter (%)`      = Reads_passed_filter_pct,
    `Q30 – raw (%)`        = Raw_Q30_pct,
    `Q30 – trimmed (%)`    = Clean_Q30_pct,
    `GC content (%)`       = Raw_GC_pct,
    `Contaminants (%)`     = Mapped_to_contaminants_pct,
    `Host retained (%)`    = Retained_unmapped_pct,
    `Uniquely mapped (%)`  = Uniquely_Mapped_pct,
    `Multi-mapped (%)`     = Multi_Mapped_pct,
    `Unmapped-short (%)`   = Unmapped_TooShort_pct,
    `Assigned (%)`         = Assigned_pct,
    `No feature (%)`       = Unassigned_NoFeatures_pct,
    `Ambiguous (%)`        = Unassigned_Ambiguity_pct
  ) %>%
  as.matrix()

rownames(hm_mat) <- short_id(hm_df$Sample)

# Column annotation for pheatmap (annotation on samples, here rows in original
# but we want columns = samples, rows = metrics — so transpose)
hm_mat_t <- t(hm_mat)   # rows = metrics, cols = samples

ann_col <- data.frame(
  Group      = as.character(hm_df$Group),
  `Time Point` = as.character(hm_df$Time_Point),
  row.names  = short_id(hm_df$Sample),
  check.names = FALSE
)

ann_colors <- list(
  Group        = col_group,
  `Time Point` = col_time
)

# Row annotation: pipeline stage grouping
stage_ann <- data.frame(
  Stage = c("Trimming", "Trimming", "Trimming", "Trimming",
            "Decontamination", "Decontamination",
            "STAR alignment", "STAR alignment", "STAR alignment",
            "Quantification", "Quantification", "Quantification"),
  row.names = rownames(hm_mat_t)
)
ann_colors$Stage <- c(
  "Trimming"         = "#3498DB",
  "Decontamination"  = "#E67E22",
  "STAR alignment"   = "#27AE60",
  "Quantification"   = "#9B59B6"
)

ph12 <- pheatmap(
  hm_mat_t,
  scale            = "row",        # z-score per metric so colour encodes deviation
  cluster_rows     = FALSE,        # keep pipeline order
  cluster_cols     = TRUE,         # cluster samples to reveal outliers
  annotation_col   = ann_col,
  annotation_row   = stage_ann,
  annotation_colors = ann_colors,
  color            = colorRampPalette(c("#2980B9", "white", "#C0392B"))(100),
  border_color     = "grey85",
  fontsize         = 10,
  fontsize_row     = 10,
  fontsize_col     = 10,
  cellwidth        = 32,
  cellheight       = 26,
  display_numbers  = round(hm_mat_t, 1),
  number_format    = "%.1f",
  number_color     = "grey20",
  fontsize_number  = 7,
  main             = "QC Metrics Summary — All Pipeline Steps\n(colour = row z-score; numbers = raw values)",
  angle_col        = 45,
  silent           = TRUE
)

save_pheatmap(ph12, "12_qc_summary_heatmap", width = 14, height = 8,
              desc = "Pheatmap of all key QC metrics across samples")

# =============================================================================
# Print manifest
# =============================================================================

cat("\n=== Plots Complete:", format(Sys.time()), "===\n\n")
cat("Output directory:", out_dir, "\n\n")
cat("Files generated:\n")
for (entry in plot_log) cat("  ", entry, "\n")
cat("\n")

# Write manifest to file
manifest <- data.frame(
  File        = sub("  .*", "", trimws(plot_log)),
  Description = sub(".*  ", "", trimws(plot_log)),
  stringsAsFactors = FALSE
)
write.csv(manifest,
          file.path(out_dir, "00_plot_manifest.csv"),
          row.names = FALSE)
cat("Plot manifest written to: 00_plot_manifest.csv\n\n")

#!/usr/bin/env Rscript

# Nano-FIRSTseq: dcDNA vs FIRST-seq Comparison
# Usage: ./2_dcDNA_FIRST.R <dcDNA_stats.tsv> <FIRST_stats.tsv> <mod_positions.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || any(args %in% c("-h", "--help"))) {
  cat("Usage:\n")
  cat("  ./2_dcDNA_FIRST.R dcDNA.STATS FIRST.STATS mod_positions.tsv\n")
  quit(save = "no", status = 1)
}

dcDNA_file    <- args[1]
firstseq_file <- args[2]
mod_file      <- args[3]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggbeeswarm)
})

output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)

order_mods <- c("m1A", "I", "m3C", "m1G", "m22G", "m1acp3Y", "m3U")
mod_positions <- read.delim(mod_file, sep = "", header = TRUE, comment.char = "#")

# ---------------------------
# Function: Process stats file
# ---------------------------
process_func <- function(stats_file_path, method) {
  stats_file <- read.delim(stats_file_path, sep = "")
  stats_file <- stats_file %>% filter(pos > 20)

  bases <- c("A", "T", "C", "G")
  mismatch <- numeric(nrow(stats_file))
  for (b in bases) {
    mismatch <- mismatch + ifelse(stats_file$ref_nuc != b, stats_file[[b]], 0)
  }

  stats_file <- stats_file %>%
    mutate(
      mismatch   = mismatch,
      mis_freq   = round(mismatch / coverage, 3),
      Method     = method,
      uniq_coord = paste(chr, pos, sep = "_")
    )

  merged <- merge(stats_file, mod_positions,
                  by.x = c("chr", "pos"), by.y = c("Reference", "Position"),
                  all = TRUE)

  merged <- merged %>%
    mutate(
      Unique = replace_na(Unique, 'Unm'),
      Mod    = ifelse(is.na(Mod), ref_nuc, Mod)
    )

  return(merged)
}

# ---------------------------
# Function: Aggregate mismatch frequencies
# ---------------------------
clean_mis_input_agg <- function(data) {
  data <- data %>% filter(ref_nuc %in% c("A", "T", "C", "G")) %>%
    mutate(sum = A + T + C + G,
           A = A / sum, T = T / sum, C = C / sum, G = G / sum)

  aggregate(data[, c("A", "T", "C", "G")],
            by = list(Mod = data$Mod, ref_base = data$ref_nuc),
            FUN = mean, na.rm = TRUE)
}

# ---------------------------
# Function: Barplot
# ---------------------------
mismatch_profile_barplot_facet <- function(data, label) {
  base_colors <- c("A" = "#1fab89", "T" = "#eb4d55", "C" = "#1e56a0", "G" = "#f0cf85")
  ref_color   <- "#888888"

  data_long <- pivot_longer(data, cols = c("A", "T", "C", "G"),
                            names_to = "base", values_to = "value")

  data_long <- data_long %>%
    filter(Mod %in% order_mods) %>%
    mutate(
      fill = ifelse(base == ref_base, "ref", base),
      fill = factor(fill, levels = c("A", "T", "C", "G", "ref"))
    ) %>%
    group_by(Mod) %>%
    arrange(fill, .by_group = TRUE) %>%
    ungroup()

  data_long$Mod <- factor(data_long$Mod, levels = order_mods)
  color_map <- c("ref" = ref_color, base_colors)

  pdf(file = file.path(output_dir, paste0(label, "_mismatch_profile_barplot.pdf")), height = 1.5, width = 6)
  print(
    ggplot(data_long, aes(x = Mod, y = value, fill = fill)) +
      geom_col(position = "stack", colour = "black") +
      scale_fill_manual(values = color_map, drop = FALSE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

# ---------------------------
# Function: Dotplot
# ---------------------------
dotplot_facet <- function(data, label) {
  data <- data %>% filter(Mod %in% order_mods)
  data$Mod <- factor(data$Mod, levels = order_mods)

  pdf(file = file.path(output_dir, paste0(label, "_Dotplot.pdf")), height = 3, width = 15)
  print(
    ggplot(data, aes(x = Method, y = mis_freq, color = Method)) +
      geom_quasirandom(varwidth = TRUE, alpha = 0.6) +
      geom_boxplot(width = 0.3, size = 0.5, alpha = 1) +
      scale_color_manual(values = c("dcDNA-seq" = "#006989", "FIRST-seq" = "#7C5591")) +
      facet_wrap(~ Mod, nrow = 1, drop = TRUE) +
      theme_classic()
  )
  dev.off()
}

# ---------------------------
# Function: RT Drop / Coverage Plot
# ---------------------------
plot_area_coverage_facet <- function(merged, name, palette) {
  for (chro in c("18s_rRNA", "25s_rRNA")) {
    subs <- subset(merged, chr == chro)
    subs$pos <- subs$pos - 14
    subs$position <- paste(subs$chr, subs$pos, sep = "_")

    pdf(file = file.path(output_dir, paste0(chro, "_", name, "_coverage_facet.pdf")), height = 3, width = 7)
    print(
      ggplot(subs, aes(x = pos, y = coverage / max(coverage, na.rm = TRUE), colour = Method)) +
        geom_line() +
        geom_area(aes(fill = Method, group = Method), alpha = 1/3, position = 'identity') +
        theme_classic() +
        scale_fill_manual(values = c("dcDNA-seq" = "#006989", "FIRST-seq" = "#7C5591")) +
        scale_colour_manual(values = c("dcDNA-seq" = "#006989", "FIRST-seq" = "#7C5591")) +
        geom_vline(data = subset(subs, position %in% c("18s_rRNA_1191", "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843")),
                   aes(xintercept = pos), linetype = "dashed")
    )
    dev.off()
  }
}

# ---------------------------
# Run All
# ---------------------------
dcDNA.processed    <- process_func(dcDNA_file,    "dcDNA-seq")
firstseq.processed <- process_func(firstseq_file, "FIRST-seq")

dcDNA.mis    <- clean_mis_input_agg(dcDNA.processed)
firstseq.mis <- clean_mis_input_agg(firstseq.processed)

mismatch_profile_barplot_facet(dcDNA.mis, "dcDNA-seq")
mismatch_profile_barplot_facet(firstseq.mis, "FIRST-seq")

dotplot_facet(rbind(dcDNA.processed, firstseq.processed), "Modifications_Mismatch")

plot_area_coverage_facet(
  rbind(dcDNA.processed, firstseq.processed),
  "Firstseq_dcDNAseq",
  c("dcDNA-seq" = "#006989", "FIRST-seq" = "#7C5591")
)

#!/usr/bin/env Rscript

# ---------------------------
# Nano-FIRSTseq: dRNA vs dcDNA Comparison
# CLI Executable R Script (Final Version)
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3 || any(args %in% c("--help", "-h"))) {
  cat("Usage: ./1_dRNA_dcDNA.R <dRNA_stats.tsv> <dcDNA_stats.tsv> <mod_positions.tsv>\n")
  cat("Example: ./1_dRNA_dcDNA.R data/stats/dRNA_YeastTotalRNA_tRNA.STATS data/stats/dcDNA_Maxima.STATS data/references/yeast_mod_positions.tsv\n")
  quit(save = "no", status = 1)
}

dRNA_input  <- args[1]
dcDNA_input <- args[2]
mod_file    <- args[3]

suppressPackageStartupMessages({
  library(reshape2)
  library(EnvStats)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(ggbeeswarm)
})

# Output directory
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)

# Color definitions
dRNA_col  <- "#E88D67"
dcDNA_col <- "#006989"

# Load mod positions
mod_positions <- read.delim(mod_file)

# Reference modification order
order_mods <- c("A", "C", "G", "T", "Am", "i6A", "m66A", "t6A", "ac4C", "Cm",
                "m5C", "Gm", "m2G", "m7G", "m5U", "Um", "Y", "m1A", "I", "m3C",
                "m1G", "m22G", "m1acp3Y", "m3U")

# ---------------------------
# Function: Process Input
# ---------------------------
process_func <- function(stats_file_path, method) {
  stats_file <- read.delim(stats_file_path, sep = "")
  stats_file <- stats_file %>%
    filter(pos > 20)

  # Create mismatch by subtracting reference base count
  bases <- c("A", "T", "C", "G")
  mismatch <- numeric(nrow(stats_file))

  for (b in bases) {
    mismatch <- mismatch + ifelse(stats_file$ref_nuc != b, stats_file[[b]], 0)
  }

  stats_file <- stats_file %>%
    mutate(
      mismatch = mismatch,
      mis_freq = round(mismatch / coverage, 3),
      uniq_coord = paste(chr, pos, sep = "_"),
      Method = method
    )

  merged <- merge(stats_file, mod_positions,
                  by.x = c("chr", "pos"), by.y = c("Reference", "Position"),
                  all = TRUE) %>%
    mutate(
      Unique = replace_na(Unique, 'Unm'),
      Mod = replace_na(Mod, 'Unm'),
      Mod = if_else(Mod == "Unm", ref_nuc, Mod)
    )

  return(merged)
}


# ---------------------------
# Function: Aggregate mismatch frequencies
# ---------------------------
clean_mis_input_agg <- function(data) {
  data <- data %>%
    filter(ref_nuc %in% c("A", "T", "C", "G")) %>%
    mutate(sum = A + T + C + G,
           A = A / sum, T = T / sum, C = C / sum, G = G / sum)

  data_count <- count(data, Mod)

  data_avg <- aggregate(data[, c("A", "T", "C", "G")],
                        by = list(Mod = data$Mod, ref_base = data$ref_nuc),
                        FUN = mean, na.rm = TRUE)

  final_data <- merge(data_avg, data_count, by = "Mod")
  return(final_data)
}

# ---------------------------
# Function: Mismatch Barplot (Ref at Bottom)
# ---------------------------
mismatch_profile_barplot_facet <- function(data, label) {
  base_colors <- c("A" = "#1fab89", "T" = "#eb4d55", "C" = "#1e56a0", "G" = "#f0cf85")
  ref_color   <- "#888888"

  data_long <- pivot_longer(data, cols = c("A", "T", "C", "G"),
                            names_to = "base", values_to = "value")

  #   Filter only ordered mods
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
  # Filter only known modifications
  data <- data %>%
    filter(Mod %in% order_mods)

  # Enforce order for facetting
  data$Mod <- factor(data$Mod, levels = order_mods)

  pdf(file = file.path(output_dir, paste0(label, "_Dotplot.pdf")), height = 3, width = 15)
  print(
    ggplot(data, aes(x = Method, y = mis_freq, color = Method)) +
      geom_quasirandom(varwidth = TRUE, alpha = 0.6) +
      geom_boxplot(width = 0.3, size = 0.5, alpha = 1) +
      scale_color_manual(values = c(dcDNA_col, dRNA_col)) +
      facet_wrap(~ Mod, nrow = 1, drop = FALSE) +
      theme_classic()
  )
  dev.off()
}

# ---------------------------
# Run Analysis
# ---------------------------
dRNA.processed <- process_func(dRNA_input, "dRNA-seq")
cDNA.processed <- process_func(dcDNA_input, "dcDNA-seq")

dRNA.mis <- clean_mis_input_agg(dRNA.processed)
cDNA.mis <- clean_mis_input_agg(cDNA.processed)

mismatch_profile_barplot_facet(dRNA.mis, "dRNA-seq")
mismatch_profile_barplot_facet(cDNA.mis, "dcDNA-seq")

dotplot_facet(rbind(dRNA.processed, cDNA.processed), "Modifications_Mismatch")

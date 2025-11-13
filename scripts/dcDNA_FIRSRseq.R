#!/usr/bin/env Rscript

# Nano-FIRSTseq: dcDNA vs FIRST-seq Comparison

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || any(args %in% c("--help", "-h"))) {
  cat("Usage: ./2_dcDNA_FIRST.R <dcDNA_stats.tsv> <FIRST_stats.tsv> <mod_positions.tsv>\n")
  quit(save = "no", status = 1)
}

dcDNA_input     <- args[1]
firstseq_input  <- args[2]
mod_file        <- args[3]

suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(reshape2)
})

# Output directory
dir.create("figures", showWarnings = FALSE)

# Load mod positions
mod_positions <- read.delim(mod_file, sep = "", header = TRUE, comment.char = "#")

mod_positions_WC <- mod_positions %>%
  filter(Mod %in% c("m1A","I","m3C","m1G","m22G","m1acp3Y","m3U")) %>%
  mutate(Position = as.numeric(Position))

dcDNA_col     <- "#006989"
firstseq_col  <- "#7C5591"
palette       <- c(dcDNA_col, firstseq_col)

mod <- data.frame(
  ModName = c("Y", "m1acp3Y", "m66A", "m66A", "m1A", "m1A", "m3U", "m3U"),
  position = c("18s_rRNA_999", "18s_rRNA_1191", "18s_rRNA_1781", "18s_rRNA_1782", 
               "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843"),
  Chr = c("18s_rRNA", "18s_rRNA", "18s_rRNA", "18s_rRNA", "25s_rRNA", "25s_rRNA", "25s_rRNA", "25s_rRNA"),
  pos_num = c(999, 1191, 1781, 1782, 645, 2142, 2634, 2843)
)

# Function: Process input stats
process_func <- function(stats_file_path, method) {
  stats_file <- read.delim(stats_file_path, sep = "")

  bases <- c("A", "T", "C", "G")
  mis_all <- vector()

  for (base in unique(stats_file$ref_nuc)) {
    subs <- subset(stats_file, ref_nuc == base)
    mis_bases <- bases[bases != base]
    subs$mismatch <- subs[, mis_bases[1]] + subs[, mis_bases[2]] + subs[, mis_bases[3]]
    mis_all <- rbind(mis_all, subs)
  }

  processed_data <- mis_all %>%
    filter(pos > 20) %>%
    mutate(Biotype = sub(".*_", "", chr), 
           Method = method, 
           uniq_coord = paste(chr, pos, sep = "_"),
           mis_freq = round(mismatch / coverage, 3)) %>%
    merge(mod_positions, by.x = c("chr", "pos"), by.y = c("Reference", "Position"), all = TRUE) %>%
    mutate(Unique = tidyr::replace_na(Unique, 'Unm'),
           Mod = tidyr::replace_na(Mod, 'Unm'),
           Mod = if_else(Mod == "Unm", ref_nuc, Mod))

  return(processed_data)
}

# Function: Clean for barplot
clean_mis_input_agg <- function(data) {
  data <- data %>%
    filter(ref_nuc %in% c("A", "T", "C", "G")) %>%
    mutate(sum = A + T + C + G,
           A = A / sum, T = T / sum, C = C / sum, G = G / sum)

  data_count <- count(data, "Mod")
  data_avg <- aggregate(data[, c("A", "T", "C", "G")],
                        by = list(data$Mod, data$ref_nuc),
                        FUN = mean, na.rm = TRUE)

  data_melted <- melt(data_avg, id.vars = c("Group.1", "Group.2"))
  final_data <- merge(data_melted, data_count, by.x = "Group.1", by.y = "Mod")
  return(final_data)
}

# Function: Barplot
mismatch_profile_barplot_facet <- function(data, label) {
  base_colors <- c("A" = "#1fab89", "T" = "#eb4d55", "C" = "#1e56a0", "G" = "#f0cf85")
  ref_color   <- "#888888"

  colnames(data)[1:2] <- c("Mod", "ref_base")

  data <- data %>%
    filter(Mod %in% unique(Mod)) %>%
    mutate(
      ref_base = as.character(ref_base),
      fill = ifelse(variable == ref_base, "ref", variable),
      fill = factor(fill, levels = c("A", "T", "C", "G", "ref"))  # ref at top
    )

  color_map <- c("ref" = ref_color, base_colors)

  pdf(file = paste0("figures/", label, "_mismatch_profile_barplot.pdf"), height = 1.5, width = 6)
  print(
    ggplot(data, aes(x = Mod, y = value, fill = fill)) +
      geom_col(position = "stack", colour = "black") +
      scale_fill_manual(values = color_map, drop = FALSE) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

# Function: RT drop plot
plot_area_coverage_facet <- function(merged, name, palette) {
  for (chro in c("18s_rRNA", "25s_rRNA")) {
    subs <- subset(merged, chr == chro)
    subs$pos <- subs$pos - 14
    subs$position <- paste(subs$chr, subs$pos, sep = "_")

    pdf(file = paste0("figures/", chro, "_", name, "_coverage_facet.pdf"), height = 6, width = 7)
    print(
      ggplot(subs, aes(x = pos, y = norm_cov, colour = Enzyme)) +
        geom_line() +
        geom_area(aes(fill = Enzyme, group = Enzyme), alpha = 1/3, position = 'identity') +
        theme_classic() +
        scale_fill_manual(values = palette) +
        scale_colour_manual(values = palette) +
        geom_vline(data = subset(subs, position %in% c("18s_rRNA_1191", "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843")),
                   aes(xintercept = pos), linetype = "dashed")
    )
    dev.off()
  }
}

# RUN
dcDNA.processed    <- process_func(dcDNA_input,    "dcDNA-seq")
firstseq.processed <- process_func(firstseq_input, "FIRST-seq")

dcDNA.mis    <- clean_mis_input_agg(dcDNA.processed)
firstseq.mis <- clean_mis_input_agg(firstseq.processed)

mismatch_profile_barplot_facet(dcDNA.mis, "dcDNA-seq")
mismatch_profile_barplot_facet(firstseq.mis, "FIRST-seq")

runs_forcoverage <- rbind(dcDNA.processed, firstseq.processed)
plot_area_coverage_facet(runs_forcoverage, "Firstseq_dcDNAseq", palette)

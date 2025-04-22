#!/usr/bin/env Rscript

# Nano-FIRSTseq: Various RT Enzymes Comparison
# Usage: ./3_Various_Enzymes.R <data_folder> <mod_positions.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || any(args %in% c("-h", "--help"))) {
  cat("Usage:\n")
  cat("  ./3_Various_Enzymes.R <data_folder> <mod_positions.tsv>\n")
  cat("  Example: ./3_Various_Enzymes.R data/ references/yeast_mod_positions.tsv\n")
  quit(save = "no", status = 1)
}

# Paths
data_dir     <- args[1]
mod_file     <- args[2]

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggbeeswarm)
  library(plyr)
  library(ggrepel)
  library(ggcorrplot)
})

# Config
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)

# Enzyme + Buffer combinations
enzymes <- c("PS2", "SS2", "SS3", "SS4", "Maxima", "TGIRT", "Induro", "Marathon")
buffers <- c("Mg", "Mn")

# Mod order
order_mods <- c("m1A", "I", "m3C", "m1G", "m22G", "m1acp3Y", "m3U")
mod_positions <- read.delim(mod_file, sep = "", header = TRUE, comment.char = "#")

# ---------------------------
# PART 1: Processing Functions
# ---------------------------

# Load & normalize chromosome-level counts
process_count <- function(input, enzyme, buffer) {
  data <- read.delim(input, header = FALSE)
  data2 <- data[, c("V1", "V3")]
  colnames(data2) <- c("Chr", "Count")
  data2$Count <- ((data2$Count) / sum(data2$Count)) * 1e6
  data2$Enzyme <- enzyme
  data2$Buffer <- buffer
  return(data2)
}

# Load & process modification mismatch tables
process_func <- function(stats_path, readends_path, enzyme, buffer) {
  stats_file <- read.delim(stats_path, sep = "")
  readends_file <- read.delim(readends_path, sep = "")

  merged <- merge(stats_file, readends_file[c("Ref", "Pos", "Ends")],
                  by.x = c("chr", "pos"), by.y = c("Ref", "Pos"))

  clean <- merged[, c("chr", "pos", "ref_nuc", "coverage", "rtstop", "Ends", "ins", "del", "A", "T", "C", "G")]

  bases <- c("A", "T", "C", "G")
  mis_all <- vector("list", length = 0)
  for (base in unique(clean$ref_nuc)) {
    subs <- subset(clean, ref_nuc == base)
    mis_bases <- bases[bases != base]
    subs$mismatch <- rowSums(subs[, mis_bases], na.rm = TRUE)
    mis_all[[base]] <- subs
  }
  mis_all <- bind_rows(mis_all) %>% filter(pos > 20)
  mis_all$Biotype <- gsub(".*_", "", mis_all$chr)
  mis_all$Enzyme <- enzyme
  mis_all$Buffer <- buffer
  mis_all$uniq_coord <- paste(mis_all$chr, mis_all$pos, sep = "_")

  merged_final <- lapply(unique(mis_all$chr), function(ref) {
    subs <- subset(mis_all, chr == ref)
    subs_3end <- subset(subs, pos > nrow(subs) / 2)
    subs$norm_cov <- round(subs$coverage / max(subs_3end$coverage), 3)
    subs
  }) %>% bind_rows()

  merged_final$rtstop <- as.numeric(merged_final$rtstop)
  merged_final$norm_rt <- round(merged_final$rtstop / merged_final$coverage, 3)
  merged_final$norm_rt_end <- round(merged_final$Ends / merged_final$coverage, 3)
  merged_final$mis_freq <- round(merged_final$mismatch / merged_final$coverage, 3)

  merged2 <- merge(merged_final, mod_positions,
                   by.x = c("chr", "pos"), by.y = c("Reference", "Position"),
                   all = TRUE)
  merged2$Unique <- replace_na(merged2$Unique, "Unm")
  merged2$Mod <- replace_na(merged2$Mod, "Unm")

  # Relabel unmodified bases with their reference
  unm <- subset(merged2, Mod == "Unm")
  unm$Mod <- unm$ref_nuc
  mod <- subset(merged2, Mod != "Unm")
  final <- rbind(unm, mod)
  return(final)
}

# ---------------------------
# PART 2: Run Analysis
# ---------------------------

all_data <- list()
for (e in enzymes) {
  for (b in buffers) {
    stat_path <- file.path(data_dir, "stats", paste0("FIRST_", b, "_", e, ".STATS"))
    readend_path <- file.path(data_dir, "readends", paste0("FIRST_", b, "_", e, ".readends.tsv"))
    if (file.exists(stat_path) && file.exists(readend_path)) {
      df <- process_func(stat_path, readend_path, enzyme = e, buffer = b)
      all_data[[paste(e, b, sep = "_")]] <- df
    }
  }
}

data_combined <- bind_rows(all_data) %>% filter(!is.na(Enzyme))

# Normalize norm_cov individually per Enzyme/Buffer
library(forcats)
data_combined <- data_combined %>%
  group_by(Enzyme, Buffer, chr) %>%
  mutate(norm_cov = norm_cov / max(norm_cov, na.rm = TRUE)) %>%
  ungroup()

# Palette
enzyme_palette <- c(
  "PS2" = "#aa6f73", "SS2" = "#eea990", "SS3" = "#f6e0b5", "SS4" = "#D193C1",
  "Maxima" = "#9C8CC3", "TGIRT" = "#5294a3", "Induro" = "#a3d9d9", "Marathon" = "#60bfae"
)

# Faceted RT coverage plot (one per Buffer)
plot_area_coverage_facet <- function(data, label, palette) {
  for (chro in c("18s_rRNA", "25s_rRNA")) {
    subs <- subset(data, chr == chro)
    subs$pos <- subs$pos - 14
    subs$position <- paste(subs$chr, subs$pos, sep = "_")
    subs$Group <- interaction(subs$Enzyme, subs$Buffer)

    pdf(file = file.path(output_dir, paste0(chro, "_", label, "_coverage_facet.pdf")), height = 10, width = 7)
    print(
      ggplot(subs, aes(x = pos, y = norm_cov, color = Group)) +
        geom_line() +
        geom_area(aes(fill = Group, group = Group), alpha = 1/3, position = "identity") +
        theme_classic() +
        scale_fill_manual(values = rep(palette, each = 2)) +
        scale_color_manual(values = rep(palette, each = 2)) +
        geom_vline(data = subset(subs, position %in% c("18s_rRNA_1191", "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843")),
                   aes(xintercept = pos), linetype = "dashed") +
        facet_wrap(~Enzyme, nrow = length(unique(subs$Enzyme)))
    )
    dev.off()
  }
}

# Generate plots separately for Mg and Mn
plot_area_coverage_facet(subset(data_combined, Buffer == "Mg"), "Various_Enzymes_Mg", enzyme_palette)
plot_area_coverage_facet(subset(data_combined, Buffer == "Mn"), "Various_Enzymes_Mn", enzyme_palette)

# ---------------------------
# PART 3: Mismatch Barplots and Dotplots
# ---------------------------

# Mismatch frequency aggregation
clean_mis_input_agg <- function(data) {
  data <- data %>% filter(ref_nuc %in% c("A", "T", "C", "G")) %>%
    mutate(sum = A + T + C + G,
           A = A / sum, T = T / sum, C = C / sum, G = G / sum)

  data_avg <- aggregate(data[, c("A", "T", "C", "G")],
                        by = list(Mod = data$Mod, ref_base = data$ref_nuc, Enzyme = data$Enzyme, Buffer = data$Buffer),
                        FUN = mean, na.rm = TRUE)

  return(data_avg)
}

# Barplot
mismatch_profile_barplot_facet <- function(data, label) {
  base_colors <- c("A" = "#1fab89", "T" = "#eb4d55", "C" = "#1e56a0", "G" = "#f0cf85")
  ref_color   <- "#888888"

  data_long <- pivot_longer(data, cols = c("A", "T", "C", "G"), names_to = "base", values_to = "value") %>%
    filter(Mod %in% order_mods) %>%
    mutate(
      fill = ifelse(base == ref_base, "ref", base),
      fill = factor(fill, levels = c("A", "T", "C", "G", "ref")),
      Mod = factor(Mod, levels = order_mods)
    )

  pdf(file = file.path(output_dir, paste0(label, "_mismatch_profile_barplot.pdf")), height = 6, width = 8)
  print(
    ggplot(data_long, aes(x = Mod, y = value, fill = fill)) +
      geom_col(position = "stack", colour = "black") +
      scale_fill_manual(values = c("ref" = ref_color, base_colors)) +
      facet_grid(Buffer ~ Enzyme) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}

# Dotplot
plot_dotplot <- function(data, label) {
  data <- data %>% filter(Mod %in% order_mods)
  data$Mod <- factor(data$Mod, levels = order_mods)

  pdf(file = file.path(output_dir, paste0(label, "_Dotplot.pdf")), height = 3, width = 15)
  print(
    ggplot(data, aes(x = interaction(Enzyme, Buffer), y = mis_freq, color = Buffer)) +
      geom_quasirandom(varwidth = TRUE, alpha = 0.6) +
      geom_boxplot(width = 0.3, size = 0.5, alpha = 1, outlier.shape = NA) +
      facet_wrap(~ Mod, nrow = 1) +
      theme_classic() +
      labs(x = "Enzyme_Buffer", y = "Mismatch Frequency")
  )
  dev.off()
}

# Compute and plot
mis_agg <- clean_mis_input_agg(data_combined)
mismatch_profile_barplot_facet(mis_agg, "Various_Enzymes")
plot_dotplot(data_combined, "Various_Enzymes")



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
# PART 2: Plotting Functions
# ---------------------------

# (plot_barplot, plot_dotplot, plot_coverage, etc. to follow)

# ---------------------------
# Run Processing
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

# Ready to call plotting functions on data_combined
print(unique(data_combined$Enzyme))
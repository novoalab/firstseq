# Nano-FIRSTseq: dRNA vs dcDNA Comparison
# R Analysis Script
# Author: Oguzhan Begik
# Description: This script analyzes mismatch frequencies from dRNA-seq and dcDNA-seq datasets.

# Load necessary libraries
library(reshape2)
library(EnvStats)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(viridis)
library(ggbeeswarm)

# Define colors for plots
dRNA_col <- "#E88D67"  # Orange for dRNA-seq
dcDNA_col <- "#006989"  # Blue for dcDNA-seq

# Load known modification positions
mod_positions <- read.delim("data/references/yeast_mod_positions.tsv")

# Define modification order lists for analysis
order_mods <- c("A", "C", "G", "T", "Am", "i6A", "m66A", "t6A", "ac4C", "Cm", 
                "m5C", "Gm", "m2G", "m7G", "m5U", "Um", "Y", "m1A", "I", "m3C", 
                "m1G", "m22G", "m1acp3Y", "m3U")
order_mods_WConly <- c("m1A", "I", "m3C", "m1G", "m22G", "m1acp3Y", "m3U")

# Function: Process input statistics file
# Reads sequencing statistics, calculates mismatch frequencies, and merges with known modification positions
process_func <- function(stats_file_path, method) {
    stats_file <- read.delim(stats_file_path, sep="")
    
    # Filter out low-quality positions (first 20 bases)
    stats_file <- stats_file %>%
        filter(pos > 20) %>%
        mutate(
            mismatch = rowSums(stats_file[, c("A", "T", "C", "G")] - stats_file$ref_nuc, na.rm = TRUE),
            mis_freq = round(mismatch / coverage, 3),
            uniq_coord = paste(chr, pos, sep="_"),
            Method = method
        )

    # Merge with known modification positions
    merged <- merge(stats_file, mod_positions, 
                    by.x = c("chr", "pos"), by.y = c("Reference", "Position"), 
                    all = TRUE)
    
    # Replace missing modification data with 'Unm' (unmodified)
    merged <- merged %>%
        mutate(Unique = replace_na(Unique, 'Unm'),
               Mod = replace_na(Mod, 'Unm')) %>%
        mutate(Mod = if_else(Mod == "Unm", ref_nuc, Mod))
    
    return(merged)
}

# Process dRNA-seq and dcDNA-seq data
dRNA.processed <- process_func("data/stats/dRNA_YeastTotalRNA_tRNA.STATS", "dRNA-seq")
cDNA.processed <- process_func("data/stats/dcDNA_Maxima.STATS", "dcDNA-seq")

# Function: Aggregate mismatch frequencies
# Normalizes mismatch counts and computes mean mismatch frequencies per modification type
clean_mis_input_agg <- function(data) {
    # Filter for valid reference nucleotides
    data <- data %>% filter(ref_nuc %in% c("A", "T", "C", "G"))
    
    # Normalize mismatch counts
    data <- data %>%
        mutate(sum = rowSums(select(data, A, T, C, G)),
               A = A / sum, T = T / sum, C = C / sum, G = G / sum)
    
    # Count occurrences of each modification type
    data_count <- count(data, Mod)
    
    # Compute mean mismatch frequencies per modification
    data_avg <- aggregate(data[, c("A", "T", "C", "G")], 
                          by = list(data$Mod, data$ref_nuc), 
                          FUN = mean, na.rm = TRUE)
    
    # Melt data for visualization
    data_melted <- melt(data_avg, id.vars = c("Group.1", "Group.2"))
    
    # Merge with count data
    final_data <- merge(data_melted, data_count, 
                        by.x = "Group.1", by.y = "Mod")
    
    return(final_data)
}

# Compute mismatch profiles
dRNA.mis <- clean_mis_input_agg(dRNA.processed)
cDNA.mis <- clean_mis_input_agg(cDNA.processed)

# Function: Generate mismatch profile barplots
# Creates and saves a barplot of mismatch frequencies for each dataset
mismatch_profile_barplot_facet <- function(data, label) {
    pdf(file = paste0("figures/", label, "_mismatch_profile_barplot.pdf"), 
        height = 1.5, width = 6)
    
    print(
        ggplot(data, aes(x = Group.1, y = value, fill = variable)) +
            scale_fill_manual(values = c("#1fab89", "#eb4d55", "#1e56a0", "#f0cf85", "#888888")) +
            geom_bar(stat = 'identity', colour = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    
    dev.off()
}

# Generate barplots for dRNA-seq and dcDNA-seq
mismatch_profile_barplot_facet(dRNA.mis, "dRNA-seq")
mismatch_profile_barplot_facet(cDNA.mis, "dcDNA-seq")

# Function: Generate dot plots of modifications
# Creates dot plots showing mismatch frequencies per modification type
dotplot_facet <- function(data, label) {
    # Filter out standard nucleotides and modifications not relevant for analysis
    data <- data %>% filter(!(Mod %in% c("A", "G", "C", "T", "a", "t", "g", "c", 
                                         "N", "xU", "mcm5s2U", "ncm5U", "mcm5U", "m1I")))
    
    # Order modifications according to predefined order
    data$Mod <- factor(data$Mod, levels = order_mods)
    
    # Save dot plot as PDF
    pdf(file = paste0("figures/", label, "_Dotplot.pdf"), height = 3, width = 15)
    
    print(
        ggplot(data, aes(x = Method, y = mis_freq, color = Method)) +
            geom_quasirandom(varwidth = TRUE, alpha = 0.6) +
            geom_boxplot(width = 0.3, size = 0.5, alpha = 1) +
            scale_color_manual(values = c(dcDNA_col, dRNA_col)) +
            facet_wrap(~ Mod, nrow = 1) +
            theme_classic()
    )
    
    dev.off()
}

# Generate dot plots for mismatch modifications
dotplot_facet(rbind(dRNA.processed, cDNA.processed), "Modifications_Mismatch")

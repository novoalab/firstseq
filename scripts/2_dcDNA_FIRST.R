# Nano-FIRSTseq: dcDNA vs FIRST-seq Comparison
# R Analysis Script
# Author: Oguzhan Begik
# Description: This script compares mismatch frequencies and RT drop patterns between dcDNA-seq and FIRST-seq.

# Load necessary libraries
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)

# Load known modification positions
mod_positions <- read.delim("data/stats/mod_positions.tsv")

# Filter for Watson-Crick base modifications
mod_positions_WC <- mod_positions %>%
    filter(Mod %in% c("m1A","I","m3C","m1G","m22G","m1acp3Y","m3U")) %>%
    mutate(Position = as.numeric(Position))

# Define colors for plots
dcDNA_col <- "#006989"      # Blue for dcDNA-seq
firstseq_col <- "#7C5591"   # Purple for FIRST-seq
palette <- c(dcDNA_col, firstseq_col)

# Create a table for modifications on rRNAs
mod <- data.frame(
    ModName = c("Y", "m1acp3Y", "m66A", "m66A", "m1A", "m1A", "m3U", "m3U"),
    position = c("18s_rRNA_999", "18s_rRNA_1191", "18s_rRNA_1781", "18s_rRNA_1782", 
                 "25s_rRNA_645", "25s_rRNA_2142", "25s_rRNA_2634", "25s_rRNA_2843"),
    Chr = c("18s_rRNA", "18s_rRNA", "18s_rRNA", "18s_rRNA", "25s_rRNA", "25s_rRNA", "25s_rRNA", "25s_rRNA"),
    pos_num = c(999, 1191, 1781, 1782, 645, 2142, 2634, 2843)
)

# Function: Process input statistics for mismatch analysis
process_func <- function(stats_file_path, method) {
    # Read input statistics file
    stats_file <- read.delim(stats_file_path, sep="")
    
    # Define nucleotide bases
    bases <- c("A", "T", "C", "G")
    
    # Calculate mismatches
    mis_all <- vector()
    for (base in unique(stats_file$ref_nuc)) {
        subs <- subset(stats_file, ref_nuc == base)
        mis_bases <- bases[bases != base]
        subs$mismatch <- subs[, mis_bases[1]] + subs[, mis_bases[2]] + subs[, mis_bases[3]]
        mis_all <- rbind(mis_all, subs)
    }

    # Process data and merge with known modifications
    processed_data <- mis_all %>%
        filter(pos > 20) %>%
        mutate(Biotype = sub(".*_", "", chr), 
               Method = method, 
               uniq_coord = paste(chr, pos, sep="_"),
               mis_freq = round(mismatch / coverage, 3)) %>%
        merge(mod_positions, by.x = c("chr", "pos"), by.y = c("Reference", "Position"), all = TRUE) %>%
        mutate(Unique = replace_na(Unique, 'Unm'),
               Mod = replace_na(Mod, 'Unm'),
               Mod = if_else(Mod == "Unm", ref_nuc, Mod))  # Relabel unmodified bases
    
    return(processed_data)
}

# Process dRNA-seq and FIRST-seq data
dcDNA.processed <- process_func("data/stats/dcDNA_Maxima.STATS", "dcDNA-seq")
firstseq.processed <- process_func("data/stats/FIRST_Maxima_Short.STATS", "FIRST-seq")

# Function: Aggregate mismatch frequencies
clean_mis_input_agg <- function(data) {
    data <- data %>% filter(ref_nuc %in% c("A", "T", "C", "G")) %>%
        mutate(sum = A + T + C + G,
               A = A / sum, T = T / sum, C = C / sum, G = G / sum)

    data_count <- count(data, "Mod")
    
    # Compute mean mismatch frequencies per modification
    data_avg <- aggregate(data[, c("A", "T", "C", "G")], by = list(data$Mod, data$ref_nuc), FUN = mean, na.rm = TRUE)
    
    # Melt data for visualization
    data_melted <- melt(data_avg, id.vars = c("Group.1", "Group.2"))
    
    # Merge with count data
    final_data <- merge(data_melted, data_count, by.x = "Group.1", by.y = "Mod")
    
    return(final_data)
}

# Compute mismatch profiles
dcDNA.mis <- clean_mis_input_agg(dcDNA.processed)
firstseq.mis <- clean_mis_input_agg(firstseq.processed)

# Function: Generate mismatch profile barplots
mismatch_profile_barplot_facet <- function(data, label) {
    pdf(file = paste0("figures/", label, "_mismatch_profile_barplot.pdf"), height = 1.5, width = 6)
    
    print(
        ggplot(data, aes(x = Group.1, y = value, fill = variable)) +
            scale_fill_manual(values = c("#1fab89", "#eb4d55", "#1e56a0", "#f0cf85", "#888888")) +
            geom_bar(stat = 'identity', colour = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    
    dev.off()
}

# Generate barplots for dcDNA-seq and FIRST-seq
mismatch_profile_barplot_facet(dcDNA.mis, "dcDNA-seq")
mismatch_profile_barplot_facet(firstseq.mis, "FIRST-seq")

# Function: Generate RT drop and coverage plots
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

# Process coverage data
runs_forcoverage <- rbind(dcDNA.processed, firstseq.processed)
palette_coverage <- c(dcDNA_col, firstseq_col)

# Generate RT drop and coverage plots
plot_area_coverage_facet(runs_forcoverage, "Firstseq_dcDNAseq", palette_coverage)

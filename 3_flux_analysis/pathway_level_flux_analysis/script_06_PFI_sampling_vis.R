# =================== R Visualization with Reaction Count Threshold ===================
# Purpose: Visualize Z-test results filtered by minimum reaction count threshold
# ====================================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(RColorBrewer)

# =================== CONFIGURATION ===================
base_dir <- "E:/Projects/revision/Code_2_upload/pathway_level_flux_analysis"


load_file <- 'sampling'  # 
threshold_rxn <- 2  # Minimum number of reactions required

# =================== SETUP DIRECTORIES ===================
output_dir <- file.path(base_dir, 
                        paste0("6_sampling_PFI_vis_", threshold_rxn),
                        load_file)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

matlab_results_dir <- file.path(base_dir, paste0("5_PFI_zTest_", load_file))

# =================== LOAD DATA ===================
if (load_file == "fba"){
  file_name <- "Permutation_test_FBA_results_Bonferroni_sorted_cleaned.xlsx"
} else if (load_file == "pfba"){
  file_name <- "Permutation_test_pFBA_results_Bonferroni_sorted_cleaned.xlsx"
} else {
  file_name <- paste0("Z_test_", load_file, "_results_Bonferroni_sorted_cleaned.xlsx")
}

z_results_all <- read_excel(file.path(matlab_results_dir, file_name))

z_results_all <- z_results_all %>%
  filter(Pathway != "Artificial reactions")

z_results_all <- z_results_all %>%
  filter(Pathway != "Exchange/demand reactions")

# =================== APPLY REACTION COUNT THRESHOLD ===================
cat("\n=================== THRESHOLD FILTERING ===================\n")
cat("Reaction count threshold:", threshold_rxn, "\n\n")

# Before filtering statistics
total_pathways <- nrow(z_results_all)
below_threshold <- sum(z_results_all$Reaction_count < threshold_rxn)
above_threshold <- sum(z_results_all$Reaction_count >= threshold_rxn)

cat("BEFORE FILTERING:\n")
cat("  Total pathways:", total_pathways, "\n")
cat("  Pathways below threshold (<", threshold_rxn, "rxns):", below_threshold, 
    sprintf("(%.1f%%)\n", 100*below_threshold/total_pathways))
cat("  Pathways at/above threshold (≥", threshold_rxn, "rxns):", above_threshold,
    sprintf("(%.1f%%)\n", 100*above_threshold/total_pathways))

# Significance in below-threshold pathways
below_thresh_data <- z_results_all %>% filter(Reaction_count < threshold_rxn)
below_thresh_sig <- sum(below_thresh_data$P_adj < 0.05, na.rm = TRUE)

cat("\n  In BELOW-threshold pathways:\n")
cat("    Significant (p_adj < 0.05):", below_thresh_sig,
    sprintf("out of %d (%.1f%%)\n", nrow(below_thresh_data), 
            100*below_thresh_sig/nrow(below_thresh_data)))

# Apply threshold filter
z_results <- z_results_all %>%
  filter(Reaction_count >= threshold_rxn)

cat("\nAFTER FILTERING (≥", threshold_rxn, "reactions):\n")
cat("  Remaining pathways:", nrow(z_results), "\n")
cat("  Significant pathways (p_adj < 0.05):", sum(z_results$P_adj < 0.05, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100*sum(z_results$P_adj < 0.05)/nrow(z_results)))
cat("  Mean |Log2FC|:", sprintf("%.3f\n", mean(abs(z_results$Log2FoldChange_flux), na.rm = TRUE)))
cat("  Mean reaction count:", sprintf("%.1f\n", mean(z_results$Reaction_count, na.rm = TRUE)))
cat("==========================================================\n\n")

# =================== EXPORT FILTERED RESULTS ===================
# Export pathways that passed the threshold
write_xlsx(z_results, 
           file.path(output_dir, paste0("Z_test_results_threshold", threshold_rxn, "_filtered.xlsx")))

# Export pathways that were excluded (below threshold)
below_threshold_data <- z_results_all %>%
  filter(Reaction_count < threshold_rxn) %>%
  arrange(desc(abs(Log2FoldChange_flux)))

write_xlsx(below_threshold_data,
           file.path(output_dir, paste0("Z_test_results_below_threshold", threshold_rxn, ".xlsx")))

# =================== DOTPLOT: SIGNIFICANT PATHWAYS ===================
sig_pathways <- z_results %>%
  filter(P_adj < 0.05) %>%
  arrange(desc(Log2FoldChange_flux)) %>%
  mutate(
    Direction = ifelse(Log2FoldChange_flux > 0, "Up in HSD", "Down in HSD"),
    Direction = factor(Direction, levels = c("Up in HSD", "Down in HSD")),
    Pathway_short = ifelse(nchar(Pathway) > 60, 
                           paste0(substr(Pathway, 1, 60), "..."), 
                           Pathway)
  )

if (nrow(sig_pathways) > 0) {
  p1 <- ggplot(sig_pathways, aes(x = Log2FoldChange_flux, 
                                 y = reorder(Pathway_short, Log2FoldChange_flux))) +
    geom_point(aes(size = log10(Reaction_count + 1), color = Direction), alpha = 0.7) +
    scale_color_manual(values = c("Up in HSD" = "#e74c3c", "Down in HSD" = "#3498db")) +
    scale_size_continuous(name = "log10(Reaction\nCount + 1)", range = c(2, 6)) +
    labs(
      title = paste0("Significant Pathways (≥", threshold_rxn, " reactions): Log2FC (HSD vs NSD)"),
      subtitle = paste0(nrow(sig_pathways), " significant pathways with ≥", threshold_rxn, " reactions"),
      x = "Log2 Fold Change (HSD vs NSD)",
      y = "Metabolic Pathways"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
  
  ggsave(file.path(output_dir, "dotplot_significant_log2FC.png"), 
         plot = p1, width = 12, height = max(8, nrow(sig_pathways)*0.3), dpi = 300)
  ggsave(file.path(output_dir, "dotplot_significant_log2FC.svg"), 
         plot = p1, width = 12, height = max(8, nrow(sig_pathways)*0.3))
  
  cat("Created dotplot for significant pathways\n")
}


# =================== DOTPLOT: ALL PATHWAYS WITH -log10(P_adj) COLOR ===================
all_pathways <- z_results %>%
  # REMOVED: filter(P_adj < 0.05) %>%
  arrange(desc(Log2FoldChange_flux)) %>%
  mutate(
    Pathway_short = ifelse(nchar(Pathway) > 60, 
                           paste0(substr(Pathway, 1, 60), "..."), 
                           Pathway),
    neg_log10_padj = -log10(P_adj),
    # Cap extremely small p-values for visualization
    neg_log10_padj_capped = pmin(neg_log10_padj, 20)
  )

if (nrow(all_pathways) > 0) {
  
  # Check p-value distribution
  cat("\n-log10(p_adj) statistics:\n")
  cat("  Range:", round(min(all_pathways$neg_log10_padj), 2), "to", 
      round(max(all_pathways$neg_log10_padj), 2), "\n")
  cat("  Mean:", round(mean(all_pathways$neg_log10_padj), 2), "\n")
  cat("  Median:", round(median(all_pathways$neg_log10_padj), 2), "\n\n")
  
  cat("Reaction count statistics:\n")
  cat("  Range:", min(all_pathways$Reaction_count), "to", 
      max(all_pathways$Reaction_count), "\n")
  cat("  Mean:", round(mean(all_pathways$Reaction_count), 1), "\n")
  cat("  Median:", median(all_pathways$Reaction_count), "\n\n")
  
  cat("Significance breakdown:\n")
  cat("  Significant (p_adj < 0.05):", sum(all_pathways$P_adj < 0.05), "\n")
  cat("  Non-significant:", sum(all_pathways$P_adj >= 0.05), "\n\n")
  
  # -log10(P_adj) as color with log10 reaction count
  p_all <- ggplot(all_pathways, aes(x = Log2FoldChange_flux, 
                                    y = reorder(Pathway_short, Log2FoldChange_flux))) +
    geom_point(aes(size = log10(Reaction_count + 1), 
                   color = neg_log10_padj_capped), 
               alpha = 0.8) +
    scale_color_viridis_c(
      name = "-log10(p_adj)",
      option = "plasma",
      direction = -1,
      limits = c(min(all_pathways$neg_log10_padj_capped), 
                 max(all_pathways$neg_log10_padj_capped)),
      breaks = pretty(c(min(all_pathways$neg_log10_padj_capped), 
                        max(all_pathways$neg_log10_padj_capped)), n = 5)
    ) +
    scale_size_continuous(
      name = "log10(Reaction\nCount + 1)", 
      range = c(3, 10),
      breaks = pretty(log10(all_pathways$Reaction_count + 1), n = 4)
    ) +
    labs(
      title = paste0("All Pathways (≥", threshold_rxn, " reactions): Colored by Significance"),
      subtitle = paste0(sum(all_pathways$P_adj < 0.05), " significant / ", 
                        nrow(all_pathways), " total pathways"),
      x = "Log2 Fold Change (HSD vs NSD)",
      y = "Metabolic Pathways"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5)
  
  ggsave(file.path(output_dir, "dotplot_all_pathways.png"), 
         plot = p_all, width = 12, height = max(8, nrow(all_pathways)*0.3), dpi = 300)
  ggsave(file.path(output_dir, "dotplot_all_pathways.svg"), 
         plot = p_all, width = 12, height = max(8, nrow(all_pathways)*0.3))
  
  cat("✅ Created dotplot with ALL pathways:\n")
  cat("  - Color: -log10(p_adj) - lighter colors = less significant\n")
  cat("  - Size: log10(Reaction Count + 1)\n\n")
}
print(p_all)


# =================== Highlight commonly affected pathways ===================
# =================== LOAD THREE-WAY COMPARISON DATA ===================
upset_file <- file.path("E:/Projects/revision/Code_2_upload/pathway_level_flux_analysis/4_compare_all", "UpSet_pathway_data.xlsx")

axis_colors <- NULL  # Initialize

if (file.exists(upset_file)) {
  # Load up-regulated pathways
  up_data <- read_excel(upset_file, sheet = "Up_Regulated")
  
  # Filter: pathways where both pFBA = 1 AND FVA_sampling = 1
  up_all_three <- up_data %>%
    filter(pFBA == 1 & FVA_sampling == 1) %>%
    pull(Pathway)
  
  # Load down-regulated pathways
  down_data <- read_excel(upset_file, sheet = "Down_Regulated")
  
  # Filter: pathways where both pFBA = 1 AND FVA_sampling = 1
  down_all_three <- down_data %>%
    filter(pFBA == 1 & FVA_sampling == 1) %>%
    pull(Pathway)
  
  # Further filter by statistical significance (p_adj < 0.05) from z_results
  significant_pathways <- z_results %>%
    filter(P_adj < 0.05) %>%
    pull(Pathway)

  up_all_three <- up_all_three[up_all_three %in% significant_pathways]
  down_all_three <- down_all_three[down_all_three %in% significant_pathways]
  
  cat("\nPathways common to pFBA and FVA_sampling:\n")
  cat("  Down-regulated:", length(down_all_three), "\n")
  cat("  Up-regulated:", length(up_all_three), "\n\n")
} else {
  cat("Warning: UpSet data file not found\n")
  down_all_three <- character(0)
  up_all_three <- character(0)
}



# =================== DOTPLOT: ALL PATHWAYS WITH COLORED Y-AXIS LABELS ===================
all_pathways <- z_results %>%
  arrange(desc(Log2FoldChange_flux)) %>%
  mutate(
    Pathway_short = ifelse(nchar(Pathway) > 60, 
                           paste0(substr(Pathway, 1, 60), "..."), 
                           Pathway),
    neg_log10_padj = -log10(P_adj),
    neg_log10_padj_capped = pmin(neg_log10_padj, 20),
    abs_zscore = abs(Z_score)
  )

if (nrow(all_pathways) > 0) {
  
  cat("\n-log10(p_adj) statistics:\n")
  cat("  Range:", round(min(all_pathways$neg_log10_padj), 2), "to", 
      round(max(all_pathways$neg_log10_padj), 2), "\n")
  
  cat("\n|Z-score| statistics:\n")
  cat("  Range:", round(min(all_pathways$abs_zscore), 2), "to", 
      round(max(all_pathways$abs_zscore), 2), "\n")
  cat("  Mean:", round(mean(all_pathways$abs_zscore), 2), "\n\n")
  
  cat("Significance breakdown:\n")
  cat("  Significant (p_adj < 0.05):", sum(all_pathways$P_adj < 0.05), "\n")
  cat("  Non-significant:", sum(all_pathways$P_adj >= 0.05), "\n\n")
  
  # Create color vector for y-axis labels based on pathway order in plot
  pathway_order <- all_pathways %>%
    arrange(Log2FoldChange_flux) %>%
    pull(Pathway)
  
  axis_colors <- sapply(pathway_order, function(p) {
    if (p %in% down_all_three) {
      "#3498db"  # Blue for down-regulated in all 3 methods
    } else if (p %in% up_all_three) {
      "#e74c3c"  # Red for up-regulated in all 3 methods
    } else {
      "black"    # Black for sampling only
    }
  })
  
  # Count colored pathways
  n_blue <- sum(axis_colors == "#3498db")
  n_red <- sum(axis_colors == "#e74c3c")
  cat("Y-axis label colors:\n")
  cat("  Blue (down in all 3):", n_blue, "\n")
  cat("  Red (up in all 3):", n_red, "\n")
  cat("  Black (sampling only):", length(axis_colors) - n_blue - n_red, "\n\n")
  
  # Create dotplot
  p_all <- ggplot(all_pathways, aes(x = Log2FoldChange_flux, 
                                    y = reorder(Pathway_short, Log2FoldChange_flux))) +
    geom_point(aes(size = log10(abs_zscore+1), 
                   color = neg_log10_padj_capped), 
               alpha = 0.8) +
    scale_color_viridis_c(
      name = "-log10(p_adj)",
      option = "plasma",
      direction = -1,
      limits = c(min(all_pathways$neg_log10_padj_capped), 
                 max(all_pathways$neg_log10_padj_capped))
    ) +
    scale_size_continuous(
      name = "log10(|Z-score|+1)",
      range = c(1, 8)
    ) +
    labs(
      title = "Pathway-level flux differences",
      x = expression(log[2](bar(X)[HSD] / bar(X)[NSD])),
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, color = axis_colors),  # Colored labels
      axis.text.x = element_text(size = 16),
      axis.title = element_text(size = 25),
      plot.title = element_text(size = 25),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20),
      legend.position = "right", # right
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5)
  
  ggsave(file.path(output_dir, "dotplot_all_pathways_zscore.png"), 
         plot = p_all, width = 20, height = max(8, nrow(all_pathways)*0.2), dpi = 300)
  ggsave(file.path(output_dir, "dotplot_all_pathways_zscore.svg"), 
         plot = p_all, width = 20, height = max(8, nrow(all_pathways)*0.2))
  ggsave(file.path(output_dir, "dotplot_all_pathways_zscore.pdf"), 
         plot = p_all, width = 20, height = max(8, nrow(all_pathways)*0.2))
  
  cat("✅ Created dotplot with ALL pathways:\n")
  cat("  - Color: -log10(p_adj) - significance level\n")
  cat("  - Size: log10(|Z-score|+1) - statistical strength\n")
  cat("  - Y-axis: Blue = down in all methods, Red = up in all methods\n\n")
}

print(p_all)

# =================== DOTPLOT: ALL PATHWAYS WITH FLUX DIFFERENCE AS SIZE ===================
all_pathways <- z_results %>%
  arrange(desc(Log2FoldChange_flux)) %>%
  mutate(
    Pathway_short = ifelse(nchar(Pathway) > 60, 
                           paste0(substr(Pathway, 1, 60), "..."), 
                           Pathway),
    neg_log10_padj = -log10(P_adj),
    neg_log10_padj = ifelse(is.infinite(neg_log10_padj), 20, neg_log10_padj),  # Fix Inf
    neg_log10_padj_capped = pmin(neg_log10_padj, 20),
    abs_flux_diff = abs(HSD_flux - NSD_flux)  # Absolute flux difference
  )

if (nrow(all_pathways) > 0) {
  
  cat("\n-log10(p_adj) statistics:\n")
  cat("  Range:", round(min(all_pathways$neg_log10_padj), 2), "to", 
      round(max(all_pathways$neg_log10_padj), 2), "\n")
  
  cat("\nAbsolute flux difference statistics:\n")
  cat("  Range:", round(min(all_pathways$abs_flux_diff), 2), "to", 
      round(max(all_pathways$abs_flux_diff), 2), "\n")
  cat("  Mean:", round(mean(all_pathways$abs_flux_diff), 2), "\n")
  cat("  Median:", round(median(all_pathways$abs_flux_diff), 2), "\n\n")
  
  cat("Significance breakdown:\n")
  cat("  Significant (p_adj < 0.05):", sum(all_pathways$P_adj < 0.05), "\n")
  cat("  Non-significant:", sum(all_pathways$P_adj >= 0.05), "\n\n")
  
  # Create color vector for y-axis labels based on pathway order in plot
  pathway_order <- all_pathways %>%
    arrange(desc(Log2FoldChange_flux)) %>%  # Changed to desc() to match the flipped axis
    pull(Pathway)
  
  axis_colors <- sapply(pathway_order, function(p) {
    if (p %in% down_all_three) {
      "#3498db"  # Blue for down-regulated in all 3 methods
    } else if (p %in% up_all_three) {
      "#e74c3c"  # Red for up-regulated in all 3 methods
    } else {
      "black"    # Black for sampling only
    }
  })
  
  # Count colored pathways
  n_blue <- sum(axis_colors == "#3498db")
  n_red <- sum(axis_colors == "#e74c3c")
  cat("Y-axis label colors:\n")
  cat("  Blue (down in all 3):", n_blue, "\n")
  cat("  Red (up in all 3):", n_red, "\n")
  cat("  Black (sampling only):", length(axis_colors) - n_blue - n_red, "\n\n")
  
  # Create dotplot
  p_all <- ggplot(all_pathways, aes(x = Log2FoldChange_flux, 
                                    y = reorder(Pathway_short, -Log2FoldChange_flux))) +
    geom_point(aes(size = log10(abs_flux_diff + 1),  # Log10 of flux difference
                   color = neg_log10_padj_capped), 
               alpha = 0.8) +
    scale_color_viridis_c(
      name = "-log10(p_adj)",
      option = "plasma",
      direction = -1,
      limits = c(min(all_pathways$neg_log10_padj_capped), 
                 max(all_pathways$neg_log10_padj_capped))
    ) +
    scale_size_continuous(
      name = "log10(|Flux diff.| + 1)",
      range = c(1, 8)
    ) +
    labs(
      title = "Pathway-level flux differences",
      x = expression(log[2](bar(X)[HSD] / bar(X)[NSD])),
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, color = axis_colors),
      axis.text.x = element_text(size = 16),
      axis.title = element_text(size = 25),
      plot.title = element_text(size = 25),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5)
  
  ggsave(file.path(output_dir, "dotplot_all_pathways_flux_diff.png"), 
         plot = p_all, width = 20, height = max(16, nrow(all_pathways)*0.2), dpi = 300)
  ggsave(file.path(output_dir, "dotplot_all_pathways_flux_diff.svg"), 
         plot = p_all, width = 20, height = max(16, nrow(all_pathways)*0.2))
  ggsave(file.path(output_dir, "dotplot_all_pathways_flux_diff.pdf"), 
         plot = p_all, width = 20, height = max(16, nrow(all_pathways)*0.2))
  
  cat("✅ Created dotplot with ALL pathways:\n")
  cat("  - Color: -log10(p_adj) - significance level\n")
  cat("  - Size: log10(|HSD_flux - NSD_flux| + 1) - magnitude of flux change\n")
  cat("  - Y-axis: Blue = down in all methods, Red = up in all methods\n\n")
}

print(p_all)


# =================== EXPORT VISUALIZATION DATA ===================
cat("\n=================== EXPORTING VISUALIZATION DATA ===================\n")

# Prepare the complete visualization dataset
viz_data <- all_pathways %>%
  arrange(desc(Log2FoldChange_flux)) %>%  # Same order as y-axis (flipped)
  mutate(
    # Add the axis color information
    Y_axis_color = sapply(Pathway, function(p) {
      if (p %in% down_all_three) {
        "Blue (Down in pFBA & Sampling)"
      } else if (p %in% up_all_three) {
        "Red (Up in pFBA & Sampling)"
      } else {
        "Black (Sampling only)"
      }
    }),
    # Add plot position for reference
    Plot_order = row_number(),
    # Size value used in plot
    Size_value = log10(abs_flux_diff + 1),
    # Color value used in plot
    Color_value = neg_log10_padj_capped
  )

# Reorder columns manually
viz_data <- viz_data[, c(
  "Plot_order",
  "Pathway",
  "Pathway_short",
  "Y_axis_color",
  "Log2FoldChange_flux",
  "P_value",
  "P_adj",
  "neg_log10_padj",
  "neg_log10_padj_capped",
  "Color_value",
  "Z_score",
  "HSD_flux",
  "NSD_flux",
  "abs_flux_diff",
  "Size_value",
  "Reaction_count",
  names(viz_data)[!names(viz_data) %in% c(
    "Plot_order", "Pathway", "Pathway_short", "Y_axis_color",
    "Log2FoldChange_flux", "P_value", "P_adj", "neg_log10_padj",
    "neg_log10_padj_capped", "Color_value", "Z_score", "HSD_flux",
    "NSD_flux", "abs_flux_diff", "Size_value", "Reaction_count"
  )]
)]

# Create summary sheet
viz_summary <- data.frame(
  Metric = c(
    "Total pathways plotted",
    "Significant pathways (p_adj < 0.05)",
    "Non-significant pathways",
    "Blue labels (Down in pFBA & Sampling)",
    "Red labels (Up in pFBA & Sampling)",
    "Black labels (Sampling only)",
    "Min Log2FC",
    "Max Log2FC",
    "Min -log10(p_adj)",
    "Max -log10(p_adj)",
    "Min |Flux difference|",
    "Max |Flux difference|",
    "Plot width (pixels)",
    "Plot height (pixels)"
  ),
  Value = c(
    nrow(viz_data),
    sum(viz_data$P_adj < 0.05),
    sum(viz_data$P_adj >= 0.05),
    sum(viz_data$Y_axis_color == "Blue (Down in pFBA & Sampling)"),
    sum(viz_data$Y_axis_color == "Red (Up in pFBA & Sampling)"),
    sum(viz_data$Y_axis_color == "Black (Sampling only)"),
    round(min(viz_data$Log2FoldChange_flux), 4),
    round(max(viz_data$Log2FoldChange_flux), 4),
    round(min(viz_data$neg_log10_padj), 4),
    round(max(viz_data$neg_log10_padj), 4),
    round(min(viz_data$abs_flux_diff), 4),
    round(max(viz_data$abs_flux_diff), 4),
    20 * 72,  # Width in pixels at 72 dpi
    max(16, nrow(viz_data)*0.2) * 72  # Height in pixels at 72 dpi
  )
)

# Create legend explanation sheet
legend_info <- data.frame(
  Plot_Element = c(
    "X-axis",
    "Y-axis",
    "Point color",
    "Point size",
    "Y-axis label color (Blue)",
    "Y-axis label color (Red)",
    "Y-axis label color (Black)",
    "Vertical dashed line"
  ),
  Represents = c(
    "Log2 Fold Change (HSD/NSD)",
    "Metabolic pathways (ordered by Log2FC, flipped)",
    "-log10(adjusted p-value), capped at 20",
    "log10(|HSD_flux - NSD_flux| + 1)",
    "Down-regulated in both pFBA and Sampling",
    "Up-regulated in both pFBA and Sampling",
    "Significant in Sampling only",
    "Zero fold change reference"
  ),
  Details = c(
    "Positive values = higher in HSD; Negative = higher in NSD",
    "Pathways with highest Log2FC at top after flipping",
    "Higher values (lighter colors) = more significant",
    "Larger points = bigger flux difference between conditions",
    paste(sum(viz_data$Y_axis_color == "Blue (Down in pFBA & Sampling)"), "pathways"),
    paste(sum(viz_data$Y_axis_color == "Red (Up in pFBA & Sampling)"), "pathways"),
    paste(sum(viz_data$Y_axis_color == "Black (Sampling only)"), "pathways"),
    "Visual reference for no change"
  )
)

# Export to Excel with multiple sheets
viz_export_list <- list(
  "Visualization_Data" = viz_data,
  "Summary_Statistics" = viz_summary,
  "Legend_Explanation" = legend_info
)

# Add a sheet with just the colored pathways if they exist
if (length(down_all_three) > 0 || length(up_all_three) > 0) {
  colored_pathways <- viz_data %>%
    filter(Y_axis_color != "Black (Sampling only)")
  
  # Select specific columns using bracket notation
  colored_pathways <- colored_pathways[, c("Plot_order", "Pathway", "Y_axis_color", 
                                           "Log2FoldChange_flux", "P_adj", "neg_log10_padj", 
                                           "abs_flux_diff", "Reaction_count")]
  
  viz_export_list[["Colored_Pathways"]] = colored_pathways
}

# =================== EXPORT HIGHLIGHTED PATHWAYS ===================

# Create summary of highlighted pathways
highlighted_pathways <- data.frame(
  Pathway = c(down_all_three, up_all_three),
  Direction = c(rep("Down-regulated (pFBA & FVA, p<0.05)", length(down_all_three)),
                rep("Up-regulated (pFBA & FVA, p<0.05)", length(up_all_three))),
  stringsAsFactors = FALSE
)

# Add detailed statistics for these pathways
if (nrow(highlighted_pathways) > 0) {
  highlighted_detail <- z_results %>%
    filter(Pathway %in% highlighted_pathways$Pathway) %>%
    select(Pathway, Z_score, P_value, P_adj, Log2FoldChange_flux, 
           NSD_flux, HSD_flux, Reaction_count) %>%
    left_join(highlighted_pathways, by = "Pathway") %>%
    arrange(Direction, desc(abs(Log2FoldChange_flux)))
  
  # Export to Excel
  library(writexl)
  write_xlsx(highlighted_detail, 
             file.path(output_dir, "highlighted_pathways_pFBA_FVA_overlap.xlsx"))
  
  cat("\n✅ Exported highlighted pathways:\n")
  cat("  File: highlighted_pathways_pFBA_FVA_overlap.xlsx\n")
  cat("  Down-regulated (blue labels):", length(down_all_three), "pathways\n")
  cat("  Up-regulated (red labels):", length(up_all_three), "pathways\n\n")
  
  # Print pathway names to console
  if (length(down_all_three) > 0) {
    cat("DOWN-REGULATED pathways (blue labels):\n")
    for (p in sort(down_all_three)) {
      cat("  -", p, "\n")
    }
    cat("\n")
  }
  
  if (length(up_all_three) > 0) {
    cat("UP-REGULATED pathways (red labels):\n")
    for (p in sort(up_all_three)) {
      cat("  -", p, "\n")
    }
    cat("\n")
  }
} else {
  cat("No pathways met both pFBA/FVA overlap and p_adj < 0.05 criteria\n")
}






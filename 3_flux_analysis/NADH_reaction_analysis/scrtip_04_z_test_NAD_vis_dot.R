# scrtip_04_z_test_NAD_vis_dot
#  NAD(H) REACTION DOTPLOT WITH Z-TEST RESULTS
# Purpose: Create dotplot visualization for NAD(H)-dependent reactions
# showing flux changes, statistical significance, and highlighting common reactions

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(forcats)
library(writexl)
library(ggtext)  # For colored y-axis labels

rm(list = ls())

# =================== SETUP ===================
# Path to Z-test results from MATLAB
current_folder = "E:/Projects/revision/Code_2_upload/NAD_reaction_analysis"
matlab_output <- file.path(current_folder,"02_z_test_NAD")
z_test_file <- file.path(matlab_output, "Z_test_NAD_reactions_for_R_visualization.xlsx")

# Path to common reactions (for highlighting)
reaction_folder <- "E:/Projects/revision/Code_2_upload/NAD_reaction_analysis/03_nad_analysis_all"
dec_rxn_file <- file.path(reaction_folder, "NAD_All_Common_Decreased_Reactions.csv")
up_rxn_file <- file.path(reaction_folder, "NAD_All_Common_Increased_Reactions.csv")

# Output directory
output_dir <- file.path(current_folder, "04_z_test_NAD_vis")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =================== LOAD DATA ===================
cat("Loading Z-test results...\n")
z_results <- read_excel(z_test_file)

cat("Loading common reaction lists...\n")
dec_rxns <- read.csv(dec_rxn_file)$Reaction_ID
up_rxns <- read.csv(up_rxn_file)$Reaction_ID

cat("Total NAD(H) reactions:", nrow(z_results), "\n")
cat("Common decreased reactions:", length(dec_rxns), "\n")
cat("Common increased reactions:", length(up_rxns), "\n\n")

# =================== PROCESS DATA ===================
df_plot <- z_results %>%
  mutate(
    # Calculate metrics
    neg_log10_padj = -log10(P_adj),
    neg_log10_padj = ifelse(is.infinite(neg_log10_padj), 20, neg_log10_padj),
    neg_log10_padj_capped = pmin(neg_log10_padj, 20),
    abs_flux_diff = abs(HSD_flux_mean - NSD_flux_mean),
    
    # Extract first two genes from Gene_Rules
    gene_parts = str_split(Gene_Rules, "\\s*;\\s*|\\s+or\\s+|\\s+and\\s+"),
    first_two_genes = sapply(gene_parts, function(x) {
      parts <- x[1:2]
      parts <- parts[!is.na(parts) & parts != ""]
      if (length(parts) == 0) return("No gene")
      paste(parts, collapse = " ; ")
    }),
    
    # Extract compartment codes from Associated_Metabolites
    compartment_codes = ifelse(is.na(Compartments) | Compartments == "", 
                               "", 
                               paste0("[", Compartments, "]")),
    
    
    # Create y-axis label: Gene (Reaction_ID) [compartment]
    y_axis_label = paste0(first_two_genes, " (", Reaction_ID, ") ", compartment_codes),
    
    # Clean up encoding issues
    y_axis_label = str_replace_all(y_axis_label, "Î±", "α"),
    y_axis_label = str_replace_all(y_axis_label, "Î²", "β"),
    y_axis_label = str_replace_all(y_axis_label, "\\bTSG101\\b", "Ldh"),
    
    # Highlight color for common reactions
    highlight_color = case_when(
      Reaction_ID %in% dec_rxns ~ "#3498db",  # Blue for decreased
      Reaction_ID %in% up_rxns ~ "#e74c3c",   # Red for increased
      TRUE ~ "black"                           # Black for others
    ),
    
    # Reorder by log2FC
    y_axis_label = fct_reorder(y_axis_label, Log2FoldChange, .desc = TRUE)
  )

# Save processed data
write_xlsx(df_plot, file.path(output_dir, "NAD_reactions_dotplot_ready.xlsx"))

# =================== SUMMARY STATISTICS ===================
cat("\n=== Data Summary ===\n")
cat("-log10(p_adj) range:", round(min(df_plot$neg_log10_padj_capped), 2), 
    "to", round(max(df_plot$neg_log10_padj_capped), 2), "\n")
cat("Absolute flux difference:\n")
cat("  Range:", round(min(df_plot$abs_flux_diff), 2), "to", 
    round(max(df_plot$abs_flux_diff), 2), "\n")
cat("  Mean:", round(mean(df_plot$abs_flux_diff), 2), "\n")
cat("  Median:", round(median(df_plot$abs_flux_diff), 2), "\n\n")

cat("Significance breakdown:\n")
cat("  Significant (p_adj < 0.05):", sum(df_plot$P_adj < 0.05, na.rm = TRUE), "\n")
cat("  p_adj < 0.01:", sum(df_plot$P_adj < 0.01, na.rm = TRUE), "\n")
cat("  p_adj < 0.001:", sum(df_plot$P_adj < 0.001, na.rm = TRUE), "\n")
cat("  Non-significant:", sum(df_plot$P_adj >= 0.05, na.rm = TRUE), "\n\n")

cat("Y-axis label colors:\n")
cat("  Blue (decreased in all 3 methods):", sum(df_plot$highlight_color == "#3498db"), "\n")
cat("  Red (increased in all 3 methods):", sum(df_plot$highlight_color == "#e74c3c"), "\n")
cat("  Black (not common to all methods):", sum(df_plot$highlight_color == "black"), "\n\n")

# =================== CREATE DOTPLOT ===================
cat("Creating dotplot...\n")

# Get ordered colors for y-axis
axis_colors <- df_plot %>%
  arrange(Log2FoldChange) %>%
  pull(highlight_color)

p_all <- ggplot(df_plot, aes(x = Log2FoldChange, y = y_axis_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(aes(size = log10(abs_flux_diff + 1),  # Log10 of flux difference
                 color = neg_log10_padj_capped), 
             alpha = 0.8) +
  scale_color_viridis_c(
    name = "-log10(p_adj)",
    option = "plasma",
    direction = -1,
    limits = c(min(df_plot$neg_log10_padj_capped), 
               max(df_plot$neg_log10_padj_capped))
  ) +
  scale_size_continuous(
    name = "log10(|Flux diff.| + 1)",
    range = c(1, 8)
  ) +
  labs(
    title = "NAD(H)-dependent reaction flux changes",
    subtitle = "HSD vs NSD (FVA-sampling)",
    x = expression(log[2](bar(X)[HSD] / bar(X)[NSD])),
    y = ""
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7, color = axis_colors, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p_all)

# =================== SAVE PLOTS ===================
cat("\nSaving plots...\n")

plot_height <- max(12, nrow(df_plot) * 0.15)

ggsave(file.path(output_dir, "dotplot_NAD_reactions_all.png"), 
       plot = p_all, width = 14, height = plot_height, dpi = 300)
ggsave(file.path(output_dir, "dotplot_NAD_reactions_all.svg"), 
       plot = p_all, width = 14, height = plot_height)
ggsave(file.path(output_dir, "dotplot_NAD_reactions_all.pdf"), 
       plot = p_all, width = 14, height = plot_height)

cat("✅ Plots saved to:", output_dir, "\n")
cat("  - Color scale: -log10(p_adj) (statistical significance)\n")
cat("  - Dot size: log10(|flux difference| + 1)\n")
cat("  - Y-axis colors: Blue = decreased (all methods), Red = increased (all methods)\n\n")






# =================== FILTERED PLOT (SIGNIFICANT ONLY) ===================
cat("Creating filtered plot (significant reactions only)...\n")

df_plot_sig <- df_plot %>%
  filter(P_adj < 0.05)


## Make log2FC filter
log2fc_threshold <- 0.5  # Adjust this value

df_plot_sig <- df_plot_sig %>%
  filter(abs(Log2FoldChange) > log2fc_threshold)
cat("After log2FC filter (|log2FC| >", log2fc_threshold, "):", nrow(df_plot_sig), "reactions\n")


if (nrow(df_plot_sig) > 0) {
  cat("Significant reactions:", nrow(df_plot_sig), "\n")
  
  # Get ordered colors for filtered plot (REVERSED ORDER - decreasing on top)
  axis_colors_sig <- df_plot_sig %>%
    arrange(desc(Log2FoldChange)) %>%  # Changed to desc for reversed order
    pull(highlight_color)
  
  p_sig <- ggplot(df_plot_sig, aes(x = Log2FoldChange, 
                                   y = fct_reorder(y_axis_label, Log2FoldChange, .desc = TRUE))) +  # Added .desc = TRUE
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(aes(size = log10(abs_flux_diff + 1),
                   color = neg_log10_padj_capped), 
               alpha = 0.8) +
    scale_color_viridis_c(
      name = "-log10(p_adj)",
      option = "plasma",
      direction = -1
    ) +
    scale_size_continuous(
      name = "log10(|Flux diff.| + 1)",
      range = c(1, 6)
    #  range = c(2, 12)
    ) +
    labs(
    #  title = "Significant NAD(H)-dependent reaction changes",
      title = "Differential fluxes of NAD(H)-dependent reactions",
    #  subtitle = "HSD vs NSD (p_adj < 0.05, FVA-sampling)",
      subtitle = "FVA-sampling (Filtered)",
      x = expression(log[2](bar(X)[HSD] / bar(X)[NSD])),
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, color = axis_colors_sig, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 18),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  print(p_sig)
  
  # ... rest of saving code stays the same
  
  plot_height_sig <- max(10, nrow(df_plot_sig) * 0.2)
  plot_height_sig <- 20
  
  # Create filename suffix based on filter used
  filter_suffix <- paste0("_log2FC_gt_", log2fc_threshold)
  
  ggsave(file.path(output_dir, paste0("dotplot_NAD_reactions_significant", filter_suffix, ".png")), 
         plot = p_sig, width = 16, height = plot_height_sig, dpi = 300)
  ggsave(file.path(output_dir, paste0("dotplot_NAD_reactions_significant", filter_suffix, ".svg")), 
         plot = p_sig, width = 16, height = plot_height_sig)
  ggsave(file.path(output_dir, paste0("dotplot_NAD_reactions_significant", filter_suffix, ".pdf")), 
         plot = p_sig, width = 16, height = plot_height_sig)
  
  cat("✅ Filtered plot saved with threshold info:", filter_suffix, "\n\n")
  
  
  cat("✅ Filtered plot saved (significant reactions only)\n\n")
} else {
  cat("⚠️  No significant reactions to plot\n\n")
}

# =================== EXPORT SUMMARY TABLE ===================
# =================== EXPORT SUMMARY TABLES ===================

# Export 1: ALL reactions with statistics
summary_all <- df_plot %>%
  mutate(
    sig_category = case_when(
      P_adj < 0.001 ~ "p < 0.001",
      P_adj < 0.01 ~ "p < 0.01",
      P_adj < 0.05 ~ "p < 0.05",
      TRUE ~ "ns"
    ),
    direction = ifelse(Log2FoldChange > 0, "Increased", "Decreased"),
    common_status = case_when(
      Reaction_ID %in% dec_rxns ~ "Common decreased",
      Reaction_ID %in% up_rxns ~ "Common increased",
      TRUE ~ "Not common"
    )
  ) %>%
  select(Reaction_ID, Subsystem, Compartments, Gene_Rules, y_axis_label,
         Log2FoldChange, P_adj, Z_score, 
         NSD_flux_mean, HSD_flux_mean, abs_flux_diff,
         sig_category, direction, common_status,
         everything())  # Add any remaining columns

# Export 2: FILTERED reactions shown in the plot
summary_filtered <- df_plot_sig %>%
  mutate(
    sig_category = case_when(
      P_adj < 0.001 ~ "p < 0.001",
      P_adj < 0.01 ~ "p < 0.01",
      P_adj < 0.05 ~ "p < 0.05",
      TRUE ~ "ns"
    ),
    direction = ifelse(Log2FoldChange > 0, "Increased", "Decreased"),
    common_status = case_when(
      Reaction_ID %in% dec_rxns ~ "Common decreased",
      Reaction_ID %in% up_rxns ~ "Common increased",
      TRUE ~ "Not common"
    )
  ) %>%
  select(Reaction_ID, Subsystem, Compartments, Gene_Rules, y_axis_label,
         Log2FoldChange, P_adj, Z_score, 
         NSD_flux_mean, HSD_flux_mean, abs_flux_diff,
         sig_category, direction, common_status,
         everything())  # Add any remaining columns

# Export both detailed tables
write_xlsx(
  list(
    "All_Reactions" = summary_all,
    "Filtered_Plot_Data" = summary_filtered
  ), 
  file.path(output_dir, paste0("NAD_reactions_detailed", filter_suffix, ".xlsx"))
)

cat("✅ Detailed reaction data saved (all + filtered plot data)\n")
cat("   Total reactions analyzed:", nrow(summary_all), "\n")
cat("   Reactions shown in plot:", nrow(summary_filtered), "\n")
cat("   Filter applied:", filter_suffix, "\n")
cat("\n=== Analysis Complete ===\n")


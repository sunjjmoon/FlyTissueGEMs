# Load necessary libraries
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggrepel)
library(scales)
library(ggpubr)
library(writexl)

# Remove all objects
rm(list=ls())

# Set file paths
filepath <- 'C:/Users/sunjj/tissueGEM_revision/v2/code_to_upload/3_Flux_analysis_f_20250513/D_flux_analysis/log/nadC'
setwd(filepath)

output_folder <- file.path(filepath, "bubble")

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Load data
filename <- 'nadhC_summary_gene_filled'
df <- read_excel(paste0(filepath, '/', filename, '.xlsx'))

# Calculate log2 Fold Change (log2FC) and Mean Difference
df <- df %>%
  mutate(log2FC = log2(abs(relChg))) %>%
  mutate(mean_diff = abs(exp_mean - ctr_mean)) %>%  # Mean Difference for bubble size
  arrange(desc(log2FC))  # Sort by log2FC

# Reorder factor levels using Rxns (unique)
df$Rxns <- factor(df$Rxns, levels = df$Rxns[order(df$log2FC)])
df$Gene <- factor(df$Gene, levels = df$Gene[order(df$log2FC)])

# Save updated data to Excel
write_xlsx(df, paste0(output_folder, "/nadhC_summary_updated.xlsx"))

# Create dot plot with Mean Difference as bubble size & Color by `ctr_mean`
dot_plot <- ggplot(df, aes(x = log2FC, y = Gene, size = mean_diff, color = log10(abs(ctr_mean + 1)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.1) +  # Add x-line at 0
  geom_point() +  # Now using color gradient
  scale_color_gradient(low = "lightgray", high = "darkred") +  # Color gradient: low = light, high = dark
  scale_size_continuous(range = c(3, 15)) +  # Adjust bubble sizes
  labs(
    x = expression(paste("log"[2], " |Fold Change|")),
    y = NULL, 
    title = "NAD(H)-dependent reactions",
    size = "|flux difference|",
    color = "log10(|control flux + 1|)"  # Legend title for color scale
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 15, face = "bold"),  # Adjust y-axis text
    axis.text.x = element_text(size = 13),  # Adjust x-axis text
    axis.title.x = element_text(size = 15, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),  # Increase legend title size
    legend.text = element_text(size = 12),  # Increase legend text size
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1)
  )

# Print the plot
print(dot_plot)

# Save visualization as SVG & PNG (NO PDFs)
ggsave(paste0(output_folder, "/dot_plot_nadhC.svg"), plot = dot_plot, width = 12, height = 8, dpi = 300, device = "svg")
ggsave(paste0(output_folder, "/dot_plot_nadhC.png"), plot = dot_plot, width = 12, height = 8, dpi = 300)


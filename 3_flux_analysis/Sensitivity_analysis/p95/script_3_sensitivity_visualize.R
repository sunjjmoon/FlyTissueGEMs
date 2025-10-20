library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

rm(list=ls())
   
# Set working directory and load data
current_path <- "E:/Projects/revision/Code_2_upload/Sensitivity_analysis/p95_f"
setwd(current_path)
file_path <- paste0(current_path,"/2_results/","pyruvate_consumption_sensitivity_results_v3.xlsx")

# Create output directory
output_dir <- "3_results_vis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

print(paste("Working directory:", getwd()))

# Load all sheets from Excel
print("Loading Excel sheets...")
detailed_summary <- read_excel(file_path, sheet = "Detailed_Summary")
raw_sample_data <- read_excel(file_path, sheet = "Raw_Sample_Data")
individual_reactions <- read_excel(file_path, sheet = "Individual_Reactions")
summary_data <- read_excel(file_path, sheet = "Summary")

# Clean column names
names(detailed_summary) <- make.names(names(detailed_summary))
names(raw_sample_data) <- make.names(names(raw_sample_data))
names(individual_reactions) <- make.names(names(individual_reactions))

print("Data loaded successfully!")
print(paste("Detailed summary:", nrow(detailed_summary), "rows"))
print(paste("Raw samples:", nrow(raw_sample_data), "rows"))
print(paste("Individual reactions:", nrow(individual_reactions), "rows"))

# Load pyruvate reaction metadata for subsystem information
pyruvate_file_path <- file.path((current_path), "files", "Reactions for metabolite MAM02819c_pyr_c.tsv")


pyruvate_metadata <- read.delim(pyruvate_file_path, stringsAsFactors = FALSE)

# Clean column names
names(pyruvate_metadata) <- make.names(names(pyruvate_metadata))

print("Pyruvate reaction metadata loaded:")
print(paste("Found", nrow(pyruvate_metadata), "pyruvate reactions with subsystem info"))
print("Available subsystems:")
print(unique(pyruvate_metadata$Subsystem))



# =============================================================================
# LOAD GLYCOLYSIS DATA AND CREATE ENHANCED ENZYME LABELS
# =============================================================================
# Load glycolysis reaction data
glycolysis_file_path <- file.path(dirname(pyruvate_file_path), "Reactions for subsystem glycolysis_gluconeogenesis.tsv")
glycolysis_data <- read.delim(glycolysis_file_path, stringsAsFactors = FALSE)

# Clean column names
names(glycolysis_data) <- make.names(names(glycolysis_data))

print("Glycolysis reaction data loaded:")
print(paste("Found", nrow(glycolysis_data), "glycolysis reactions"))
print("Available columns:")
print(colnames(glycolysis_data))

# Function to create enhanced enzyme labels
create_enzyme_label <- function(reaction_id) {
  # Find this reaction in glycolysis_data
  match_row <- glycolysis_data[glycolysis_data$Reaction.ID == reaction_id, ]
  
  if (nrow(match_row) > 0 && !is.na(match_row$Genes) && match_row$Genes != "") {
    genes_str <- match_row$Genes
    
    # Split genes by ';' and take first two
    genes_split <- unlist(strsplit(genes_str, ";"))
    genes_split <- trimws(genes_split)  # Remove whitespace
    
    if (length(genes_split) >= 2) {
      first_two <- paste(genes_split[1:2], collapse = " ")
    } else if (length(genes_split) == 1) {
      first_two <- genes_split[1]
    } else {
      first_two <- ""
    }
    
    # Apply custom gene name mapping (based on your MATLAB code)
    first_two <- switch(first_two,
                        "Ald1" = "Aldo",
                        "Gapdh2 Gapdh1" = "Gapdh1/2",
                        "CG7059 Pglym78" = "Pgam1", 
                        "CG30486 Pgk" = "Pgk",
                        "Hex-A" = "Hex-a/b",
                        "Ih PyK" = "Pyk",
                        first_two)  # default: keep original
    
    # Create combined label
    if (first_two != "") {
      return(paste0(first_two, " (", reaction_id, ")"))
    } else {
      return(reaction_id)
    }
  } else {
    # Fallback: return just the reaction ID
    return(reaction_id)
  }
}

# =============================================================================
# 5. ALL PYRUVATE REACTIONS ANALYSIS (UPDATED)
# =============================================================================

# Top 10 most sensitive enzymes
top_sensitive <- detailed_summary %>%
  arrange(desc(Sensitivity)) %>%
  slice_head(n = 10)

# Load pyruvate reaction annotation filea

pyruvate_file_path <- file.path((current_path), "files", "Reactions for metabolite MAM02819c_pyr_c.tsv")
pyruvate_annotations <- read.delim(pyruvate_file_path, stringsAsFactors = FALSE)

# Clean column names for annotations
names(pyruvate_annotations) <- make.names(names(pyruvate_annotations))

print("Pyruvate reaction annotations loaded:")
print(paste("Found", nrow(pyruvate_annotations), "pyruvate-associated reactions"))


# Create mapping function for gene names
get_gene_label <- function(reaction_id) {
  match_row <- pyruvate_annotations[pyruvate_annotations$Reaction.ID == reaction_id, ]
  if (nrow(match_row) > 0 && !is.na(match_row$Genes) && match_row$Genes != "") {
    # Take first gene if multiple genes separated by semicolon
    first_gene <- strsplit(match_row$Genes, ";")[[1]][1]
    return(paste0(first_gene, " (", reaction_id, ")"))
  } else {
    return(reaction_id)
  }
}

# Create all_reactions data frame
all_reactions <- individual_reactions %>%
  group_by(Pyruvate_Reaction) %>%
  summarise(
    Avg_Baseline = mean(abs(Baseline_Median), na.rm = TRUE),
    Avg_Perturbed = mean(abs(Perturbed_Median), na.rm = TRUE),
    Avg_Change = mean(Flux_Change, na.rm = TRUE),
    N_Observations = n(),
    .groups = "drop"
  ) %>%
  # Add gene labels
  mutate(Gene_Label = sapply(Pyruvate_Reaction, get_gene_label)) %>%
  arrange(desc(Avg_Baseline))

print(paste("Total unique pyruvate reactions found:", nrow(all_reactions)))

print(paste("Total unique pyruvate reactions found:", nrow(all_reactions)))

p5 <- all_reactions %>%
  pivot_longer(cols = c(Avg_Baseline, Avg_Perturbed),
               names_to = "Condition", values_to = "Flux") %>%
  mutate(Condition = ifelse(Condition == "Avg_Baseline", "Baseline", "Perturbed")) %>%
  ggplot(aes(x = reorder(Gene_Label, Flux), y = Flux, fill = Condition)) +
  geom_col(position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = "#2E86AB", "Perturbed" = "#A23B72")) +
  labs(
    title = "5. All Pyruvate-Associated Reactions",
    x = "Pyruvate Reaction (Gene)",
    y = "Average Flux (absolute)",
    fill = "Condition",
    caption = "All reactions shown, ordered by baseline flux magnitude"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )

print(p5)
ggsave(file.path(output_dir, "5_top_pyruvate_reactions.png"), p5, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "5_top_pyruvate_reactions.pdf"), p5, width = 12, height = 8)

# =============================================================================
# 6. CLEAN BUBBLE PLOT: BASELINE VS RELATIVE CHANGE
# =============================================================================
# Create enhanced reaction data
enhanced_reaction_data <- individual_reactions %>%
  # Add gene labels
  mutate(Gene_Label = sapply(Pyruvate_Reaction, get_gene_label)) %>%
  # Calculate baseline flux absolute values
  mutate(
    Baseline_Flux_Abs = abs(Baseline_Median),
    Perturbed_Flux_Abs = abs(Perturbed_Median),
    # Calculate relative flux change
    Relative_Flux_Change = ifelse(abs(Baseline_Median) > 1e-10, 
                                  Flux_Change / abs(Baseline_Median), 
                                  NA)
  ) %>%
  # Filter out very small baseline fluxes
  filter(Baseline_Flux_Abs > 0.01)

# =============================================================================
# 7. STATISTICAL SUMMARY AND INSIGHTS
# =============================================================================

# Calculate key statistics
# Calculate key statistics
stats_summary <- list(
  Most_Sensitive_Enzyme = detailed_summary$Enzyme[which.max(detailed_summary$Sensitivity)],
  Max_Sensitivity = max(detailed_summary$Sensitivity, na.rm = TRUE),
  Max_Sensitivity_SE = detailed_summary$Sensitivity_SE[which.max(detailed_summary$Sensitivity)],
  
  Largest_Response_Enzyme = detailed_summary$Enzyme[which.max(abs(detailed_summary$System_Response))],
  Max_System_Response = max(abs(detailed_summary$System_Response), na.rm = TRUE),
  
  Most_Active_Pyruvate_Reaction = all_reactions$Pyruvate_Reaction[1],
  Max_Pyruvate_Flux = max(all_reactions$Avg_Baseline),
  
  N_Enzymes_Analyzed = nrow(detailed_summary),
  N_Pyruvate_Reactions = length(unique(individual_reactions$Pyruvate_Reaction)),
  N_Total_Samples = nrow(raw_sample_data)
)

# Print summary
cat("\n=== PYRUVATE CONSUMPTION SENSITIVITY ANALYSIS COMPLETE ===\n")
cat("Generated plots in:", output_dir, "\n")
cat("Key Results:\n")
cat("- Most sensitive enzyme:", stats_summary$Most_Sensitive_Enzyme, "\n")
cat("- Max sensitivity:", round(stats_summary$Max_Sensitivity, 4), "±", round(stats_summary$Max_Sensitivity_SE, 4), "\n")
cat("- Largest system response:", round(stats_summary$Max_System_Response, 4), "\n")
cat("- Most active pyruvate reaction:", stats_summary$Most_Active_Pyruvate_Reaction, "\n")
cat("- Total enzymes analyzed:", stats_summary$N_Enzymes_Analyzed, "\n")
cat("- Total pyruvate reactions:", stats_summary$N_Pyruvate_Reactions, "\n")
cat("- Total sample points:", stats_summary$N_Total_Samples, "\n")

cat("\nFiles generated:\n")
cat("- 6 publication-ready plots (.png and .pdf)\n")
cat("- 4 CSV files with processed data\n")
cat("- analysis_statistics.csv with key metrics\n")





# =============================================================================
# 8. INDIVIDUAL PERTURBATION EFFECTS ON PYRUVATE REACTIONS
# =============================================================================

# First, let's check what columns we have in individual_reactions
print("Columns in individual_reactions:")
print(colnames(individual_reactions))

# Prepare data for individual perturbation analysis
individual_perturbation_data <- individual_reactions %>%
  # Add gene labels
  mutate(Gene_Label = sapply(Pyruvate_Reaction, get_gene_label)) %>%
  # Calculate absolute flux change
  mutate(Abs_Flux_Change = abs(Flux_Change)) %>%
  # Filter for meaningful changes and reactions
  filter(abs(Baseline_Median) > 0.1 | abs(Perturbed_Median) > 0.1) %>%
  # Get top reactions by baseline activity
  group_by(Gene_Label) %>%
  mutate(Max_Baseline = max(abs(Baseline_Median), na.rm = TRUE)) %>%
  ungroup() %>%
  # Focus on top 15 most active pyruvate reactions
  filter(Gene_Label %in% (
    individual_reactions %>%
      mutate(Gene_Label = sapply(Pyruvate_Reaction, get_gene_label)) %>%
      group_by(Gene_Label) %>%
      summarise(Max_Baseline = max(abs(Baseline_Median), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Max_Baseline)) %>%
      slice_head(n = 15) %>%
      pull(Gene_Label)
  ))


# Add subsystem information to individual_perturbation_data
individual_perturbation_data <- individual_perturbation_data %>%
  left_join(
    pyruvate_metadata %>% select(Reaction.ID, Subsystem, Equation, C.P),
    by = c("Pyruvate_Reaction" = "Reaction.ID")
  ) %>%
  # Handle missing subsystem info
  mutate(
    Subsystem = ifelse(is.na(Subsystem), "Unknown", Subsystem),
    Subsystem = ifelse(Subsystem == "", "Unknown", Subsystem)
  )

# Check the subsystems we have in our data
print("Subsystems in our perturbation data:")
print(table(individual_perturbation_data$Subsystem))



# Add enhanced enzyme labels to individual_perturbation_data
individual_perturbation_data <- individual_perturbation_data %>%
  mutate(
    Enzyme_Label = sapply(Enzyme, create_enzyme_label)
  )

print("Enhanced enzyme labels created!")
print("Sample of new enzyme labels:")
print(head(unique(individual_perturbation_data[c("Enzyme", "Enzyme_Label")]), 10))


# Option B: Grouped bar chart showing top perturbations for each reaction
p8b <- individual_perturbation_data %>%
  # For each pyruvate reaction, get top 3 perturbations
  group_by(Gene_Label) %>%
  arrange(desc(abs(Flux_Change))) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  # Focus on top 8 most responsive reactions
  filter(Gene_Label %in% (
    individual_perturbation_data %>%
      group_by(Gene_Label) %>%
      summarise(Max_Change = max(abs(Flux_Change), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Max_Change)) %>%
      slice_head(n = 8) %>%
      pull(Gene_Label)
  )) %>%
  ggplot(aes(x = reorder(Gene_Label, abs(Flux_Change)), y = Flux_Change, fill = Enzyme)) +
  geom_col(position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Perturbed\nEnzyme") +
  labs(
    title = "8B. Top Perturbation Effects by Pyruvate Reaction",
    subtitle = "Top 3 enzyme perturbations for each of the 8 most responsive reactions",
    x = "Pyruvate Reaction (Gene)",
    y = "Flux Change",
    caption = "Shows which enzyme perturbations have the strongest effects on each reaction"
  ) +
  labs(
    title = "8B. Top Perturbation Effects by Pyruvate Reaction",
    subtitle = "Top 3 enzyme perturbations for each of the 8 most responsive reactions",
    x = "Pyruvate Reaction (Gene)",
    y = "Flux Change",
    caption = "Shows which enzyme perturbations have the strongest effects on each reaction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p8b)
ggsave(file.path(output_dir, "8b_top_perturbations_per_reaction.png"), p8b, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "8b_top_perturbations_per_reaction.pdf"), p8b, width = 12, height = 10)


cat("\n=== FIGURE 8 ANALYSIS COMPLETE ===\n")
cat("Generated 4 different visualizations of individual perturbation effects:\n")
cat("- 8A: Heatmap showing perturbation effects across reactions\n")
cat("- 8B: Grouped bar chart of top perturbations per reaction\n")
cat("- 8C: Dot plot of top 30 individual perturbation effects\n")
cat("- 8D: Network visualization of strongest perturbation-reaction pairs\n")
cat("\nChoose the visualization that best fits your publication needs!\n")


# Add relative flux change column with proper sign handling
individual_perturbation_data <- individual_perturbation_data %>%
  mutate(
    Rel_Flux_Change = case_when(
      abs(Baseline_Median) < 1e-10 ~ NA_real_,  # Avoid division by very small numbers
      Baseline_Median < 0 ~ Flux_Change / abs(Baseline_Median),  # For negative baseline
      Baseline_Median > 0 ~ Flux_Change / Baseline_Median,       # For positive baseline
      TRUE ~ NA_real_
    )
  )


# =============================================================================
# 10. ALL PYRUVATE REACTIONS - RELATIVE FLUX CHANGES
# =============================================================================

# Calculate everything first
n_genes <- length(unique(individual_perturbation_data$Gene_Label))

p10c <- individual_perturbation_data %>%
  filter(!is.na(Flux_Change)) %>%
  group_by(Gene_Label) %>%
  mutate(Avg_Flux_Change = mean(Flux_Change, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Gene_Label_Ordered = reorder(Gene_Label, Avg_Flux_Change)) %>%
  ggplot(aes(x = Gene_Label_Ordered, y = Flux_Change, fill = Enzyme_Label)) +
  geom_col(position = "dodge", alpha = 1) +
  coord_flip() +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Perturbed\nreactions") +
  labs(
    x = "Pyruvate consumption reactions",
    y = "Flux change (v_pert - v_base)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 8, face = "bold"),   # shrink legend text
    legend.title = element_text(size = 10, face = "bold"),  # shrink legend title
    legend.key.size = unit(0.4, "cm"),
    
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey90", size = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = seq(1.5, n_genes - 0.5, 1),  # Use geom_vline instead
             color = "grey70", size = 0.3, alpha = 0.7)

print(p10c)


ggsave(file.path(output_dir, "10c_all_pyruvate_reactions_flux_changes.png"), p10c, width = 7, height = 6, dpi = 300)
ggsave(file.path(output_dir, "10c_all_pyruvate_reactions_flux_changes.pdf"), p10c, width = 7, height = 6)
ggsave(file.path(output_dir, "10c_all_pyruvate_reactions_flux_changes.svg"), p10c, width = 7, height = 6)





# =============================================================================
# 12. PERTURBATION EFFECTS WITH SUBSYSTEM IN LABELS
# =============================================================================
p12b <- individual_perturbation_data %>%
  filter(!is.na(Flux_Change)) %>%
  arrange(Subsystem, Gene_Label) %>%
  mutate(Gene_Label_Ordered = factor(Gene_Label, levels = unique(Gene_Label))) %>%
  ggplot(aes(x = Enzyme_Label, y = Gene_Label_Ordered)) +  # Use Enzyme_Label for readable names
  geom_point(aes(size = abs(Flux_Change), 
                 color = ifelse(Flux_Change > 0, "Positive", "Negative")),
             alpha = 0.7) +
  scale_size_continuous(name = "Flux\nChange", 
                        range = c(3, 12),  # Increased from (1,8) to (3,12)
                        breaks = c(50, 150, 300),  # Specify actual break values
                        labels = c("Small", "Medium", "Large")) +
  scale_color_manual(values = c("Positive" = "#A23B72", "Negative" = "#2E86AB"),
                     name = "Effect\nDirection") +
  labs(
    title = "Perturbation Effects on Pyruvate Reactions (Absolute Flux Changes)",
    subtitle = "Absolute flux changes by enzyme perturbation",
    x = "Perturbed Enzyme",
    y = "Pyruvate Gene (rxn)",
    caption = "Point size shows absolute flux change magnitude; reactions grouped by subsystem"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),  # Increased from 12
    plot.subtitle = element_text(size = 12),  # Added subtitle size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Increased from 8
    axis.text.y = element_text(size = 10),  # Increased from 7
    axis.title.x = element_text(size = 12),  # Added axis title sizes
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),  # Added legend text size
    legend.title = element_text(size = 11)  # Added legend title size
  )

print(p12b)
ggsave(file.path(output_dir, "12b_perturbation_effects_absolute_flux.png"), p12b, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "12b_perturbation_effects_absolute_flux.pdf"), p12b, width = 8, height = 6)
ggsave(file.path(output_dir, "12b_perturbation_effects_absolute_flux.svg"), p12b, width = 8, height = 6)


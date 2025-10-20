# Script: Create Three-Way Venn Diagrams for FBA vs pFBA vs FVA_Sampling

# Load required libraries
library(readxl)
library(VennDiagram)
library(ggvenn)
library(ggplot2)
library(grid)
library(gridExtra)

rm(list = ls())

# Set working directory and file paths
excel_file <- "E:/Projects/revision/Code_2_upload/pFBA_FBA_FVA/04_sampling_comp/venn_diagram_data.xlsx"
output_dir <- "E:/Projects/revision/Code_2_upload/pFBA_FBA_FVA/06_venn"

# Check if file exists
if (!file.exists(excel_file)) {
  stop("Excel file not found. Please check the path.")
}

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Reading data from:", excel_file, "\n")

# Read the summary counts for manual Venn diagram creation
summary_counts <- read_excel(excel_file, sheet = "Summary_Counts")
print(summary_counts)

# Read the binary matrix for automatic Venn diagram creation
binary_data <- read_excel(excel_file, sheet = "R_Binary_Matrix")

# Create lists for increased reactions
fba_inc <- binary_data$Reaction_ID[binary_data$FBA_Inc == 1]
pfba_inc <- binary_data$Reaction_ID[binary_data$pFBA_Inc == 1]
sampling_inc <- binary_data$Reaction_ID[binary_data$Sampling_Inc == 1]

# Create lists for decreased reactions
fba_dec <- binary_data$Reaction_ID[binary_data$FBA_Dec == 1]
pfba_dec <- binary_data$Reaction_ID[binary_data$pFBA_Dec == 1]
sampling_dec <- binary_data$Reaction_ID[binary_data$Sampling_Dec == 1]

# Remove any NA values
fba_inc <- fba_inc[!is.na(fba_inc)]
pfba_inc <- pfba_inc[!is.na(pfba_inc)]
sampling_inc <- sampling_inc[!is.na(sampling_inc)]
fba_dec <- fba_dec[!is.na(fba_dec)]
pfba_dec <- pfba_dec[!is.na(pfba_dec)]
sampling_dec <- sampling_dec[!is.na(sampling_dec)]

cat("Data loaded successfully:\n")
cat("Increased reactions - FBA:", length(fba_inc), "pFBA:", length(pfba_inc), "Sampling:", length(sampling_inc), "\n")
cat("Decreased reactions - FBA:", length(fba_dec), "pFBA:", length(pfba_dec), "Sampling:", length(sampling_dec), "\n")



# Create ggvenn versions (cleaner looking)
# Prepare data for ggvenn
inc_data <- list(
  FBA = fba_inc,
  pFBA = pfba_inc,
  FVA_sampling = sampling_inc
)

dec_data <- list(
  FBA = fba_dec,
  pFBA = pfba_dec,
  FVA_sampling = sampling_dec
)

# Create increased reactions Venn diagram with ggvenn
p1 <- ggvenn(inc_data, 
             fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
             stroke_size = 1.5,
             set_name_size = 5,
             text_size = 4) +
  ggtitle("Overlap of increased reaction rates (HSD > NSD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

ggsave(file.path(output_dir, "1_Venn_Increased_Reactions_ggvenn.png"), p1, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "1_Venn_Increased_Reactions_ggvenn.svg"), p1, width = 10, height = 8, dpi = 300, bg = "white")

# Create decreased reactions Venn diagram with ggvenn
p2 <- ggvenn(dec_data,
             fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
             stroke_size = 1.5,
             set_name_size = 5,
             text_size = 4) +
  ggtitle("Overlap of decreased reaction rates (HSD < NSD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

ggsave(file.path(output_dir, "1_Venn_Decreased_Reactions_ggvenn.png"), p2, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "1_Venn_Decreased_Reactions_ggvenn.svg"), p2, width = 10, height = 8, dpi = 300, bg = "white")



#################### Euler Venn ##################
# Install if needed: install.packages("eulerr")
library(eulerr)

# Create Euler diagrams (proportional to actual set sizes)
# Prepare data for eulerr
inc_data_euler <- list(
  FBA = fba_inc,
  pFBA = pfba_inc,
  FVA_sampling = sampling_inc
)

dec_data_euler <- list(
  FBA = fba_dec,
  pFBA = pfba_dec,
  FVA_sampling = sampling_dec
)

# Create Euler diagram for increased reactions
inc_euler <- euler(inc_data_euler)

png(file.path(output_dir, "2_Euler_Increased_Reactions.png"), width = 10, height = 8, res = 300, units = "in")
plot(inc_euler, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of increased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()


svg(file.path(output_dir, "2_Euler_Increased_Reactions.svg"),
    width = 10, height = 8)   # no res/units needed, since SVG is vector
plot(inc_euler, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of increased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()



# Create Euler diagram for decreased reactions
dec_euler <- euler(dec_data_euler)

png(file.path(output_dir, "2_Euler_Decreased_Reactions.png"), width = 10, height = 8, res = 300, units = "in")
plot(dec_euler, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of decreased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()

svg(file.path(output_dir, "2_Euler_Decreased_Reactions.svg"),
    width = 10, height = 8)   # no res/units needed, since SVG is vector
plot(dec_euler, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of decreased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()


###############################################################################################################
# Calculate and print overlap statistics
cat("\n=== Venn Diagram Overlap Analysis ===\n")

# For increased reactions
inc_fba_pfba <- length(intersect(fba_inc, pfba_inc))
inc_fba_sampling <- length(intersect(fba_inc, sampling_inc))
inc_pfba_sampling <- length(intersect(pfba_inc, sampling_inc))
inc_all_three <- length(intersect(intersect(fba_inc, pfba_inc), sampling_inc))

cat("\nIncreased reactions overlaps:\n")
cat("FBA only:", length(fba_inc) - inc_fba_pfba - inc_fba_sampling + inc_all_three, "\n")
cat("pFBA only:", length(pfba_inc) - inc_fba_pfba - inc_pfba_sampling + inc_all_three, "\n")
cat("Sampling only:", length(sampling_inc) - inc_fba_sampling - inc_pfba_sampling + inc_all_three, "\n")
cat("FBA + pFBA only:", inc_fba_pfba - inc_all_three, "\n")
cat("FBA + Sampling only:", inc_fba_sampling - inc_all_three, "\n")
cat("pFBA + Sampling only:", inc_pfba_sampling - inc_all_three, "\n")
cat("All three methods:", inc_all_three, "\n")

# For decreased reactions
dec_fba_pfba <- length(intersect(fba_dec, pfba_dec))
dec_fba_sampling <- length(intersect(fba_dec, sampling_dec))
dec_pfba_sampling <- length(intersect(pfba_dec, sampling_dec))
dec_all_three <- length(intersect(intersect(fba_dec, pfba_dec), sampling_dec))

cat("\nDecreased reactions overlaps:\n")
cat("FBA only:", length(fba_dec) - dec_fba_pfba - dec_fba_sampling + dec_all_three, "\n")
cat("pFBA only:", length(pfba_dec) - dec_fba_pfba - dec_pfba_sampling + dec_all_three, "\n")
cat("Sampling only:", length(sampling_dec) - dec_fba_sampling - dec_pfba_sampling + dec_all_three, "\n")
cat("FBA + pFBA only:", dec_fba_pfba - dec_all_three, "\n")
cat("FBA + Sampling only:", dec_fba_sampling - dec_all_three, "\n")
cat("pFBA + Sampling only:", dec_pfba_sampling - dec_all_three, "\n")
cat("All three methods:", dec_all_three, "\n")

# Calculate agreement percentages
total_reactions <- nrow(binary_data)
inc_agreement_percent <- (inc_all_three / max(length(fba_inc), length(pfba_inc), length(sampling_inc))) * 100
dec_agreement_percent <- (dec_all_three / max(length(fba_dec), length(pfba_dec), length(sampling_dec))) * 100

cat("\nAgreement Analysis:\n")
cat("Total reactions analyzed:", total_reactions, "\n")
cat("Three-method agreement on increased reactions:", round(inc_agreement_percent, 1), "%\n")
cat("Three-method agreement on decreased reactions:", round(dec_agreement_percent, 1), "%\n")

# Create summary table for export
overlap_summary <- data.frame(
  Category = c("Increased_FBA_only", "Increased_pFBA_only", "Increased_Sampling_only",
               "Increased_FBA_pFBA", "Increased_FBA_Sampling", "Increased_pFBA_Sampling", "Increased_All_three",
               "Decreased_FBA_only", "Decreased_pFBA_only", "Decreased_Sampling_only",
               "Decreased_FBA_pFBA", "Decreased_FBA_Sampling", "Decreased_pFBA_Sampling", "Decreased_All_three"),
  Count = c(length(fba_inc) - inc_fba_pfba - inc_fba_sampling + inc_all_three,
            length(pfba_inc) - inc_fba_pfba - inc_pfba_sampling + inc_all_three,
            length(sampling_inc) - inc_fba_sampling - inc_pfba_sampling + inc_all_three,
            inc_fba_pfba - inc_all_three, inc_fba_sampling - inc_all_three, inc_pfba_sampling - inc_all_three, inc_all_three,
            length(fba_dec) - dec_fba_pfba - dec_fba_sampling + dec_all_three,
            length(pfba_dec) - dec_fba_pfba - dec_pfba_sampling + dec_all_three,
            length(sampling_dec) - dec_fba_sampling - dec_pfba_sampling + dec_all_three,
            dec_fba_pfba - dec_all_three, dec_fba_sampling - dec_all_three, dec_pfba_sampling - dec_all_three, dec_all_three)
)

# Save summary
write.csv(overlap_summary, file.path(output_dir, "3_Venn_Overlap_Summary.csv"), row.names = FALSE)

cat("\nFiles created:\n")
cat("- Venn_Increased_Reactions_VennDiagram.png\n")
cat("- Venn_Decreased_Reactions_VennDiagram.png\n")
cat("- Venn_Increased_Reactions_ggvenn.png/.svg\n")
cat("- Venn_Decreased_Reactions_ggvenn.png/.svg\n")
cat("- Venn_Overlap_Summary.csv\n")

cat("\nVenn diagram analysis complete!\n")



###############################################################################################################

# Define detail columns for exports
detail_cols <- c("Reaction_ID", "Subsystem", "Reaction_Formula", "Associated_Metabolites")

# Check if additional columns exist and add them if available
available_cols <- colnames(binary_data)
if ("Metabolite_IDs" %in% available_cols) {
  detail_cols <- c(detail_cols, "Metabolite_IDs")
} else if ("Metabolite_ID" %in% available_cols) {
  detail_cols <- c(detail_cols, "Metabolite_ID")
}

# Add gene columns if they exist                    <- ADD THIS SECTION
if ("Associated_Genes" %in% available_cols) {
  detail_cols <- c(detail_cols, "Associated_Genes")
}
if ("Gene_Rules" %in% available_cols) {
  detail_cols <- c(detail_cols, "Gene_Rules")
}



# Add flux difference columns if they exist
flux_diff_cols <- c("FBA_Difference", "pFBA_Difference", "Sampling_Difference")
existing_flux_cols <- flux_diff_cols[flux_diff_cols %in% available_cols]
detail_cols <- c(detail_cols, existing_flux_cols)

cat("Detail columns defined:", paste(detail_cols, collapse = ", "), "\n")


# Get the actual reaction IDs for each intersection
inc_all_three_rxns <- intersect(intersect(fba_inc, pfba_inc), sampling_inc)
dec_all_three_rxns <- intersect(intersect(fba_dec, pfba_dec), sampling_dec)

# Subsystem Analysis for Common Reactions
cat("\n=== Subsystem Analysis for Common Reactions ===\n")

# Analyze subsystems for reactions agreed upon by all three methods
if (length(inc_all_three_rxns) > 0) {
  inc_common_data <- binary_data[binary_data$Reaction_ID %in% inc_all_three_rxns, ]
  inc_subsystem_counts <- table(inc_common_data$Subsystem)
  inc_subsystem_df <- data.frame(
    Subsystem = names(inc_subsystem_counts),
    Count = as.numeric(inc_subsystem_counts),
    Direction = "Increased"
  )
  
  cat("Subsystems for commonly INCREASED reactions:\n")
  print(inc_subsystem_df[order(inc_subsystem_df$Count, decreasing = TRUE), ])
} else {
  inc_subsystem_df <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly increased reactions found.\n")
}

if (length(dec_all_three_rxns) > 0) {
  dec_common_data <- binary_data[binary_data$Reaction_ID %in% dec_all_three_rxns, ]
  dec_subsystem_counts <- table(dec_common_data$Subsystem)
  dec_subsystem_df <- data.frame(
    Subsystem = names(dec_subsystem_counts),
    Count = as.numeric(dec_subsystem_counts),
    Direction = "Decreased"
  )
  
  cat("\nSubsystems for commonly DECREASED reactions:\n")
  print(dec_subsystem_df[order(dec_subsystem_df$Count, decreasing = TRUE), ])
} else {
  dec_subsystem_df <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly decreased reactions found.\n")
}

# Combine subsystem data for visualization
if (nrow(inc_subsystem_df) > 0 || nrow(dec_subsystem_df) > 0) {
  combined_subsystem_df <- rbind(inc_subsystem_df, dec_subsystem_df)
  
  # Create subsystem visualization
  if (nrow(combined_subsystem_df) > 0) {
    # Create bar plot of subsystems
    p_subsystem <- ggplot(combined_subsystem_df, aes(x = reorder(Subsystem, Count), y = Count, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Increased" = "#e74c3c", "Decreased" = "#3498db")) +
      coord_flip() +
      labs(title = "Commonly affected subsystems",
           subtitle = "identified by all three methods (FBA, pFBA, FVA_sampling)",
           x = "",
           y = "Number of reactions",
           fill = "Direction") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 10),
            legend.position = "bottom")
    
    ggsave(file.path(output_dir, "4_Common_Reactions_Subsystem_Analysis.png"), p_subsystem, 
           width = 12, height = 8, dpi = 300, bg = "white")
    
    ggsave(file.path(output_dir, "4_Common_Reactions_Subsystem_Analysis.svg"), p_subsystem, 
           width = 12, height = 8, dpi = 300, bg = "white")
    
    # Export subsystem analysis to CSV
    write.csv(combined_subsystem_df, file.path(output_dir, "4_Common_Reactions_Subsystem_Counts.csv"), row.names = FALSE)
    
    # Create detailed tables for each direction
    # Create detailed tables for each direction
    if (length(inc_all_three_rxns) > 0) {
      inc_detailed <- inc_common_data[, detail_cols]
      write.csv(inc_detailed, file.path(output_dir, "4_Common_Increased_Reactions_Detailed.csv"), row.names = FALSE)
    }
    
    if (length(dec_all_three_rxns) > 0) {
      dec_detailed <- dec_common_data[, detail_cols]
      write.csv(dec_detailed, file.path(output_dir, "4_Common_Decreased_Reactions_Detailed.csv"), row.names = FALSE)
    }
    
    cat("\nSubsystem analysis files created:\n")
    cat("- Common_Reactions_Subsystem_Analysis.png\n")
    cat("- Common_Reactions_Subsystem_Counts.csv\n")
    if (length(inc_all_three_rxns) > 0) cat("- Common_Increased_Reactions_Detailed.csv\n")
    if (length(dec_all_three_rxns) > 0) cat("- Common_Decreased_Reactions_Detailed.csv\n")
  }
}

cat("\nSubsystem analysis for common reactions complete!\n")





# Create separate subsystem visualizations for increased and decreased
if (nrow(inc_subsystem_df) > 0 || nrow(dec_subsystem_df) > 0) {
  
  # Plot for INCREASED reactions
  if (nrow(inc_subsystem_df) > 0) {
    p_inc_subsystem <- ggplot(inc_subsystem_df, aes(x = reorder(Subsystem, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#e74c3c") +
      coord_flip() +
      labs(title = "Commonly INCREASED subsystems",
           subtitle = "identified by all three methods (FBA, pFBA, FVA_sampling)",
           x = "",
           y = "Number of reactions") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 10))
    
    ggsave(file.path(output_dir, "4a_Increased_Subsystem_Analysis.png"), p_inc_subsystem, 
           width = 12, height = max(6, nrow(inc_subsystem_df) * 0.4), dpi = 300, bg = "white")
    
    ggsave(file.path(output_dir, "4a_Increased_Subsystem_Analysis.svg"), p_inc_subsystem, 
           width = 12, height = max(6, nrow(inc_subsystem_df) * 0.4), bg = "white")
    
    cat("Increased subsystem plot created: 4a_Increased_Subsystem_Analysis.png\n")
  }
  
  # Plot for DECREASED reactions
  if (nrow(dec_subsystem_df) > 0) {
    p_dec_subsystem <- ggplot(dec_subsystem_df, aes(x = reorder(Subsystem, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db") +
      coord_flip() +
      labs(title = "Commonly DECREASED subsystems",
           subtitle = "identified by all three methods (FBA, pFBA, FVA_sampling)",
           x = "",
           y = "Number of reactions") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 10))
    
    ggsave(file.path(output_dir, "4b_Decreased_Subsystem_Analysis.png"), p_dec_subsystem, 
           width = 12, height = max(6, nrow(dec_subsystem_df) * 0.4), dpi = 300, bg = "white")
    
    ggsave(file.path(output_dir, "4b_Decreased_Subsystem_Analysis.svg"), p_dec_subsystem, 
           width = 12, height = max(6, nrow(dec_subsystem_df) * 0.4), bg = "white")
    
    cat("Decreased subsystem plot created: 4b_Decreased_Subsystem_Analysis.png\n")
  }
  
  # Export subsystem counts
  combined_subsystem_df <- rbind(inc_subsystem_df, dec_subsystem_df)
  write.csv(combined_subsystem_df, file.path(output_dir, "4_Subsystem_Counts.csv"), row.names = FALSE)
  
  # Create detailed tables for each direction
  if (length(inc_all_three_rxns) > 0) {
    inc_detailed <- inc_common_data[, detail_cols]
    write.csv(inc_detailed, file.path(output_dir, "4_Common_Increased_Reactions_Detailed.csv"), row.names = FALSE)
  }
  
  if (length(dec_all_three_rxns) > 0) {
    dec_detailed <- dec_common_data[, detail_cols]
    write.csv(dec_detailed, file.path(output_dir, "4_Common_Decreased_Reactions_Detailed.csv"), row.names = FALSE)
  }
  
  cat("\nSubsystem analysis files created:\n")
  if (nrow(inc_subsystem_df) > 0) cat("- 4a_Increased_Subsystem_Analysis.png/svg\n")
  if (nrow(dec_subsystem_df) > 0) cat("- 4b_Decreased_Subsystem_Analysis.png/svg\n")
  cat("- 4_Subsystem_Counts.csv\n")
  if (length(inc_all_three_rxns) > 0) cat("- 4_Common_Increased_Reactions_Detailed.csv\n")
  if (length(dec_all_three_rxns) > 0) cat("- 4_Common_Decreased_Reactions_Detailed.csv\n")
}

cat("\nSubsystem analysis for common reactions complete!\n")



########################################################################################################################
## 5.Simple gene-level visualization for commonly affected reactions
# Simple gene-level visualization for commonly affected reactions

# Function to extract genes and count them
count_genes <- function(common_data, direction) {
  if (nrow(common_data) == 0 || !"Associated_Genes" %in% colnames(common_data)) {
    return(data.frame())
  }
  
  all_genes <- c()
  for (i in 1:nrow(common_data)) {
    gene_string <- common_data$Associated_Genes[i]
    if (!is.na(gene_string) && gene_string != "No genes associated") {
      genes <- strsplit(gene_string, "; ")[[1]]
      all_genes <- c(all_genes, genes)
    }
  }
  
  if (length(all_genes) == 0) return(data.frame())
  
  gene_counts <- table(all_genes)
  data.frame(
    Gene = names(gene_counts),
    Count = as.numeric(gene_counts),
    Direction = direction,
    stringsAsFactors = FALSE
  )
}

# Count genes for each direction
inc_genes <- count_genes(inc_common_data, "Increased")
dec_genes <- count_genes(dec_common_data, "Decreased")

# Combine all genes
all_genes <- rbind(inc_genes, dec_genes)

if (nrow(all_genes) > 0) {
  # Create separate plots for increased and decreased
  if (nrow(inc_genes) > 0) {
    p_inc_genes <- ggplot(inc_genes, aes(x = reorder(Gene, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#e74c3c") +
      coord_flip() +
      labs(title = "Genes in commonly INCREASED reactions",
           x = "Gene",
           y = "Number of reactions") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "5a_Increased_Gene_Analysis.png"), 
           p_inc_genes, 
           width = 10, height = min(max(6, nrow(inc_genes) * 0.3), 20), dpi = 300, bg = "white")
    cat("Increased gene plot created: 5a_Increased_Gene_Analysis.png\n")
  }
  
  if (nrow(dec_genes) > 0) {
    p_dec_genes <- ggplot(dec_genes, aes(x = reorder(Gene, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db") +
      coord_flip() +
      labs(title = "Genes in commonly DECREASED reactions",
           x = "Gene", 
           y = "Number of reactions") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "5b_Decreased_Gene_Analysis.png"), p_dec_genes, 
           width = 10, height = min(max(6, nrow(inc_genes) * 0.3), 20), dpi = 300, bg = "white")
    cat("Decreased gene plot created: 5b_Decreased_Gene_Analysis.png\n")
  }
  
  # Export gene counts
  write.csv(all_genes, file.path(output_dir, "5_Gene_Counts.csv"), row.names = FALSE)
  cat("Gene data exported: 5_Gene_Counts.csv\n")
  
  # Print gene summary
  cat("\nGene summary:\n")
  if (nrow(inc_genes) > 0) {
    cat("Increased reactions genes:\n")
    print(inc_genes[order(inc_genes$Count, decreasing = TRUE), ])
  }
  if (nrow(dec_genes) > 0) {
    cat("\nDecreased reactions genes:\n")
    print(dec_genes[order(dec_genes$Count, decreasing = TRUE), ])
  }
} else {
  cat("No gene data available.\n")
}


#################################################################################################################
##################### Filtered ###############
# Create filtered Venn diagrams (excluding exchange/demand reactions)
cat("\n=== Creating Filtered Venn Diagrams (No Exchange/Demand) ===\n")

# Filter out exchange/demand reactions based on subsystem
cat("Filtering out exchange and demand reactions...\n")
cat("Total reactions before filtering:", nrow(binary_data), "\n")

# Filter out reactions where Subsystem is "Exchange/demand reactions"
binary_data_filtered <- binary_data[binary_data$Subsystem != "Exchange/demand reactions", ]

cat("Total reactions after filtering:", nrow(binary_data_filtered), "\n")
cat("Removed", nrow(binary_data) - nrow(binary_data_filtered), "exchange/demand reactions\n")

# Create new filtered lists for increased reactions
fba_inc_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$FBA_Inc == 1]
pfba_inc_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$pFBA_Inc == 1]
sampling_inc_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$Sampling_Inc == 1]

# Create new filtered lists for decreased reactions
fba_dec_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$FBA_Dec == 1]
pfba_dec_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$pFBA_Dec == 1]
sampling_dec_filtered <- binary_data_filtered$Reaction_ID[binary_data_filtered$Sampling_Dec == 1]

# Remove any NA values
fba_inc_filtered <- fba_inc_filtered[!is.na(fba_inc_filtered)]
pfba_inc_filtered <- pfba_inc_filtered[!is.na(pfba_inc_filtered)]
sampling_inc_filtered <- sampling_inc_filtered[!is.na(sampling_inc_filtered)]
fba_dec_filtered <- fba_dec_filtered[!is.na(fba_dec_filtered)]
pfba_dec_filtered <- pfba_dec_filtered[!is.na(pfba_dec_filtered)]
sampling_dec_filtered <- sampling_dec_filtered[!is.na(sampling_dec_filtered)]

cat("Filtered data loaded successfully:\n")
cat("Increased reactions - FBA:", length(fba_inc_filtered), "pFBA:", length(pfba_inc_filtered), "Sampling:", length(sampling_inc_filtered), "\n")
cat("Decreased reactions - FBA:", length(fba_dec_filtered), "pFBA:", length(pfba_dec_filtered), "Sampling:", length(sampling_dec_filtered), "\n")

# Create filtered Venn diagrams for INCREASED reactions
inc_data_filtered <- list(
  FBA = fba_inc_filtered,
  pFBA = pfba_inc_filtered,
  FVA_sampling = sampling_inc_filtered
)

p1_filtered <- ggvenn(inc_data_filtered, 
                      fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
                      stroke_size = 1.5,
                      set_name_size = 5,
                      text_size = 4) +
  ggtitle("Overlap of increased reaction rates (HSD > NSD)") +

  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "6_Venn_Increased_Reactions_Filtered.png"), p1_filtered, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "6_Venn_Increased_Reactions_Filtered.svg"), p1_filtered, width = 10, height = 8, dpi = 300, bg = "white")

# Create filtered Venn diagrams for DECREASED reactions
dec_data_filtered <- list(
  FBA = fba_dec_filtered,
  pFBA = pfba_dec_filtered,
  FVA_sampling = sampling_dec_filtered
)

p2_filtered <- ggvenn(dec_data_filtered,
                      fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
                      stroke_size = 1.5,
                      set_name_size = 5,
                      text_size = 4) +
  ggtitle("Overlap of decreased reaction rates (HSD > NSD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "6_Venn_Decreased_Reactions_Filtered.png"), p2_filtered, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "6_Venn_Decreased_Reactions_Filtered.svg"), p2_filtered, width = 10, height = 8, dpi = 300, bg = "white")

# Calculate filtered overlap statistics
cat("\n=== Filtered Venn Diagram Overlap Analysis ===\n")

# For filtered increased reactions
inc_fba_pfba_filtered <- length(intersect(fba_inc_filtered, pfba_inc_filtered))
inc_fba_sampling_filtered <- length(intersect(fba_inc_filtered, sampling_inc_filtered))
inc_pfba_sampling_filtered <- length(intersect(pfba_inc_filtered, sampling_inc_filtered))
inc_all_three_filtered <- length(intersect(intersect(fba_inc_filtered, pfba_inc_filtered), sampling_inc_filtered))

cat("\nFiltered Increased reactions overlaps:\n")
cat("FBA only:", length(fba_inc_filtered) - inc_fba_pfba_filtered - inc_fba_sampling_filtered + inc_all_three_filtered, "\n")
cat("pFBA only:", length(pfba_inc_filtered) - inc_fba_pfba_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered, "\n")
cat("Sampling only:", length(sampling_inc_filtered) - inc_fba_sampling_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered, "\n")
cat("FBA + pFBA only:", inc_fba_pfba_filtered - inc_all_three_filtered, "\n")
cat("FBA + Sampling only:", inc_fba_sampling_filtered - inc_all_three_filtered, "\n")
cat("pFBA + Sampling only:", inc_pfba_sampling_filtered - inc_all_three_filtered, "\n")
cat("All three methods:", inc_all_three_filtered, "\n")

# For filtered decreased reactions
dec_fba_pfba_filtered <- length(intersect(fba_dec_filtered, pfba_dec_filtered))
dec_fba_sampling_filtered <- length(intersect(fba_dec_filtered, sampling_dec_filtered))
dec_pfba_sampling_filtered <- length(intersect(pfba_dec_filtered, sampling_dec_filtered))
dec_all_three_filtered <- length(intersect(intersect(fba_dec_filtered, pfba_dec_filtered), sampling_dec_filtered))

cat("\nFiltered Decreased reactions overlaps:\n")
cat("FBA only:", length(fba_dec_filtered) - dec_fba_pfba_filtered - dec_fba_sampling_filtered + dec_all_three_filtered, "\n")
cat("pFBA only:", length(pfba_dec_filtered) - dec_fba_pfba_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered, "\n")
cat("Sampling only:", length(sampling_dec_filtered) - dec_fba_sampling_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered, "\n")
cat("FBA + pFBA only:", dec_fba_pfba_filtered - dec_all_three_filtered, "\n")
cat("FBA + Sampling only:", dec_fba_sampling_filtered - dec_all_three_filtered, "\n")
cat("pFBA + Sampling only:", dec_pfba_sampling_filtered - dec_all_three_filtered, "\n")
cat("All three methods:", dec_all_three_filtered, "\n")

# Export filtered summary
filtered_overlap_summary <- data.frame(
  Category = c("Filtered_Increased_FBA_only", "Filtered_Increased_pFBA_only", "Filtered_Increased_Sampling_only",
               "Filtered_Increased_FBA_pFBA", "Filtered_Increased_FBA_Sampling", "Filtered_Increased_pFBA_Sampling", "Filtered_Increased_All_three",
               "Filtered_Decreased_FBA_only", "Filtered_Decreased_pFBA_only", "Filtered_Decreased_Sampling_only",
               "Filtered_Decreased_FBA_pFBA", "Filtered_Decreased_FBA_Sampling", "Filtered_Decreased_pFBA_Sampling", "Filtered_Decreased_All_three"),
  Count = c(length(fba_inc_filtered) - inc_fba_pfba_filtered - inc_fba_sampling_filtered + inc_all_three_filtered,
            length(pfba_inc_filtered) - inc_fba_pfba_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered,
            length(sampling_inc_filtered) - inc_fba_sampling_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered,
            inc_fba_pfba_filtered - inc_all_three_filtered, inc_fba_sampling_filtered - inc_all_three_filtered, inc_pfba_sampling_filtered - inc_all_three_filtered, inc_all_three_filtered,
            length(fba_dec_filtered) - dec_fba_pfba_filtered - dec_fba_sampling_filtered + dec_all_three_filtered,
            length(pfba_dec_filtered) - dec_fba_pfba_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered,
            length(sampling_dec_filtered) - dec_fba_sampling_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered,
            dec_fba_pfba_filtered - dec_all_three_filtered, dec_fba_sampling_filtered - dec_all_three_filtered, dec_pfba_sampling_filtered - dec_all_three_filtered, dec_all_three_filtered)
)

write.csv(filtered_overlap_summary, file.path(output_dir, "6_Venn_Overlap_Summary_Filtered.csv"), row.names = FALSE)

cat("\nFiltered Venn diagram files created:\n")
cat("- Venn_Increased_Reactions_Filtered.png/.svg\n")
cat("- Venn_Decreased_Reactions_Filtered.png/.svg\n")
cat("- Venn_Overlap_Summary_Filtered.csv\n")

cat("\nFiltered Venn diagram analysis complete!\n")


#################### Euler Venn ##################
# Install if needed: install.packages("eulerr")
#library(eulerr)

# Create Euler diagrams (proportional to actual set sizes)
# Prepare FILTERED data for eulerr
inc_data_euler_filtered <- list(
  FBA = fba_inc_filtered,
  pFBA = pfba_inc_filtered,
  FVA_sampling = sampling_inc_filtered
)

dec_data_euler_filtered <- list(
  FBA = fba_dec_filtered,
  pFBA = pfba_dec_filtered,
  FVA_sampling = sampling_dec_filtered
)

# Create Euler diagram for increased reactions
inc_euler_filtered <- euler(inc_data_euler_filtered)

png(file.path(output_dir, "7_Euler_Increased_Reactions_Filtered.png"), width = 10, height = 8, res = 300, units = "in")
plot(inc_euler_filtered, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of increased reaction rates (HSD > NSD) - Filtered",
     main.cex = 1.5)
dev.off()

svg(file.path(output_dir, "7_Euler_Increased_Reactions_Filtered.svg"),
    width = 10, height = 8)   # no res/units needed, since SVG is vector
plot(inc_euler_filtered, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of increased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()



# Create Euler diagram for decreased reactions
dec_euler_filtered <- euler(dec_data_euler_filtered)

png(file.path(output_dir, "7_Euler_Decreased_Reactions_filtered.png"), width = 10, height = 8, res = 300, units = "in")
plot(dec_euler_filtered, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of decreased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()

svg(file.path(output_dir, "7_Euler_Decreased_Reactions_filtered.svg"),
    width = 10, height = 8)   # no res/units needed, since SVG is vector
plot(dec_euler_filtered, 
     fills = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
     alpha = 0.7,
     labels = list(font = 2, cex = 1.2),
     edges = list(col = "black", lwd = 2),
     main = "Overlap of decreased reaction rates (HSD > NSD)",
     main.cex = 1.5)
dev.off()


# Calculate and print FILTERED overlap statistics
cat("\n=== Filtered Venn Diagram Overlap Analysis ===\n")

# For FILTERED increased reactions
inc_fba_pfba_filtered <- length(intersect(fba_inc_filtered, pfba_inc_filtered))
inc_fba_sampling_filtered <- length(intersect(fba_inc_filtered, sampling_inc_filtered))
inc_pfba_sampling_filtered <- length(intersect(pfba_inc_filtered, sampling_inc_filtered))
inc_all_three_filtered <- length(intersect(intersect(fba_inc_filtered, pfba_inc_filtered), sampling_inc_filtered))

cat("\nFiltered Increased reactions overlaps:\n")
cat("FBA only:", length(fba_inc_filtered) - inc_fba_pfba_filtered - inc_fba_sampling_filtered + inc_all_three_filtered, "\n")
cat("pFBA only:", length(pfba_inc_filtered) - inc_fba_pfba_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered, "\n")
cat("Sampling only:", length(sampling_inc_filtered) - inc_fba_sampling_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered, "\n")
cat("FBA + pFBA only:", inc_fba_pfba_filtered - inc_all_three_filtered, "\n")
cat("FBA + Sampling only:", inc_fba_sampling_filtered - inc_all_three_filtered, "\n")
cat("pFBA + Sampling only:", inc_pfba_sampling_filtered - inc_all_three_filtered, "\n")
cat("All three methods:", inc_all_three_filtered, "\n")

# For FILTERED decreased reactions
dec_fba_pfba_filtered <- length(intersect(fba_dec_filtered, pfba_dec_filtered))
dec_fba_sampling_filtered <- length(intersect(fba_dec_filtered, sampling_dec_filtered))
dec_pfba_sampling_filtered <- length(intersect(pfba_dec_filtered, sampling_dec_filtered))
dec_all_three_filtered <- length(intersect(intersect(fba_dec_filtered, pfba_dec_filtered), sampling_dec_filtered))

cat("\nFiltered Decreased reactions overlaps:\n")
cat("FBA only:", length(fba_dec_filtered) - dec_fba_pfba_filtered - dec_fba_sampling_filtered + dec_all_three_filtered, "\n")
cat("pFBA only:", length(pfba_dec_filtered) - dec_fba_pfba_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered, "\n")
cat("Sampling only:", length(sampling_dec_filtered) - dec_fba_sampling_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered, "\n")
cat("FBA + pFBA only:", dec_fba_pfba_filtered - dec_all_three_filtered, "\n")
cat("FBA + Sampling only:", dec_fba_sampling_filtered - dec_all_three_filtered, "\n")
cat("pFBA + Sampling only:", dec_pfba_sampling_filtered - dec_all_three_filtered, "\n")
cat("All three methods:", dec_all_three_filtered, "\n")

# Calculate FILTERED agreement percentages
total_reactions_filtered <- nrow(binary_data_filtered)
inc_agreement_percent_filtered <- (inc_all_three_filtered / max(length(fba_inc_filtered), length(pfba_inc_filtered), length(sampling_inc_filtered))) * 100
dec_agreement_percent_filtered <- (dec_all_three_filtered / max(length(fba_dec_filtered), length(pfba_dec_filtered), length(sampling_dec_filtered))) * 100

cat("\nFiltered Agreement Analysis:\n")
cat("Total reactions analyzed (filtered):", total_reactions_filtered, "\n")
cat("Three-method agreement on increased reactions:", round(inc_agreement_percent_filtered, 1), "%\n")
cat("Three-method agreement on decreased reactions:", round(dec_agreement_percent_filtered, 1), "%\n")

# Create FILTERED summary table for export
overlap_summary_filtered <- data.frame(
  Category = c("Filtered_Increased_FBA_only", "Filtered_Increased_pFBA_only", "Filtered_Increased_Sampling_only",
               "Filtered_Increased_FBA_pFBA", "Filtered_Increased_FBA_Sampling", "Filtered_Increased_pFBA_Sampling", "Filtered_Increased_All_three",
               "Filtered_Decreased_FBA_only", "Filtered_Decreased_pFBA_only", "Filtered_Decreased_Sampling_only",
               "Filtered_Decreased_FBA_pFBA", "Filtered_Decreased_FBA_Sampling", "Filtered_Decreased_pFBA_Sampling", "Filtered_Decreased_All_three"),
  Count = c(length(fba_inc_filtered) - inc_fba_pfba_filtered - inc_fba_sampling_filtered + inc_all_three_filtered,
            length(pfba_inc_filtered) - inc_fba_pfba_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered,
            length(sampling_inc_filtered) - inc_fba_sampling_filtered - inc_pfba_sampling_filtered + inc_all_three_filtered,
            inc_fba_pfba_filtered - inc_all_three_filtered, inc_fba_sampling_filtered - inc_all_three_filtered, inc_pfba_sampling_filtered - inc_all_three_filtered, inc_all_three_filtered,
            length(fba_dec_filtered) - dec_fba_pfba_filtered - dec_fba_sampling_filtered + dec_all_three_filtered,
            length(pfba_dec_filtered) - dec_fba_pfba_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered,
            length(sampling_dec_filtered) - dec_fba_sampling_filtered - dec_pfba_sampling_filtered + dec_all_three_filtered,
            dec_fba_pfba_filtered - dec_all_three_filtered, dec_fba_sampling_filtered - dec_all_three_filtered, dec_pfba_sampling_filtered - dec_all_three_filtered, dec_all_three_filtered)
)

# Save FILTERED summary
write.csv(overlap_summary_filtered, file.path(output_dir, "8_Venn_Overlap_Summary_Filtered.csv"), row.names = FALSE)

cat("\nFiltered files created:\n")
cat("- 6_Venn_Increased_Reactions_Filtered.png/.svg\n")
cat("- 6_Venn_Decreased_Reactions_Filtered.png/.svg\n")
cat("- 7_Euler_Increased_Reactions_Filtered.png/.svg\n")
cat("- 7_Euler_Decreased_Reactions_Filtered.png/.svg\n")
cat("- 8_Venn_Overlap_Summary_Filtered.csv\n")

cat("\nFiltered Venn diagram analysis complete!\n")


# Get the actual reaction IDs for each FILTERED intersection
inc_all_three_rxns_filtered <- intersect(intersect(fba_inc_filtered, pfba_inc_filtered), sampling_inc_filtered)
dec_all_three_rxns_filtered <- intersect(intersect(fba_dec_filtered, pfba_dec_filtered), sampling_dec_filtered)

# Subsystem Analysis for Common FILTERED Reactions
cat("\n=== Subsystem Analysis for Common FILTERED Reactions ===\n")

# Analyze subsystems for FILTERED reactions agreed upon by all three methods
if (length(inc_all_three_rxns_filtered) > 0) {
  inc_common_data_filtered <- binary_data_filtered[binary_data_filtered$Reaction_ID %in% inc_all_three_rxns_filtered, ]
  inc_subsystem_counts_filtered <- table(inc_common_data_filtered$Subsystem)
  inc_subsystem_df_filtered <- data.frame(
    Subsystem = names(inc_subsystem_counts_filtered),
    Count = as.numeric(inc_subsystem_counts_filtered),
    Direction = "Increased"
  )
  
  cat("Subsystems for commonly INCREASED reactions (filtered):\n")
  print(inc_subsystem_df_filtered[order(inc_subsystem_df_filtered$Count, decreasing = TRUE), ])
} else {
  inc_subsystem_df_filtered <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly increased reactions found (filtered).\n")
}

if (length(dec_all_three_rxns_filtered) > 0) {
  dec_common_data_filtered <- binary_data_filtered[binary_data_filtered$Reaction_ID %in% dec_all_three_rxns_filtered, ]
  dec_subsystem_counts_filtered <- table(dec_common_data_filtered$Subsystem)
  dec_subsystem_df_filtered <- data.frame(
    Subsystem = names(dec_subsystem_counts_filtered),
    Count = as.numeric(dec_subsystem_counts_filtered),
    Direction = "Decreased"
  )
  
  cat("\nSubsystems for commonly DECREASED reactions (filtered):\n")
  print(dec_subsystem_df_filtered[order(dec_subsystem_df_filtered$Count, decreasing = TRUE), ])
} else {
  dec_subsystem_df_filtered <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly decreased reactions found (filtered).\n")
}

# Create separate FILTERED subsystem visualizations for increased and decreased
if (nrow(inc_subsystem_df_filtered) > 0 || nrow(dec_subsystem_df_filtered) > 0) {
  
  # Plot for INCREASED FILTERED reactions
  if (nrow(inc_subsystem_df_filtered) > 0) {
    p_inc_subsystem_filtered <- ggplot(inc_subsystem_df_filtered, aes(x = reorder(Subsystem, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#e74c3c") +
      coord_flip() +
      labs(title = "Commonly affected subsystems",
           subtitle = "identified by all three methods",
           x = "",
           y = "Number of reactions") +
      
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size = 15),     # axis labe
          axis.text.y = element_text(size = 15))
    
    ggsave(file.path(output_dir, "9a_Increased_Subsystem_Analysis_Filtered.png"), p_inc_subsystem_filtered, 
           width = 12, height = max(6, nrow(inc_subsystem_df_filtered) * 0.4), dpi = 300, bg = "white")
    
    ggsave(file.path(output_dir, "9a_Increased_Subsystem_Analysis_Filtered.svg"), p_inc_subsystem_filtered, 
           width = 10, height = max(6, nrow(inc_subsystem_df_filtered) * 0.4), bg = "white")
    
    ggsave(file.path(output_dir, "9a_Increased_Subsystem_Analysis_Filtered.pdf"), p_inc_subsystem_filtered, 
           width = 10, height = max(6, nrow(inc_subsystem_df_filtered) * 0.4))
    
    cat("Filtered increased subsystem plot created: 9a_Increased_Subsystem_Analysis_Filtered.png\n")
  }
  
  # Plot for DECREASED FILTERED reactions
  if (nrow(dec_subsystem_df_filtered) > 0) {
    p_dec_subsystem_filtered <- ggplot(dec_subsystem_df_filtered, aes(x = reorder(Subsystem, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db") +
      coord_flip() +
      labs(title = "Commonly affected subsystems",
           subtitle = "identified by all three methods",
           x = "",
           y = "Number of reactions") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 10))
    
    ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered.png"), p_dec_subsystem_filtered, 
           width = 12, height = max(6, nrow(dec_subsystem_df_filtered) * 0.4), dpi = 300, bg = "white")
    
    ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered.pdf"), p_dec_subsystem_filtered, 
           width = 10, height = max(6, nrow(inc_subsystem_df_filtered) * 0.4), bg = "white")
    
    ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered.svg"), p_dec_subsystem_filtered, 
           width = 10, height = max(6, nrow(dec_subsystem_df_filtered) * 0.4), bg = "white")
    
    cat("Filtered decreased subsystem plot created: 9b_Decreased_Subsystem_Analysis_Filtered.png\n")
  }
  
  # Export FILTERED subsystem counts
  combined_subsystem_df_filtered <- rbind(inc_subsystem_df_filtered, dec_subsystem_df_filtered)
  write.csv(combined_subsystem_df_filtered, file.path(output_dir, "9_Subsystem_Counts_Filtered.csv"), row.names = FALSE)
  
  # Create detailed FILTERED tables for each direction
  if (length(inc_all_three_rxns_filtered) > 0) {
    inc_detailed_filtered <- inc_common_data_filtered[, detail_cols]
    write.csv(inc_detailed_filtered, file.path(output_dir, "9_Common_Increased_Reactions_Detailed_Filtered.csv"), row.names = FALSE)
  }
  
  if (length(dec_all_three_rxns_filtered) > 0) {
    dec_detailed_filtered <- dec_common_data_filtered[, detail_cols]
    write.csv(dec_detailed_filtered, file.path(output_dir, "9_Common_Decreased_Reactions_Detailed_Filtered.csv"), row.names = FALSE)
  }
  
  cat("\nFiltered subsystem analysis files created:\n")
  if (nrow(inc_subsystem_df_filtered) > 0) cat("- 9a_Increased_Subsystem_Analysis_Filtered.png/svg\n")
  if (nrow(dec_subsystem_df_filtered) > 0) cat("- 9b_Decreased_Subsystem_Analysis_Filtered.png/svg\n")
  cat("- 9_Subsystem_Counts_Filtered.csv\n")
  if (length(inc_all_three_rxns_filtered) > 0) cat("- 9_Common_Increased_Reactions_Detailed_Filtered.csv\n")
  if (length(dec_all_three_rxns_filtered) > 0) cat("- 9_Common_Decreased_Reactions_Detailed_Filtered.csv\n")
}

cat("\nFiltered subsystem analysis for common reactions complete!\n")



# Plot for DECREASED FILTERED reactions
if (nrow(dec_subsystem_df_filtered) > 0) {
  # Keep only top 10 subsystems
  dec_subsystem_df_top10 <- head(dec_subsystem_df_filtered[order(dec_subsystem_df_filtered$Count, decreasing = TRUE), ], 10)
  
  p_dec_subsystem_filtered <- ggplot(dec_subsystem_df_top10, aes(x = reorder(Subsystem, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "#3498db") +
    coord_flip() +
    labs(title = "Top 10 commonly affected subsystems",
         subtitle = "identified by all three methods",
         x = "",
         y = "Number of reactions") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size = 15),     # axis labe
          axis.text.y = element_text(size = 15))
  
  
  ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered_top10.png"), p_dec_subsystem_filtered, 
         width = 12, height = max(6, nrow(dec_subsystem_df_filtered) * 0.4), dpi = 300, bg = "white")
  
  ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered_top10.pdf"), p_dec_subsystem_filtered, 
         width = 10, height = max(6, nrow(inc_subsystem_df_filtered) * 0.4), bg = "white")
  
  ggsave(file.path(output_dir, "9b_Decreased_Subsystem_Analysis_Filtered_top10.svg"), p_dec_subsystem_filtered, 
         width = 10, height = max(6, nrow(dec_subsystem_df_filtered) * 0.4), bg = "white")
}


# Simple gene-level visualization for commonly affected FILTERED reactions

# Count genes for each direction (FILTERED data)
inc_genes_filtered <- count_genes(inc_common_data_filtered, "Increased")
dec_genes_filtered <- count_genes(dec_common_data_filtered, "Decreased")

# Combine all FILTERED genes
all_genes_filtered <- rbind(inc_genes_filtered, dec_genes_filtered)

if (nrow(all_genes_filtered) > 0) {
  # Create separate plots for increased and decreased FILTERED reactions
  if (nrow(inc_genes_filtered) > 0) {
    p_inc_genes_filtered <- ggplot(inc_genes_filtered, aes(x = reorder(Gene, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#e74c3c") +
      coord_flip() +
      labs(title = "Genes in commonly INCREASED reactions (Filtered)",
           subtitle = "no exchange/demand reactions",
           x = "Gene",
           y = "Number of reactions") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "10a_Increased_Gene_Analysis_Filtered.png"), p_inc_genes_filtered, 
           width = 10, height = max(6, nrow(inc_genes_filtered) * 0.3), dpi = 300, bg = "white")
    cat("Filtered increased gene plot created: 10a_Increased_Gene_Analysis_Filtered.png\n")
  }
  
  if (nrow(dec_genes_filtered) > 0) {
    p_dec_genes_filtered <- ggplot(dec_genes_filtered, aes(x = reorder(Gene, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db") +
      coord_flip() +
      labs(title = "Genes in commonly DECREASED reactions (Filtered)",
           subtitle = "no exchange/demand reactions",
           x = "Gene", 
           y = "Number of reactions") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "10b_Decreased_Gene_Analysis_Filtered.png"), p_dec_genes_filtered, 
           width = 10, height = max(6, nrow(dec_genes_filtered) * 0.3), dpi = 300, bg = "white")
    cat("Filtered decreased gene plot created: 10b_Decreased_Gene_Analysis_Filtered.png\n")
  }
  
  # Export FILTERED gene counts
  write.csv(all_genes_filtered, file.path(output_dir, "10_Gene_Counts_Filtered.csv"), row.names = FALSE)
  cat("Filtered gene data exported: 10_Gene_Counts_Filtered.csv\n")
  
  # Print FILTERED gene summary
  cat("\nFiltered gene summary:\n")
  if (nrow(inc_genes_filtered) > 0) {
    cat("Increased reactions genes (filtered):\n")
    print(inc_genes_filtered[order(inc_genes_filtered$Count, decreasing = TRUE), ])
  }
  if (nrow(dec_genes_filtered) > 0) {
    cat("\nDecreased reactions genes (filtered):\n")
    print(dec_genes_filtered[order(dec_genes_filtered$Count, decreasing = TRUE), ])
  }
} else {
  cat("No filtered gene data available.\n")
}

# Get the actual reaction IDs for each intersection (FILTERED DATA)
inc_all_three_rxns_filtered <- intersect(intersect(fba_inc_filtered, pfba_inc_filtered), sampling_inc_filtered)
dec_all_three_rxns_filtered <- intersect(intersect(fba_dec_filtered, pfba_dec_filtered), sampling_dec_filtered)

# Subsystem Analysis for Common Reactions (FILTERED - No Exchange/Demand)
cat("\n=== Subsystem Analysis for Common Metabolic Reactions ===\n")

# Analyze subsystems for reactions agreed upon by all three methods (filtered data)
if (length(inc_all_three_rxns_filtered) > 0) {
  inc_common_data_filtered <- binary_data_filtered[binary_data_filtered$Reaction_ID %in% inc_all_three_rxns_filtered, ]
  inc_subsystem_counts_filtered <- table(inc_common_data_filtered$Subsystem)
  inc_subsystem_df_filtered <- data.frame(
    Subsystem = names(inc_subsystem_counts_filtered),
    Count = as.numeric(inc_subsystem_counts_filtered),
    Direction = "Increased"
  )
  
  cat("Subsystems for commonly INCREASED metabolic reactions:\n")
  print(inc_subsystem_df_filtered[order(inc_subsystem_df_filtered$Count, decreasing = TRUE), ])
} else {
  inc_subsystem_df_filtered <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly increased metabolic reactions found.\n")
}

if (length(dec_all_three_rxns_filtered) > 0) {
  dec_common_data_filtered <- binary_data_filtered[binary_data_filtered$Reaction_ID %in% dec_all_three_rxns_filtered, ]
  dec_subsystem_counts_filtered <- table(dec_common_data_filtered$Subsystem)
  dec_subsystem_df_filtered <- data.frame(
    Subsystem = names(dec_subsystem_counts_filtered),
    Count = as.numeric(dec_subsystem_counts_filtered),
    Direction = "Decreased"
  )
  
  cat("\nSubsystems for commonly DECREASED metabolic reactions:\n")
  print(dec_subsystem_df_filtered[order(dec_subsystem_df_filtered$Count, decreasing = TRUE), ])
} else {
  dec_subsystem_df_filtered <- data.frame(Subsystem = character(0), Count = numeric(0), Direction = character(0))
  cat("No commonly decreased metabolic reactions found.\n")
}




# Combine filtered subsystem data for visualization
if (nrow(inc_subsystem_df_filtered) > 0 || nrow(dec_subsystem_df_filtered) > 0) {
  combined_subsystem_df_filtered <- rbind(inc_subsystem_df_filtered, dec_subsystem_df_filtered)
  
  # Create subsystem visualization (filtered)
  if (nrow(combined_subsystem_df_filtered) > 0) {
    # Create bar plot of subsystems (filtered)
    p_subsystem_filtered <- ggplot(combined_subsystem_df_filtered, aes(x = reorder(Subsystem, Count), y = Count, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("Increased" = "#e74c3c", "Decreased" = "#3498db")) +
      coord_flip() +
      labs(title = "Subsystem Distribution of Commonly Altered Metabolic Reactions",
           subtitle = "Reactions identified by all three methods (FBA, pFBA, Sampling) - Exchange/Demand Excluded",
           x = "Subsystem",
           y = "Number of Reactions",
           fill = "Direction") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 10),
            legend.position = "bottom")
    
    ggsave(file.path(output_dir, "9_Common_Metabolic_Reactions_Subsystem_Analysis.png"), p_subsystem_filtered, 
           width = 12, height = 8, dpi = 300, bg = "white")
    
    # Export filtered subsystem analysis to CSV
    write.csv(combined_subsystem_df_filtered, file.path(output_dir, "9_Common_Metabolic_Reactions_Subsystem_Counts.csv"), row.names = FALSE)
    
    # Create detailed tables for each direction (filtered)
    if (length(inc_all_three_rxns_filtered) > 0) {
      inc_detailed_filtered <- inc_common_data_filtered[, detail_cols]
      write.csv(inc_detailed_filtered, file.path(output_dir, "9_Common_Increased_Metabolic_Reactions_Detailed.csv"), row.names = FALSE)
    }
    
    if (length(dec_all_three_rxns_filtered) > 0) {
      dec_detailed_filtered <- dec_common_data_filtered[, detail_cols]
      write.csv(dec_detailed_filtered, file.path(output_dir, "9_Common_Decreased_Metabolic_Reactions_Detailed.csv"), row.names = FALSE)
    }
    
    cat("\nFiltered subsystem analysis files created:\n")
    cat("- Common_Metabolic_Reactions_Subsystem_Analysis.png\n")
    cat("- Common_Metabolic_Reactions_Subsystem_Counts.csv\n")
    if (length(inc_all_three_rxns_filtered) > 0) cat("- Common_Increased_Metabolic_Reactions_Detailed.csv\n")
    if (length(dec_all_three_rxns_filtered) > 0) cat("- Common_Decreased_Metabolic_Reactions_Detailed.csv\n")
  }
}

cat("\nFiltered subsystem analysis for common metabolic reactions complete!\n")


# Install if needed: install.packages("UpSetR")
# UpSet plots for better intersection visualization
library(UpSetR)

# Create UpSet plot for increased reactions
increased_reactions_all <- unique(c(fba_inc, pfba_inc, sampling_inc))

upset_matrix_inc <- data.frame(
  FBA = as.numeric(increased_reactions_all %in% fba_inc),
  pFBA = as.numeric(increased_reactions_all %in% pfba_inc),
  FVA_sampling = as.numeric(increased_reactions_all %in% sampling_inc)
)

png(file.path(output_dir, "10_UpSet_Increased_Reactions.png"), width = 10, height = 6, res = 300, units = "in")
upset(upset_matrix_inc, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

svg(file.path(output_dir, "10_UpSet_Increased_Reactions.svg"), width = 10, height = 6)

upset(upset_matrix_inc, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

# Create UpSet plot for decreased reactions
decreased_reactions_all <- unique(c(fba_dec, pfba_dec, sampling_dec))

upset_matrix_dec <- data.frame(
  FBA = as.numeric(decreased_reactions_all %in% fba_dec),
  pFBA = as.numeric(decreased_reactions_all %in% pfba_dec),
  FVA_sampling = as.numeric(decreased_reactions_all %in% sampling_dec)
)

png(file.path(output_dir, "10_UpSet_Decreased_Reactions.png"), width = 10, height = 6, res = 300, units = "in")
upset(upset_matrix_dec, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

svg(file.path(output_dir, "10_UpSet_Decreased_Reactions.svg"), width = 10, height = 6)
upset(upset_matrix_dec, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

# Also create UpSet plots for filtered data (metabolic reactions only)
increased_reactions_filtered <- unique(c(fba_inc_filtered, pfba_inc_filtered, sampling_inc_filtered))

upset_matrix_inc_filtered <- data.frame(
  FBA = as.numeric(increased_reactions_filtered %in% fba_inc_filtered),
  pFBA = as.numeric(increased_reactions_filtered %in% pfba_inc_filtered),
  FVA_sampling = as.numeric(increased_reactions_filtered %in% sampling_inc_filtered)
)

png(file.path(output_dir, "10_UpSet_Increased_Reactions_Filtered.png"), width = 10, height = 6, res = 300, units = "in")
upset(upset_matrix_inc_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

pdf(file.path(output_dir, "10_UpSet_Increased_Reactions_Filtered.pdf"), width = 8, height = 6)
upset(upset_matrix_inc_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      line.size = 1)
dev.off()

svg(file.path(output_dir, "10_UpSet_Increased_Reactions_Filtered.svg"), width = 10, height = 6)
upset(upset_matrix_inc_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()


decreased_reactions_filtered <- unique(c(fba_dec_filtered, pfba_dec_filtered, sampling_dec_filtered))

upset_matrix_dec_filtered <- data.frame(
  FBA = as.numeric(decreased_reactions_filtered %in% fba_dec_filtered),
  pFBA = as.numeric(decreased_reactions_filtered %in% pfba_dec_filtered),
  FVA_sampling = as.numeric(decreased_reactions_filtered %in% sampling_dec_filtered)
)

png(file.path(output_dir, "10_UpSet_Decreased_Reactions_Filtered.png"), width = 8, height = 6, res = 300, units = "in")
upset(upset_matrix_dec_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      scale.sets = "identity",
     # set_size.scale_max = 10,
      
      line.size = 1)
dev.off()

pdf(file.path(output_dir, "10_UpSet_Decreased_Reactions_Filtered_pdf.pdf"),
    width = 8, height = 6)   # no res needed for vector formats
upset(upset_matrix_dec_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      scale.sets = "identity",
      line.size = 1)
dev.off()



pdf(file.path(output_dir, "10_UpSet_Decreased_Reactions_Filtered_grid.pdf"),
    width = 8, height = 6)
grid::grid.newpage()
upset(upset_matrix_dec_filtered, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      scale.sets = "identity",
      line.size = 1)
dev.off()


# Panel D Options for Both Increased and Decreased Reactions

# Function to save plots in both PDF and SVG
save_plot_both <- function(plot_obj, filename, width = 6, height = 6) {
  # PDF
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot_obj, 
         width = width, height = height, bg = "white")
  # SVG
  ggsave(file.path(output_dir, paste0(filename, ".svg")), plot_obj, 
         width = width, height = height, bg = "white")
  
  ggsave(file.path(output_dir, paste0(filename, ".png")), plot_obj, 
         width = width, height = height, bg = "white")
}

# Function to save non-ggplot objects (tables, networks)
save_other_both <- function(plot_code, filename, width = 4, height = 4) {
  # PDF
  pdf(file.path(output_dir, paste0(filename, ".pdf")), width = width, height = height)
  eval(plot_code)
  dev.off()
  # SVG
  svg(file.path(output_dir, paste0(filename, ".svg")), width = width, height = height)
  eval(plot_code)
  dev.off()
}

###############################################################################
# DECREASED REACTIONS ANALYSIS
###############################################################################

if (exists("dec_subsystem_df_filtered") && exists("dec_common_data_filtered")) {
  
  # Get top 3 subsystems for decreased
  top_subsystems_dec <- head(dec_subsystem_df_filtered[order(dec_subsystem_df_filtered$Count, decreasing = TRUE), ], 3)
  
  # Extract genes for top subsystems
  top_subsystem_genes_dec <- list()
  for (subsys in top_subsystems_dec$Subsystem) {
    subsys_reactions <- dec_common_data_filtered[dec_common_data_filtered$Subsystem == subsys, ]
    all_genes <- c()
    for (i in 1:nrow(subsys_reactions)) {
      gene_string <- subsys_reactions$Associated_Genes[i]
      if (!is.na(gene_string) && gene_string != "No genes associated") {
        genes <- strsplit(gene_string, "; ")[[1]]
        all_genes <- c(all_genes, genes)
      }
    }
    gene_counts <- table(all_genes)
    gene_df <- data.frame(
      Gene = names(gene_counts),
      Count = as.numeric(gene_counts),
      Subsystem = subsys
    )
    top_subsystem_genes_dec[[subsys]] <- head(gene_df[order(gene_df$Count, decreasing = TRUE), ], 5)
  }
  
  combined_gene_data_dec <- do.call(rbind, top_subsystem_genes_dec)
  
  # OPTION 1: Network for decreased
  library(igraph)
  library(ggraph)
  
  edges_dec <- combined_gene_data_dec[, c("Subsystem", "Gene", "Count")]
  colnames(edges_dec) <- c("from", "to", "weight")
  
  subsystem_nodes_dec <- data.frame(
    name = top_subsystems_dec$Subsystem,
    type = "subsystem",
    size = top_subsystems_dec$Count
  )
  gene_nodes_dec <- data.frame(
    name = unique(combined_gene_data_dec$Gene),
    type = "gene", 
    size = 3
  )
  nodes_dec <- rbind(subsystem_nodes_dec, gene_nodes_dec)
  
  g_dec <- graph_from_data_frame(edges_dec, vertices = nodes_dec, directed = TRUE)
  
  p_network_dec <- ggraph(g_dec, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray60") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("subsystem" = "#3498db", "gene" = "#A6CEE3"),
                       breaks = c("subsystem", "gene"), 
                       name = 'Level') +  # Both blue tones
#    scale_size_continuous(range = c(3, 8)) +
    scale_size_continuous(range = c(3, 8), breaks = c(9, 15, 21), name = "# of rxns")+
 #   scale_edge_width_continuous(name = "# Reactions", range = c(0.5, 3)) +
    scale_edge_width_continuous(name = "# of genes", range = c(0.5, 3), 
                                breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
    theme_void() +
  #  labs(title = "Key genes in decreased subsystems") +
 #   theme(legend.position = "right")
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),     # shrink legend keys
      legend.text = element_text(size = 7),  # smaller legend text
      legend.title = element_text(size = 8)  # smaller legend title
    )
  
  save_plot_both(p_network_dec, "PanelD_Network_Decreased")
  
  p_network_dec_v2 <- ggraph(g_dec, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray60") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("subsystem" = "#3498db", "gene" = "#A6CEE3")) +  # Both blue tones
    #    scale_size_continuous(range = c(3, 8)) +
    scale_size_continuous(range = c(3, 8), breaks = c(9, 15, 21), name = "# of rxns")+
    #   scale_edge_width_continuous(name = "# Reactions", range = c(0.5, 3)) +
    scale_edge_width_continuous(name = "# of genes", range = c(0.5, 3), 
                                breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
    theme_void() +
    #  labs(title = "Key genes in decreased subsystems") +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("PanelD_Network_Decreased", "_v2.svg")), p_network_dec_v2, 
         width = 6, height = 6, bg = "white")
  
  
  # OPTION 2: Table for decreased
  table_data_dec <- data.frame()
  for (subsys in top_subsystems_dec$Subsystem) {
    genes <- top_subsystem_genes_dec[[subsys]]
    gene_list <- paste(head(genes$Gene, 5), collapse = ", ")
    table_data_dec <- rbind(table_data_dec, data.frame(
      Subsystem = subsys,
      Reactions = top_subsystems_dec$Count[top_subsystems_dec$Subsystem == subsys],
      Key_Genes = gene_list
    ))
  }
  
  table_data_dec$Subsystem <- gsub("(.{30})", "\\1\n", table_data_dec$Subsystem)
  
  library(gridExtra)
  library(grid)
  
  p_table_dec <- tableGrob(table_data_dec, rows = NULL, theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 10)),
    colhead = list(fg_params = list(fontsize = 12, fontface = "bold"))
  ))
  
  save_other_both(quote(grid.draw(p_table_dec)), "PanelD_Table_Decreased", width = 10, height = 4)
  
  # OPTION 3: Focused bar for decreased
  gene_subsystem_count_dec <- combined_gene_data_dec %>%
    group_by(Gene) %>%
    summarise(
      Subsystem_count = n_distinct(Subsystem),
      Total_reactions = sum(Count)
    ) %>%
    filter(Subsystem_count > 1 | Total_reactions > 2) %>%
    arrange(desc(Total_reactions))
  
  if (nrow(gene_subsystem_count_dec) > 0) {
    p_focused_bar_dec <- ggplot(head(gene_subsystem_count_dec, 10), 
                                aes(x = reorder(Gene, Total_reactions), y = Total_reactions)) +
      geom_bar(stat = "identity", fill = "#3498db") +
      coord_flip() +
      labs(title = "Multi-pathway genes in decreased reactions",
           x = "Gene", 
           y = "Total reactions") +
      theme_minimal()
    
    save_plot_both(p_focused_bar_dec, "PanelD_FocusedBar_Decreased")
  }
}

###############################################################################
# INCREASED REACTIONS ANALYSIS
###############################################################################

if (exists("inc_subsystem_df_filtered") && exists("inc_common_data_filtered")) {
  
  # Get top 3 subsystems for increased
  top_subsystems_inc <- head(inc_subsystem_df_filtered[order(inc_subsystem_df_filtered$Count, decreasing = TRUE), ], 3)
  
  # Extract genes for top subsystems
  top_subsystem_genes_inc <- list()
  for (subsys in top_subsystems_inc$Subsystem) {
    subsys_reactions <- inc_common_data_filtered[inc_common_data_filtered$Subsystem == subsys, ]
    all_genes <- c()
    for (i in 1:nrow(subsys_reactions)) {
      gene_string <- subsys_reactions$Associated_Genes[i]
      if (!is.na(gene_string) && gene_string != "No genes associated") {
        genes <- strsplit(gene_string, "; ")[[1]]
        all_genes <- c(all_genes, genes)
      }
    }
    gene_counts <- table(all_genes)
    gene_df <- data.frame(
      Gene = names(gene_counts),
      Count = as.numeric(gene_counts),
      Subsystem = subsys
    )
    top_subsystem_genes_inc[[subsys]] <- head(gene_df[order(gene_df$Count, decreasing = TRUE), ], 5)
  }
  
  combined_gene_data_inc <- do.call(rbind, top_subsystem_genes_inc)
  
  
  
  # OPTION 1: Network for increased
  edges_inc <- combined_gene_data_inc[, c("Subsystem", "Gene", "Count")]
  colnames(edges_inc) <- c("from", "to", "weight")
  
  subsystem_nodes_inc <- data.frame(
    name = top_subsystems_inc$Subsystem,
    type = "subsystem",
    size = top_subsystems_inc$Count
  )
  gene_nodes_inc <- data.frame(
    name = unique(combined_gene_data_inc$Gene),
    type = "gene", 
    size = 3
  )
  nodes_inc <- rbind(subsystem_nodes_inc, gene_nodes_inc)
  
  g_inc <- graph_from_data_frame(edges_inc, vertices = nodes_inc, directed = TRUE)
  
  p_network_inc <- ggraph(g_inc, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray60") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("subsystem" = "#e74c3c", "gene" = "#f5b7b1"),
                       breaks = c("subsystem", "gene"), 
                       name = 'Level') +  # Both blue tones) +  # Both red tones
    scale_size_continuous(range = c(3, 8), breaks = c(3, 5, 7),name="# of rxns") +
   # scale_edge_width_continuous(name = "# Reactions", range = c(0.5, 3)) +
    scale_edge_width_continuous(name = "# of genes", range = c(0.5, 3), 
                                breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
    theme_void() +
 #   labs(title = "Key genes in increased subsystems") +
  #  theme(legend.position = "right")
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),     # shrink legend keys
      legend.text = element_text(size = 7),  # smaller legend text
      legend.title = element_text(size = 8)  # smaller legend title
    )
  
  save_plot_both(p_network_inc, "PanelD_Network_Increased")
  
    # SVG - Separately, use the legend
  p_network_inc_v2 <- ggraph(g_inc, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray60") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("subsystem" = "#e74c3c", "gene" = "#f5b7b1")) +  # Both red tones
    scale_size_continuous(range = c(3, 8), breaks = c(3, 5, 7)) +
    # scale_edge_width_continuous(name = "# Reactions", range = c(0.5, 3)) +
    scale_edge_width_continuous(name = "# Reactions", range = c(0.5, 3), 
                                breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
    theme_void() +
    #   labs(title = "Key genes in increased subsystems") +
    theme(legend.position = "right")
  

    ggsave(file.path(output_dir, paste0("PanelD_Network_Increased", "_v2.svg")), p_network_inc_v2, 
           width = 6, height = 6, bg = "white")
  
  # OPTION 2: Table for increased
  table_data_inc <- data.frame()
  for (subsys in top_subsystems_inc$Subsystem) {
    genes <- top_subsystem_genes_inc[[subsys]]
    gene_list <- paste(head(genes$Gene, 5), collapse = ", ")
    table_data_inc <- rbind(table_data_inc, data.frame(
      Subsystem = subsys,
      Reactions = top_subsystems_inc$Count[top_subsystems_inc$Subsystem == subsys],
      Key_Genes = gene_list
    ))
  }
  
  table_data_inc$Subsystem <- gsub("(.{30})", "\\1\n", table_data_inc$Subsystem)
  
  p_table_inc <- tableGrob(table_data_inc, rows = NULL, theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 10)),
    colhead = list(fg_params = list(fontsize = 12, fontface = "bold"))
  ))
  
  save_other_both(quote(grid.draw(p_table_inc)), "PanelD_Table_Increased", width = 10, height = 4)
  
  # OPTION 3: Focused bar for increased
  gene_subsystem_count_inc <- combined_gene_data_inc %>%
    group_by(Gene) %>%
    summarise(
      Subsystem_count = n_distinct(Subsystem),
      Total_reactions = sum(Count)
    ) %>%
    filter(Subsystem_count > 1 | Total_reactions > 2) %>%
    arrange(desc(Total_reactions))
  
  if (nrow(gene_subsystem_count_inc) > 0) {
    p_focused_bar_inc <- ggplot(head(gene_subsystem_count_inc, 10), 
                                aes(x = reorder(Gene, Total_reactions), y = Total_reactions)) +
      geom_bar(stat = "identity", fill = "#e74c3c") +
      coord_flip() +
      labs(title = "Multi-pathway genes in increased reactions",
           x = "Gene", 
           y = "Total reactions") +
      theme_minimal()
    
    save_plot_both(p_focused_bar_inc, "PanelD_FocusedBar_Increased")
  }
}

# Export data for both
if (exists("combined_gene_data_dec")) {
  write.csv(table_data_dec, file.path(output_dir, "PanelD_Summary_Table_Decreased.csv"), row.names = FALSE)
  write.csv(combined_gene_data_dec, file.path(output_dir, "PanelD_All_Gene_Data_Decreased.csv"), row.names = FALSE)
}

if (exists("combined_gene_data_inc")) {
  write.csv(table_data_inc, file.path(output_dir, "PanelD_Summary_Table_Increased.csv"), row.names = FALSE)
  write.csv(combined_gene_data_inc, file.path(output_dir, "PanelD_All_Gene_Data_Increased.csv"), row.names = FALSE)
}

cat("Panel D options created for both increased and decreased (PDF + SVG):\n")
cat("DECREASED:\n")
cat("- PanelD_Network_Decreased.pdf/.svg\n")
cat("- PanelD_Table_Decreased.pdf/.svg\n")
if (exists("p_focused_bar_dec")) cat("- PanelD_FocusedBar_Decreased.pdf/.svg\n")
cat("INCREASED:\n")
cat("- PanelD_Network_Increased.pdf/.svg\n")
cat("- PanelD_Table_Increased.pdf/.svg\n")
if (exists("p_focused_bar_inc")) cat("- PanelD_FocusedBar_Increased.pdf/.svg\n")
cat("CSV exports created for both directions\n")

# Simple donut plot for decreased subsystems
if (exists("dec_subsystem_df") && nrow(dec_subsystem_df) > 0) {
  # Calculate percentages
  dec_subsystem_df$Percentage <- round(dec_subsystem_df$Count / sum(dec_subsystem_df$Count) * 100, 1)
  
  # Calculate positions for labels
  dec_subsystem_df$ypos <- cumsum(dec_subsystem_df$Count) - 0.5 * dec_subsystem_df$Count
  
  # Create donut plot
  p_dec_donut <- ggplot(dec_subsystem_df, aes(x = 2, y = Count, fill = Subsystem)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Creates the donut hole
    scale_fill_viridis_d(option = "mako", alpha = 0.8) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray50"),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold")
    ) +
    labs(
      title = "Subsystem Distribution: Decreased Reactions",
      subtitle = paste0("Total: ", sum(dec_subsystem_df$Count), " reactions"),
      fill = "Subsystem"
    ) +
    geom_text(
      aes(y = ypos, label = ifelse(Percentage >= 3, paste0(Count, "\n(", Percentage, "%)"), "")), 
      x = 2, size = 3, fontface = "bold", color = "white"
    ) +
    annotate("text", x = 0.5, y = 0, 
             label = paste0(sum(dec_subsystem_df$Count), "\nReactions"), 
             size = 5, fontface = "bold", color = "gray40")
  
  # Save the plot
  ggsave(file.path(output_dir, "Decreased_Subsystem_Donut.png"), p_dec_donut, 
         width = 10, height = 8, dpi = 300, bg = "white")
  
  cat("Decreased subsystem donut plot created: Decreased_Subsystem_Donut.png\n")
} else {
  cat("No decreased subsystem data available for donut plot.\n")
}
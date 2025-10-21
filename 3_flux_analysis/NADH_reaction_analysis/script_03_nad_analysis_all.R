# script_03_nad_analysis_all
#All NAD-Associated Reactions Analysis
# Extract all reactions containing NAD-related metabolites across all compartments
# Update such that we include all the nad reactions, including abs(diff) < 1, to get subtle changes

# Load required libraries
library(readxl)
library(VennDiagram)
library(ggvenn)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

rm(list = ls())

# Set working directory and file paths
current_dir <- "E:/Projects/revision/Code_2_upload/NAD_reaction_analysis"
base_dir <- file.path(current_dir,"files")
output_dir <- file.path(current_dir, "03_nad_analysis_all")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the main analysis data
excel_file <- file.path(base_dir, "venn_diagram_data_th0_forNAD.xlsx")

cat("Loading main data from:", excel_file, "\n")

# Check if file exists
if (!file.exists(excel_file)) {
  stop("Main Excel file not found. Please check the path.")
}

# Load the binary matrix data
binary_data <- read_excel(excel_file, sheet = "R_Binary_Matrix")
cat("Loaded main data with", nrow(binary_data), "reactions\n")



# Function to check if a reaction contains any NAD metabolites
# Function to check if a reaction contains any NAD metabolites
contains_nad_metabolite <- function(metabolite_string) {
  if (is.na(metabolite_string) || metabolite_string == "") {
    return(FALSE)
  }
  
  # Search for NAD+ and NADH in human-readable format
  # Focus only on NAD+/NADH (not NADP+/NADPH)
  nad_patterns <- c("NAD\\+ \\[", "NADH \\[")
  
  # Check if any NAD pattern appears in the metabolite string
  any(sapply(nad_patterns, function(pattern) {
    grepl(pattern, metabolite_string, ignore.case = FALSE)
  }))
}
#

# Extract reactions containing NAD metabolites
if ("Associated_Metabolites" %in% colnames(binary_data)) {
  nad_reactions_mask <- sapply(binary_data$Associated_Metabolites, contains_nad_metabolite)
} else {
  # If no Associated_Metabolites column, look for any mention of NAD in reaction names/formulas
  nad_reactions_mask <- grepl("NAD|nad", binary_data$Reaction_Formula, ignore.case = TRUE)
  cat("Warning: Using reaction formula to identify NAD reactions (no metabolite IDs available)\n")
}


# Filter for NAD-associated reactions
nad_all_data <- binary_data[nad_reactions_mask, ]
cat("Found", nrow(nad_all_data), "reactions containing NAD-related metabolites\n")

# Remove exchange/demand reactions
nad_all_data_filtered <- nad_all_data[nad_all_data$Subsystem != "Exchange/demand reactions", ]
cat("After removing exchange/demand:", nrow(nad_all_data_filtered), "NAD reactions remain\n")

# Analyze compartment distribution
if ("Associated_Metabolites" %in% colnames(nad_all_data_filtered)) {
  compartment_analysis <- data.frame(
    Compartment = c("Cytosol", "Mitochondria", "Extracellular", "Peroxisome",
                    "Lysosome", "Endoplasmic reticulum", "Golgi", "Nucleus", "Intermediate"),
    Code = c("c", "m", "e", "p", "l", "r", "g", "n", "i"),
    Count = 0
  )
  
  for (i in 1:nrow(compartment_analysis)) {
    comp_code <- compartment_analysis$Code[i]
    comp_pattern <- paste0("NAD.*\\[", compartment_analysis$Compartment[i], "\\]")
    compartment_analysis$Count[i] <- sum(grepl(comp_pattern, nad_all_data_filtered$Associated_Metabolites))
  }
  
  cat("\nCompartment distribution of NAD reactions:\n")
  print(compartment_analysis[compartment_analysis$Count > 0, ])
  
  # Create compartment visualization
  comp_plot_data <- compartment_analysis[compartment_analysis$Count > 0, ]
  if (nrow(comp_plot_data) > 0) {
    p_compartment <- ggplot(comp_plot_data, aes(x = reorder(Compartment, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db", alpha = 0.8) +
      coord_flip() +
      labs(title = "NAD-Associated Reactions by Compartment",
           x = "Compartment",
           y = "Number of Reactions") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
      geom_text(aes(label = Count), hjust = -0.3)
    
    ggsave(file.path(output_dir, "NAD_Reactions_Compartment_Distribution.png"), p_compartment, 
           width = 10, height = 6, dpi = 300, bg = "white")
  }
  
  # Export compartment analysis
  write.csv(compartment_analysis, file.path(output_dir, "NAD_Compartment_Analysis.csv"), row.names = FALSE)
}

# Create pie chart for original compartment analysis (all NAD reactions)
comp_plot_data_original <- compartment_analysis[compartment_analysis$Count > 0, ]

if (nrow(comp_plot_data_original) > 0) {
  # Calculate percentages
  comp_plot_data_original$percentage <- comp_plot_data_original$Count / sum(comp_plot_data_original$Count) * 100
  comp_plot_data_original$label <- paste0(comp_plot_data_original$Compartment, "\n", 
                                          comp_plot_data_original$Count, " (", 
                                          round(comp_plot_data_original$percentage, 1), "%)")
  
  # Create pie chart
  library(ggrepel)
  
  p_compartment_pie_all <- ggplot(comp_plot_data_original, aes(x = "", y = Count, fill = Compartment)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text_repel(aes(label = label), 
                    position = position_stack(vjust = 0.5),
                    size = 3, 
                    force = 1, # push labels further apart
                #    box.padding = 0.5, # add space around text
                    max.overlaps = 10) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    theme_void() +
    theme(legend.position = "right") 
    #labs(title = "All NAD Reactions by Compartment")
  
  # Save pie chart in multiple formats
  ggsave(file.path(output_dir, "NAD_Compartment_All_Pie.png"), p_compartment_pie_all, 
         width = 5, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "NAD_Compartment_All_Pie.pdf"), p_compartment_pie_all, 
         width = 5, height = 5, bg = "white")
  ggsave(file.path(output_dir, "NAD_Compartment_All_Pie.svg"), p_compartment_pie_all, 
         width = 5, height = 5, bg = "white")
  
  cat("Pie chart saved: NAD_Compartment_All_Pie (png/pdf/svg)\n")
}

# Export the original compartment analysis
write.csv(compartment_analysis, file.path(output_dir, "NAD_Compartment_Analysis_All.csv"), row.names = FALSE)


# Create lists for all NAD increased reactions
fba_nad_all_inc <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$FBA_Inc == 1]
pfba_nad_all_inc <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$pFBA_Inc == 1]
sampling_nad_all_inc <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$Sampling_Inc == 1]

# Create lists for all NAD decreased reactions
fba_nad_all_dec <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$FBA_Dec == 1]
pfba_nad_all_dec <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$pFBA_Dec == 1]
sampling_nad_all_dec <- nad_all_data_filtered$Reaction_ID[nad_all_data_filtered$Sampling_Dec == 1]

# Remove any NA values
fba_nad_all_inc <- fba_nad_all_inc[!is.na(fba_nad_all_inc)]
pfba_nad_all_inc <- pfba_nad_all_inc[!is.na(pfba_nad_all_inc)]
sampling_nad_all_inc <- sampling_nad_all_inc[!is.na(sampling_nad_all_inc)]
fba_nad_all_dec <- fba_nad_all_dec[!is.na(fba_nad_all_dec)]
pfba_nad_all_dec <- pfba_nad_all_dec[!is.na(pfba_nad_all_dec)]
sampling_nad_all_dec <- sampling_nad_all_dec[!is.na(sampling_nad_all_dec)]

cat("\n=== All NAD-Associated Reactions Analysis Summary ===\n")
cat("NAD Increased reactions - FBA:", length(fba_nad_all_inc), "pFBA:", length(pfba_nad_all_inc), "Sampling:", length(sampling_nad_all_inc), "\n")
cat("NAD Decreased reactions - FBA:", length(fba_nad_all_dec), "pFBA:", length(pfba_nad_all_dec), "Sampling:", length(sampling_nad_all_dec), "\n")


########## Extract the compartment analysis for those that have at least reactions having diff. greater 1 from any of the analysis
# Filter for reactions that have at least one change (any 1 in the binary columns)
has_change <- rowSums(nad_all_data_filtered[, c("FBA_Inc", "pFBA_Inc", "Sampling_Inc", 
                                                "FBA_Dec", "pFBA_Dec", "Sampling_Dec")]) > 0

nad_all_data_filtered_with_changes <- nad_all_data_filtered[has_change, ]

cat("Reactions with at least one change:", nrow(nad_all_data_filtered_with_changes), 
    "out of", nrow(nad_all_data_filtered), "total NAD reactions\n")

# Create compartment analysis for filtered data only
compartment_analysis_filtered <- data.frame(
  Compartment = c("Cytosol", "Mitochondria", "Extracellular", "Peroxisome",
                  "Lysosome", "Endoplasmic reticulum", "Golgi", "Nucleus", "Intermediate"),
  Code = c("c", "m", "e", "p", "l", "r", "g", "n", "i"),
  Count = 0
)

for (i in 1:nrow(compartment_analysis_filtered)) {
  comp_pattern <- paste0("NAD.*\\[", compartment_analysis_filtered$Compartment[i], "\\]")
  compartment_analysis_filtered$Count[i] <- sum(grepl(comp_pattern, nad_all_data_filtered_with_changes$Associated_Metabolites))
}

cat("\nCompartment distribution of NAD reactions WITH CHANGES:\n")
print(compartment_analysis_filtered[compartment_analysis_filtered$Count > 0, ])

# Update your filtered data for further analysis
nad_all_data_filtered <- nad_all_data_filtered_with_changes


# Create pie chart for compartment distribution of NAD reactions with changes
comp_plot_data_filtered <- compartment_analysis_filtered[compartment_analysis_filtered$Count > 0, ]

if (nrow(comp_plot_data_filtered) > 0) {
  # Calculate percentages
  comp_plot_data_filtered$percentage <- comp_plot_data_filtered$Count / sum(comp_plot_data_filtered$Count) * 100
  comp_plot_data_filtered$label <- paste0(comp_plot_data_filtered$Compartment, "\n", 
                                          comp_plot_data_filtered$Count, " (", 
                                          round(comp_plot_data_filtered$percentage, 1), "%)")
  
  # Create pie chart
  p_compartment_pie <- ggplot(comp_plot_data_filtered, aes(x = "", y = Count, fill = Compartment)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    theme_void() +
    theme(legend.position = "right")
   # labs(title = "NAD Reactions with Changes by Compartment")
  
  # Save pie chart in multiple formats
  ggsave(file.path(output_dir, "NAD_Compartment_Changes_Pie.png"), p_compartment_pie, 
         width = 6, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "NAD_Compartment_Changes_Pie.pdf"), p_compartment_pie, 
         width = 6, height = 6, bg = "white")
  ggsave(file.path(output_dir, "NAD_Compartment_Changes_Pie.svg"), p_compartment_pie, 
         width = 6, height = 6, bg = "white")
  
  cat("Pie chart saved: NAD_Compartment_Changes_Pie (png/pdf/svg)\n")
}

# Export the filtered compartment analysis
write.csv(compartment_analysis_filtered, file.path(output_dir, "NAD_Compartment_Analysis_Filtered.csv"), row.names = FALSE)




# Fix binary indicators - remove the abs(1) threshold
# Any non-zero difference should be counted as increased/decreased

library(openxlsx)

# Create a copy for modification
nad_all_data_filtered_v2 <- nad_all_data_filtered

# Fix FBA indicators
nad_all_data_filtered_v2$FBA_Inc <- ifelse(nad_all_data_filtered_v2$FBA_Difference > 0, 1, 0)
nad_all_data_filtered_v2$FBA_Dec <- ifelse(nad_all_data_filtered_v2$FBA_Difference < 0, 1, 0)

# Fix pFBA indicators  
nad_all_data_filtered_v2$pFBA_Inc <- ifelse(nad_all_data_filtered_v2$pFBA_Difference > 0, 1, 0)
nad_all_data_filtered_v2$pFBA_Dec <- ifelse(nad_all_data_filtered_v2$pFBA_Difference < 0, 1, 0)

# Fix Sampling indicators
nad_all_data_filtered_v2$Sampling_Inc <- ifelse(nad_all_data_filtered_v2$Sampling_Difference > 0, 1, 0)
nad_all_data_filtered_v2$Sampling_Dec <- ifelse(nad_all_data_filtered_v2$Sampling_Difference < 0, 1, 0)

# Check the changes
cat("=== BEFORE (original) vs AFTER (v2) Comparison ===\n")
cat("FBA_Inc: original =", sum(nad_all_data_filtered$FBA_Inc), ", v2 =", sum(nad_all_data_filtered_v2$FBA_Inc), "\n")
cat("FBA_Dec: original =", sum(nad_all_data_filtered$FBA_Dec), ", v2 =", sum(nad_all_data_filtered_v2$FBA_Dec), "\n")
cat("pFBA_Inc: original =", sum(nad_all_data_filtered$pFBA_Inc), ", v2 =", sum(nad_all_data_filtered_v2$pFBA_Inc), "\n")
cat("pFBA_Dec: original =", sum(nad_all_data_filtered$pFBA_Dec), ", v2 =", sum(nad_all_data_filtered_v2$pFBA_Dec), "\n")
cat("Sampling_Inc: original =", sum(nad_all_data_filtered$Sampling_Inc), ", v2 =", sum(nad_all_data_filtered_v2$Sampling_Inc), "\n")
cat("Sampling_Dec: original =", sum(nad_all_data_filtered$Sampling_Dec), ", v2 =", sum(nad_all_data_filtered_v2$Sampling_Dec), "\n")

# Show examples of changes
cat("\n=== Examples of Fixed Reactions ===\n")
example_rows <- which(nad_all_data_filtered_v2$Sampling_Inc != nad_all_data_filtered$Sampling_Inc | 
                        nad_all_data_filtered_v2$Sampling_Dec != nad_all_data_filtered$Sampling_Dec)[1:5]

for (i in example_rows) {
  if (!is.na(i)) {
    cat(sprintf("Row %d: Reaction %s, Sampling_Diff = %.3f\n", 
                i, nad_all_data_filtered_v2$Reaction_ID[i], nad_all_data_filtered_v2$Sampling_Difference[i]))
    cat(sprintf("  Original: Inc=%d, Dec=%d\n", nad_all_data_filtered$Sampling_Inc[i], nad_all_data_filtered$Sampling_Dec[i]))
    cat(sprintf("  Fixed:    Inc=%d, Dec=%d\n", nad_all_data_filtered_v2$Sampling_Inc[i], nad_all_data_filtered_v2$Sampling_Dec[i]))
  }
}

# Save the updated data to new Excel file
excel_file_v2 <- file.path(base_dir, "venn_diagram_data_v2_noDiff1.xlsx")

# Create workbook and add the updated data
wb <- createWorkbook()
addWorksheet(wb, "R_Binary_Matrix")
writeData(wb, "R_Binary_Matrix", nad_all_data_filtered_v2)

# Save the workbook
saveWorkbook(wb, excel_file_v2, overwrite = TRUE)

cat(sprintf("\n=== RESULTS SAVED ===\n"))
cat("Updated data saved to:", excel_file_v2, "\n")
cat("Variable 'nad_all_data_filtered_v2' contains the corrected data\n")

# Now create the corrected reaction lists for analysis
fba_nad_all_inc_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$FBA_Inc == 1]
pfba_nad_all_inc_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$pFBA_Inc == 1]
sampling_nad_all_inc_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$Sampling_Inc == 1]

fba_nad_all_dec_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$FBA_Dec == 1]
pfba_nad_all_dec_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$pFBA_Dec == 1]
sampling_nad_all_dec_v2 <- nad_all_data_filtered_v2$Reaction_ID[nad_all_data_filtered_v2$Sampling_Dec == 1]

# Remove any NA values
fba_nad_all_inc_v2 <- fba_nad_all_inc_v2[!is.na(fba_nad_all_inc_v2)]
pfba_nad_all_inc_v2 <- pfba_nad_all_inc_v2[!is.na(pfba_nad_all_inc_v2)]
sampling_nad_all_inc_v2 <- sampling_nad_all_inc_v2[!is.na(sampling_nad_all_inc_v2)]
fba_nad_all_dec_v2 <- fba_nad_all_dec_v2[!is.na(fba_nad_all_dec_v2)]
pfba_nad_all_dec_v2 <- pfba_nad_all_dec_v2[!is.na(pfba_nad_all_dec_v2)]
sampling_nad_all_dec_v2 <- sampling_nad_all_dec_v2[!is.na(sampling_nad_all_dec_v2)]

cat("\n=== V2 NAD Reactions Summary ===\n")
cat("NAD Increased reactions - FBA:", length(fba_nad_all_inc_v2), "pFBA:", length(pfba_nad_all_inc_v2), "Sampling:", length(sampling_nad_all_inc_v2), "\n")
cat("NAD Decreased reactions - FBA:", length(fba_nad_all_dec_v2), "pFBA:", length(pfba_nad_all_dec_v2), "Sampling:", length(sampling_nad_all_dec_v2), "\n")

# Calculate common reactions with new data
nad_all_common_inc_v2 <- intersect(intersect(fba_nad_all_inc_v2, pfba_nad_all_inc_v2), sampling_nad_all_inc_v2)
nad_all_common_dec_v2 <- intersect(intersect(fba_nad_all_dec_v2, pfba_nad_all_dec_v2), sampling_nad_all_dec_v2)

cat("Common increased (all 3 methods):", length(nad_all_common_inc_v2), "\n")
cat("Common decreased (all 3 methods):", length(nad_all_common_dec_v2), "\n")

# Simple save to output_dir
write.csv(nad_all_data_filtered_v2, file.path(output_dir, "nad_all_data_filtered_v2.csv"), row.names = FALSE)

# Save the Excel file to output_dir too
excel_file_v2 <- file.path(output_dir, "venn_diagram_data_v2_noDiff1.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "R_Binary_Matrix")
writeData(wb, "R_Binary_Matrix", nad_all_data_filtered_v2)
saveWorkbook(wb, excel_file_v2, overwrite = TRUE)

cat("Files saved to:", output_dir, "\n")


# Overwrite the original variables (remove _v2 suffix)
nad_all_data_filtered <- nad_all_data_filtered_v2

# Overwrite the reaction lists too
fba_nad_all_inc <- fba_nad_all_inc_v2
pfba_nad_all_inc <- pfba_nad_all_inc_v2
sampling_nad_all_inc <- sampling_nad_all_inc_v2
fba_nad_all_dec <- fba_nad_all_dec_v2
pfba_nad_all_dec <- pfba_nad_all_dec_v2
sampling_nad_all_dec <- sampling_nad_all_dec_v2

# Overwrite common reactions
nad_all_common_inc <- nad_all_common_inc_v2
nad_all_common_dec <- nad_all_common_dec_v2

# Save with original filenames
write.csv(nad_all_data_filtered, file.path(output_dir, "nad_all_data_filtered.csv"), row.names = FALSE)

excel_file <- file.path(output_dir, "venn_diagram_data_v2_noDiff1.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "R_Binary_Matrix")
writeData(wb, "R_Binary_Matrix", nad_all_data_filtered)
saveWorkbook(wb, excel_file, overwrite = TRUE)

cat("Variables overwritten and saved to:", output_dir, "\n")


# Create Venn diagrams for all NAD INCREASED reactions
nad_all_inc_data <- list(
  FBA = fba_nad_all_inc,
  pFBA = pfba_nad_all_inc,
  Sampling = sampling_nad_all_inc
)

p1_nad_all <- ggvenn(nad_all_inc_data, 
                     fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
                     stroke_size = 1.5,
                     set_name_size = 5,
                     text_size = 4) +
  ggtitle("All NAD-Associated Reactions: Increased in Disease (HSD > NSD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "NAD_All_Venn_Increased_Reactions.png"), p1_nad_all, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "NAD_All_Venn_Increased_Reactions.svg"), p1_nad_all, width = 10, height = 8, dpi = 300, bg = "white")

# Create Venn diagrams for all NAD DECREASED reactions
nad_all_dec_data <- list(
  FBA = fba_nad_all_dec,
  pFBA = pfba_nad_all_dec,
  Sampling = sampling_nad_all_dec
)

p2_nad_all <- ggvenn(nad_all_dec_data,
                     fill_color = c("#4ECDC4", "#FF6B6B", "#45B7D1"),
                     stroke_size = 1.5,
                     set_name_size = 5,
                     text_size = 4) +
  ggtitle("All NAD-Associated Reactions: Decreased in Disease (HSD < NSD)") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(output_dir, "NAD_All_Venn_Decreased_Reactions.png"), p2_nad_all, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "NAD_All_Venn_Decreased_Reactions.svg"), p2_nad_all, width = 10, height = 8, dpi = 300, bg = "white")

# Calculate overlap statistics for all NAD reactions
cat("\n=== All NAD Reactions Venn Diagram Overlap Analysis ===\n")

# For all NAD increased reactions
nad_all_inc_fba_pfba <- length(intersect(fba_nad_all_inc, pfba_nad_all_inc))
nad_all_inc_fba_sampling <- length(intersect(fba_nad_all_inc, sampling_nad_all_inc))
nad_all_inc_pfba_sampling <- length(intersect(pfba_nad_all_inc, sampling_nad_all_inc))
nad_all_inc_all_three <- length(intersect(intersect(fba_nad_all_inc, pfba_nad_all_inc), sampling_nad_all_inc))

cat("\nAll NAD Increased reactions overlaps:\n")
cat("FBA only:", length(fba_nad_all_inc) - nad_all_inc_fba_pfba - nad_all_inc_fba_sampling + nad_all_inc_all_three, "\n")
cat("pFBA only:", length(pfba_nad_all_inc) - nad_all_inc_fba_pfba - nad_all_inc_pfba_sampling + nad_all_inc_all_three, "\n")
cat("Sampling only:", length(sampling_nad_all_inc) - nad_all_inc_fba_sampling - nad_all_inc_pfba_sampling + nad_all_inc_all_three, "\n")
cat("FBA + pFBA only:", nad_all_inc_fba_pfba - nad_all_inc_all_three, "\n")
cat("FBA + Sampling only:", nad_all_inc_fba_sampling - nad_all_inc_all_three, "\n")
cat("pFBA + Sampling only:", nad_all_inc_pfba_sampling - nad_all_inc_all_three, "\n")
cat("All three methods:", nad_all_inc_all_three, "\n")

# For all NAD decreased reactions
nad_all_dec_fba_pfba <- length(intersect(fba_nad_all_dec, pfba_nad_all_dec))
nad_all_dec_fba_sampling <- length(intersect(fba_nad_all_dec, sampling_nad_all_dec))
nad_all_dec_pfba_sampling <- length(intersect(pfba_nad_all_dec, sampling_nad_all_dec))
nad_all_dec_all_three <- length(intersect(intersect(fba_nad_all_dec, pfba_nad_all_dec), sampling_nad_all_dec))

cat("\nAll NAD Decreased reactions overlaps:\n")
cat("FBA only:", length(fba_nad_all_dec) - nad_all_dec_fba_pfba - nad_all_dec_fba_sampling + nad_all_dec_all_three, "\n")
cat("pFBA only:", length(pfba_nad_all_dec) - nad_all_dec_fba_pfba - nad_all_dec_pfba_sampling + nad_all_dec_all_three, "\n")
cat("Sampling only:", length(sampling_nad_all_dec) - nad_all_dec_fba_sampling - nad_all_dec_pfba_sampling + nad_all_dec_all_three, "\n")
cat("FBA + pFBA only:", nad_all_dec_fba_pfba - nad_all_dec_all_three, "\n")
cat("FBA + Sampling only:", nad_all_dec_fba_sampling - nad_all_dec_all_three, "\n")
cat("pFBA + Sampling only:", nad_all_dec_pfba_sampling - nad_all_dec_all_three, "\n")
cat("All three methods:", nad_all_dec_all_three, "\n")

# Export common reactions (all three methods agree)
nad_all_common_inc <- intersect(intersect(fba_nad_all_inc, pfba_nad_all_inc), sampling_nad_all_inc)
nad_all_common_dec <- intersect(intersect(fba_nad_all_dec, pfba_nad_all_dec), sampling_nad_all_dec)

if (length(nad_all_common_inc) > 0) {
  common_inc_data <- nad_all_data_filtered[nad_all_data_filtered$Reaction_ID %in% nad_all_common_inc, ]
  write.csv(common_inc_data, file.path(output_dir, "NAD_All_Common_Increased_Reactions.csv"), row.names = FALSE)
  cat("\nExported", length(nad_all_common_inc), "commonly increased NAD reactions\n")
}

if (length(nad_all_common_dec) > 0) {
  common_dec_data <- nad_all_data_filtered[nad_all_data_filtered$Reaction_ID %in% nad_all_common_dec, ]
  write.csv(common_dec_data, file.path(output_dir, "NAD_All_Common_Decreased_Reactions.csv"), row.names = FALSE)
  cat("Exported", length(nad_all_common_dec), "commonly decreased NAD reactions\n")
}

# Export overlap summary
nad_all_overlap_summary <- data.frame(
  Category = c("NAD_All_Increased_FBA_only", "NAD_All_Increased_pFBA_only", "NAD_All_Increased_Sampling_only",
               "NAD_All_Increased_FBA_pFBA", "NAD_All_Increased_FBA_Sampling", "NAD_All_Increased_pFBA_Sampling", "NAD_All_Increased_All_three",
               "NAD_All_Decreased_FBA_only", "NAD_All_Decreased_pFBA_only", "NAD_All_Decreased_Sampling_only",
               "NAD_All_Decreased_FBA_pFBA", "NAD_All_Decreased_FBA_Sampling", "NAD_All_Decreased_pFBA_Sampling", "NAD_All_Decreased_All_three"),
  Count = c(length(fba_nad_all_inc) - nad_all_inc_fba_pfba - nad_all_inc_fba_sampling + nad_all_inc_all_three,
            length(pfba_nad_all_inc) - nad_all_inc_fba_pfba - nad_all_inc_pfba_sampling + nad_all_inc_all_three,
            length(sampling_nad_all_inc) - nad_all_inc_fba_sampling - nad_all_inc_pfba_sampling + nad_all_inc_all_three,
            nad_all_inc_fba_pfba - nad_all_inc_all_three, nad_all_inc_fba_sampling - nad_all_inc_all_three, nad_all_inc_pfba_sampling - nad_all_inc_all_three, nad_all_inc_all_three,
            length(fba_nad_all_dec) - nad_all_dec_fba_pfba - nad_all_dec_fba_sampling + nad_all_dec_all_three,
            length(pfba_nad_all_dec) - nad_all_dec_fba_pfba - nad_all_dec_pfba_sampling + nad_all_dec_all_three,
            length(sampling_nad_all_dec) - nad_all_dec_fba_sampling - nad_all_dec_pfba_sampling + nad_all_dec_all_three,
            nad_all_dec_fba_pfba - nad_all_dec_all_three, nad_all_dec_fba_sampling - nad_all_dec_all_three, nad_all_dec_pfba_sampling - nad_all_dec_all_three, nad_all_dec_all_three)
)

write.csv(nad_all_overlap_summary, file.path(output_dir, "NAD_All_Venn_Overlap_Summary.csv"), row.names = FALSE)

# Export the complete NAD data for further analysis
write.csv(nad_all_data_filtered, file.path(output_dir, "NAD_All_Reactions_Analysis_Data.csv"), row.names = FALSE)

cat("\nAll NAD Analysis files created in:", output_dir, "\n")
cat("- NAD_All_Venn_Increased_Reactions.png/.svg\n")
cat("- NAD_All_Venn_Decreased_Reactions.png/.svg\n")
cat("- NAD_All_Venn_Overlap_Summary.csv\n")
cat("- NAD_All_Reactions_Analysis_Data.csv\n")
if (exists("comp_plot_data") && nrow(comp_plot_data) > 0) {
  cat("- NAD_Reactions_Compartment_Distribution.png\n")
  cat("- NAD_Compartment_Analysis.csv\n")
}
if (length(nad_all_common_inc) > 0) cat("- NAD_All_Common_Increased_Reactions.csv\n")
if (length(nad_all_common_dec) > 0) cat("- NAD_All_Common_Decreased_Reactions.csv\n")

cat("\nAll NAD-associated reactions analysis complete!\n")



# Detailed analysis of common NAD reactions
cat("\n=== Detailed Analysis of Common NAD Reactions ===\n")

# Get the common reactions (all three methods agree)
nad_all_common_inc <- intersect(intersect(fba_nad_all_inc, pfba_nad_all_inc), sampling_nad_all_inc)
nad_all_common_dec <- intersect(intersect(fba_nad_all_dec, pfba_nad_all_dec), sampling_nad_all_dec)

cat("Common increased NAD reactions:", length(nad_all_common_inc), "\n")
cat("Common decreased NAD reactions:", length(nad_all_common_dec), "\n")

# Analyze subsystems for common reactions
if (length(nad_all_common_inc) > 0 || length(nad_all_common_dec) > 0) {
  
  # Combine all common reactions
  all_common_nad <- c(nad_all_common_inc, nad_all_common_dec)
  common_nad_data <- nad_all_data_filtered[nad_all_data_filtered$Reaction_ID %in% all_common_nad, ]
  
  # Add direction information
  common_nad_data$Direction <- ifelse(common_nad_data$Reaction_ID %in% nad_all_common_inc, "Increased", "Decreased")
  
  # Subsystem analysis
  subsystem_counts <- table(common_nad_data$Subsystem)
  subsystem_df <- data.frame(
    Subsystem = names(subsystem_counts),
    Count = as.numeric(subsystem_counts)
  )
  subsystem_df <- subsystem_df[order(subsystem_df$Count, decreasing = TRUE), ]
  
  cat("\nSubsystem distribution of common NAD reactions:\n")
  print(subsystem_df)
  
  # Create subsystem bar chart
  p_subsystem_bar <- ggplot(subsystem_df, aes(x = reorder(Subsystem, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "#e74c3c", alpha = 0.8) +
    coord_flip() +
    labs(title = "Common NAD Reactions by Subsystem",
         subtitle = "Reactions identified by all three methods",
         x = "Subsystem",
         y = "Number of Reactions") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12)) +
    geom_text(aes(label = Count), hjust = -0.3)
  
  ggsave(file.path(output_dir, "Common_NAD_Subsystems_Bar.png"), p_subsystem_bar, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  # Create subsystem pie chart
  p_subsystem_pie <- ggplot(subsystem_df, aes(x = "", y = Count, fill = Subsystem)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = "Common NAD Reactions by Subsystem (Proportion)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "right") +
    geom_text(aes(label = paste0(Count, "\n(", round(100*Count/sum(Count), 1), "%)")), 
              position = position_stack(vjust = 0.5))
  
  ggsave(file.path(output_dir, "Common_NAD_Subsystems_Pie.png"), p_subsystem_pie, 
         width = 12, height = 8, dpi = 300, bg = "white")
  
  # Compartment analysis for common reactions
  # Extract compartment information from metabolite strings
  extract_nad_compartments <- function(metabolite_string) {
    if (is.na(metabolite_string) || metabolite_string == "") return(character(0))
    
    # Find NAD+ and NADH with their compartments
    nad_matches <- regmatches(metabolite_string, gregexpr("NAD[H+]* \\[[^\\]]+\\]", metabolite_string))[[1]]
    
    if (length(nad_matches) > 0) {
      # Extract just the compartment names
      compartments <- regmatches(nad_matches, gregexpr("\\[([^\\]]+)\\]", nad_matches))
      compartments <- unique(unlist(lapply(compartments, function(x) gsub("\\[|\\]", "", x))))
      return(compartments)
    }
    return(character(0))
  }
  
  # Get compartments for each common reaction
  compartment_list <- lapply(common_nad_data$Associated_Metabolites, extract_nad_compartments)
  all_compartments <- unlist(compartment_list)
  
  if (length(all_compartments) > 0) {
    compartment_counts <- table(all_compartments)
    compartment_df <- data.frame(
      Compartment = names(compartment_counts),
      Count = as.numeric(compartment_counts)
    )
    compartment_df <- compartment_df[order(compartment_df$Count, decreasing = TRUE), ]
    
    cat("\nCompartment distribution of common NAD reactions:\n")
    print(compartment_df)
    
    # Create compartment bar chart
    p_compartment_bar <- ggplot(compartment_df, aes(x = reorder(Compartment, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#3498db", alpha = 0.8) +
      coord_flip() +
      labs(title = "Common NAD Reactions by Compartment",
           subtitle = "Reactions identified by all three methods",
           x = "Compartment",
           y = "Number of Reactions") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)) +
      geom_text(aes(label = Count), hjust = -0.3)
    
    ggsave(file.path(output_dir, "Common_NAD_Compartments_Bar.png"), p_compartment_bar, 
           width = 10, height = 6, dpi = 300, bg = "white")
    
    # Create compartment pie chart
    p_compartment_pie <- ggplot(compartment_df, aes(x = "", y = Count, fill = Compartment)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      labs(title = "Common NAD Reactions by Compartment (Proportion)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right") +
      geom_text(aes(label = paste0(Count, "\n(", round(100*Count/sum(Count), 1), "%)")), 
                position = position_stack(vjust = 0.5))
    
    ggsave(file.path(output_dir, "Common_NAD_Compartments_Pie.png"), p_compartment_pie, 
           width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  # Export detailed summary tables
  write.csv(subsystem_df, file.path(output_dir, "Common_NAD_Subsystem_Summary.csv"), row.names = FALSE)
  if (exists("compartment_df")) {
    write.csv(compartment_df, file.path(output_dir, "Common_NAD_Compartment_Summary.csv"), row.names = FALSE)
  }
  
  cat("\nDetailed analysis files created:\n")
  cat("- Common_NAD_Subsystems_Bar.png\n")
  cat("- Common_NAD_Subsystems_Pie.png\n")
  cat("- Common_NAD_Compartments_Bar.png\n")
  cat("- Common_NAD_Compartments_Pie.png\n")
  cat("- Common_NAD_Subsystem_Summary.csv\n")
  cat("- Common_NAD_Compartment_Summary.csv\n")
  
} else {
  cat("No common NAD reactions found for detailed analysis.\n")
}


# Install if needed: install.packages("UpSetR")
# UpSet plot for all NAD reactions
library(UpSetR)

# Create binary matrix directly from the reaction lists instead
# For decreased all NAD reactions
decreased_reactions_all <- unique(c(fba_nad_all_dec, pfba_nad_all_dec, sampling_nad_all_dec))

# Create a proper binary matrix
upset_matrix_dec_all <- data.frame(
  FBA = as.numeric(decreased_reactions_all %in% fba_nad_all_dec),
  pFBA = as.numeric(decreased_reactions_all %in% pfba_nad_all_dec),
  FVA_sampling = as.numeric(decreased_reactions_all %in% sampling_nad_all_dec)
)

# Create UpSet plot for decreased all NAD reactions
png(file.path(output_dir, "NAD_All_Decreased_UpSet.png"), width = 10, height = 6, res = 300, units = "in")
upset(upset_matrix_dec_all, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()

pdf(file.path(output_dir, "NAD_All_Decreased_UpSet.pdf"),
    width = 8, height = 6)   # no res needed for vector formats
upset(upset_matrix_dec_all, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#3498db",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      scale.sets = "identity",
      line.size = 1)
dev.off()

# For increased all NAD reactions
increased_reactions_all <- unique(c(fba_nad_all_inc, pfba_nad_all_inc, sampling_nad_all_inc))

upset_matrix_inc_all <- data.frame(
  FBA = as.numeric(increased_reactions_all %in% fba_nad_all_inc),
  pFBA = as.numeric(increased_reactions_all %in% pfba_nad_all_inc),
  FVA_sampling = as.numeric(increased_reactions_all %in% sampling_nad_all_inc)
)

png(file.path(output_dir, "NAD_All_Increased_UpSet.png"), width = 10, height = 6, res = 300, units = "in")
upset(upset_matrix_inc_all, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 1.2,
      point.size = 3,
      line.size = 1)
dev.off()


pdf(file.path(output_dir, "NAD_All_Increased_UpSet.pdf"), width = 8, height = 6)
upset(upset_matrix_inc_all, 
      sets = c("FBA", "pFBA", "FVA_sampling"),
      order.by = "freq",
      main.bar.color = "#e74c3c",
      sets.bar.color = "#2c3e50",
      text.scale = 2,
      point.size = 3,
      line.size = 1)
dev.off()

# NAD Network Analysis - Gene-Subsystem Networks

# Function to save plots in both PDF and SVG
save_plot_both <- function(plot_obj, filename, width = 8, height = 6) {
  # PDF
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot_obj, 
         width = width, height = height, bg = "white")
  # SVG
  ggsave(file.path(output_dir, paste0(filename, ".svg")), plot_obj, 
         width = width, height = height, bg = "white")
}

# Function to extract genes and count them
count_genes_nad <- function(common_data, direction) {
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

###############################################################################
# NAD DECREASED REACTIONS NETWORK ANALYSIS
###############################################################################
library(stringr)

if (length(nad_all_common_dec) > 0) {
  
  # Get data for commonly decreased NAD reactions
  nad_common_data_dec <- nad_all_data_filtered[nad_all_data_filtered$Reaction_ID %in% nad_all_common_dec, ]

  nad_common_data_dec = nad_common_data_dec %>%
    mutate(
      Associated_Genes = str_replace(Associated_Genes, "TSG101", "LDH"),
      Gene_Rules = str_replace(Gene_Rules, "TSG101", "LDH")
    )
  
  # Get subsystem distribution
  nad_subsystem_counts_dec <- table(nad_common_data_dec$Subsystem)
  nad_subsystem_df_dec <- data.frame(
    Subsystem = names(nad_subsystem_counts_dec),
    Count = as.numeric(nad_subsystem_counts_dec),
    Direction = "Decreased"
  )
  nad_subsystem_df_dec <- nad_subsystem_df_dec[order(nad_subsystem_df_dec$Count, decreasing = TRUE), ]
  
  # Get top 3 subsystems for network
  top_nad_subsystems_dec <- head(nad_subsystem_df_dec, 3)
  
  # Extract genes for top NAD subsystems
  top_nad_subsystem_genes_dec <- list()
  for (subsys in top_nad_subsystems_dec$Subsystem) {
    subsys_reactions <- nad_common_data_dec[nad_common_data_dec$Subsystem == subsys, ]
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
    top_nad_subsystem_genes_dec[[subsys]] <- head(gene_df[order(gene_df$Count, decreasing = TRUE), ], 5)
  }
  
  combined_nad_gene_data_dec <- do.call(rbind, top_nad_subsystem_genes_dec)
  
  # Create network for NAD decreased
  library(igraph)
  library(ggraph)
  
  edges_nad_dec <- combined_nad_gene_data_dec[, c("Subsystem", "Gene", "Count")]
  colnames(edges_nad_dec) <- c("from", "to", "weight")
  
  subsystem_nodes_nad_dec <- data.frame(
    name = top_nad_subsystems_dec$Subsystem,
    type = "subsystem",
    size = top_nad_subsystems_dec$Count
  )
  gene_nodes_nad_dec <- data.frame(
    name = unique(combined_nad_gene_data_dec$Gene),
    type = "gene", 
    size = 3
  )
  nodes_nad_dec <- rbind(subsystem_nodes_nad_dec, gene_nodes_nad_dec)
  
  # Convert type to factor with specific order
  nodes_nad_dec$type <- factor(nodes_nad_dec$type, levels = c("subsystem", "gene"))
  
  g_nad_dec <- graph_from_data_frame(edges_nad_dec, vertices = nodes_nad_dec, directed = TRUE)
  
  p_network_nad_dec <- ggraph(g_nad_dec, layout = "fr") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "gray60") +
    geom_node_point(aes(size = size, color = type)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    scale_color_manual(values = c("subsystem" = "#3498db", "gene" = "#2980b9"), 
                       breaks = c("subsystem", "gene"), 
                       name = "") +
    scale_size_continuous(range = c(3, 8), 
                          breaks = function(x) {
                            range_vals <- range(top_nad_subsystems_dec$Count)
                            seq(ceiling(range_vals[1]), floor(range_vals[2]), by = max(1, floor(diff(range_vals)/3)))
                          }, 
                          name = "Reactions") +
    scale_edge_width_continuous(name = "# of rxns", range = c(0.5, 3), 
                                breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
    theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.4, "cm"),
          legend.margin = margin(l = 5))
  
  save_plot_both(p_network_nad_dec, "NAD_Network_Decreased")
  
  cat("NAD decreased network created: NAD_Network_Decreased.pdf/.svg\n")
  
  # Export NAD gene data
  write.csv(combined_nad_gene_data_dec, file.path(output_dir, "NAD_Gene_Data_Decreased.csv"), row.names = FALSE)
  
} else {
  cat("No commonly decreased NAD reactions found for network analysis.\n")
}

# Load necessary libraries
library(readr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(colorRamp2)

#Change to your working directory
setwd("C:/Users/argendall/OneDrive - LA TROBE UNIVERSITY/Pumpkin_Etch_Project/Sequencing/DE_TEST")

#load(".RData") # only to re-load if required

#Define  colors for the three different scales, column group vector to fix the order of plotted columns, must have same number as group i.e. here 3.
column_group <- c("Unaffected", "Edge", "Etched")


custom_colors <- colorRampPalette(c("white", brewer.pal(n = 9, name = "YlGnBu")))(256)
#custom_colors <- colorRampPalette(brewer.pal(1,"RdYlBu"))(256) # for log2 data alternative

custom_colors_REDS <- colorRampPalette(c("white", brewer.pal(n = 9, name = "Reds")))(256) #TPM display

########################################################################################################
# Filter the MERGED_1 dataset based on the symbols in linput_file
# MERGED_1 is a large RNASeq data set from the Pumpkin Etch project. 

# VARIABLES ARE:
# group - Unaffected, Edge, Etched
# symbol - LOC (i.e NCBI gene code)
# sample - id/name , here they are called  1001 through 1009 (Unaffected: 1001-2-3, Edge: 1004-5-6, Etched: 1007-8-9)
# replicate - here 1,2 or 3
# tpm - transcripts per million/expression level
#
# all other columns are not used/relevant here

#if necessary - import MERGED_1 data
MERGED_1 <- read_csv("MERGED_1.csv")


########################################################################################################
# From here, "process_data", line 48,  you can run it once, then call each file as specified in the 
##                    INPUTS                  ##  
## section, lines 198 onwards, at the end.
########################################################################################################
#PLOTTING MEAN data across replicates
# Define the function that will generate 6 heatmaps for each input file, first as means, then individual replicates
process_data <- function(input_file, ncol = 3, fin_width = 35, fin_height = 10) {
  
  # Extract the prefix of the input file name
  file_prefix <- tools::file_path_sans_ext(basename(input_file))
  
  # Read the CSV file
  input_data <- read_csv(input_file)
  
  # Filter the MERGED_1 dataset based on the symbols in input_data
  filtered_data <- MERGED_1 %>%
    filter(symbol %in% input_data$symbol)
  
  # Aggregate data
  aggregated_data <- aggregate(tpm ~ group + symbol, data = filtered_data, FUN = mean)
  
  # Reshape the data for the heatmap
  heatmap_data <- dcast(aggregated_data, symbol ~ group, value.var = "tpm")
  rownames(heatmap_data) <- heatmap_data$symbol
  heatmap_data <- heatmap_data[, -1]  # Remove the 'symbol' column
  
  # Ensure row order matches the input file order - useful for plotting pathways etc. If you don't want to do this, you remove thse two lines, and edit the heatmap scripts to cluster_rows=TRUE
  row_order <- input_data$symbol
  heatmap_data <- heatmap_data[match(row_order, rownames(heatmap_data)), ]
  
  
  heatmap_data_log2 <- log2(heatmap_data + 1) # Create log2 (TPM+1) transformation for alternate display. 
  
  #column_group <- c("Unaffected", "Edge", "Etched")  # Example column group vector to fix the order of plotted columns
  ordered_columns <- match(column_group, colnames(heatmap_data))
  heatmap_data_2 <- heatmap_data[, ordered_columns]
  heatmap_data_log2 <- heatmap_data_log2[, ordered_columns]
  
  ####################################################################
  # Generate heatmaps
  #HEATMAP 1
  png(paste0(file_prefix, "_mean_TPM.png"))
  mean_TPM <- pheatmap(heatmap_data_2, 
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE, 
                       display_numbers = FALSE, 
                       color = custom_colors_REDS,
                       cellwidth = 10, 
                       cellheight = 10, 
                       main = paste(file_prefix, "(TPM)"))
  dev.off()
  
  #HEATMAP 2
  png(paste0(file_prefix, "_mean_TPM_row_scaled.png"))
  mean_Row_scale <- pheatmap(heatmap_data_2, 
                             cluster_rows = FALSE, 
                             cluster_cols = FALSE, 
                             display_numbers = FALSE, 
                             scale = "row",
                             cellwidth = 10, 
                             cellheight = 10, 
                             main = paste(file_prefix, "row scaled (TPM)"))
  dev.off()
  
  #HEATMAP 3
  png(paste0(file_prefix,"_mean_log2.png"))
  mean_log2 <- pheatmap(heatmap_data_log2, 
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE, 
                        display_numbers = FALSE, 
                        color = custom_colors,
                        cellwidth = 10, 
                        cellheight = 10, 
                        main = paste(file_prefix, " Mean log2(TPM+1)"))
  dev.off()
  
  # Arrange the "means" plots side by side
  png(paste0(file_prefix, "_Means.png"), width = 3500, height = 1000, res = 300)
  combined_means <- grid.arrange(mean_TPM$gtable, mean_log2$gtable, mean_Row_scale$gtable, ncol = 3)
                                 
  dev.off()



################################################
### REPLICATES #################################
################################################
  # Aggregate data
  aggregated_data <- filtered_data %>%
    group_by(symbol, group, replicate, sample) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE), .groups = 'drop')
  
  # Reshape the data for the heatmap
  heatmap_data <- dcast(aggregated_data, symbol ~ sample, value.var = "tpm")
  rownames(heatmap_data) <- heatmap_data$symbol
  heatmap_data <- heatmap_data[, -1]  # Remove the 'symbol' column
  
  # Ensure row order matches the input file order
  row_order <- input_data$symbol
  heatmap_data <- heatmap_data[match(row_order, rownames(heatmap_data)), ]
  
  # Generate heatmaps
  png(paste0(file_prefix, "_reps_(TPM).png"))
  reps_TPM <- pheatmap(heatmap_data, 
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE, 
                       display_numbers = FALSE, 
                       color = custom_colors_REDS,
                       cellwidth = 10, 
                       cellheight = 10, 
                       main = paste(file_prefix, "reps (TPM)"))
  dev.off()
  
  png(paste0(file_prefix, "_reps_scaled_TPM.png"))
  reps_scaled_TPM <- pheatmap(heatmap_data, 
                              cluster_rows = FALSE, 
                              cluster_cols = FALSE, 
                              display_numbers = FALSE, 
                              scale = "row",
                              cellwidth = 10, 
                              cellheight = 10, 
                              main = paste(file_prefix, " reps row scaled (TPM)"))
  dev.off()
  
  # LOG2 transformed
  heatmap_data_log2 <- log2(heatmap_data + 1)
  
  png(paste0(file_prefix, "_reps_log2_TPM.png"))
  #custom_colors <- brewer.pal(n = 9, name = "YlGnBu")
  reps_log2_TPM <- pheatmap(heatmap_data_log2, 
                            cluster_rows = FALSE, 
                            cluster_cols = FALSE, 
                            display_numbers = FALSE, 
                            legend = TRUE,
                            color = custom_colors,
                            cellwidth = 10, 
                            cellheight = 10, 
                            main = paste(file_prefix, " reps log2(TPM+1)"))
  dev.off()
  
  # Arrange the replicates plots side by side
  png(paste0(file_prefix, "_Replicates.png"), width = 3500, height = 1000, res = 300)
  combined_replicates <- grid.arrange(reps_TPM$gtable, reps_log2_TPM$gtable, reps_scaled_TPM$gtable, ncol = 3)
  dev.off()

  
  # Arrange the mean and replicate plots
  png(paste0(file_prefix, "_combined_heatmaps.png"), width = fin_width, height = fin_height, res = 300)
  All_HMs <- grid.arrange(combined_replicates, combined_means, nrow = 2, top = paste(file_prefix,"Heatmaps"))
  dev.off()
}

####################################################################################################################################
##     INPUTS                                                                                                                     ##  
####################################################################################################################################
#This specifies the input file, and final width and height of teh 6 heatmaps - adjust so they fit in a single page or as desired.
process_data("MERGED_1_patho.csv", fin_width = 4000, fin_height = 8000)
process_data("Lignin_TF_anno_2.csv", fin_width = 4000, fin_height = 4000)


#process_data("CADs_v2.csv", fin_width = 3000, fin_height = 2000)
#process_data("LIPOX.csv", fin_width = 3000, fin_height = 3000)
#process_data("DESEQ_MERGED_1_TOP50_LOC.csv", fin_width = 3000, fin_height = 10000)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tools)



########################################################################################################
# Filter the MERGED_1 dataset based on the symbols in linput_file
# MERGED_1 is a large RNASeq data set from the Pumpkin Etch project. 

# VARIABLES ARE:
# group - Unaffected, Edge, Etched (sample/treatment groups) 
# symbol - LOC (i.e NCBI gene code)
# sample - id/name , here they are called  1001 through 1009 (Unaffected: 1001-2-3, Edge: 1004-5-6, Etched: 1007-8-9)
# replicate - here 1,2 or 3
# tpm - transcripts per million/expression level
# name -gene nanme/description, used to annotated the facet headings
# all other columns are not used/relevant here


########################################################################################################


# Set working directory and load environment (assuming .RData exists)
setwd("C:/Users/argendall/OneDrive - LA TROBE UNIVERSITY/Pumpkin_Etch_Project/Sequencing/DE_TEST")
#load(".RData")

#Import MERGED_1 data -replace wiht you actual data file
MERGED_1 <- read_csv("MERGED_1.csv")

# Define custom colors
BOX_COLORS <- c("Unaffected" = "#F5F5F5", "Edge" = "#DFC27D", "Etched" = "#8C510A")

# Read input data outside the function (assuming it's in CSV format)
input_data <- read_csv("MERGED_1.csv")


########################################################################################################
# From here, "process_data", line 46,  you can run it once, then call each file as specified in the 
##                    INPUTS                  ##  
## section, lines 198 onwards, at the end.
########################################################################################################

# Function to process data and generate plot
process_data <- function(input_file, ncol = 3, plot_width = 800, plot_height = 600) {
  # Extract file prefix
  file_prefix <- tools::file_path_sans_ext(basename(input_file))
  
  # Input file validation
  if (!file.exists(input_file)) stop(paste0("Error: File '", input_file, "' not found."))
  
  # Read the input file
  input_symbols <- read_csv(input_file)
  
  # Filter data based on input symbols
  filtered_data <- input_data %>%
    filter(symbol %in% input_symbols$symbol)
  
  # Set the order of symbols as in input_data
  symbol_order <- unique(input_symbols$symbol)
  filtered_data$symbol <- factor(filtered_data$symbol, levels = symbol_order)
  
  # Arrange filtered data to maintain order
  filtered_data <- filtered_data %>%
    arrange(factor(symbol, levels = symbol_order))
  
  # Y lim 0 Generate plot
  png(paste0(file_prefix, "_BOXPLOTS.png"), width = plot_width, height = plot_height)
  
  plot <- ggplot(filtered_data, aes(x = factor(group, levels = c("Unaffected", "Edge", "Etched")), y = tpm, fill = factor(group, levels = c("Unaffected", "Edge", "Etched")))) +
    geom_boxplot() +
    labs(title = paste(file_prefix, " BOXPLOTS"),
         x = NULL,  # Remove x-axis label
         y = "TPM") +
    theme(axis.text.x = element_blank(),  # Remove x-axis labels
          axis.title.x = element_blank()) +
    facet_wrap(~ factor(symbol, levels = symbol_order) + name, ncol = ncol, scales = "free_y") +
    scale_fill_manual(values = BOX_COLORS, labels = c("Unaffected", "Edge", "Etched"), name = "Tissue", guide = guide_legend(title.position = "top")) +  # Customize legend title and order of labels
    ylim(0, NA)  # Set the minimum y-axis limit to 0
  
  print(plot)
  dev.off()
}

####################################################################################################################################
##     INPUTS                                                                                                                     ##  
####################################################################################################################################
#This specifies the input file, and final width and height of teh bxoplot s - adjust so they fit in a single page or as desired.
# Example usage with ncol, width, and height specified
process_data("Lignin_TF_anno_2.csv", ncol = 3, plot_width = 1000, plot_height = 1000)

#process_data("CADs_v2.csv", ncol = 5, plot_width = 1500, plot_height = 600)

#process_data("MERGED_1_patho.csv", ncol = 5, plot_width = 1500, plot_height = 3000)
#process_data("DESEQ_MERGED_1_TOP50_LOC.csv", ncol = 8, plot_width = 2000, plot_height = 3000)






** RNASeq
Scripts for RNA seq analysis**

A folder for sharing various scripts for RNAseq analysis

**TPM_HEATMAP_GoI_2.R:**                      an R script for automating the production of 6 versions of heamtmaps to quicukly examine RNAseq expression patterns for a  list of genes of interest. It uses:
  MERGED_1.csv (here as a .zip file):     sample data, of genes and tpm/exxression values
  Lignin_TF_anno_2.csv:                   a list of genes (here as LOC codes)  to filter MERGED_1
  MERGED_1_patho.csv:                     another list of genes (here as LOC codes)  to filter MERGED_1
  -----------------
  Lignin_TF_anno_2_combined_heatmaps.png  Output of TPM_HEATMAP_GoI_2.R with  Lignin_TF_anno_2.csv, and individual images are also produced.


  

BOXPLOTS_3.R                               an R script for automating the production ofboxplots to quickly examine RNAseq expression patterns for a  list of genes of interest. It uses:
  MERGED_1.csv (here as a .zip file):     sample data, of genes and tpm/exxression values
   Lignin_TF_anno_2.csv:                  a list of genes (here as LOC codes)  to filter MERGED_1
    -----------------
  Lignin_TF_anno_2_BOXPLOTS.png                    Output of BOXPLOTS_3.R  with  Lignin_TF_anno_2.csv, and individual images are also produced.
                  

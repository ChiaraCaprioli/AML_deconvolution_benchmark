# Prepare input to CIBERSORTx
# Chiara Caprioli
# 11th Dec 2023

## Set up
# Load libraries
library(tidyverse)
library(Seurat)
library(SummarizedExperiment)
source("~/Desktop/working/bulk_ct_deconvolution/utils.R")

# sc reference to build signature matrix
sc_reference <- readRDS("~/Desktop/working/bulk_ct_deconvolution/sc_reference/sc_reference_sample.rds")
sig_matrix <- SigMatrixCibersortx(seurat = sc_reference, label = "aggregated_lineage2")
write.table(
  sig_matrix, 
  row.names = F, 
  col.names = T,
  quote = F,
  file = "~/Desktop/working/bulk_ct_deconvolution/deconvolution/CIBERSORTx/sig_matrix.txt",
  sep="\t"
) 

# pb data
pb_data_sample <- readRDS("~/Desktop/working/bulk_ct_deconvolution/sc_reference/pb_data_sample.rds")
mix_matrix <- MixMatrixPB(pb_counts = pb_data_sample$bulk.counts)

write.table(
  mix_matrix, 
  row.names = F, 
  quote = F,
  sep="\t",
  file = "~/Desktop/working/bulk_ct_deconvolution/deconvolution/CIBERSORTx/mix_matrix_summed.txt"
)


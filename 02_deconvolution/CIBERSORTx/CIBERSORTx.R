---
title: "CIBERSORTx - Cell type deconvolution"
author: "Chiara Caprioli"
date: "11 Dec 2023"
---

# Libraries
library(tidyverse)
library(Seurat)

# Utils
wdir <- "~/Desktop/working/bulk_ct_deconvolution/"
source(paste0(wdir, "utils.R"))
outdir <- paste0(wdir, "02_deconvolution/CIBERSORTx/") 
if (!dir.exists(outdir)) {
  dir.create(outdir)
}
outdir2 <- paste0(wdir, "02_deconvolution/results")

### Prepare input to CIBERSORTX ###
# sc reference to build signature matrix
sc_reference <- readRDS(paste0(wdir, "01_prepare_input/sc_reference_sample.rds"))

sig_matrix <- SigMatrixCibersortx(seurat = sc_reference, label = "aggregated_lineage2")

write.table(
  sig_matrix, 
  row.names = F, 
  col.names = T,
  quote = F,
  file = paste0(outdir, "sig_matrix.txt"),
  sep="\t"
) 

# pb data
pb_data_sample <- readRDS(paste0(wdir, "01_prepare_input/pb_data_sample.rds"))

L_bulk <- list("summed" = pb_data_sample$bulk.counts, "tpm" = pb_data_sample$bulk.counts.tpm)

for (b in names(L_bulk)) {
  
  mix_matrix <- MixMatrixPB(pb_counts = L_bulk[[b]])
  write.table(
    mix_matrix, 
    row.names = F, 
    quote = F,
    sep="\t",
    file = paste0(outdir, "mix_matrix_", b, ".txt") 
  )
  
}

### Run CIBERSORTX on website (https://cibersortx.stanford.edu) ###
# Summed counts S-mode
# Summed counts no batch correction
# TPM counts S-mode
# TPM counts no batch correction

### Process CIBERSORTX output ###
res <- WrapCibersortx(path = outdir)
saveRDS(res, paste0(outdir2, "/CIBERSORTx_proportions.rds"))

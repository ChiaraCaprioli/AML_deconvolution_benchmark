---
title: "DWLS - Cell type deconvolution"
author: "Chiara Caprioli"
date: "1 Dec 2023"
---
  
# Libraries
# install.packages("DWLS")
library(DWLS)
library(Seurat)
library(Biobase)
library(ggpubr)
library(tidyverse)

# Utils
wdir <- "~/Desktop/working/bulk_ct_deconvolution/"
outdir <- paste0(wdir, "benchmarking/DWLS/") 
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Data preparation
## Bulk
pb_summed <- readRDS(paste0(wdir, "sc_reference/pb_data_sample_es_summed.rds"))
bulk.mtx.summed = exprs(pb_summed)

pb_tpm <- readRDS(paste0(wdir, "sc_reference/pb_data_sample_es_tpm.rds"))
bulk.mtx.tpm = exprs(pb_tpm)

## Single-cell
seurat <- readRDS(paste0(wdir, "sc_reference/sc_reference_sample.rds"))

## Signature from single-cell data
# !!! time consuming

signature_mast <- buildSignatureMatrixUsingSeurat( ### positive markers by Seurat::FindMarkers (MAST method) 
  scdata = as.matrix(seurat@assays$RNA@counts),
  id = seurat$aggregated_lineage2,
  path = outdir,
  diff.cutoff = 0.5,
  pval.cutoff = 0.01
)

# Deconvolution
est.prop <- apply(bulk.mtx.summed, 2, function(x){
  b = setNames(x, rownames(bulk.mtx.summed))
  tr <- trimData(signature_mast, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(est.prop) <- colnames(tr$sig)
est.prop = apply(est.prop,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
est.prop = apply(est.prop,2,function(x) x/sum(x)) #explicit STO constraint
est.prop <- t(est.prop)
saveRDS(est.prop, paste0(outdir, "/DWLS_proportions.rds"))

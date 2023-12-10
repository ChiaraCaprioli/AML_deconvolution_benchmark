---
title: "BisqueRNA - Cell type deconvolution"
author: "Chiara Caprioli"
date: "10 Dec 2023"
---
  
# Libraries
# install.packages("BisqueRNA")
library(BisqueRNA)
library(Seurat)
library(Biobase)
library(tidyverse)

# Utils
wdir <- "~/Desktop/working/bulk_ct_deconvolution/"
source(paste0(wdir, "utils.R"))
outdir <- paste0(wdir, "deconvolution/BisqueRNA/") 
if (!dir.exists(outdir)) {
  dir.create(outdir)
}
outdir2 <- paste0(wdir, "deconvolution/results")

# Data preparation
## Bulk
pb_summed <- readRDS(paste0(wdir, "sc_reference/pb_data_sample_es_summed.rds"))
pb_tpm <- readRDS(paste0(wdir, "sc_reference/pb_data_sample_es_tpm.rds"))

## List of bulk matrices
list_bulk <- list("summed" = pb_summed, "tpm" = pb_tpm)

## Single-cell
sc_reference_sample_es <- readRDS(paste0(wdir, "sc_reference/sc_reference_sample_es.rds"))

## Markers
sc_ct_markers <- read_csv(paste0(wdir, "sc_reference/sc_ct_markers.csv"))
selected_markers <- sc_ct_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 2) 

list_markers <- list("no_markers" = NULL, "with_markers" = unique(selected_markers$gene))

# Deconvolution

## Reference-based decomposition
refbased_est_prop <- BenchmarkBisqueRefBased(
  L_bulk = list_bulk,
  sc.eset = sc_reference_sample_es, 
  L_markers = list_markers
)
saveRDS(refbased_est_prop, paste0(outdir2, "/Bisque.refbased_proportions.rds"))


# AML_deconvolution_benchmark

## Contents
- [Overview](#Overview)
- [Input](#Input)
- [Output](#Output)
- [Folder structure](#Folder-structure)
- [Required packages](#Required-packages)
- [References](#References)

## Overview
*Cell type deconvolution* refers to the task of inferring proportions of individual cell types in bulk RNA-sequenced samples, using a reference of either predetermined cell type-specific gene expression signatures or single-cell RNA-sequencing (scRNA-seq) data from the same tissue type. The use of a scRNA-seq reference is deemed advantageous in terms of lack of bias and resolution of cell type information.

Many different deconvolution approaches exist and have been extensively benchmarked; the most relevant factors that have been found to affect deconvolution results include (i) the type of data normalization, (ii) the reference matrix (which should include all cell types being part of the mixtures), iii) the strategy for marker selection.

Our primary objective is estimating cell type proportions from bulk RNA-seq data of AML samples (e.g., TCGA and BeatAML cohorts), using our scRNA-seq dataset of AML samples as a reference.

To pick the deconvolution method with the highest performance when using our specific scRNA-seq reference, we perform benchmarking by creating pseudo-bulk data matching the scRNA-seq reference and running different deconvolution approaches. (The preferred strategy would be comparing scRNA-seq data with matched RNA-seq data from the same samples, which we currently lack).

For now, the following methods are implemented:
- ```MuSiC```
- ```DWLS```
- ```Bisque``` 
- ```CIBERSORTx``` 

Pseudo-bulk data are created from raw scRNA-seq counts treated in two distinct ways, which are subject to benchmarking as well: 
- sum by sample
- TPM transformation by sample

## Input
A ```Seurat``` object. Metadata should include *sample information* and *cell type labels* from previous annotation with any method of choice.

## Output
Output for each deconvolution method consists of a list of matrices with estimated cell type proportions by sample for each mode (e.g., summed/tpm counts, with/without markers, etc).

Estimated vs real cell type proportions are then compared across modalities using different metrics:
- root-mean-squared deviation (RMSD)
- mean absolute difference (mAD)
- Pearson’s correlation (R)
  
## Folder structure
```
AML_deconvolution_benchmark/
├── 01_prepare_input
│   └── 01_prepare_input.Rmd
├── 02_deconvolution
│   ├── BisqueRNA
│   │   └── BisqueRNA.R
│   ├── CIBERSORTx
│   │   └── CIBERSORTx.R
│   ├── DWLS
│   │   └── DWLS.R
│   ├── MuSiC
│   │   └── MuSiC.Rmd
│   └── results
├── 03_benchmarking
│   ├── benchmarking.Rmd
│   └── results
└── utils.R
```
## Required packages 
```
> sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] compiler_4.2.2    fastmap_1.1.1     cli_3.6.1         htmltools_0.5.5   tools_4.2.2       rstudioapi_0.15.0 yaml_2.3.7       
 [8] rmarkdown_2.25    knitr_1.45        xfun_0.41         digest_0.6.33     rlang_1.1.2       evaluate_0.23 
```

```
Biobase
BisqueRNA
ComplexHeatmap
DWLS
ggdendro
ggpubr
MuSiC #devtools::install_github('xuranw/MuSiC') 
rsample
scuttle
Seurat
SingleCellExperiment
tidyverse
```

## References 
- Avila Cobos F et al., Benchmarking of cell type deconvolution pipelines for transcriptomics data. <https://doi.org/10.1038/s41467-020-19015-1>
- Avila Cobos F et al., Effective methods for bulk RNA-seq deconvolution using scnRNA-seq transcriptomes. <https://doi.org/10.1186/s13059-023-03016-6>
- Zeng AGX et al., A cellular hierarchy framework for understanding heterogeneity and predicting drug response in acute myeloid leukemia. https://doi.org/10.1038/s41591-022-01819-x
- Wang X et al., Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. <https://doi.org/10.1038/s41467-018-08023-x>
- Tsoucas D et al., Accurate estimation of cell-type composition from gene expression data. https://doi.org/10.1038/s41467-019-10802-z
- Jew B et al., Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. https://doi.org/10.1038/s41467-020-15816-6
- Steen CB et al., Profiling cell type abundance and expression in bulk tissues with CIBERSORTx. https://cibersortx.stanford.edu/about.php

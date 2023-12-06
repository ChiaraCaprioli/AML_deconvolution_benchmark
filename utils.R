###################### Preparation of single-cell reference ######################
###################### Get pseudo-bulk counts from single-cell reference ######################
# The function returns a list including: 
# [[1]] matrix of summed counts (bulk.counts)  
# [[2]] matrix of real cell type counts (num.real) 
# [[3]] matrix of TPM counts (bulk.counts.tpm)  
# [[4]] matrix of real cell type proportions (prop.real) 

GetPsedoBulkCounts <- function(seurat, cluster_var, sample_var) { 
  sce <- as.SingleCellExperiment(seurat) 
  
  ## aggregate counts and cell type counts 
  pb_counts <- MuSiC::bulk_construct(
    sce, 
    clusters = cluster_var, 
    samples = sample_var
  )  
  
  ## TPM counts 
  pb_counts$bulk.counts.tpm <- calculateTPM(
    pb_counts$bulk.counts, 
    lengths = NULL # no length required for UMI-based methods
  ) 
  
  ## calculate cell type proportions by sample  
  pb_counts$prop.real = relative.ab(pb_counts$num.real, by.col = F) 
  
  return(pb_counts) 
  
}

###################### MuSiC ######################
###################### Benchmark MuSiC direct mode ######################

# requires a named list
# option can be one of "both", "no_markers" and "with_markers"

BenchmarkMusicDirect <- function(list_bulk, sc.sce, clusters, markers, samples, option) {
  
  L_prop_est <- list()
  for (i in names(list_bulk)) {
    
    if (option == "both" | option == "no_markers") {
      message(
        glue::glue(
          "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with direct mode, {i} counts and no markers"
        )
      )
      direct_est_prop_no_markers <- music_prop(
        bulk.mtx = list_bulk[[i]], 
        sc.sce = sce, 
        markers = NULL,
        clusters = clusters,
        samples = samples,
        verbose = F
      )
      L_prop_est[[i]][["no_markers"]] <- direct_est_prop_no_markers
    }
    
    if (option == "both" | option == "with_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with direct mode, {i} counts and markers"
          )
        )
      direct_est_prop_markers <- music_prop(
        bulk.mtx = list_bulk[[i]], 
        sc.sce = sce, 
        markers = markers,
        clusters = clusters,
        samples = samples,
        verbose = F
      )
      L_prop_est[[i]][["with_markers"]] <- direct_est_prop_markers
      }
  }
  return(L_prop_est)
  }

###################### Get cell type hierarchical clusters ######################
GetCellTypeHierarchies <- function(x, clusters, markers, samples, option) {
  
  L_plots <- list()
  if (option == "both" | option == "no_markers") {
    message(
      glue::glue(
        "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Getting hierarchies with no markers"
      )
    )
    basis <- music_basis(
      x = sce, 
      markers = NULL,
      clusters = clusters,
      samples = samples, 
      verbose = F
    )
  
    d <- dist(t(log(basis$Disgn.mtx + 1e-6)), method = "euclidean")
    hc <- hclust(d, method = "complete")
    p <- ggdendrogram(hc) +
      theme_bw() +
      labs(
        title = "No markers",
        subtitle = "Cluster log(Design Matrix)",
        y = "Height"
      ) +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2),
        plot.subtitle = element_text(hjust = 0.5, vjust = 1) 
      )
    
    L_plots[["no_markers"]] <- p
    
  }
  
  if (option == "both" | option == "with_markers") {
    message(
      glue::glue(
        "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Getting hierarchies with markers"
      )
    )
    basis <- music_basis(
      x = sce, 
      markers = markers,
      clusters = clusters,
      samples = samples, 
      verbose = F
    )
    
    d <- dist(t(log(basis$Disgn.mtx + 1e-6)), method = "euclidean")
    hc <- hclust(d, method = "complete")
    p <- ggdendrogram(hc) +
      theme_bw() +
      labs(
        title = "With markers",
        subtitle = "Cluster log(Design Matrix)",
        y = "Height"
      ) +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, vjust = 2),
        plot.subtitle = element_text(hjust = 0.5, vjust = 1) 
        )
    
    L_plots[["with_markers"]] <- p
  }
  return(L_plots)
}

###################### CIBERSORTx ######################
SigMatrixCibersortx <- function(seurat, label) {
  sig_matrix <- as.matrix(GetAssayData(seurat, slot = "count"))
  colnames(sig_matrix) <- as.character(seurat@meta.data[[label]])
  sig_matrix <- cbind(rownames(sig_matrix), sig_matrix)
  colnames(sig_matrix)[1] <- "gene_name"
  return(sig_matrix)
}

MixMatrixPB <- function(pb_counts) {
  mix_matrix <- pb_counts
  mix_matrix <- cbind(rownames(mix_matrix), mix_matrix)
  colnames(mix_matrix)[1] <- "gene_name"
  return(mix_matrix)
} 


 


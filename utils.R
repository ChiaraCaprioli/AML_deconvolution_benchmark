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

BenchmarkMusicDirect <- function(L_bulk, L_markers, sc.sce, clusters, samples) {
  
  L_final <- list()

  for (i in names(L_bulk)) {
    
    L_prop <- list()
    for (m in names(L_markers)) {
      if (m == "no_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with direct mode, {i} counts and no markers"
          )
        )
      }
      if (m == "with_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with direct mode, {i} counts and markers"
          )
        )
      }
      
      est.prop <- music_prop(
        bulk.mtx = L_bulk[[i]], 
        sc.sce = sce, 
        markers = L_markers[[m]],
        clusters = clusters, 
        samples = samples,
        verbose = F
      )
      
      names(est.prop)[1] <- paste(i,m,sep = "_")
      L_prop[[m]] <- est.prop[1]
    }
    L_final[[i]] <- L_prop
  }
  L_save <- do.call(c, unlist(L_final, recursive = F, use.names = T))
  L_names <- list()
  for (i in names(L_save)) {
    x = paste(word(i, 1, sep = fixed('.')), word(i, 2, sep = fixed('.')), sep = ".")
    L_names[[i]] <- x
  }
  names(L_save) <- unlist(L_names)
  
  return(L_save)
  
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

###################### Benchmark MuSiC recursive mode ######################

BenchmarkMusicRecursive <- function(L_bulk, L_markers, L_clusters.type, sc.sce, clusters, samples) {
  
  L_final <- list()
  
  for (i in names(L_bulk)) {
    
    L_prop <- list()
    
    for (m in names(L_markers)) {
      if (m == "no_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with recursive mode, {i} counts and no markers"
          )
        )
      }
      if (m == "with_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running MuSiC with recursive mode, {i} counts and markers"
          )
        )
      }
      
      groups = paste0("clusterType_", m)
      
      est.prop <- music_prop.cluster(
        bulk.mtx = L_bulk[[i]], 
        sc.sce = sce, 
        clusters = clusters, 
        samples = samples,
        group.markers = L_markers[[m]], 
        groups = paste0("clusterType_", m), 
        clusters.type = L_clusters.type[[m]],
        verbose = F
      )
      
      names(est.prop)[1] <- paste(i,m,sep = "_")
      L_prop[[m]] <- est.prop[1]
    }
    L_final[[i]] <- L_prop
  }
  
  L_save <- do.call(c, unlist(L_final, recursive = F, use.names = T))
  L_names <- list()
  for (i in names(L_save)) {
    x = paste(word(i, 1, sep = fixed('.')), word(i, 2, sep = fixed('.')), sep = ".")
    L_names[[i]] <- x
    }
  names(L_save) <- unlist(L_names)
  
  return(L_save)

}

###################### BisqueRNA ######################
###################### Benchmark BisqueRNA (reference-based decomposition) ######################
BenchmarkBisqueRefBased <- function(L_bulk, sc.eset, L_markers) {
  
  L_final <- list()
  
  for (i in names(L_bulk)) {
    
    L_prop <- list()
    for (m in names(L_markers)) {
      if (m == "no_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running Bisque with {i} counts and no markers"
          )
        )
      }
      if (m == "with_markers") {
        message(
          glue::glue(
            "{format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')} Running Bisque with {i} counts and markers"
          )
        )
      }
      
      est.prop <- ReferenceBasedDecomposition(
        bulk.eset = L_bulk[[i]],
        sc.eset = sc.eset,
        markers = L_markers[[m]],
        cell.types = "cellType",
        subject.names = "SubjectName",
        use.overlap = F,
        verbose = F,
        old.cpm = TRUE
      )
    
      est.prop[1] <- list(t(est.prop[[1]])) # set samples as rows to be consistent with other tools
      names(est.prop)[1] <- paste(i,m,sep = "_")
      L_prop[[m]] <- est.prop[1]
    }
    L_final[[i]] <- L_prop
  }
  L_save <- do.call(c, unlist(L_final, recursive = F, use.names = T))
  L_names <- list()
  for (i in names(L_save)) {
    x = paste(word(i, 1, sep = fixed('.')), word(i, 2, sep = fixed('.')), sep = ".")
    L_names[[i]] <- x
  }
  names(L_save) <- unlist(L_names)
  
  return(L_save)
  
}

###################### CIBERSORTx ######################
###################### Prepare signature matrix ######################
SigMatrixCibersortx <- function(seurat, label) {
  sig_matrix <- as.matrix(GetAssayData(seurat, slot = "count"))
  colnames(sig_matrix) <- as.character(seurat@meta.data[[label]])
  sig_matrix <- sig_matrix[,order(seurat@meta.data[[label]],colnames(sig_matrix))] # reorder columns according to labels/levels
  sig_matrix <- cbind(rownames(sig_matrix), sig_matrix)
  colnames(sig_matrix)[1] <- "gene_name"
  return(sig_matrix)
}

###################### Prepare mixture matrix ######################
MixMatrixPB <- function(pb_counts) {
  mix_matrix <- pb_counts
  mix_matrix <- cbind(rownames(mix_matrix), mix_matrix)
  colnames(mix_matrix)[1] <- "gene_name"
  return(mix_matrix)
} 

###################### Benchmarking ######################
###################### Load results from different methods ######################
LoadDeconvolutionResults <- function(path) {
  files <-  paste0(path, list.files(path, pattern = "proportions"))
  L_prop <- lapply(files, function(x){
    res <- readRDS(x)
  })
  names(L_prop) <- gsub("_.*$", "", basename(files))
  return(L_prop)
}

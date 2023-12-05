


list_bulk <- list(bulk.mtx.summed, bulk.mtx.tpm)

sce = sce

markers = selected_markers
clusters = 'aggregated_lineage2'
samples = 'sample'

BenchmarkMusicDirect <- function(list_bulk, sc.sce, markers, clusters, samples, mode) {
  
  c(
    "all", 
    "summed_no_markers",
    "summed_with_markers",
    "tpm_no_markers",
    "tpm_with_markers"
    )
  
  message("Running MuSiC on summed counts & direct mode & without markers")
  
}

music_prop(
  bulk.mtx = bulk.mtx, 
  sc.sce = sce, 
  markers = NULL,
  clusters = clusters,
  samples = samples,
  verbose = F
)
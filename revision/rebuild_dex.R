# Recompute the manuscript-dependent chunks of the Dex notebook from raw counts.
# The app-only STRING network is independent of the paper's local enrichment.
source("revision/notebook_chunks.R")
root <- getwd()
notebook <- "analysis/03_chondrocytes_dexamethasone.Rmd"
env <- new.env(parent = globalenv())
chunks <- read_notebook_chunks(notebook)
setwd("analysis")
eval(chunks$bodyclocks_setup, env)
setwd(root)
annotation <- readRDS("revision/inputs/annotation/dexamethasone_gene_info.rds")
env$connect_biomart <- function(species_dataset) {
  stopifnot(species_dataset == "mmusculus_gene_ensembl")
  species_dataset
}
env$get_biomart_annotation <- function(ids, mart, attributes, filters, cache_key) {
  stopifnot(identical(attributes, annotation$attributes), identical(filters, annotation$filters),
            all(annotation$data$ensembl_gene_id %in% ids), cache_key == "dexamethasone_gene_info")
  annotation$data
}
for (label in c("libraries", "load_data", "deseq2_normalization", "run_rain_annotate",
                "prepare_heatmap", "prepare_timeseries_data", "compute_amplitude_display_table")) {
  message("Dex: ", label)
  eval(chunks[[label]], env)
}
stopifnot(identical(as.integer(env$measure_sequence), c(rep(3L,6),2L,rep(3L,5))),
          ncol(env$norm_counts) == 35L, sum(env$merged$BH.Q < 0.05) == 2043L)
manifest <- read.csv("revision/inputs/inputs.csv")
stopifnot(digest::digest(file = "results/tables/dex_circ.csv", algo = "sha256") ==
            manifest$sha256[manifest$dataset == "dex"])


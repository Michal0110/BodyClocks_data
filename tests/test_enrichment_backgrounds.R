#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path("R", "io_utils.R"))
source(file.path("R", "string_utils.R"))
source(file.path("R", "validation_utils.R"))

calls <- new.env(parent = emptyenv())
calls$backgrounds <- list()

fake_stringdb <- function(background = NULL, fail_custom = FALSE) {
  force(background)
  force(fail_custom)
  list(
    map = function(genes, column, removeUnmappedRows = TRUE) {
      data.frame(
        gene = genes[[column]],
        STRING_id = paste0("10090.", genes[[column]]),
        stringsAsFactors = FALSE
      )
    },
    get_proteins = function() {
      data.frame(
        protein_external_id = paste0("10090.", LETTERS[1:4]),
        stringsAsFactors = FALSE
      )
    },
    get_enrichment = function(string_ids, category = "All") {
      calls$backgrounds[length(calls$backgrounds) + 1L] <- list(background)
      stopifnot(identical(category, "All"))
      if (fail_custom && length(background)) {
        stop("deliberate custom-background failure")
      }
      data.frame(
        category = "KEGG",
        term = "mmu00010",
        description = "Example pathway",
        fdr = 0.01,
        preferredNames = I(list(sub("^10090[.]", "", string_ids))),
        stringsAsFactors = FALSE
      )
    }
  )
}

real_init_stringdb <- init_stringdb
init_stringdb <- function(species_id, version = "12.0", score_threshold = 700,
                          input_directory = tempdir(), backgroundV = NULL) {
  fake_stringdb(backgroundV)
}

dual <- get_string_enrichment_backgrounds(
  foreground_genes = c("A", "B"),
  tested_genes = c("A", "B", "C"),
  species_id = 10090
)

stopifnot(
  length(calls$backgrounds) == 2L,
  length(calls$backgrounds[[1L]]) == 0L,
  identical(calls$backgrounds[[2L]], paste0("10090.", c("A", "B", "C"))),
  dual$backgrounds$dataset_tested$metadata$mapped_foreground_count == 2L,
  dual$backgrounds$dataset_tested$metadata$mapped_background_count == 3L,
  dual$backgrounds$all_string$metadata$mapped_background_count == 4L
)

nodes <- data.frame(
  symbol = c("A", "B"), phase = c(4, 8), stringsAsFactors = FALSE
)
gene_info <- data.frame(
  preferred_symbol = c("A", "B"),
  description = c("Gene A", "Gene B"),
  stringsAsFactors = FALSE
)
bundle <- build_string_enrichment_background_bundle(
  dual,
  nodes_data = nodes,
  gene_info = gene_info,
  gene_mapping = list(A = "A", B = "B"),
  phase_colors = c("4" = "red", "8" = "blue")
)
stopifnot(validate_enrichment_background_bundle(bundle))

empty_dual <- empty_string_enrichment_backgrounds(
  species_id = 10090,
  all_string_background_count = 4L
)
empty_bundle <- build_string_enrichment_background_bundle(
  empty_dual,
  nodes_data = data.frame(
    symbol = character(), phase = numeric(), stringsAsFactors = FALSE
  ),
  gene_info = data.frame(
    preferred_symbol = character(), description = character(),
    stringsAsFactors = FALSE
  ),
  gene_mapping = list(),
  phase_colors = character()
)
stopifnot(
  validate_enrichment_background_bundle(empty_bundle),
  nrow(empty_bundle$backgrounds$dataset_tested$node_titles) == 0L,
  nrow(empty_bundle$backgrounds$all_string$cluster_descriptions) == 0L
)

go_terms <- data.frame(
  cluster = paste0("GO:", 1:6),
  cluster_description = paste("Term", 1:6),
  category = "Process",
  fdr = c(0.04, 0.01, 0.03, 0.02, 0.005, 0.025),
  stringsAsFactors = FALSE
)
tooltip <- format_string_enrichment_tooltip(go_terms)
stopifnot(
  grepl("Term 5, Term 2, Term 4, Term 6, Term 3 (and 1 more)",
        tooltip, fixed = TRUE),
  !grepl("Term 1", tooltip, fixed = TRUE)
)

calls$backgrounds <- list()
init_stringdb <- function(species_id, version = "12.0", score_threshold = 700,
                          input_directory = tempdir(), backgroundV = NULL) {
  fake_stringdb(backgroundV, fail_custom = TRUE)
}
custom_error <- try(
  get_string_enrichment_backgrounds(
    foreground_genes = "A", tested_genes = c("A", "B"), species_id = 10090
  ),
  silent = TRUE
)
stopifnot(
  inherits(custom_error, "try-error"),
  grepl("deliberate custom-background failure", custom_error, fixed = TRUE)
)

init_stringdb <- real_init_stringdb
cat("enrichment background data tests passed\n")

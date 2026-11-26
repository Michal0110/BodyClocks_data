# Optional reconstruction of the compact reference from the four pinned upstream files.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
params <- list(pathway_dir = normalizePath(args[1]))
refs <- read.delim("revision/inputs/pathway_references.tsv", header = FALSE)
for (i in seq_len(nrow(refs))) {
  stopifnot(digest::digest(file = file.path(params$pathway_dir, refs[i,1]), algo = "sha256") == refs[i,2])
}
source("revision/notebook_chunks.R")
chunks <- read_notebook_chunks("analysis/comparison_analysis_paper.rmd")
eval(chunks$load_libraries)
eval(chunks$helper_functions)
source("R/paper_pathway_annotations.R")
snapshot <- read_paper_pathway_annotations()
saveRDS(snapshot, "revision/inputs/pathway_annotations.rds", compress = "xz")


# Run from the repository root, after prepare_inputs.py and (for RAIN) rebuild_dex.R.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(0L, 2L))

# No arguments renders corrected RAIN. For CircaCompare pass the calls directory
# followed by a separate figure output directory.
cc <- length(args) == 2L
out <- if (cc) args[2] else "results/publication_figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)


# Both runs use the frozen pathway map. Export every eligible term, including zero-hit
# terms, so the comparison can distinguish a lost discovery from a missing test.
rmarkdown::render("analysis/comparison_analysis_paper.rmd",
  output_file = if (cc) "CircaCompare.html" else "corrected_RAIN.html",
  output_dir = normalizePath(out), intermediates_dir = tempdir(),
  output_options = list(code_folding = "hide"),
  params = list(circacompare_dir = if (cc) args[1] else "", expression_dir = "results/tables",
                pathway_snapshot = "revision/inputs/pathway_annotations.rds",
                figure_output = out, export_all_terms = TRUE))


# Split the paper's combined table by background for inspection.
overlap <- read.csv(file.path(out, "pathway_overlap_background_comparison.csv"))
for (bg in unique(overlap$Background_Strategy))
  write.csv(overlap[overlap$Background_Strategy == bg, ],
            file.path(out, paste0("pathway_overlap_", bg, ".csv")), row.names=FALSE)

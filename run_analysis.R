#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) flag %in% args
arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[[1]]) else NULL
}

dry_run <- has_flag("--dry-run")
fail_fast <- has_flag("--fail-fast")
subset_arg <- arg_value("--scripts")

if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is required. Activate the conda env and run renv::restore().",
       call. = FALSE)
}

status <- renv::status()
if (is.list(status) && isFALSE(status$synchronized)) {
  stop("renv is out of sync. Run renv::restore() first.", call. = FALSE)
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("rmarkdown is required. Run renv::restore().", call. = FALSE)
}

source(file.path("R", "io_utils.R"))
source(file.path("R", "validation_utils.R"))

scripts <- c(
  "01" = "analysis/01_cartilage.Rmd",
  "02" = "analysis/02_tendon.Rmd",
  "03" = "analysis/03_chondrocytes_dexamethasone.Rmd",
  "04" = "analysis/04_chondrocytes_heat_shock.Rmd",
  "05" = "analysis/05_chondrocytes_osmotic_stress.Rmd",
  "06" = "analysis/06_glomeruli.Rmd",
  "07" = "analysis/07_mammary_gland.Rmd",
  "08" = "analysis/08_nih3t3.Rmd",
  "09" = "analysis/09_podocytes.Rmd",
  "10" = "analysis/10_xiphoid.Rmd",
  "11" = "analysis/11_mouse_atlas.Rmd",
  "12" = "analysis/12_baboon_atlas.Rmd",
  "13" = "analysis/13_liver_rnaseq.Rmd",
  "paper" = "analysis/comparison_analysis_paper.rmd",
  "pairwise" = "analysis/pairwise_similarity_analysis.R"
)

script_dependencies <- list(
  # The manuscript comparison includes Figures 1E and 1F from the pairwise
  # similarity analysis, which requires outputs from every mouse-processing
  # script. Script 12 is the baboon atlas and is intentionally omitted.
  paper = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
            "11", "13")
)

script_post_dependencies <- list(
  paper = "pairwise"
)

expand_script_dependencies <- function(requested_ids) {
  expanded <- character(0)
  add_with_deps <- function(id) {
    deps <- script_dependencies[[id]]
    if (length(deps)) {
      for (dep in deps) add_with_deps(dep)
    }
    expanded <<- c(expanded, id)
  }
  for (id in requested_ids) {
    add_with_deps(id)
    post_deps <- script_post_dependencies[[id]]
    if (length(post_deps)) {
      for (post_dep in post_deps) add_with_deps(post_dep)
    }
  }
  unique(expanded)
}

if (!is.null(subset_arg)) {
  requested <- trimws(strsplit(subset_arg, ",", fixed = TRUE)[[1]])
  missing <- setdiff(requested, names(scripts))
  if (length(missing)) {
    stop("Unknown script id(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  scripts <- scripts[expand_script_dependencies(requested)]
}

ensure_dir(file.path("results", "reports"))
ensure_dir(file.path("results", "tables"))
ensure_dir(file.path("results", "shiny_data"))

log_path <- file.path("results", "run_analysis.log")
cat("Run started:", format(Sys.time()), "\n", file = log_path, append = TRUE)

if (dry_run) {
  missing <- scripts[!file.exists(scripts)]
  if (length(missing)) {
    stop("Missing analysis script(s):\n", paste(missing, collapse = "\n"),
         call. = FALSE)
  }
  validate_output_manifest()
  message("Dry run OK: scripts and output manifest are present.")
  quit(status = 0)
}

failures <- character()
run_script <- function(script) {
  extension <- tolower(tools::file_ext(script))
  if (identical(extension, "rmd")) {
    rmarkdown::render(script, output_dir = file.path("results", "reports"),
                      quiet = FALSE)
  } else if (identical(extension, "r")) {
    source(script, local = new.env(parent = globalenv()))
  } else {
    stop("Unsupported analysis script type: ", script, call. = FALSE)
  }
}

for (id in names(scripts)) {
  script <- scripts[[id]]
  message("[", Sys.time(), "] Running: ", script)
  result <- tryCatch(
    run_script(script),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    msg <- paste0(script, " FAILED: ", conditionMessage(result))
    failures <- c(failures, msg)
    cat(format(Sys.time()), msg, "\n", file = log_path, append = TRUE)
    if (fail_fast) stop(msg, call. = FALSE)
  } else {
    cat(format(Sys.time()), script, "OK\n", file = log_path, append = TRUE)
  }
}

if (length(failures)) {
  cat("Failures:\n", paste(failures, collapse = "\n"), "\n",
      file = log_path, append = TRUE)
  quit(status = 1)
}

cat("Run completed:", format(Sys.time()), "\n", file = log_path, append = TRUE)
message("Analysis complete.")

validate_output_manifest <- function(path = file.path("config", "output_manifest.csv")) {
  if (!file.exists(path)) {
    stop("Missing output manifest: ", path, call. = FALSE)
  }
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("dataset_id", "analysis_script", "species", "app_folder",
                "prefix", "output_type", "expected_file", "required", "notes")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Output manifest missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  manifest
}

validate_rds_exists <- function(manifest = validate_output_manifest(),
                                root = file.path("results", "shiny_data")) {
  required <- manifest[as.logical(manifest$required), , drop = FALSE]
  paths <- file.path(root, sub("^data/", "", required$app_folder),
                     required$expected_file)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing required RDS output(s):\n", paste(missing, collapse = "\n"),
         call. = FALSE)
  }
  invisible(paths)
}

validate_display_columns <- function(x) {
  required <- c("symbol", "description", "pVal", "BH.Q")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Display table missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  TRUE
}

validate_plot_data_columns <- function(x) {
  required <- c("Gene", "time_point", "mean", "BH.Q")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Plot data missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  TRUE
}

validate_network_columns <- function(nodes, edges) {
  if (!nrow(nodes) || !nrow(edges)) {
    return(TRUE)
  }
  if (!"id" %in% names(nodes)) {
    stop("Network nodes missing id column.", call. = FALSE)
  }
  edge_required <- c("from", "to")
  missing_edges <- setdiff(edge_required, names(edges))
  if (length(missing_edges)) {
    stop("Network edges missing columns: ", paste(missing_edges, collapse = ", "),
         call. = FALSE)
  }
  TRUE
}


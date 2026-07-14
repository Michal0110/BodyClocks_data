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
  network_datasets <- unique(manifest$dataset_id[
    manifest$output_type == "nodes"
  ])
  bundle_rows <- manifest[
    manifest$output_type == "enrichment_backgrounds", , drop = FALSE
  ]
  bundle_counts <- table(bundle_rows$dataset_id)
  missing_bundles <- setdiff(network_datasets, names(bundle_counts))
  duplicate_bundles <- names(bundle_counts)[bundle_counts != 1L]
  if (length(missing_bundles) || length(duplicate_bundles)) {
    details <- c(
      if (length(missing_bundles)) paste0(
        "missing for: ", paste(missing_bundles, collapse = ", ")
      ),
      if (length(duplicate_bundles)) paste0(
        "not unique for: ", paste(duplicate_bundles, collapse = ", ")
      )
    )
    stop("Each network dataset must have one enrichment-background bundle (",
         paste(details, collapse = "; "), ").", call. = FALSE)
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
  bundle_paths <- paths[required$output_type == "enrichment_backgrounds"]
  for (path in bundle_paths) {
    validate_enrichment_background_bundle(readRDS(path), path = path)
  }
  invisible(paths)
}

validate_enrichment_background_bundle <- function(x, path = NULL) {
  context <- if (is.null(path)) "Enrichment-background bundle" else path
  if (!is.list(x) || !identical(as.integer(x$schema_version), 1L)) {
    stop(context, " must use schema_version 1.", call. = FALSE)
  }
  if (!is.list(x$metadata)) {
    stop(context, " metadata must be a list.", call. = FALSE)
  }

  required_backgrounds <- c("dataset_tested", "all_string")
  required_metadata <- c(
    "species_id", "string_version", "categories", "fdr_cutoff",
    "foreground_input_count", "tested_input_count", "backgrounds"
  )
  missing_metadata <- setdiff(required_metadata, names(x$metadata))
  if (length(missing_metadata)) {
    stop(context, " metadata is missing: ",
         paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  expected_categories <- c("Process", "KEGG", "RCTM", "WikiPathways")
  if (length(x$metadata$species_id) != 1L ||
      is.na(x$metadata$species_id) ||
      length(x$metadata$string_version) != 1L ||
      is.na(x$metadata$string_version) ||
      !nzchar(as.character(x$metadata$string_version)) ||
      !setequal(as.character(x$metadata$categories), expected_categories) ||
      length(x$metadata$fdr_cutoff) != 1L ||
      is.na(x$metadata$fdr_cutoff) ||
      as.numeric(x$metadata$fdr_cutoff) != 0.05) {
    stop(context, " has invalid STRING enrichment provenance metadata.",
         call. = FALSE)
  }
  input_counts <- unlist(x$metadata[c(
    "foreground_input_count", "tested_input_count"
  )], use.names = FALSE)
  if (length(input_counts) != 2L || anyNA(input_counts) ||
      any(!is.finite(input_counts)) || any(input_counts < 0)) {
    stop(context, " has invalid input gene counts.", call. = FALSE)
  }
  if (!is.list(x$metadata$backgrounds)) {
    stop(context, " background metadata must be a list.", call. = FALSE)
  }
  missing_background_metadata <- setdiff(
    required_backgrounds, names(x$metadata$backgrounds)
  )
  if (length(missing_background_metadata)) {
    stop(context, " metadata is missing background(s): ",
         paste(missing_background_metadata, collapse = ", "), call. = FALSE)
  }
  if (!is.list(x$backgrounds)) {
    stop(context, " backgrounds must be a list.", call. = FALSE)
  }
  missing_backgrounds <- setdiff(required_backgrounds, names(x$backgrounds))
  if (length(missing_backgrounds)) {
    stop(context, " is missing background(s): ",
         paste(missing_backgrounds, collapse = ", "), call. = FALSE)
  }

  expected_columns <- list(
    node_titles = c("id", "title"),
    cluster_descriptions = c("category", "cluster", "description", "fdr"),
    cluster_to_genes = c("category", "cluster", "genes")
  )
  for (background in required_backgrounds) {
    background_metadata <- x$metadata$backgrounds[[background]]
    required_background_metadata <- c(
      "background", "mapped_foreground_count", "mapped_background_count"
    )
    if (!is.list(background_metadata) ||
        length(setdiff(required_background_metadata,
                       names(background_metadata)))) {
      stop(context, " has incomplete metadata for background '", background,
           "'.", call. = FALSE)
    }
    if (!identical(as.character(background_metadata$background), background)) {
      stop(context, " has mismatched metadata for background '", background,
           "'.", call. = FALSE)
    }
    mapped_foreground_count <- background_metadata$mapped_foreground_count
    mapped_background_count <- background_metadata$mapped_background_count
    if (length(mapped_foreground_count) != 1L ||
        is.na(mapped_foreground_count) ||
        !is.finite(mapped_foreground_count) || mapped_foreground_count < 0 ||
        length(mapped_background_count) != 1L ||
        is.na(mapped_background_count) ||
        !is.finite(mapped_background_count) || mapped_background_count < 0 ||
        mapped_foreground_count > mapped_background_count) {
      stop(context, " has invalid mapped counts for background '", background,
           "'.", call. = FALSE)
    }

    variant <- x$backgrounds[[background]]
    if (!is.list(variant)) {
      stop(context, " background '", background, "' must be a list.",
           call. = FALSE)
    }
    for (object_name in names(expected_columns)) {
      object <- variant[[object_name]]
      if (!is.data.frame(object) ||
          !identical(names(object), expected_columns[[object_name]])) {
        stop(context, " background '", background, "' has invalid ",
             object_name, " columns.", call. = FALSE)
      }
    }

    if (!is.character(variant$node_titles$id) ||
        !is.character(variant$node_titles$title)) {
      stop(context, " background '", background,
           "' node_titles columns must be character.", call. = FALSE)
    }
    if (anyNA(variant$node_titles$id) ||
        any(!nzchar(variant$node_titles$id)) ||
        anyDuplicated(variant$node_titles$id)) {
      stop(context, " background '", background,
           "' node IDs must be unique and nonempty.", call. = FALSE)
    }
    descriptions <- variant$cluster_descriptions
    if (!is.character(descriptions$category) ||
        !is.character(descriptions$cluster) ||
        !is.character(descriptions$description) ||
        !is.numeric(descriptions$fdr) ||
        anyNA(descriptions$category) ||
        any(!descriptions$category %in% expected_categories) ||
        anyNA(descriptions$cluster) || any(!nzchar(descriptions$cluster)) ||
        anyNA(descriptions$description) || anyNA(descriptions$fdr) ||
        any(!is.finite(descriptions$fdr)) ||
        any(descriptions$fdr < 0 | descriptions$fdr > x$metadata$fdr_cutoff)) {
      stop(context, " background '", background,
           "' cluster_descriptions has invalid column types.", call. = FALSE)
    }
    memberships <- variant$cluster_to_genes
    if (!is.character(memberships$category) ||
        !is.character(memberships$cluster) ||
        !is.list(memberships$genes) ||
        any(!vapply(memberships$genes, is.character, logical(1))) ||
        anyNA(memberships$category) ||
        any(!memberships$category %in% expected_categories) ||
        anyNA(memberships$cluster) || any(!nzchar(memberships$cluster))) {
      stop(context, " background '", background,
           "' cluster_to_genes has invalid column types.", call. = FALSE)
    }

    description_keys <- paste(
      descriptions$category, descriptions$cluster, sep = "\r"
    )
    membership_keys <- paste(
      memberships$category, memberships$cluster, sep = "\r"
    )
    if (anyDuplicated(description_keys) || anyDuplicated(membership_keys)) {
      stop(context, " background '", background,
           "' contains duplicate category/term keys.", call. = FALSE)
    }
    if (!setequal(description_keys, membership_keys)) {
      stop(context, " background '", background,
           "' term descriptions and memberships do not have matching keys.",
           call. = FALSE)
    }
    invalid_genes <- vapply(memberships$genes, function(genes) {
      !length(genes) || anyNA(genes) || any(!nzchar(genes)) ||
        anyDuplicated(genes) > 0L
    }, logical(1))
    if (any(invalid_genes)) {
      stop(context, " background '", background,
           "' gene memberships must be unique and nonempty.", call. = FALSE)
    }
  }
  TRUE
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

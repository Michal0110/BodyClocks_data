source_if_exists <- function(path) {
  if (file.exists(path)) source(path)
}

source_if_exists(file.path("R", "io_utils.R"))

init_stringdb <- function(species_id, version = "12.0", score_threshold = 700,
                          input_directory = file.path("results", "cache", "stringdb"),
                          backgroundV = NULL) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb is required for init_stringdb().", call. = FALSE)
  }
  ensure_dir(input_directory)
  tryCatch(
    if (is.null(backgroundV)) {
      STRINGdb::STRINGdb$new(
        version = version,
        species = species_id,
        score_threshold = score_threshold,
        input_directory = input_directory
      )
    } else {
      STRINGdb::STRINGdb$new(
        version = version,
        species = species_id,
        score_threshold = score_threshold,
        input_directory = input_directory,
        backgroundV = unique(as.character(backgroundV))
      )
    },
    error = function(e) {
      stop("Unable to initialise STRINGdb for species ", species_id,
           " version ", version, ": ", conditionMessage(e), call. = FALSE)
    }
  )
}

STRING_SCORE_THRESHOLDS <- c(700, 900)
DEFAULT_STRING_SCORE_THRESHOLD <- 700
STRING_ENRICHMENT_CATEGORIES <- c("Process", "KEGG", "RCTM", "WikiPathways")
STRING_ENRICHMENT_BACKGROUNDS <- c("dataset_tested", "all_string")
STRING_ENRICHMENT_SCHEMA_VERSION <- 1L

create_gene_mapping <- function(gene_ids = NULL, symbols = NULL, species_id = 10090) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("httr and jsonlite are required for create_gene_mapping().", call. = FALSE)
  }

  if (is.null(symbols)) {
    symbols <- gene_ids
  }
  symbols <- unique(stats::na.omit(as.character(symbols)))
  symbols <- symbols[nzchar(symbols)]

  mapping <- stats::setNames(as.list(symbols), symbols)
  if (!length(symbols)) {
    return(mapping)
  }

  symbol_mapping <- tryCatch({
    response <- httr::POST(
      "https://string-db.org/api/json/get_string_ids",
      body = list(
        identifiers = paste(symbols, collapse = "\n"),
        species = species_id,
        limit = 1,
        caller_identity = "bodyclocks"
      ),
      encode = "form"
    )
    if (httr::status_code(response) >= 300) {
      stop("STRINGdb mapping status ", httr::status_code(response))
    }
    jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  }, error = function(e) {
    stop("STRINGdb gene mapping failed: ", conditionMessage(e), call. = FALSE)
  })

  if (!is.null(symbol_mapping) && nrow(symbol_mapping) > 0) {
    for (i in seq_len(nrow(symbol_mapping))) {
      query <- symbol_mapping$queryItem[i]
      preferred <- symbol_mapping$preferredName[i]
      if (!is.na(query) && nzchar(query) &&
          !is.na(preferred) && nzchar(preferred)) {
        mapping[[preferred]] <- query
      }
    }
  }

  mapping
}

map_to_symbol <- function(name, mapping) {
  if (!is.null(mapping) && name %in% names(mapping)) {
    return(mapping[[name]])
  }
  name
}

map_to_mgi <- map_to_symbol

resolve_string_nodes_id_column <- function(nodes_data) {
  candidates <- c("Row.names", "ENSEMBL", "CycID", "gene", "ensembl_gene_id")
  hit <- candidates[candidates %in% names(nodes_data)]
  if (!length(hit)) {
    stop("Could not infer Ensembl/id column in nodes_data. Available columns: ",
         paste(names(nodes_data), collapse = ", "), call. = FALSE)
  }
  hit[[1]]
}

get_string_network <- function(gene_list, gene_mapping, string_db, nodes_data = NULL,
                               id_col = NULL, symbol_col = "symbol") {
  if (is.null(nodes_data) && exists("nodes_data", envir = parent.frame())) {
    nodes_data <- get("nodes_data", envir = parent.frame())
  }
  if (is.null(nodes_data)) {
    stop("nodes_data is required for get_string_network().", call. = FALSE)
  }
  if (!symbol_col %in% names(nodes_data)) {
    stop("nodes_data is missing symbol column '", symbol_col, "'.", call. = FALSE)
  }
  if (is.null(id_col)) {
    id_col <- resolve_string_nodes_id_column(nodes_data)
  }

  gene_list <- unique(stats::na.omit(as.character(gene_list)))
  if (!length(gene_list)) {
    return(data.frame(
      from = character(),
      to = character(),
      combined_score = numeric()
    ))
  }

  tryCatch({
    gene_df <- data.frame(gene = gene_list, stringsAsFactors = FALSE)
    cat("Mapping", length(gene_list), "genes to STRING database...\n")
    mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

    cat("Retrieving interactions for", nrow(mapped_genes), "mapped genes...\n")
    interactions <- string_db$get_interactions(mapped_genes$STRING_id)

    id_to_ensembl <- stats::setNames(mapped_genes$gene, mapped_genes$STRING_id)
    ensembl_to_symbol <- stats::setNames(nodes_data[[symbol_col]], nodes_data[[id_col]])

    resolve <- function(string_id) {
      ensembl_id <- id_to_ensembl[[string_id]]
      if (!is.null(ensembl_id) && ensembl_id %in% names(ensembl_to_symbol)) {
        return(ensembl_to_symbol[[ensembl_id]])
      }
      string_id
    }

    edges <- data.frame(
      from = vapply(interactions$from, resolve, character(1)),
      to = vapply(interactions$to, resolve, character(1)),
      combined_score = as.numeric(interactions$combined_score),
      stringsAsFactors = FALSE
    )

    unresolved <- unique(c(edges$from, edges$to))
    unresolved <- unresolved[grepl("^[0-9]+[.]ENS", unresolved)]
    if (length(unresolved) > 0) {
      stop(
        "STRINGdb network retrieval returned unresolved STRING protein IDs. ",
        "This usually means gene_list identifiers do not match the node IDs. ",
        "Use the same identifier type for gene_list and nodes_data[[id_col]], ",
        "or pass id_col explicitly. Example unresolved ID: ",
        unresolved[[1]],
        call. = FALSE
      )
    }

    cat("Retrieved", nrow(edges), "interactions between",
        length(unique(c(edges$from, edges$to))), "proteins\n")
    edges
  }, error = function(e) {
    stop("STRINGdb network retrieval failed: ", conditionMessage(e), call. = FALSE)
  })
}

get_string_network_bundle <- function(gene_list, gene_mapping, species_id,
                                      nodes_data = NULL,
                                      thresholds = STRING_SCORE_THRESHOLDS,
                                      default_threshold = DEFAULT_STRING_SCORE_THRESHOLD,
                                      version = "12.0",
                                      input_directory = file.path("results", "cache", "stringdb"),
                                      id_col = NULL,
                                      symbol_col = "symbol",
                                      include_enrichment = TRUE) {
  thresholds <- unique(as.integer(thresholds))
  default_threshold <- as.integer(default_threshold)
  if (!default_threshold %in% thresholds) {
    thresholds <- c(default_threshold, thresholds)
  }

  edges_by_threshold <- list()
  default_string_db <- NULL

  for (threshold in thresholds) {
    cat("Building STRING network at score_threshold =", threshold, "\n")
    string_db <- init_stringdb(
      species_id = species_id,
      version = version,
      score_threshold = threshold,
      input_directory = input_directory
    )
    edges_by_threshold[[as.character(threshold)]] <- get_string_network(
      gene_list = gene_list,
      gene_mapping = gene_mapping,
      string_db = string_db,
      nodes_data = nodes_data,
      id_col = id_col,
      symbol_col = symbol_col
    )
    if (threshold == default_threshold) {
      default_string_db <- string_db
    }
  }

  if (is.null(default_string_db)) {
    default_string_db <- init_stringdb(
      species_id = species_id,
      version = version,
      score_threshold = default_threshold,
      input_directory = input_directory
    )
  }

  list(
    edges_by_threshold = edges_by_threshold,
    edges = edges_by_threshold[[as.character(default_threshold)]],
    cluster_data = if (isTRUE(include_enrichment)) {
      get_string_enrichment(gene_list, default_string_db)
    } else {
      empty_string_enrichment()
    }
  )
}

empty_string_edges_by_threshold <- function(thresholds = STRING_SCORE_THRESHOLDS) {
  thresholds <- unique(as.integer(thresholds))
  stats::setNames(
    lapply(thresholds, function(x) data.frame(
      from = character(),
      to = character(),
      combined_score = numeric()
    )),
    as.character(thresholds)
  )
}

empty_string_enrichment <- function() {
  data.frame(
    category = character(),
    term = character(),
    description = character(),
    fdr = numeric(),
    preferredNames = I(list()),
    stringsAsFactors = FALSE
  )
}

normalise_string_enrichment <- function(result,
                                        categories = STRING_ENRICHMENT_CATEGORIES,
                                        fdr_cutoff = 0.05) {
  if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
    return(empty_string_enrichment())
  }

  required <- c("category", "term", "description", "fdr", "preferredNames")
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop("STRING enrichment result is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  result$category <- as.character(result$category)
  result$term <- as.character(result$term)
  result$description <- as.character(result$description)
  result$fdr <- as.numeric(result$fdr)
  keep <- !is.na(result$fdr) & result$fdr <= fdr_cutoff &
    result$category %in% categories
  result <- result[keep, , drop = FALSE]
  if (!nrow(result)) {
    return(empty_string_enrichment())
  }

  result <- result[order(result$fdr, result$category, result$description,
                         result$term, na.last = TRUE), , drop = FALSE]
  result <- result[!duplicated(paste(result$category, result$term, sep = "\r")),
                   , drop = FALSE]
  rownames(result) <- NULL
  result
}

map_genes_to_string_ids <- function(gene_list, string_db, label = "genes") {
  gene_list <- unique(stats::na.omit(trimws(as.character(gene_list))))
  gene_list <- gene_list[nzchar(gene_list)]
  if (!length(gene_list)) {
    return(character())
  }

  mapped <- tryCatch(
    string_db$map(
      data.frame(gene = gene_list, stringsAsFactors = FALSE),
      "gene",
      removeUnmappedRows = TRUE
    ),
    error = function(e) {
      stop("Unable to map ", label, " to STRING: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  if (!"STRING_id" %in% names(mapped)) {
    stop("STRING mapping for ", label, " did not return STRING_id.",
         call. = FALSE)
  }
  unique(stats::na.omit(as.character(mapped$STRING_id)))
}

query_string_enrichment <- function(string_ids, string_db,
                                    categories = STRING_ENRICHMENT_CATEGORIES,
                                    fdr_cutoff = 0.05,
                                    background_label = "all_string") {
  string_ids <- unique(stats::na.omit(as.character(string_ids)))
  if (!length(string_ids)) {
    return(empty_string_enrichment())
  }

  result <- tryCatch(
    string_db$get_enrichment(string_ids, category = "All"),
    error = function(e) {
      stop("STRING enrichment failed for background '", background_label,
           "': ", conditionMessage(e), call. = FALSE)
    }
  )
  normalise_string_enrichment(result, categories, fdr_cutoff)
}

get_string_enrichment <- function(gene_list, string_db,
                                  categories = STRING_ENRICHMENT_CATEGORIES,
                                  fdr_cutoff = 0.05) {
  gene_list <- unique(stats::na.omit(as.character(gene_list)))
  if (!length(gene_list)) {
    return(empty_string_enrichment())
  }

  tryCatch({
    cat("Mapping", length(gene_list), "genes to STRING database for enrichment...\n")
    mapped_ids <- map_genes_to_string_ids(gene_list, string_db,
                                           label = "enrichment genes")

    cat("Performing enrichment analysis...\n")
    query_string_enrichment(
      mapped_ids,
      string_db,
      categories = categories,
      fdr_cutoff = fdr_cutoff,
      background_label = "all_string"
    )
  }, error = function(e) {
    warning("STRINGdb enrichment failed: ", conditionMessage(e), call. = FALSE)
    empty_string_enrichment()
  })
}

get_all_string_protein_ids <- function(string_db) {
  proteins <- tryCatch(
    string_db$get_proteins(),
    error = function(e) {
      stop("Unable to obtain the STRING protein universe: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  id_candidates <- c("STRING_id", "protein_external_id", "string_protein_id")
  id_col <- id_candidates[id_candidates %in% names(proteins)]
  if (!length(id_col)) {
    stop("STRING protein table has no recognised protein identifier column.",
         call. = FALSE)
  }
  ids <- unique(stats::na.omit(as.character(proteins[[id_col[[1]]]])))
  if (!length(ids)) {
    stop("STRING returned an empty protein universe.", call. = FALSE)
  }
  ids
}

empty_string_enrichment_backgrounds <- function(
    species_id,
    all_string_background_count,
    dataset_tested_background_count = 0L,
    version = "12.0",
    categories = STRING_ENRICHMENT_CATEGORIES,
    fdr_cutoff = 0.05,
    foreground_input_count = 0L,
    tested_input_count = 0L) {
  background_counts <- list(
    dataset_tested = dataset_tested_background_count,
    all_string = all_string_background_count
  )
  if (any(lengths(background_counts) != 1L) ||
      any(vapply(background_counts, is.na, logical(1))) ||
      any(!vapply(background_counts, is.finite, logical(1))) ||
      any(vapply(background_counts, function(x) x < 0, logical(1)))) {
    stop("Mapped background counts must be finite non-negative values.",
         call. = FALSE)
  }
  empty_variant <- function(background) {
    list(
      cluster_data = empty_string_enrichment(),
      metadata = list(
        background = background,
        mapped_foreground_count = 0L,
        mapped_background_count = as.integer(background_counts[[background]])
      )
    )
  }
  list(
    schema_version = STRING_ENRICHMENT_SCHEMA_VERSION,
    metadata = list(
      species_id = as.integer(species_id),
      string_version = as.character(version),
      categories = as.character(categories),
      fdr_cutoff = as.numeric(fdr_cutoff),
      foreground_input_count = as.integer(foreground_input_count),
      tested_input_count = as.integer(tested_input_count)
    ),
    backgrounds = stats::setNames(
      lapply(STRING_ENRICHMENT_BACKGROUNDS, empty_variant),
      STRING_ENRICHMENT_BACKGROUNDS
    )
  )
}

get_string_enrichment_backgrounds <- function(
    foreground_genes,
    tested_genes,
    species_id,
    version = "12.0",
    input_directory = file.path("results", "cache", "stringdb"),
    score_threshold = DEFAULT_STRING_SCORE_THRESHOLD,
    categories = STRING_ENRICHMENT_CATEGORIES,
    fdr_cutoff = 0.05) {
  foreground_genes <- unique(stats::na.omit(trimws(as.character(foreground_genes))))
  foreground_genes <- foreground_genes[nzchar(foreground_genes)]
  tested_genes <- unique(stats::na.omit(trimws(as.character(tested_genes))))
  tested_genes <- tested_genes[nzchar(tested_genes)]

  categories <- unique(as.character(categories))
  unexpected_categories <- setdiff(categories, STRING_ENRICHMENT_CATEGORIES)
  if (length(unexpected_categories)) {
    stop("Unsupported STRING enrichment categories: ",
         paste(unexpected_categories, collapse = ", "), call. = FALSE)
  }
  if (!length(tested_genes)) {
    stop("tested_genes must contain the genes eligible for rhythmicity testing.",
         call. = FALSE)
  }

  # This object is enrichment-only. It is deliberately not reused by the PPI
  # network code, whose 700/900 score threshold objects remain independent.
  all_string_db <- init_stringdb(
    species_id = species_id,
    version = version,
    score_threshold = score_threshold,
    input_directory = input_directory
  )
  mapped_foreground <- map_genes_to_string_ids(
    foreground_genes, all_string_db, label = "rhythmic foreground genes"
  )
  mapped_tested <- map_genes_to_string_ids(
    tested_genes, all_string_db, label = "tested-gene background"
  )
  if (!length(mapped_tested)) {
    stop("No tested background genes mapped to STRING for species ", species_id,
         ".", call. = FALSE)
  }

  all_string_ids <- get_all_string_protein_ids(all_string_db)
  all_string_foreground <- intersect(mapped_foreground, all_string_ids)
  dataset_foreground <- intersect(mapped_foreground, mapped_tested)

  all_string_clusters <- query_string_enrichment(
    all_string_foreground,
    all_string_db,
    categories = categories,
    fdr_cutoff = fdr_cutoff,
    background_label = "all_string"
  )

  dataset_clusters <- empty_string_enrichment()
  if (length(dataset_foreground)) {
    # Construct a second STRINGdb instance with an explicit background. Errors
    # are fatal: silently substituting the species-wide background would make
    # the two UI choices analytically indistinguishable.
    dataset_string_db <- init_stringdb(
      species_id = species_id,
      version = version,
      score_threshold = score_threshold,
      input_directory = input_directory,
      backgroundV = mapped_tested
    )
    dataset_clusters <- query_string_enrichment(
      dataset_foreground,
      dataset_string_db,
      categories = categories,
      fdr_cutoff = fdr_cutoff,
      background_label = "dataset_tested"
    )
  }

  list(
    schema_version = STRING_ENRICHMENT_SCHEMA_VERSION,
    metadata = list(
      species_id = as.integer(species_id),
      string_version = as.character(version),
      categories = categories,
      fdr_cutoff = as.numeric(fdr_cutoff),
      foreground_input_count = as.integer(length(foreground_genes)),
      tested_input_count = as.integer(length(tested_genes))
    ),
    backgrounds = list(
      dataset_tested = list(
        cluster_data = dataset_clusters,
        metadata = list(
          background = "dataset_tested",
          mapped_foreground_count = as.integer(length(dataset_foreground)),
          mapped_background_count = as.integer(length(mapped_tested))
        )
      ),
      all_string = list(
        cluster_data = all_string_clusters,
        metadata = list(
          background = "all_string",
          mapped_foreground_count = as.integer(length(all_string_foreground)),
          mapped_background_count = as.integer(length(all_string_ids))
        )
      )
    )
  )
}

map_enrichment_genes <- function(cluster_data, gene_mapping) {
  if (!nrow(cluster_data) || !"preferredNames" %in% names(cluster_data)) {
    return(cluster_data)
  }

  mapped_cluster_data <- cluster_data
  for (i in seq_len(nrow(mapped_cluster_data))) {
    genes <- mapped_cluster_data$preferredNames[[i]]
    if (is.null(genes)) {
      next
    }
    if (length(genes) == 1 && grepl(",", genes)) {
      genes <- trimws(strsplit(genes, ",")[[1]])
    }
    mapped_cluster_data$preferredNames[[i]] <- paste(
      vapply(genes, map_to_symbol, character(1), mapping = gene_mapping),
      collapse = ","
    )
  }

  mapped_genes <- unique(unlist(strsplit(as.character(mapped_cluster_data$preferredNames), ",")))
  cat("Mapped", length(mapped_genes), "unique genes in enrichment results\n")
  mapped_cluster_data
}

split_preferred_names <- function(preferred_names) {
  lapply(preferred_names, function(gene_str) {
    if (is.character(gene_str) && length(gene_str) == 1) {
      return(trimws(strsplit(gene_str, ",")[[1]]))
    }
    gene_str
  })
}

string_enrichment_category_label <- function(category) {
  labels <- c(
    Process = "GO Biological Process",
    KEGG = "KEGG Pathway",
    RCTM = "Reactome",
    WikiPathways = "WikiPathways"
  )
  if (category %in% names(labels)) unname(labels[[category]]) else category
}

format_string_enrichment_tooltip <- function(cluster_info, go_term_limit = 5L) {
  if (is.null(cluster_info) || !is.data.frame(cluster_info) || !nrow(cluster_info)) {
    return("")
  }
  required <- c("cluster", "cluster_description", "category", "fdr")
  missing <- setdiff(required, names(cluster_info))
  if (length(missing)) {
    stop("Tooltip enrichment data is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  cluster_info <- cluster_info[
    !duplicated(paste(cluster_info$category, cluster_info$cluster, sep = "\r")),
    , drop = FALSE
  ]
  cluster_info <- cluster_info[
    order(is.na(cluster_info$fdr), cluster_info$fdr,
          cluster_info$cluster_description, cluster_info$cluster),
    , drop = FALSE
  ]
  categories <- c(
    STRING_ENRICHMENT_CATEGORIES,
    setdiff(unique(as.character(cluster_info$category)),
            STRING_ENRICHMENT_CATEGORIES)
  )
  categories <- categories[categories %in% cluster_info$category]

  sections <- vapply(categories, function(category) {
    category_terms <- cluster_info[cluster_info$category == category,
                                   , drop = FALSE]
    n_more <- 0L
    if (category == "Process" && nrow(category_terms) > go_term_limit) {
      n_more <- nrow(category_terms) - as.integer(go_term_limit)
      category_terms <- category_terms[seq_len(go_term_limit), , drop = FALSE]
    }
    term_text <- paste(category_terms$cluster_description, collapse = ", ")
    if (n_more > 0L) {
      term_text <- paste0(term_text, " (and ", n_more, " more)")
    }
    paste0("<b>", string_enrichment_category_label(category), ":</b> ",
           term_text)
  }, character(1))
  paste(sections[nzchar(sections)], collapse = "<br>")
}

build_string_network_outputs <- function(cluster_data, nodes_data, gene_info,
                                         gene_mapping, phase_colors,
                                         gene_info_symbol_col = "preferred_symbol",
                                         phase_shift = 0) {
  required_node_cols <- c("symbol", "phase")
  missing_node_cols <- setdiff(required_node_cols, names(nodes_data))
  if (length(missing_node_cols)) {
    stop("nodes_data is missing required columns: ",
         paste(missing_node_cols, collapse = ", "), call. = FALSE)
  }
  if (!gene_info_symbol_col %in% names(gene_info)) {
    stop("gene_info is missing join column '", gene_info_symbol_col, "'.",
         call. = FALSE)
  }
  if (!"description" %in% names(gene_info)) {
    stop("gene_info is missing required column 'description'.", call. = FALSE)
  }
  required_cluster_cols <- c("term", "description", "category", "preferredNames")
  missing_cluster_cols <- setdiff(required_cluster_cols, names(cluster_data))
  if (length(missing_cluster_cols)) {
    stop("cluster_data is missing required columns: ",
         paste(missing_cluster_cols, collapse = ", "), call. = FALSE)
  }
  if (!"fdr" %in% names(cluster_data)) {
    cluster_data$fdr <- NA_real_
  }

  cluster_data <- map_enrichment_genes(cluster_data, gene_mapping)
  cluster_data$preferredNames <- split_preferred_names(cluster_data$preferredNames)

  cluster_descriptions <- cluster_data %>%
    dplyr::select(cluster = term, cluster_description = description, category, fdr) %>%
    dplyr::distinct(category, cluster, .keep_all = TRUE)

  if (nrow(cluster_data)) {
    reformatted_clusters <- cluster_data %>%
      dplyr::select(cluster = term, cluster_description = description,
                    id = preferredNames, category, fdr) %>%
      tidyr::unnest(id) %>%
      dplyr::filter(!is.na(id), nzchar(id)) %>%
      dplyr::distinct(category, cluster, id, .keep_all = TRUE)
  } else {
    reformatted_clusters <- data.frame(
      cluster = character(),
      cluster_description = character(),
      id = character(),
      category = character(),
      fdr = numeric(),
      stringsAsFactors = FALSE
    )
  }

  nodes_data2 <- data.frame(
    id = nodes_data$symbol,
    label = nodes_data$symbol,
    phase = nodes_data$phase,
    stringsAsFactors = FALSE
  )
  if (phase_shift != 0) {
    nodes_data2$phase <- (nodes_data2$phase + phase_shift) %% 24
  }

  nodes <- nodes_data2 %>%
    dplyr::mutate(color = phase_colors[as.character(phase)]) %>%
    dplyr::select(-phase) %>%
    dplyr::left_join(reformatted_clusters, by = "id") %>%
    dplyr::group_by(id, label, color) %>%
    dplyr::summarize(
      clusters = list(unique(cluster)),
      cluster_description = list(unique(cluster_description)),
      .groups = "drop"
    )

  merged_nodes <- merge(nodes, gene_info, by.x = "label",
                        by.y = gene_info_symbol_col, all.x = TRUE)
  merged_nodes <- merged_nodes %>%
    dplyr::filter(!duplicated(label))

  node_tooltips <- data.frame(
    id = unique(as.character(nodes_data2$id)),
    stringsAsFactors = FALSE
  )
  node_tooltips$formatted_tooltip <- vapply(node_tooltips$id, function(node_id) {
    cluster_info <- reformatted_clusters[
      reformatted_clusters$id == node_id,
      c("cluster", "cluster_description", "category", "fdr"),
      drop = FALSE
    ]
    format_string_enrichment_tooltip(cluster_info, go_term_limit = 5L)
  }, character(1))

  nodes_with_tooltips <- merged_nodes %>%
    dplyr::left_join(node_tooltips, by = "id") %>%
    dplyr::mutate(
      formatted_tooltip = dplyr::coalesce(formatted_tooltip, ""),
      description = dplyr::coalesce(as.character(description), ""),
      title = paste0(
        "<b>", label, "</b> - ", description,
        ifelse(nchar(formatted_tooltip) > 0,
               paste0("<br><br>", formatted_tooltip),
               "")
      )
    ) %>%
    dplyr::select(-c(cluster_description, formatted_tooltip, clusters))

  list(
    cluster_data = cluster_data,
    nodes = nodes %>% dplyr::filter(!is.na(id)),
    nodes_with_tooltips = nodes_with_tooltips %>% dplyr::filter(!is.na(id)),
    cluster_for_plot = cluster_data %>%
      dplyr::select(cluster = term, genes = preferredNames),
    cluster_descriptions = cluster_descriptions,
    cluster_to_genes = cluster_data %>%
      dplyr::select(category, cluster = term, genes = preferredNames) %>%
      dplyr::distinct(category, cluster, .keep_all = TRUE)
  )
}

empty_string_enrichment_variant <- function() {
  list(
    node_titles = data.frame(
      id = character(), title = character(), stringsAsFactors = FALSE
    ),
    cluster_descriptions = data.frame(
      category = character(), cluster = character(), description = character(),
      fdr = numeric(), stringsAsFactors = FALSE
    ),
    cluster_to_genes = data.frame(
      category = character(), cluster = character(), genes = I(list()),
      stringsAsFactors = FALSE
    )
  )
}

build_string_enrichment_background_bundle <- function(
    enrichment_backgrounds,
    nodes_data,
    gene_info,
    gene_mapping,
    phase_colors,
    gene_info_symbol_col = "preferred_symbol",
    phase_shift = 0) {
  if (!is.list(enrichment_backgrounds) ||
      !identical(as.integer(enrichment_backgrounds$schema_version),
                 STRING_ENRICHMENT_SCHEMA_VERSION)) {
    stop("enrichment_backgrounds has an unsupported schema version.",
         call. = FALSE)
  }
  missing_backgrounds <- setdiff(
    STRING_ENRICHMENT_BACKGROUNDS,
    names(enrichment_backgrounds$backgrounds)
  )
  if (length(missing_backgrounds)) {
    stop("enrichment_backgrounds is missing: ",
         paste(missing_backgrounds, collapse = ", "), call. = FALSE)
  }

  variants <- stats::setNames(lapply(STRING_ENRICHMENT_BACKGROUNDS, function(background) {
    cluster_data <- enrichment_backgrounds$backgrounds[[background]]$cluster_data
    if (is.null(cluster_data)) {
      cluster_data <- empty_string_enrichment()
    }
    outputs <- build_string_network_outputs(
      cluster_data = cluster_data,
      nodes_data = nodes_data,
      gene_info = gene_info,
      gene_mapping = gene_mapping,
      phase_colors = phase_colors,
      gene_info_symbol_col = gene_info_symbol_col,
      phase_shift = phase_shift
    )

    node_titles <- outputs$nodes_with_tooltips %>%
      dplyr::transmute(id = as.character(id), title = as.character(title)) %>%
      dplyr::distinct(id, .keep_all = TRUE)
    cluster_descriptions <- outputs$cluster_descriptions %>%
      dplyr::transmute(
        category = as.character(category),
        cluster = as.character(cluster),
        description = as.character(cluster_description),
        fdr = as.numeric(fdr)
      ) %>%
      dplyr::distinct(category, cluster, .keep_all = TRUE) %>%
      dplyr::arrange(fdr, description, cluster)
    cluster_to_genes <- outputs$cluster_to_genes %>%
      dplyr::mutate(
        category = as.character(category),
        cluster = as.character(cluster),
        genes = lapply(genes, function(x) {
          unique(stats::na.omit(as.character(x[nzchar(as.character(x))])))
        })
      ) %>%
      dplyr::distinct(category, cluster, .keep_all = TRUE)

    list(
      node_titles = as.data.frame(node_titles, stringsAsFactors = FALSE),
      cluster_descriptions = as.data.frame(
        cluster_descriptions, stringsAsFactors = FALSE
      ),
      cluster_to_genes = as.data.frame(
        cluster_to_genes, stringsAsFactors = FALSE
      )
    )
  }), STRING_ENRICHMENT_BACKGROUNDS)

  background_metadata <- stats::setNames(
    lapply(STRING_ENRICHMENT_BACKGROUNDS, function(background) {
      enrichment_backgrounds$backgrounds[[background]]$metadata
    }),
    STRING_ENRICHMENT_BACKGROUNDS
  )
  metadata <- enrichment_backgrounds$metadata
  metadata$backgrounds <- background_metadata

  bundle <- list(
    schema_version = STRING_ENRICHMENT_SCHEMA_VERSION,
    metadata = metadata,
    backgrounds = variants
  )
  validate_string_enrichment_background_bundle(bundle)
  bundle
}

validate_string_enrichment_background_bundle <- function(bundle) {
  if (!is.list(bundle) ||
      !identical(as.integer(bundle$schema_version),
                 STRING_ENRICHMENT_SCHEMA_VERSION)) {
    stop("Enrichment-background bundle has an unsupported schema version.",
         call. = FALSE)
  }
  if (!is.list(bundle$metadata)) {
    stop("Enrichment-background bundle metadata must be a list.", call. = FALSE)
  }
  required_metadata <- c(
    "species_id", "string_version", "categories", "fdr_cutoff",
    "foreground_input_count", "tested_input_count", "backgrounds"
  )
  missing_metadata <- setdiff(required_metadata, names(bundle$metadata))
  if (length(missing_metadata)) {
    stop("Enrichment-background bundle metadata is missing: ",
         paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (!setequal(as.character(bundle$metadata$categories),
                STRING_ENRICHMENT_CATEGORIES) ||
      length(bundle$metadata$fdr_cutoff) != 1L ||
      is.na(bundle$metadata$fdr_cutoff) ||
      bundle$metadata$fdr_cutoff < 0 || bundle$metadata$fdr_cutoff > 1 ||
      !is.list(bundle$metadata$backgrounds)) {
    stop("Enrichment-background bundle has invalid provenance metadata.",
         call. = FALSE)
  }
  missing_backgrounds <- setdiff(
    STRING_ENRICHMENT_BACKGROUNDS, names(bundle$backgrounds)
  )
  if (length(missing_backgrounds)) {
    stop("Enrichment-background bundle is missing: ",
         paste(missing_backgrounds, collapse = ", "), call. = FALSE)
  }

  expected_columns <- list(
    node_titles = c("id", "title"),
    cluster_descriptions = c("category", "cluster", "description", "fdr"),
    cluster_to_genes = c("category", "cluster", "genes")
  )
  for (background in STRING_ENRICHMENT_BACKGROUNDS) {
    background_metadata <- bundle$metadata$backgrounds[[background]]
    required_background_metadata <- c(
      "background", "mapped_foreground_count", "mapped_background_count"
    )
    if (!is.list(background_metadata) ||
        length(setdiff(required_background_metadata,
                       names(background_metadata))) ||
        !identical(as.character(background_metadata$background), background)) {
      stop("Bundle background '", background,
           "' has invalid provenance metadata.", call. = FALSE)
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
      stop("Bundle background '", background,
           "' has invalid mapped gene counts.", call. = FALSE)
    }

    variant <- bundle$backgrounds[[background]]
    for (object_name in names(expected_columns)) {
      object <- variant[[object_name]]
      if (!is.data.frame(object) ||
          !identical(names(object), expected_columns[[object_name]])) {
        stop("Bundle background '", background, "' has invalid ",
             object_name, " columns.", call. = FALSE)
      }
    }

    if (!is.character(variant$node_titles$id) ||
        !is.character(variant$node_titles$title) ||
        anyNA(variant$node_titles$id) ||
        any(!nzchar(variant$node_titles$id)) ||
        anyNA(variant$node_titles$title) ||
        anyDuplicated(variant$node_titles$id)) {
      stop("node_titles columns must be character vectors for background '",
           background, "'.", call. = FALSE)
    }
    descriptions <- variant$cluster_descriptions
    if (!is.character(descriptions$category) ||
        !is.character(descriptions$cluster) ||
        !is.character(descriptions$description) ||
        !is.numeric(descriptions$fdr) ||
        anyNA(descriptions$category) ||
        any(!descriptions$category %in% STRING_ENRICHMENT_CATEGORIES) ||
        anyNA(descriptions$cluster) || any(!nzchar(descriptions$cluster)) ||
        anyNA(descriptions$description) || anyNA(descriptions$fdr) ||
        any(!is.finite(descriptions$fdr)) ||
        any(descriptions$fdr < 0 | descriptions$fdr > bundle$metadata$fdr_cutoff)) {
      stop("cluster_descriptions has invalid column types for background '",
           background, "'.", call. = FALSE)
    }
    memberships <- variant$cluster_to_genes
    if (!is.character(memberships$category) ||
        !is.character(memberships$cluster) ||
        !is.list(memberships$genes) ||
        any(!vapply(memberships$genes, is.character, logical(1))) ||
        anyNA(memberships$category) ||
        any(!memberships$category %in% STRING_ENRICHMENT_CATEGORIES) ||
        anyNA(memberships$cluster) || any(!nzchar(memberships$cluster))) {
      stop("cluster_to_genes has invalid column types for background '",
           background, "'.", call. = FALSE)
    }
    description_keys <- paste(
      descriptions$category, descriptions$cluster, sep = "\r"
    )
    membership_keys <- paste(
      memberships$category, memberships$cluster, sep = "\r"
    )
    if (anyDuplicated(description_keys) || anyDuplicated(membership_keys)) {
      stop("Term keys must be unique within background '", background, "'.",
           call. = FALSE)
    }
    if (!setequal(description_keys, membership_keys)) {
      stop("Term descriptions and memberships must use the same keys within ",
           "background '", background, "'.", call. = FALSE)
    }
    invalid_memberships <- vapply(memberships$genes, function(genes) {
      !length(genes) || anyNA(genes) || any(!nzchar(genes)) ||
        anyDuplicated(genes) > 0L
    }, logical(1))
    if (any(invalid_memberships)) {
      stop("Gene memberships must be unique and nonempty within background '",
           background, "'.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

string_enrichment_backgrounds_path <- function(output_dir, prefix) {
  file.path(output_dir, paste0(prefix, "_enrichment_backgrounds.rds"))
}

save_string_enrichment_background_bundle <- function(bundle, output_dir, prefix) {
  validate_string_enrichment_background_bundle(bundle)
  ensure_dir(output_dir)
  path <- string_enrichment_backgrounds_path(output_dir, prefix)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = output_dir
  )
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  saveRDS(bundle, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Failed to atomically save enrichment-background bundle: ", path,
         call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Enrichment-background bundle was not written: ", path,
         call. = FALSE)
  }
  invisible(path)
}

save_string_edges_by_threshold <- function(edges, output_dir, prefix,
                                           default_threshold = DEFAULT_STRING_SCORE_THRESHOLD) {
  ensure_dir(output_dir)
  if (is.list(edges) && !is.data.frame(edges)) {
    for (threshold in names(edges)) {
      save_rds_checked(
        edges[[threshold]],
        file = file.path(output_dir, paste0(prefix, "_edges_string", threshold, ".rds"))
      )
    }
    default_name <- as.character(default_threshold)
    if (!default_name %in% names(edges)) {
      stop("Default STRING threshold ", default_threshold,
           " is missing from edges list.", call. = FALSE)
    }
    save_rds_checked(
      edges[[default_name]],
      file = file.path(output_dir, paste0(prefix, "_edges.rds"))
    )
  } else {
    save_rds_checked(edges, file = file.path(output_dir, paste0(prefix, "_edges.rds")))
    save_rds_checked(
      edges,
      file = file.path(output_dir, paste0(prefix, "_edges_string", default_threshold, ".rds"))
    )
  }
}

save_string_network_outputs <- function(outputs, edges, output_dir, prefix,
                                        default_threshold = DEFAULT_STRING_SCORE_THRESHOLD) {
  ensure_dir(output_dir)
  save_rds_checked(outputs$nodes, file = file.path(output_dir, paste0(prefix, "_nodes.rds")))
  save_string_edges_by_threshold(edges, output_dir, prefix, default_threshold)
  save_rds_checked(outputs$nodes_with_tooltips,
          file = file.path(output_dir, paste0(prefix, "_nodes_with_tooltips.rds")))
  save_rds_checked(outputs$cluster_for_plot,
          file = file.path(output_dir, paste0(prefix, "_cluster_for_plot.rds")))
  save_rds_checked(outputs$cluster_descriptions,
          file = file.path(output_dir, paste0(prefix, "_cluster_descriptions.rds")))
}

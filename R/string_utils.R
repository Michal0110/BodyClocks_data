source_if_exists <- function(path) {
  if (file.exists(path)) source(path)
}

source_if_exists(file.path("R", "io_utils.R"))

init_stringdb <- function(species_id, version = "12.0", score_threshold = 700,
                          input_directory = file.path("results", "cache", "stringdb")) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("STRINGdb is required for init_stringdb().", call. = FALSE)
  }
  ensure_dir(input_directory)
  tryCatch(
    STRINGdb::STRINGdb$new(
      version = version,
      species = species_id,
      score_threshold = score_threshold,
      input_directory = input_directory
    ),
    error = function(e) {
      stop("Unable to initialise STRINGdb for species ", species_id,
           " version ", version, ": ", conditionMessage(e), call. = FALSE)
    }
  )
}

STRING_SCORE_THRESHOLDS <- c(700, 900)
DEFAULT_STRING_SCORE_THRESHOLD <- 700

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
    return(data.frame(from = character(), to = character()))
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
                                      symbol_col = "symbol") {
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
    cluster_data = get_string_enrichment(gene_list, default_string_db)
  )
}

empty_string_edges_by_threshold <- function(thresholds = STRING_SCORE_THRESHOLDS) {
  thresholds <- unique(as.integer(thresholds))
  stats::setNames(
    lapply(thresholds, function(x) data.frame(from = character(), to = character())),
    as.character(thresholds)
  )
}

get_string_enrichment <- function(gene_list, string_db,
                                  categories = c("KEGG", "Process", "RCTM", "WikiPathways"),
                                  fdr_cutoff = 0.05) {
  gene_list <- unique(stats::na.omit(as.character(gene_list)))
  if (!length(gene_list)) {
    return(data.frame())
  }

  tryCatch({
    gene_df <- data.frame(gene = gene_list, stringsAsFactors = FALSE)
    cat("Mapping", length(gene_list), "genes to STRING database for enrichment...\n")
    mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

    cat("Performing enrichment analysis...\n")
    get_category <- function(category) {
      result <- string_db$get_enrichment(mapped_genes$STRING_id, category = category)
      result <- result[result$fdr <= fdr_cutoff, , drop = FALSE]
      if (nrow(result) > 0) {
        result$category <- category
      }
      result
    }

    pieces <- lapply(categories, get_category)
    pieces <- Filter(function(x) nrow(x) > 0, pieces)
    if (!length(pieces)) {
      return(data.frame())
    }
    do.call(rbind, pieces)
  }, error = function(e) {
    warning("STRINGdb enrichment failed: ", conditionMessage(e), call. = FALSE)
    data.frame()
  })
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

  cluster_data <- map_enrichment_genes(cluster_data, gene_mapping)
  cluster_data$preferredNames <- split_preferred_names(cluster_data$preferredNames)

  cluster_descriptions <- cluster_data %>%
    dplyr::select(cluster = term, cluster_description = description, category)

  reformatted_clusters <- cluster_data %>%
    dplyr::select(cluster = term, cluster_description = description,
                  id = preferredNames, category) %>%
    tidyr::unnest(id)

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

  node_tooltips <- nodes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      cluster_info = list(
        cluster_descriptions %>%
          dplyr::filter(cluster %in% unlist(clusters)) %>%
          dplyr::select(cluster, cluster_description, category)
      )
    ) %>%
    dplyr::mutate(
      enrichment_sections = list(
        if (length(cluster_info) > 0) {
          categories <- unique(cluster_info$category)
          sapply(categories, function(cat) {
            cat_name <- dplyr::case_when(
              cat == "Process" ~ "GO Biological Process",
              cat == "KEGG" ~ "KEGG Pathway",
              cat == "RCTM" ~ "Reactome",
              cat == "WikiPathways" ~ "WikiPathways",
              TRUE ~ cat
            )
            cat_terms <- cluster_info %>%
              dplyr::filter(category == cat) %>%
              dplyr::pull(cluster_description)
            if (length(cat_terms) > 0) {
              paste0("<b>", cat_name, ":</b> ", paste(cat_terms, collapse = ", "))
            } else {
              character(0)
            }
          })
        } else {
          character(0)
        }
      ),
      formatted_tooltip = paste(stats::na.omit(unlist(enrichment_sections)),
                                collapse = "<br>")
    ) %>%
    dplyr::select(id, formatted_tooltip)

  nodes_with_tooltips <- merged_nodes %>%
    dplyr::left_join(node_tooltips, by = "id") %>%
    dplyr::mutate(
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
    cluster_descriptions = cluster_descriptions
  )
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

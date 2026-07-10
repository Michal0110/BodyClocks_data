source_if_exists <- function(path) {
  if (file.exists(path)) source(path)
}

source_if_exists(file.path("R", "io_utils.R"))

retry_call <- function(expr, max_retries = 3, label = "operation") {
  errors <- character()
  for (attempt in seq_len(max_retries)) {
    result <- tryCatch(
      withCallingHandlers(
        expr,
        warning = function(w) {
          message(label, " warning on attempt ", attempt, ": ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      return(result)
    }
    errors <- c(errors, paste0("attempt ", attempt, ": ", conditionMessage(result)))
    if (attempt < max_retries) {
      Sys.sleep(2^(attempt - 1))
    }
  }
  stop(label, " failed after ", max_retries, " attempts:\n",
       paste(errors, collapse = "\n"), call. = FALSE)
}

connect_biomart <- function(species_dataset, ensembl_version = NULL,
                            max_retries = 3) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt is required for connect_biomart().", call. = FALSE)
  }

  attempts <- list()
  if (!is.null(ensembl_version)) {
    attempts[[paste0("useEnsembl archive v", ensembl_version)]] <- function() {
      biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = species_dataset,
        version = ensembl_version
      )
    }
  }
  attempts[["useEnsembl current"]] <- function() {
    biomaRt::useEnsembl(biomart = "ensembl", dataset = species_dataset)
  }
  attempts[["useMart current"]] <- function() {
    biomaRt::useMart("ensembl", dataset = species_dataset)
  }
  attempts[["useast mirror"]] <- function() {
    biomaRt::useMart(
      "ENSEMBL_MART_ENSEMBL",
      host = "https://useast.ensembl.org",
      dataset = species_dataset
    )
  }
  attempts[["uswest mirror"]] <- function() {
    biomaRt::useMart(
      "ENSEMBL_MART_ENSEMBL",
      host = "https://uswest.ensembl.org",
      dataset = species_dataset
    )
  }

  errors <- character()
  for (name in names(attempts)) {
    result <- tryCatch(
      retry_call(attempts[[name]](), max_retries, paste("BioMart", name)),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      attr(result, "bodyclocks_connection_method") <- name
      return(result)
    }
    errors <- c(errors, paste0(name, ": ", conditionMessage(result)))
  }

  stop(
    "Unable to connect to Ensembl BioMart for dataset ", species_dataset, ".\n",
    paste(errors, collapse = "\n"),
    call. = FALSE
  )
}

get_biomart_annotation <- function(ids, mart, attributes, filters,
                                   max_retries = 3, cache_key = NULL,
                                   use_cache = TRUE) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("biomaRt is required for get_biomart_annotation().", call. = FALSE)
  }
  ids <- unique(stats::na.omit(ids))
  if (!length(ids)) {
    return(data.frame())
  }

  cache_path <- NULL
  if (!is.null(cache_key)) {
    cache_path <- file.path("results", "cache", "annotation", paste0(cache_key, ".rds"))
    if (use_cache && file.exists(cache_path)) {
      cached <- readRDS(cache_path)
      if (is.list(cached) && identical(cached$attributes, attributes) &&
          identical(cached$filters, filters)) {
        return(cached$data)
      }
    }
  }

  result <- retry_call(
    biomaRt::getBM(
      attributes = attributes,
      filters = filters,
      values = ids,
      mart = mart
    ),
    max_retries = max_retries,
    label = "BioMart getBM"
  )

  if (!is.null(cache_path)) {
    save_rds_checked(
      list(
        data = result,
        attributes = attributes,
        filters = filters,
        created_at = as.character(Sys.time()),
        connection_method = attr(mart, "bodyclocks_connection_method")
      ),
      cache_path
    )
  }
  result
}

batch_ensembl_rest <- function(ids, size = 500, max_retries = 3) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("httr and jsonlite are required for batch_ensembl_rest().", call. = FALSE)
  }
  ids <- unique(stats::na.omit(ids))
  chunks <- split(ids, ceiling(seq_along(ids) / size))
  results <- lapply(chunks, function(chunk) {
    retry_call({
      response <- httr::POST(
        "https://rest.ensembl.org/lookup/id",
        httr::add_headers("Content-Type" = "application/json",
                          "Accept" = "application/json"),
        body = jsonlite::toJSON(list(ids = unname(chunk)), auto_unbox = TRUE),
        encode = "raw"
      )
      if (httr::status_code(response) >= 300) {
        stop("Ensembl REST status ", httr::status_code(response))
      }
      jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"),
                         flatten = TRUE)
    }, max_retries = max_retries, label = "Ensembl REST lookup")
  })
  do.call(c, unname(results))
}


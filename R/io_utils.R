ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

write_csv_checked <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(x, path, row.names = row.names)
  if (!file.exists(path)) {
    stop("Failed to write CSV: ", path, call. = FALSE)
  }
  invisible(path)
}

save_rds_checked <- function(x, path = NULL, file = NULL) {
  if (is.null(path)) {
    path <- file
  }
  if (is.null(path)) {
    stop("save_rds_checked() requires a path or file argument.", call. = FALSE)
  }
  ensure_dir(dirname(path))
  saveRDS(x, path)
  if (!file.exists(path)) {
    stop("Failed to write RDS: ", path, call. = FALSE)
  }
  invisible(path)
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  result <- tryCatch(system2("sha256sum", path, stdout = TRUE), error = function(e) NA_character_)
  if (!length(result) || is.na(result[[1]])) {
    return(NA_character_)
  }
  strsplit(result[[1]], "[[:space:]]+")[[1]][[1]]
}

log_session_info <- function(path) {
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
  invisible(path)
}

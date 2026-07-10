clean_description <- function(description) {
  description <- gsub("\\[.*?\\]", "", description)
  trimws(description)
}

fix_excel_genes <- function(symbol) {
  symbol <- as.character(symbol)
  symbol <- sub("^Mar-([0-9]+)$", "MARCHF\\1", symbol, ignore.case = TRUE)
  symbol <- sub("^Sep-([0-9]+)$", "SEPTIN\\1", symbol, ignore.case = TRUE)
  symbol <- sub("^Dec-([0-9]+)$", "DELEC\\1", symbol, ignore.case = TRUE)
  symbol
}

annotate_with_db <- function(keys, db, keytype, column = "SYMBOL",
                             multiVals = "first") {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("AnnotationDbi is required for annotate_with_db().", call. = FALSE)
  }
  AnnotationDbi::mapIds(
    x = db,
    keys = keys,
    keytype = keytype,
    column = column,
    multiVals = multiVals
  )
}


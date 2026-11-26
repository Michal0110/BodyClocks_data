# Evaluate named chunks from the versioned paper notebooks.
read_notebook_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^```\\{r ", lines)
  labels <- sub("^```\\{r ([^,} ]+).*", "\\1", lines[starts])
  stopifnot(!anyDuplicated(labels))
  setNames(lapply(starts, function(i) {
    after <- tail(lines, -i)
    parse(text = head(after, which(after == "```")[1] - 1L))
  }), labels)
}

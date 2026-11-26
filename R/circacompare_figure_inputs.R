# Map the fitted genes to the manuscript's table schema. Expression for Fig. 1C
# comes from the exact tables used to fit CircaCompare, in the original units.
read_circacompare_figure_inputs <- function(fit_dir, expression_dir) {
  inputs <- read.csv(file.path(fit_dir, "inputs.csv"))
  layout <- read.csv(file.path(fit_dir, "sample_hours.csv"))
  ids <- c(Osmo = "osmo", Dex = "dex", HS = "hs", Cartilage = "cartilage_wt")
  stopifnot(setequal(inputs$dataset, ids), !anyDuplicated(inputs$dataset))
  expression <- list()
  tables <- list()
  for (label in names(ids)) {
    d <- ids[[label]]
    info <- inputs[inputs$dataset == d, ]
    path <- file.path(expression_dir, info$file)
    stopifnot(digest::digest(file = path, algo = "sha256") == info$sha256)
    x <- read.csv(path, check.names = FALSE)
    gene_ids <- if (d == "cartilage_wt") x$CycID else x$Row.names
    fits <- read.csv(file.path(fit_dir, paste0("single_", d, ".csv")))
    stopifnot(!anyNA(gene_ids), !anyDuplicated(gene_ids), !anyDuplicated(fits$ensembl_gene_id),
              setequal(gene_ids, fits$ensembl_gene_id))
    fits <- fits[match(gene_ids, fits$ensembl_gene_id), ]
    stopifnot(all(is.finite(fits$peak_time_hours)), all(is.finite(fits$rhythmic_q)))
    tables[[label]] <- data.frame(Row.names = gene_ids, mgi_symbol = fits$symbol,
                                  phase = fits$peak_time_hours %% 24, BH.Q = fits$rhythmic_q)
    s <- layout[layout$dataset == d, ]
    stopifnot(nrow(s) > 0, !anyDuplicated(s$sample), all(s$sample %in% names(x)),
              all(is.finite(s$elapsed_or_ZT_hour)))
    expression[[label]] <- as.matrix(x[s$sample])
    rownames(expression[[label]]) <- gene_ids
  }

  four_way <- Reduce(intersect, lapply(tables, function(t) t$Row.names[t$BH.Q < 0.05]))
  symbols <- tables$Cartilage$mgi_symbol[match(four_way, tables$Cartilage$Row.names)]
  stopifnot(!anyNA(symbols), all(nzchar(symbols)), !anyDuplicated(symbols))
  plots <- lapply(names(ids), function(label) {
    hours <- layout$elapsed_or_ZT_hour[layout$dataset == ids[[label]]]
    y <- expression[[label]][four_way, , drop = FALSE]
    do.call(rbind, lapply(sort(unique(hours)), function(h) {
      values <- y[, hours == h, drop = FALSE]
      data.frame(Gene = symbols, time_point = h, mean = rowMeans(values),
                 sd = apply(values, 1, sd))
    }))
  })
  names(plots) <- names(ids)
  list(tables = tables, plots = plots)
}

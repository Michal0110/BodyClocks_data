# Run in the project R environment from the repository root after revision/prepare_inputs.py.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
tables <- "results/tables"
out <- args[1]
workers <- if (length(args) == 2L) as.integer(args[2]) else 1L
stopifnot(is.finite(workers), workers >= 1L)
local_lib <- "results/circacompare/library"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
library(circacompare)
stopifnot(packageDescription("circacompare")$Version == "0.2.0.9000")
if (dir.exists(out) || !dir.create(out, recursive = TRUE)) stop("Provide a new output directory.")
save_csv <- function(x, file) write.csv(x, file.path(out, paste0(file, ".csv")), row.names = FALSE)
datasets <- c("osmo", "dex", "hs", "cartilage_wt")
files <- c(osmo = "osmo_circ.csv", dex = "dex_circ.csv", hs = "hs_circ.csv",
           cartilage_wt = "wt_cartilage_circ.csv")
cutoff <- 0.05
hashes <- read.csv("revision/inputs/inputs.csv")
layout <- read.csv("revision/inputs/sample_hours.csv")
read_panel <- function(d) {
  path <- file.path(tables, files[[d]])
  hash <- digest::digest(file = path, algo = "sha256")
  stopifnot(identical(hash, hashes$sha256[match(d, hashes$dataset)]))
  x <- read.csv(path, check.names = FALSE)
  s <- layout[layout$dataset == d, setdiff(names(layout), "dataset")]
  if (d == "cartilage_wt") {
    clocks <- c("9AM" = 2, "1PM" = 6, "5PM" = 10, "9PM" = 14, "1AM" = 18, "5AM" = 22)
    samples <- grep("^WT_", names(x), value = TRUE)
    h <- setNames(unname(clocks[sub("^WT_(.+)_[0-9]+$", "\\1", samples)]), samples)
    id <- x$CycID
  } else {
    meta <- read.csv(file.path("data/meta", paste0(d, "_meta.csv")))
    if (d == "osmo") {
      meta <- subset(meta, timepoint != "0h")
      h <- setNames(as.integer(sub("h$", "", meta$timepoint)), paste0("X", meta$symbol))
    } else {
      h <- setNames(as.integer(sub("h$", "", meta$time)), meta$sample)
    }
    id <- x$Row.names
  }
  stopifnot(!anyDuplicated(s$sample), setequal(names(h), s$sample),
            identical(as.numeric(h[s$sample]), as.numeric(s$elapsed_or_ZT_hour)),
            all(s$sample %in% names(x)), !anyNA(id), !anyDuplicated(id))
  y <- as.matrix(x[, s$sample])
  rownames(y) <- id
  stopifnot(all(is.finite(y)))
  list(y = y, hours = s$elapsed_or_ZT_hour, layout = data.frame(dataset = d, s),
       genes = data.frame(ensembl_gene_id = id, symbol = x$preferred_symbol,
                          rain_q = x$BH.Q), hash = hash)
}
panel <- setNames(lapply(datasets, read_panel), datasets)
save_csv(do.call(rbind, lapply(panel, "[[", "layout")), "sample_hours")
save_csv(data.frame(dataset = datasets, file = unname(files[datasets]),
                    sha256 = vapply(panel, "[[", "", "hash")), "inputs")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

# Retain the original sampling axis on the first attempt. If random starts time out,
# use the same model on a rotated clock and convert reported peak times back to input hours.
# A common rotation leaves every between-group difference invariant.
run_fit <- function(fit, shift) {
  first <- tryCatch(fit(0, 250L), error = identity)
  if (!inherits(first, "error")) return(c(first, clock_shift_h = 0))
  if (!grepl("timeout|timing out", conditionMessage(first))) stop(first)
  c(fit(shift, 10000L), clock_shift_h = shift)
}

# Each index receives the same starting RNG state regardless of worker count.
fit_rows <- function(n, seed, columns, fitter) {
  one <- function(i) {
    set.seed(seed + i)
    tryCatch({
      vals <- fitter(i)
      stopifnot(identical(names(vals), columns), all(is.finite(vals)))
      list(values = vals, error = "")
    }, error = function(e) list(values = rep(NA_real_, length(columns)), error = conditionMessage(e)))
  }
  fits <- parallel::mclapply(seq_len(n), one, mc.cores = workers, mc.set.seed = FALSE)
  values <- do.call(rbind, lapply(fits, "[[", "values"))
  colnames(values) <- columns
  data.frame(values, error = vapply(fits, "[[", "", "error"), row.names = NULL)
}
single <- list()
for (d in datasets) {
  message("Single-dataset fits: ", d)
  p <- panel[[d]]
  fields <- c("rhythmic_p", "mesor", "amplitude", "peak_time_hours")
  fits <- fit_rows(nrow(p$y), 100000L * match(d, datasets), c(fields, "clock_shift_h"), function(i) {
    run_fit(function(shift, attempts) {
      z <- circa_single(data.frame(h = p$hours + shift, y = unname(p$y[i, ])), "h", "y",
                        period = 24, return_figure = FALSE, suppress_all = TRUE, timeout_n = attempts)$summary
      v <- setNames(z$value, z$parameter)[fields]
      v["peak_time_hours"] <- (v["peak_time_hours"] - shift) %% 24
      v
    }, shift = 6)
  })
  # Same 24h cosine waveform. Joint harmonic F test also accounts for fitting the phase under the alternative.
  angle <- 2 * pi * p$hours / 24
  design <- cbind(1, cos(angle), sin(angle))
  ols <- lm.fit(design, t(p$y))
  rss <- colSums(ols$residuals^2)
  harmonic_ss <- colSums((design %*% ols$coefficients - rep(colMeans(t(p$y)), each = nrow(design)))^2)
  df <- ols$df.residual
  harmonic_p <- pf((harmonic_ss / 2) / (rss / df), 2, df, lower.tail = FALSE)
  results <- data.frame(p$genes, fits, rhythmic_q = p.adjust(fits$rhythmic_p, "BH", n = nrow(fits)),
                        harmonic_F_p = harmonic_p, harmonic_F_q = p.adjust(harmonic_p, "BH"))
  single[[d]] <- results
  save_csv(results, paste0("single_", d))
}
source("alternative_processing/circacompare/pairs.R", local = TRUE)


# Fit one 24 h rhythm per gene and dataset, then compare rhythms between datasets.
# From the repository root, after revision/prepare_inputs.py, run
# Rscript alternative_processing/circacompare/fit.R NEW_OUTPUT_DIR [WORKERS]
# The single_*.csv tables contain the calls used for the manuscript comparisons.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
tables <- "results/tables"
out <- args[1]
workers <- if (length(args) == 2L) as.integer(args[2]) else 1L
stopifnot(is.finite(workers), workers >= 1L)


# Prefer the local installation of the package version used for the revision.
local_lib <- "results/circacompare/library"
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
library(circacompare)
stopifnot(packageDescription("circacompare")$Version == "0.2.0.9000")
if (dir.exists(out) || !dir.create(out, recursive = TRUE)) stop("Provide a new output directory.")
save_csv <- function(x, file) write.csv(x, file.path(out, paste0(file, ".csv")), row.names = FALSE)


# Keep the dataset order fixed. It sets the order of each contrast and the seeds.
datasets <- c("osmo", "dex", "hs", "cartilage_wt")
files <- c(osmo = "osmo_circ.csv", dex = "dex_circ.csv", hs = "hs_circ.csv",
           cartilage_wt = "wt_cartilage_circ.csv")
cutoff <- 0.05  # In pairs.R this gates native fits on raw p and defines the BH subset on q.
hashes <- read.csv("revision/inputs/inputs.csv")
layout <- read.csv("revision/inputs/sample_hours.csv")


# Build a gene-by-sample matrix in the same sample order as the hours passed to the fit.
# The checksum also distinguishes the corrected Dex table from the submitted one.
read_panel <- function(dataset) {
  path <- file.path(tables, files[[dataset]])
  hash <- digest::digest(file = path, algo = "sha256")
  stopifnot(identical(hash, hashes$sha256[match(dataset, hashes$dataset)]))
  rain_table <- read.csv(path, check.names = FALSE)
  samples <- layout[layout$dataset == dataset, setdiff(names(layout), "dataset")]

  # Check the saved sample hours against each study's own naming and time conventions.
  if (dataset == "cartilage_wt") {
    # WT cartilage uses Zeitgeber time. In this study 9 AM corresponds to ZT2.
    clocks <- c("9AM" = 2, "1PM" = 6, "5PM" = 10, "9PM" = 14, "1AM" = 18, "5AM" = 22)
    wt_columns <- grep("^WT_", names(rain_table), value = TRUE)
    hour_by_sample <- setNames(unname(clocks[sub("^WT_(.+)_[0-9]+$", "\\1", wt_columns)]), wt_columns)
    id <- rain_table$CycID
  } else {
    meta <- read.csv(file.path("data/meta", paste0(dataset, "_meta.csv")))
    if (dataset == "osmo") {
      # The Osmo RAIN input starts at 4 h. Its exported column names have an X prefix.
      meta <- subset(meta, timepoint != "0h")
      hour_by_sample <- setNames(as.integer(sub("h$", "", meta$timepoint)), paste0("X", meta$symbol))
    } else {
      # Match by sample name. In Dex, all three retained 24 h replicates belong at 24 h.
      hour_by_sample <- setNames(as.integer(sub("h$", "", meta$time)), meta$sample)
    }
    id <- rain_table$Row.names
  }
  stopifnot(!anyDuplicated(samples$sample), setequal(names(hour_by_sample), samples$sample),
            identical(as.numeric(hour_by_sample[samples$sample]), as.numeric(samples$elapsed_or_ZT_hour)),
            all(samples$sample %in% names(rain_table)), !anyNA(id), !anyDuplicated(id))
  y <- as.matrix(rain_table[, samples$sample])
  rownames(y) <- id
  stopifnot(all(is.finite(y)))
  list(y = y, hours = samples$elapsed_or_ZT_hour, layout = data.frame(dataset = dataset, samples),
       genes = data.frame(ensembl_gene_id = id, symbol = rain_table$preferred_symbol,
                          rain_q = rain_table$BH.Q), hash = hash)
}


# Record the inputs and clock used for each dataset alongside its fitted values.
panel <- setNames(lapply(datasets, read_panel), datasets)
save_csv(do.call(rbind, lapply(panel, "[[", "layout")), "sample_hours")
save_csv(data.frame(dataset = datasets, file = unname(files[datasets]),
                    sha256 = vapply(panel, "[[", "", "hash")), "inputs")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))


# Retain the original sampling axis on the first attempt. If random starts time out,
# use the same model on a rotated clock and convert reported peak times back to input hours.
# A common rotation leaves every between-group difference invariant.
run_fit <- function(fit, shift) {
  # timeout_n counts fitting attempts. The short first pass uses the recorded hours.
  first <- tryCatch(fit(0, 250L), error = identity)
  if (!inherits(first, "error")) return(c(first, clock_shift_h = 0))
  if (!grepl("timeout|timing out", conditionMessage(first))) stop(first)
  # On timeout, allow up to 10,000 attempts after rotating the clock by six hours.
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
  # A failed gene keeps its row, NA estimates and error message for the per-gene export.
  fits <- parallel::mclapply(seq_len(n), one, mc.cores = workers, mc.set.seed = FALSE)
  values <- do.call(rbind, lapply(fits, "[[", "values"))
  colnames(values) <- columns
  data.frame(values, error = vapply(fits, "[[", "", "error"), row.names = NULL)
}


# Fit a separate cosine to each gene. Mesor is the fitted mean, amplitude is half the
# peak-to-trough range, and peak_time_hours is measured from this dataset's time origin.
single <- list()
for (dataset in datasets) {
  message("Single-dataset fits: ", dataset)
  series <- panel[[dataset]]
  fields <- c("rhythmic_p", "mesor", "amplitude", "peak_time_hours")
  fit_single_gene <- function(i) {
    run_fit(function(shift, attempts) {
      fit_summary <- circa_single(
        data.frame(h = series$hours + shift, y = unname(series$y[i, ])), "h", "y",
        period = 24,           # Fixed circadian period, in hours, for all four datasets.
        return_figure = FALSE, # Only the fitted parameters are needed for gene-level tables.
        suppress_all = TRUE,
        timeout_n = attempts
      )$summary
      v <- setNames(fit_summary$value, fit_summary$parameter)[fields]
      # Undo any common time shift before reporting the peak on the original 24 h clock.
      v["peak_time_hours"] <- (v["peak_time_hours"] - shift) %% 24
      v
    }, shift = 6)
  }

  # Use a separate seed range per dataset so a gene's fit is reproducible with any worker count.
  fits <- fit_rows(
    n = nrow(series$y), seed = 100000L * match(dataset, datasets),
    columns = c(fields, "clock_shift_h"), fitter = fit_single_gene
  )


  # Sensitivity check on the same waveform. With no rhythm the peak time is undefined.
  # Compare an intercept-only fit to one with both sine and cosine coefficients (2 df).
  angle <- 2 * pi * series$hours / 24
  design <- cbind(1, cos(angle), sin(angle))
  ols <- lm.fit(design, t(series$y))
  rss <- colSums(ols$residuals^2)
  harmonic_ss <- colSums((design %*% ols$coefficients - rep(colMeans(t(series$y)), each = nrow(design)))^2)
  df <- ols$df.residual
  harmonic_p <- pf((harmonic_ss / 2) / (rss / df), 2, df, lower.tail = FALSE)

  # The primary calls use the package's amplitude-coefficient p value (rhythmic_p).
  # Adjust across genes within this dataset, separately for each rhythmicity test.
  results <- data.frame(
    series$genes, fits,
    rhythmic_q = p.adjust(fits$rhythmic_p, "BH", n = nrow(fits)),
    harmonic_F_p = harmonic_p,
    harmonic_F_q = p.adjust(harmonic_p, "BH")
  )
  single[[dataset]] <- results
  save_csv(results, paste0("single_", dataset))
}


# Use the single-series fits to select genes and fit the six two-group contrasts.
source("alternative_processing/circacompare/pairs.R", local = TRUE)

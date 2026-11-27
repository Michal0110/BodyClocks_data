# Summarise the single-series fits saved by fit.R. Run from the repository root with
# Rscript alternative_processing/circacompare/summarize_fits.R CALLS_DIR [RAIN_TABLE_DIR]
# Supplying the corrected RAIN tables also enables a comparison on identical genes.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
out <- args[1]
save_csv <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
datasets <- c("osmo", "dex", "hs", "cartilage_wt")
single <- setNames(lapply(datasets, function(d) read.csv(file.path(out, paste0("single_", d, ".csv")))), datasets)


# These columns hold within-dataset BH values. harmonic_F is the joint sine/cosine
# sensitivity test. The clock panel and Bmal1 ID follow the paper.
methods <- c(CircaCompare = "rhythmic_q", RAIN = "rain_q", harmonic_F = "harmonic_F_q")
core <- c("Bmal1", "Clock", "Per1", "Per2", "Per3", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp", "Tef", "Npas2")
anchor <- "ENSMUSG00000055116"


# Count fitted genes and discoveries at q < 0.05. Retain the raw-p count to describe
# how many genes could enter the package's two-group fits.
by_dataset <- lapply(datasets, function(d) {
  x <- single[[d]]
  cc <- which(x$rhythmic_q < 0.05)
  data.frame(dataset = d, tested = nrow(x), fitted = sum(is.finite(x$rhythmic_p)),
             nominal_p05 = sum(x$rhythmic_p < 0.05, na.rm = TRUE), circacompare = length(cc),
             rain = sum(x$rain_q < 0.05),
             harmonic_F = sum(x$harmonic_F_q < 0.05), cc_and_rain = sum(x$rain_q[cc] < 0.05))
})
save_csv(do.call(rbind, by_dataset), "dataset_summary")


# Same ID-based intersections as the manuscript. Include full method membership.
ids <- Reduce(intersect, lapply(single, "[[", "ensembl_gene_id"))
membership <- data.frame(ensembl_gene_id = ids, symbol = single$cartilage_wt$symbol[match(ids, single$cartilage_wt$ensembl_gene_id)])
overlaps <- list()
for (method in names(methods)) {
  field <- methods[[method]]
  sets <- lapply(single, function(x) x$ensembl_gene_id[which(x[[field]] < 0.05)])
  membership[[method]] <- ids %in% Reduce(intersect, sets)
  for (groups in c(combn(datasets, 2, simplify = FALSE), list(datasets[1:3], datasets))) {

    # Restrict each list to genes measured in every dataset in this comparison.
    bg <- Reduce(intersect, lapply(single[groups], "[[", "ensembl_gene_id"))
    s <- lapply(sets[groups], intersect, y = bg)
    shared <- Reduce(intersect, s)
    is_pair <- length(s) == 2L
    expected <- if (is_pair) prod(lengths(s)) / length(bg) else NA_real_
    p <- if (is_pair) phyper(length(shared) - 1, length(s[[1]]), length(bg) - length(s[[1]]),
                             length(s[[2]]), lower.tail = FALSE) else NA_real_
    overlaps[[length(overlaps) + 1L]] <- data.frame(method, groups = paste(groups, collapse = "+"),
      jointly_tested = length(bg), shared_rhythmic = length(shared), expected, fold = length(shared) / expected,
      jaccard = length(shared) / length(Reduce(union, s)), p)
  }
}


# Adjust the six pairwise upper-tail overlap tests separately for each method.
ov <- do.call(rbind, overlaps)
ov$pair_overlap_q <- ave(ov$p, ov$method, FUN = function(p) p.adjust(p, "BH", n = sum(is.finite(p))))
save_csv(ov, "overlaps")
membership$core_clock <- membership$symbol %in% core
save_csv(membership[rowSums(membership[names(methods)]) > 0, ], "four_way_genes")


# Paper's three-stress, Bmal1-relative phase concordance endpoint (<=4 h for every pair).
stress <- single[datasets[1:3]]
ref <- do.call(rbind, lapply(single, function(x) x[x$ensembl_gene_id == anchor, ]))
stopifnot(nrow(ref) == 4L)
save_csv(data.frame(dataset = rownames(ref), ref, row.names = NULL), "bmal1")
shared <- Reduce(intersect, lapply(stress, function(x) x$ensembl_gene_id[which(x$rhythmic_q < 0.05)]))

# Wrap the signed difference into [−12, 12) before applying the four-hour criterion.
near <- function(x, y) abs((x - y + 12) %% 24 - 12) <= 4
concordant <- function(a, b, c) near(a, b) & near(a, c) & near(b, c)
if (all(ref$rhythmic_q[1:3] < 0.05) && length(shared) > 1L) {
  phase <- vapply(names(stress), function(d) {
    x <- stress[[d]]
    (x$peak_time_hours[match(shared, x$ensembl_gene_id)] - ref[d, "peak_time_hours"]) %% 24
  }, numeric(length(shared)))
  observed <- concordant(phase[, 1], phase[, 2], phase[, 3])
  save_csv(data.frame(ensembl_gene_id = shared, symbol = stress$osmo$symbol[match(shared, stress$osmo$ensembl_gene_id)],
                      phase, concordant = observed), "clock_relative_three_way_genes")


  # Keep Osmo phases fixed and independently shuffle the Dex and HS gene labels.
  # This preserves the phase distribution of each culture while breaking gene correspondence.
  permute <- function(keep, seed) {
    set.seed(seed)
    z <- phase[keep, , drop = FALSE]
    n <- nrow(z)
    trials <- 10000L  # Match the paper's three-culture permutation test.
    null <- replicate(trials, sum(concordant(z[, 1], z[sample.int(n), 2], z[sample.int(n), 3])))
    count <- sum(observed[keep])
    data.frame(n, concordant = count, fraction = count / n, null_mean_fraction = mean(null) / n,
               trials, perm_p = (1 + sum(null >= count)) / (1 + trials))
  }

  # Bmal1 has zero relative phase in every culture by construction. Repeat without counting
  # that gene, still using its peak as the origin. Use seed 42 for each test.
  crp <- rbind(data.frame(anchor_included = TRUE, permute(rep(TRUE, length(shared)), 42L)),
               data.frame(anchor_included = FALSE, permute(shared != anchor, 42L)))
  save_csv(crp, "clock_relative_three_way_summary")
}


# On identical genes, distinguish changes in the estimated phases from changes in gene selection.
if (length(args) == 2L) {
  inputs <- read.csv(file.path(out, "inputs.csv"))
  common_sets <- lapply(stress, function(x) x$ensembl_gene_id[which(x$rhythmic_q < 0.05 & x$rain_q < 0.05)])
  ids <- Reduce(intersect, common_sets)
  stopifnot(anchor %in% ids)
  matched <- data.frame(ensembl_gene_id = ids, symbol = stress$osmo$symbol[match(ids, stress$osmo$ensembl_gene_id)])
  for (method in c("CircaCompare", "RAIN")) {
    phases <- vapply(names(stress), function(d) {
      x <- stress[[d]]
      if (method == "CircaCompare") {
        peaks <- setNames(x$peak_time_hours, x$ensembl_gene_id)
      } else {
        info <- inputs[inputs$dataset == d, ]
        path <- file.path(args[2], info$file)
        stopifnot(digest::digest(file = path, algo = "sha256") == info$sha256)
        tab <- read.csv(path)
        peaks <- setNames(tab$phase, tab$Row.names)
      }
      unname((peaks[ids] - peaks[anchor]) %% 24)
    }, numeric(length(ids)))
    matched[[method]] <- concordant(phases[, 1], phases[, 2], phases[, 3])
  }
  save_csv(matched, "clock_relative_matched_genes")
}

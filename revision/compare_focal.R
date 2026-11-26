# One comparison entry point for fixed-24h detectors.
# Rscript revision/compare_focal.R RAIN|eJTK|CircaCompare CALLS_DIR OUTPUT_DIR
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, args[1] %in% c("RAIN", "eJTK", "CircaCompare"))
method <- args[1]
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("revision/notebook_chunks.R")
chunks <- read_notebook_chunks("analysis/comparison_analysis_paper.rmd")
suppressPackageStartupMessages(library(dplyr))
eval(chunks$helper_functions)
save_table <- function(x, name) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
ds <- c(Osmo = "osmo", Dex = "dex", HS = "hs", Cartilage = "cartilage_wt")
manifest <- read.csv("revision/inputs/inputs.csv")
tables <- lapply(ds, function(d) {
  info <- manifest[manifest$dataset == d, ]
  path <- file.path("results/tables", info$file)
  stopifnot(digest::digest(file = path, algo = "sha256") == info$sha256)
  x <- read.csv(path)
  id <- if (d == "cartilage_wt") x$CycID else x$Row.names
  t <- data.frame(id, symbol = x$preferred_symbol, rain_q = x$BH.Q, q = x$BH.Q, phase = x$phase)
  if (method != "RAIN") {
    if (method == "eJTK") {
      f <- read.csv(file.path(args[2], d, "gene_comparison.csv"))
      stopifnot(all(f$Period == 24))
      q <- f$GammaBH
      peak <- f$Phase
      ref <- f$RAIN_BH
    } else {
      f <- read.csv(file.path(args[2], paste0("single_", d, ".csv")))
      q <- f$rhythmic_q
      peak <- f$peak_time_hours
      ref <- f$rain_q
    }
    ix <- match(id, f$ensembl_gene_id)
    stopifnot(setequal(id, f$ensembl_gene_id), !anyDuplicated(f$ensembl_gene_id),
              isTRUE(all.equal(ref[ix], t$rain_q)))
    t$q <- q[ix]
    t$phase <- peak[ix]
  }
  stopifnot(!anyDuplicated(id), all(is.finite(t$q)), all(is.finite(t$phase)))
  t
})
# The paper defines this list in crp_setup. Evaluate that single assignment.
core_assignment <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                            identical(x[[2]], as.name("core_clock_genes")), as.list(chunks$crp_setup))
stopifnot(length(core_assignment) == 1L)
eval(core_assignment[[1]])
clock <- unique(unlist(lapply(tables, function(t) t$id[t$symbol %in% core_clock_genes])))
anchor <- "ENSMUSG00000055116"
symbols <- do.call(rbind, tables)
symbol <- setNames(symbols$symbol, symbols$id)
call <- function(q, cutoff) if (cutoff == 0.05) q < cutoff else q <= cutoff
pairs <- combn(names(ds), 2, simplify = FALSE)
common <- lapply(tables, "[[", "id")
overlaps <- list()
genes <- list()
detect <- list()
within_dataset <- list()
for (threshold in c(0.05, 0.01)) {
  sets <- lapply(tables, function(x) x$id[call(x$q, threshold)])
  for (d in names(ds)) {
    x <- tables[[d]]
    a <- call(x$q, threshold)
    b <- call(x$rain_q, threshold)
    for (drop_clock in c(FALSE, TRUE)) {
      bg <- if (drop_clock) setdiff(x$id, clock) else x$id
      test <- hypergeometric_test_pairwise_universe(x$id[a], x$id[b], bg)
      shared_clock <- intersect(intersect(x$id[a & b], bg), clock)
      within_dataset[[length(within_dataset)+1L]] <- data.frame(dataset = d, threshold, drop_clock,
        as.data.frame(test), n_clock = length(shared_clock),
        clock_genes = paste(sort(symbol[shared_clock]), collapse = "|"))
    }
    detect[[length(detect)+1L]] <- data.frame(dataset = d, threshold, tested = nrow(x), rhythmic = sum(a),
      rain = sum(b), shared = sum(a & b), method_only = sum(a & !b), rain_only = sum(b & !a),
      fraction_method_in_rain = if (any(a)) sum(a & b)/sum(a) else NA_real_,
      fraction_rain_in_method = if (any(b)) sum(a & b)/sum(b) else NA_real_)
  }
  for (drop_clock in c(FALSE, TRUE)) {
    rows <- lapply(pairs, function(g) {
      bg <- Reduce(intersect, common[g])
      if (drop_clock) bg <- setdiff(bg, clock)
      test <- hypergeometric_test_pairwise_universe(sets[[g[1]]], sets[[g[2]]], bg)
      ids <- intersect(Reduce(intersect, sets[g]), bg)
      label <- paste(g, collapse = "_vs_")
      if (length(ids)) genes[[length(genes)+1L]] <<- data.frame(pair = label, threshold, drop_clock,
        ensembl_gene_id = ids, symbol = unname(symbol[ids]), clock_panel = ids %in% clock)
      data.frame(pair = label, threshold, drop_clock, as.data.frame(test), n_clock = sum(ids %in% clock),
        clock_genes = paste(sort(symbol[intersect(ids, clock)]), collapse = "|"),
        fraction_clock = if (length(ids)) mean(ids %in% clock) else NA_real_)
    })
    family <- do.call(rbind, rows)
    family$q <- p.adjust(family$p_value, "BH")
    overlaps[[length(overlaps)+1L]] <- family
  }
}
within <- bind_rows(within_dataset) %>% group_by(threshold, drop_clock) %>%
  mutate(q = p.adjust(p_value, "BH")) %>% ungroup()
save_table(within, "within_dataset_overlap_with_RAIN")
save_table(do.call(rbind, detect), "detection")
save_table(do.call(rbind, overlaps), "overlap_tests")
save_table(bind_rows(genes), "overlap_genes")
sets <- lapply(tables, function(x) x$id[x$q < .05])
membership <- data.frame(ensembl_gene_id = unique(unlist(sets, use.names = FALSE)))
membership$symbol <- unname(symbol[membership$ensembl_gene_id])
membership$clock_panel <- membership$ensembl_gene_id %in% clock
for (d in names(ds)) membership[[d]] <- membership$ensembl_gene_id %in% sets[[d]]
save_table(membership, "rhythmic_membership")
regions <- lapply(list(names(ds)[1:3], names(ds)), function(g) {
  mask <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(g))))
  mask <- mask[rowSums(mask) > 0,,drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(mask)), function(i) {
    inside <- Reduce(intersect, sets[g[mask[i,]]])
    ids <- setdiff(inside, unlist(sets[g[!mask[i,]]], use.names = FALSE))
    data.frame(panel = paste(g, collapse = "+"), region = paste(g[mask[i,]], collapse = "+"),
               n = length(ids), clock = sum(ids %in% clock))
  }))
})
save_table(bind_rows(regions), "venn_regions")
anchors <- do.call(rbind, lapply(tables, function(t) t[t$id == anchor, ]))
stopifnot(nrow(anchors) == 4L)
save_table(data.frame(dataset = rownames(anchors), anchors), "bmal1")
relative <- lapply(names(tables), function(d) {
  t <- compute_clock_relative_phase(tables[[d]], anchors[d,"phase"])
  setNames(t$clock_relative_phase, t$id)
})
names(relative) <- names(ds)
phase_genes <- list()
phase_rows <- list()
for (rule in c("paper_anchor", "exclude_anchor_gene", "use_tested_anchor", "clock_panel_only")) {
  set.seed(42)
  rows <- lapply(c(pairs, list(names(ds)[1:3])), function(g) {
    ids <- Reduce(intersect, sets[g])
    if (rule == "exclude_anchor_gene") ids <- setdiff(ids, anchor)
    if (rule == "clock_panel_only") ids <- intersect(ids, clock)
    label <- paste(g, collapse = "_vs_")
    status <- if (!length(ids)) "empty_overlap" else if (rule != "use_tested_anchor" &&
      any(anchors[g,"q"] >= .05)) "anchor_above_cutoff" else if (length(ids) < 2) "insufficient_genes" else "complete"
    row <- data.frame(comparison = label, rule, shared = length(ids), concordant = NA_integer_,
      fraction = NA_real_, null_mean = NA_real_, permutations = 0L, exceedances = NA_integer_, p = NA_real_, status)
    if (status != "complete") return(row)
    phases <- lapply(relative[g], "[", ids)
    is_pair <- length(g) == 2L
    near <- function(a,b) circ_dist_24(a,b) <= 4
    three <- function(a,b,c) near(a,b) & near(a,c) & near(b,c)
    aligned <- if (is_pair) do.call(near, unname(phases)) else do.call(three, unname(phases))
    row$concordant <- sum(aligned)
    row$fraction <- mean(aligned)
    phase_genes[[length(phase_genes)+1L]] <<- data.frame(comparison = label, rule, ensembl_gene_id = ids,
      symbol = unname(symbol[ids]), as.data.frame(phases), concordant = unname(aligned))
    if (!is_pair) set.seed(42)
    n <- if (is_pair) 1000L else 10000L
    null <- replicate(n, if (is_pair) sum(near(phases[[1]],sample(phases[[2]]))) else
      sum(three(phases[[1]],sample(phases[[2]]),sample(phases[[3]]))))
    row$null_mean <- mean(null)/length(ids)
    row$permutations <- n
    row$exceedances <- sum(null >= row$concordant)
    row$p <- if (is_pair) row$exceedances/n else (1+row$exceedances)/(n+1)
    row
  })
  family <- do.call(rbind, rows)
  family$q <- c(p.adjust(family$p[1:6], "BH", n=6), NA_real_)
  phase_rows[[rule]] <- family
}
save_table(bind_rows(phase_rows), "phase_tests")
save_table(bind_rows(phase_genes), "phase_genes")


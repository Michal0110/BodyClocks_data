# Called by analysis/run_circacompare.R after fitting all four datasets.
# Use explicit factor levels: a positive difference means group 2 minus group 1.
pair_summaries <- list()
for (groups in combn(datasets, 2, simplify = FALSE)) {
  first <- groups[1]
  second <- groups[2]
  label <- paste(groups, collapse = "_vs_")
  a <- single[[first]]
  b <- single[[second]]
  common <- intersect(a$ensembl_gene_id, b$ensembl_gene_id)
  a <- a[match(common, a$ensembl_gene_id), ]
  b <- b[match(common, b$ensembl_gene_id), ]
  nominal <- which(a$rhythmic_p < cutoff & b$rhythmic_p < cutoff)
  ids <- common[nominal]
  message("Pairwise fits: ", label, " (", length(ids), " genes passing the package's nominal gate)")
  s1 <- panel[[first]]
  s2 <- panel[[second]]
  hours <- c(s1$hours, s2$hours)
  group <- factor(rep(groups, c(length(s1$hours), length(s2$hours))), levels = groups)
  fields <- c(mesor_delta = "Mesor difference estimate", mesor_p = "P-value for mesor difference",
              amplitude_delta = "Amplitude difference estimate", amplitude_p = "P-value for amplitude difference",
              phase_delta_h = "Phase difference estimate", phase_p = "P-value for difference in phase")
  fits <- fit_rows(length(ids), 1000000L + length(pair_summaries) * 100000L, c(names(fields), "clock_shift_h"), function(i) {
    values <- c(s1$y[ids[i], ], s2$y[ids[i], ])
    run_fit(function(shift, attempts) {
      z <- circacompare(data.frame(h = hours + shift, group = group, y = unname(values)), "h", "group", "y",
                        period = 24, alpha_threshold = cutoff, suppress_all = TRUE, timeout_n = attempts)$summary
      setNames(z$value[match(fields, z$parameter)], names(fields))
    }, shift = 6)
  })
  both_q <- a$rhythmic_q[nominal] < cutoff & b$rhythmic_q[nominal] < cutoff
  x <- data.frame(ensembl_gene_id = ids, symbol = a$symbol[nominal],
                  rhythmic_q_1 = a$rhythmic_q[nominal], rhythmic_q_2 = b$rhythmic_q[nominal],
                  both_rhythmic_q05 = both_q, fits)
  for (endpoint in c("mesor", "amplitude", "phase")) {
    p <- x[[paste0(endpoint, "_p")]]
    x[[paste0(endpoint, "_q_nominal_set")]] <- p.adjust(p, "BH", n = length(p))
    key <- paste0(endpoint, "_q_both_rhythmic_set")
    x[[key]] <- NA_real_
    x[[key]][both_q] <- p.adjust(p[both_q], "BH", n = sum(both_q))
  }
  save_csv(x, paste0("pair_", label))
  pair_summaries[[label]] <- do.call(rbind, lapply(c("both_rhythmic_set", "nominal_set"), function(family) {
    use <- if (family == "both_rhythmic_set") both_q else rep(TRUE, nrow(x))
    data.frame(first, second, family, shared_tested = length(common),
               eligible = sum(use), fitted = sum(use & x$error == ""),
               mesor_q05 = sum(x[[paste0("mesor_q_", family)]] < cutoff, na.rm = TRUE),
               amplitude_q05 = sum(x[[paste0("amplitude_q_", family)]] < cutoff, na.rm = TRUE),
               phase_q05 = sum(x[[paste0("phase_q_", family)]] < cutoff, na.rm = TRUE))
  }))
}
save_csv(do.call(rbind, pair_summaries), "pair_summary")


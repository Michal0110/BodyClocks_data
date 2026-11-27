# Sourced by fit.R with the single-series fits, expression matrices and fit helpers in scope.
# Each output row describes one gene in one dataset pair. A gene can differ in several
# endpoints. These are comparisons of fitted curves on the studies' expression scales.
pair_summaries <- list()
for (groups in combn(datasets, 2, simplify = FALSE)) {
  first <- groups[1]
  second <- groups[2]
  label <- paste(groups, collapse = "_vs_")
  first_calls <- single[[first]]
  second_calls <- single[[second]]

  # Align the two single-series tables by tested Ensembl ID before selecting eligible genes.
  common <- intersect(first_calls$ensembl_gene_id, second_calls$ensembl_gene_id)
  first_calls <- first_calls[match(common, first_calls$ensembl_gene_id), ]
  second_calls <- second_calls[match(common, second_calls$ensembl_gene_id), ]

  # CircaCompare admits a gene when both single-series raw rhythmicity p values are < 0.05.
  # Keep these fits so readers can inspect every gene admitted by the package.
  eligible_rows <- which(first_calls$rhythmic_p < cutoff & second_calls$rhythmic_p < cutoff)
  ids <- common[eligible_rows]
  message("Pairwise fits: ", label, " (", length(ids), " genes passing the package's nominal gate)")


  # Stack the original observations. The factor level order makes the signed contrast
  # second minus first. Phase uses the sampling clocks, before any Bmal1 alignment.
  first_series <- panel[[first]]
  second_series <- panel[[second]]
  hours <- c(first_series$hours, second_series$hours)
  group <- factor(rep(groups, c(length(first_series$hours), length(second_series$hours))), levels = groups)
  fields <- c(
    mesor_delta = "Mesor difference estimate",
    mesor_p = "P-value for mesor difference",
    amplitude_delta = "Amplitude difference estimate",
    amplitude_p = "P-value for amplitude difference",
    phase_delta_h = "Phase difference estimate",
    phase_p = "P-value for difference in phase"
  )

  # Fit a shared 24 h period with group-specific mesor, amplitude and phase.
  # A six-hour retry shifts both groups equally, preserving their phase difference.
  fit_pair_gene <- function(i) {
    values <- c(first_series$y[ids[i], ], second_series$y[ids[i], ])
    run_fit(function(shift, attempts) {
      fit_summary <- circacompare(
        data.frame(h = hours + shift, group = group, y = unname(values)),
        "h", "group", "y",
        period = 24,
        alpha_threshold = cutoff,  # Native single-group rhythmicity gate.
        suppress_all = TRUE, timeout_n = attempts
      )$summary
      setNames(fit_summary$value[match(fields, fit_summary$parameter)], names(fields))
    }, shift = 6)
  }

  # Reserve a separate seed range for each pair, in the dataset order chosen by fit.R.
  fits <- fit_rows(
    n = length(ids), seed = 1000000L + length(pair_summaries) * 100000L,
    columns = c(names(fields), "clock_shift_h"), fitter = fit_pair_gene
  )


  # The reviewer summary asks about genes called rhythmic after BH in both datasets.
  # Identify that subset among the genes that passed the raw-p gate.
  rhythmic_in_both <- first_calls$rhythmic_q[eligible_rows] < cutoff &
    second_calls$rhythmic_q[eligible_rows] < cutoff
  pair_results <- data.frame(
    ensembl_gene_id = ids, symbol = first_calls$symbol[eligible_rows],
    rhythmic_q_1 = first_calls$rhythmic_q[eligible_rows],
    rhythmic_q_2 = second_calls$rhythmic_q[eligible_rows],
    both_rhythmic_q05 = rhythmic_in_both, fits
  )

  # Adjust each endpoint within this pair. Recompute BH on the smaller subset as well.
  # Subsetting q values from the larger family would answer a different question.
  for (endpoint in c("mesor", "amplitude", "phase")) {
    endpoint_p <- pair_results[[paste0(endpoint, "_p")]]
    pair_results[[paste0(endpoint, "_q_nominal_set")]] <- p.adjust(endpoint_p, "BH", n = length(endpoint_p))
    key <- paste0(endpoint, "_q_both_rhythmic_set")
    pair_results[[key]] <- NA_real_
    pair_results[[key]][rhythmic_in_both] <- p.adjust(endpoint_p[rhythmic_in_both], "BH", n = sum(rhythmic_in_both))
  }

  # Count eligible and successful fits in each family before tallying significant differences.
  save_csv(pair_results, paste0("pair_", label))
  pair_summaries[[label]] <- do.call(rbind, lapply(c("both_rhythmic_set", "nominal_set"), function(family) {
    use <- if (family == "both_rhythmic_set") rhythmic_in_both else rep(TRUE, nrow(pair_results))
    data.frame(first, second, family, shared_tested = length(common),
               eligible = sum(use), fitted = sum(use & pair_results$error == ""),
               mesor_q05 = sum(pair_results[[paste0("mesor_q_", family)]] < cutoff, na.rm = TRUE),
               amplitude_q05 = sum(pair_results[[paste0("amplitude_q_", family)]] < cutoff, na.rm = TRUE),
               phase_q05 = sum(pair_results[[paste0("phase_q_", family)]] < cutoff, na.rm = TRUE))
  }))
}
save_csv(do.call(rbind, pair_summaries), "pair_summary")

# Plot Bmal1-relative phase differences for the fixed-24 h eJTK analysis.
# Run after revision/compare_focal.R has written the gene-level phase comparisons.
library(ggplot2)

method_dir <- "alternative_processing/ejtk"
results_dir <- file.path(method_dir, "results")
concordance_hours <- 4  # The paper includes genes exactly four hours apart.


# The primary analysis requires rhythmic Bmal1 in both datasets. Osmo Bmal1
# misses the q cutoff and cartilage has no eJTK calls, leaving only Dex versus HS.
phase <- read.csv(file.path(results_dir, "phase_genes.csv"))
phase <- subset(phase, rule == "paper_anchor")
stopifnot(identical(unique(phase$comparison), "Dex_vs_HS"))

membership <- read.csv(file.path(results_dir, "rhythmic_membership.csv"))
clock_ids <- membership$ensembl_gene_id[membership$clock_panel]


# Dex and HS already contain peak hours relative to each culture's Bmal1 peak.
# Use the paper's circular distance: peaks at +11 and -11 hours are two hours apart.
source("revision/notebook_chunks.R")
eval(read_notebook_chunks("analysis/comparison_analysis_paper.rmd")$helper_functions)
phase$distance <- circ_dist_24(phase$Dex, phase$HS)
stopifnot(identical(phase$distance <= concordance_hours, phase$concordant))

phase$panel <- ifelse(phase$ensembl_gene_id %in% clock_ids, "Clock panel", "Other genes")
phase$symbol <- reorder(phase$symbol, -phase$distance)


# Each row represents one shared rhythmic gene. Bmal1 is included and has zero
# distance by construction. The dashed line marks the same threshold used in the table.
fig <- ggplot(phase, aes(distance, symbol)) +
  geom_segment(aes(x = 0, xend = distance, yend = symbol), colour = "grey80") +
  geom_vline(xintercept = concordance_hours, linetype = "dashed", colour = "grey35") +
  geom_point(aes(shape = panel, colour = panel), size = 3) +
  scale_colour_manual(values = c("#176c71", "#9c5b2a")) +
  scale_x_continuous(breaks = seq(0, 12, 2), limits = c(0, 12)) +
  labs(
    x = "Absolute circular difference in Bmal1-relative peak time (h)", y = NULL,
    title = "Dexamethasone and heat shock",
    subtitle = sprintf("%d of %d shared eJTK calls fall within four hours",
                       sum(phase$concordant), nrow(phase)),
    shape = NULL, colour = NULL
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

for (ext in c("pdf", "png")) {
  output <- file.path(method_dir, "figures", paste0("clock_relative_phase.", ext))
  ggsave(output, fig, width = 6.8, height = 5.3, dpi = 300)
}

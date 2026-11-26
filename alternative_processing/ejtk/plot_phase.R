# Plot all genes in the one pair that passed the paper's anchor and overlap gates.
library(ggplot2)
p <- "alternative_processing/ejtk"
phase <- read.csv(file.path(p, "results/phase_genes.csv"))
phase <- subset(phase, rule == "paper_anchor")
stopifnot(identical(unique(phase$comparison), "Dex_vs_HS"))
core <- subset(read.csv(file.path(p, "results/rhythmic_membership.csv")), clock_panel)$ensembl_gene_id
source("revision/notebook_chunks.R")
eval(read_notebook_chunks("analysis/comparison_analysis_paper.rmd")$helper_functions)
phase$distance <- circ_dist_24(phase$Dex, phase$HS)
stopifnot(identical(phase$distance <= 4, phase$concordant))
phase$panel <- ifelse(phase$ensembl_gene_id %in% core, "Clock panel", "Other genes")
phase$symbol <- reorder(phase$symbol, -phase$distance)
fig <- ggplot(phase, aes(distance, symbol)) +
  geom_segment(aes(x=0, xend=distance, yend=symbol), colour="grey80") +
  geom_vline(xintercept=4, linetype="dashed", colour="grey35") +
  geom_point(aes(shape=panel, colour=panel), size=3) +
  scale_colour_manual(values=c("#176c71", "#9c5b2a")) +
  scale_x_continuous(breaks=seq(0,12,2), limits=c(0,12)) +
  labs(x="Absolute circular difference in Bmal1-relative peak time (h)", y=NULL,
       title="Dexamethasone and heat shock", subtitle=sprintf("%d of %d shared eJTK calls fall within four hours", sum(phase$concordant), nrow(phase)),
       shape=NULL, colour=NULL) + theme_classic() + theme(legend.position="bottom")
for (ext in c("pdf","png"))
  ggsave(file.path(p, "figures", paste0("clock_relative_phase.",ext)), fig, width=6.8, height=5.3, dpi=300)


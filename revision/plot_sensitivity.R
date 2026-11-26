# Rscript revision/plot_sensitivity.R METHOD RESULT_DIRECTORY FIGURE_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
library(ggplot2)
library(patchwork)
get_table <- function(n) read.csv(file.path(args[2], paste0(n, ".csv")))


# Within each dataset, split primary calls (q < 0.05) into shared and method-specific sets.
d <- subset(get_table("detection"), threshold == 0.05)
bars <- tidyr::pivot_longer(d, c(shared, method_only, rain_only), names_to="membership")
bars$membership <- factor(bars$membership, levels=c("shared", "method_only", "rain_only"),
                          labels=c("Both methods", paste(args[1], "only"), "RAIN only"))
p <- ggplot(bars, aes(dataset, value, fill=membership)) + geom_col() +
  scale_fill_manual(values=c("#167d83","#e69732","#667796")) +
  labs(x=NULL, y="Rhythmic genes", fill=NULL) + theme_classic()


# Between datasets, compare observed overlap to the hypergeometric expectation.
# The second panel excludes clock genes from both lists and the universe.
tests <- subset(get_table("overlap_tests"), threshold == .05)
tests$pair <- gsub("_vs_", " / ", tests$pair)
tests$panel <- ifelse(tests$drop_clock, "Clock panel removed", "All tested genes")
t <- tidyr::pivot_longer(tests, c(overlap, expected), names_to="estimate")
q <- ggplot(t, aes(value, reorder(pair, value), shape=estimate, colour=estimate)) +
  geom_point(size=2.8) + facet_wrap(~panel) + theme_classic() +
  scale_colour_manual(values=c("#626262","#157e85")) +
  labs(x="Pairwise overlap (genes)", y=NULL, shape=NULL, colour=NULL)

# Write both panels together for the reviewer summary.
fig <- (p / q) + plot_annotation(title=paste(args[1], "and corrected RAIN"))
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
for (ext in c("pdf", "png"))
  ggsave(file.path(args[3],paste0("detection_and_overlap.",ext)), fig, width=10, height=7, dpi=300)

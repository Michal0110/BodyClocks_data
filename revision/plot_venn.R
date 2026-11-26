# Rscript revision/plot_venn.R COMPARISON_RESULTS FIGURE_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)

# Read exclusive three-culture counts at the primary rhythmicity cutoff.
regions <- subset(read.csv(file.path(args[1], "venn_regions.csv")), panel == "Osmo+Dex+HS")
counts <- setNames(regions$n, regions$region)


# Unscaled circles give small intersections enough room for readable counts.
centres <- rbind(c(-0.57,0.30), c(0.57,0.30), c(0,-0.61))
locations <- rbind(Osmo=c(-1.05,0.49), Dex=c(1.05,0.49), HS=c(0,-1.17),
                   "Osmo+Dex"=c(0,0.75), "Osmo+HS"=c(-0.57,-0.45),
                   "Dex+HS"=c(0.57,-0.45), "Osmo+Dex+HS"=c(0,0))
stopifnot(setequal(names(counts), rownames(locations)))
draw <- function() {
  par(mar=rep(0.3,4))
  plot.new()
  plot.window(c(-1.85,1.85),c(-1.95,1.85),asp=1)
  colours <- c("#32858A","#AD6F38","#6C6CA1")
  symbols(centres, circles=rep(1.08,3), inches=FALSE, add=TRUE,
          bg=adjustcolor(colours, alpha.f=0.19), fg=colours, lwd=2)
  text(locations, labels=counts[rownames(locations)],cex=1.4)
  text(c(-1.0,1.0,0), c(1.60,1.60,-1.87), c("Osmotic stress","Dexamethasone","Heat shock"),font=2)
}


# Export the same labelled diagram as vector artwork and a raster image.
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
stem <- file.path(args[2],"Fig2D_three_culture_overlap")
pdf(paste0(stem,".pdf"), width=5.3, height=5.2, useDingbats=FALSE)
draw()
dev.off()
png(paste0(stem,".png"), width=5.3, height=5.2, units="in", res=300)
draw()
dev.off()

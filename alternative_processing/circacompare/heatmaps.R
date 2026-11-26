# Evaluate the manuscript heatmap construction with CircaCompare q values and phases.
source("revision/notebook_chunks.R")
library(ComplexHeatmap)
p <- "alternative_processing/circacompare"
notebooks <- c(osmo="05_chondrocytes_osmotic_stress", dex="03_chondrocytes_dexamethasone",
               hs="04_chondrocytes_heat_shock")


# For each culture, keep the paper's expression values, labels and sample order.
# Replace its q and phase columns with the matched single-series fits.
for (d in names(notebooks)) {
  merged <- read.csv(file.path("results/tables",paste0(d,"_circ.csv")),check.names=FALSE)
  fit <- read.csv(file.path(p, "calls", paste0("single_", d, ".csv")))
  ix <- match(merged$Row.names, fit$ensembl_gene_id)
  stopifnot(!anyNA(ix), !anyDuplicated(ix), length(ix)==nrow(fit))
  merged$BH.Q <- fit$rhythmic_q[ix]
  merged$phase <- fit$peak_time_hours[ix]

  # Use the original row scaling, phase sorting and heatmap colours.
  chunk <- read_notebook_chunks(file.path("analysis",paste0(notebooks[d],".Rmd")))$prepare_heatmap
  # The named chunk constructs chm1 before opening its first device.
  first_device_setup <- which(vapply(
    chunk,
    function(x) is.call(x) && identical(x[[1]], as.name("ensure_dir")),
    logical(1)
  ))[1]
  eval(chunk[seq_len(first_device_setup - 1L)])

  # Open devices here so the replay writes to the CircaCompare figure directory.
  file <- file.path(p, "figures",paste0("Fig1B_heatmap_",d))
  pdf(paste0(file,".pdf"),width=10,height=8)
  draw(chm1)
  dev.off()
  png(paste0(file,".png"),width=10,height=8,units="in",res=600)
  draw(chm1)
  dev.off()
}

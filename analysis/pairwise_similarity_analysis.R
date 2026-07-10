# pairwise_similarity_analysis.R
#
# Standalone script: compute all-pairs pairwise similarity between circadian
# datasets within each species, then produce clustered heatmaps.
#
# Three metrics per pair:
#   1. Jaccard index           (normalises for list size differences)
#   2. Fold enrichment         (overlap / expected by chance)
#   3. -log10(BH.Q p-value)    (statistical significance of overlap)
#
# Significance annotations in heatmap cells:
#   ***  BH.Q < 0.001
#   **   BH.Q < 0.01
#   *    BH.Q < 0.05
#   (blank) otherwise
#
# Usage: Rscript analysis/pairwise_similarity_analysis.R
#   Can be run from the BodyClocks_data project root or from analysis/ after
#   generating display_dt RDS files with run_analysis.R.

# ============================================================
# CONFIG  - edit these as needed
# ============================================================
find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "config", "output_manifest.csv"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find project root containing config/output_manifest.csv", call. = FALSE)
    }
    current <- parent
  }
}

PROJECT_ROOT <- find_project_root()
OUTPUT_DIR        <- file.path("results", "pairwise_similarity")
PAPER_FIGURE_DIR <- file.path("results", "publication_figures")

# Per-species BH.Q cutoffs for calling a gene rhythmic.
# Baboon atlas uses a less stringent threshold (lower statistical power).
BHQ_CUTOFF <- c(
  mouse  = 0.05,
  baboon = 0.10
)

# Which species to analyse. Set to c("mouse", "baboon") to include both.
SPECIES_TO_RUN <- c("mouse")

# Exclude dataset ids that should not enter the pairwise analysis.
EXCLUDE_IDS <- c(
  "xiphoid",        # Very low number of rhythmic genes; skews fold enrichment heatmap
  "ko_cartilage",  # Hip articular Bmal1 KO
  "liver_AL_KO",   # Liver AL KO
  "liver_TRF_KO"   # Liver TRF KO
)

# pheatmap cell sizes (points); baboon gets smaller cells due to more tissues
CELL_SIZE_DEFAULT <- 20
CELL_SIZE_LARGE   <- 11     # used when N >= 30 tissues

# PNG resolution
PNG_DPI <- 600

# ============================================================
# PACKAGES
# ============================================================
required_pkgs <- c("pheatmap", "RColorBrewer")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Install them into the project environment and snapshot renv before rerunning.",
    call. = FALSE
  )
}

library(pheatmap)
library(RColorBrewer)

# viridis is optional; fall back to RColorBrewer if unavailable
has_viridis <- requireNamespace("viridis", quietly = TRUE)
if (has_viridis) library(viridis)

# ============================================================
# DATASET REGISTRY
# ============================================================
setwd(PROJECT_ROOT)

manifest_path <- file.path("config", "output_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Missing output manifest: ", manifest_path, call. = FALSE)
}

dataset_label <- function(ds_id) {
  labels <- c(
    cartilage = "Hip articular WT",
    ko_cartilage = "Hip articular Bmal1 KO",
    xiphoid = "Xiphoid cartilage",
    osmo = "Primary Chondrocytes - Osmotic stress",
    dex = "Primary Chondrocytes - Dexamethasone",
    hs = "Primary Chondrocytes - Heat shock",
    nih3t3 = "NIH3T3 Forskolin",
    mammary = "Mammary gland",
    tendon = "Tail tendon",
    pods = "Podocytes dexamethasone",
    gloms = "Glomeruli"
  )
  if (ds_id %in% names(labels)) return(labels[[ds_id]])

  label <- sub("^baboon_", "", ds_id)
  label <- sub("^liver_", "liver ", label)
  tools::toTitleCase(gsub("_", " ", label))
}

app_folder_to_result_folder <- function(app_folder) {
  sub("^data/", "results/shiny_data/", app_folder)
}

build_dataset_registry <- function(manifest) {
  display_rows <- manifest[
    manifest$output_type %in% c("display", "paired_display") &
      !is.na(manifest$dataset_id) &
      nzchar(manifest$dataset_id),
  ]

  registry <- lapply(seq_len(nrow(display_rows)), function(i) {
    row <- display_rows[i, ]
    ds_id <- row$dataset_id
    col_type <- if (row$output_type == "paired_display") {
      if (grepl("^ko_", ds_id)) "ko_primary" else "wt_primary"
    } else {
      "single"
    }

    list(
      dataset_label = dataset_label(ds_id),
      folder = app_folder_to_result_folder(row$app_folder),
      prefix = row$prefix,
      species = row$species,
      col_type = col_type,
      display_file = file.path(app_folder_to_result_folder(row$app_folder), row$expected_file)
    )
  })

  names(registry) <- display_rows$dataset_id
  registry
}

CIRC_DATASETS <- build_dataset_registry(read.csv(manifest_path, stringsAsFactors = FALSE))

# ============================================================
# LOCAL ANALYSIS HELPERS
# ============================================================
hypergeometric_test <- function(list1, list2, genome_size) {
  overlap <- length(intersect(list1, list2))
  size1   <- length(list1)
  size2   <- length(list2)

  p_value <- phyper(overlap - 1, size1, genome_size - size1, size2,
                    lower.tail = FALSE)

  expected <- (size1 * size2) / genome_size
  fold_enrichment <- if (expected > 0) overlap / expected else NA_real_
  jaccard <- if (length(union(list1, list2)) > 0)
    overlap / length(union(list1, list2))
  else 0

  list(
    overlap         = overlap,
    size1           = size1,
    size2           = size2,
    expected        = expected,
    fold_enrichment = fold_enrichment,
    jaccard         = jaccard,
    p_value         = p_value
  )
}

extract_gene_data <- function(table_data, ds_id) {
  if (is.null(table_data)) return(NULL)
  if (!"symbol" %in% names(table_data)) return(NULL)

  table_data <- table_data[!is.na(table_data$symbol) & nzchar(table_data$symbol), ]
  if (nrow(table_data) == 0) return(NULL)

  meta <- CIRC_DATASETS[[ds_id]]
  if (is.null(meta)) return(NULL)

  switch(meta$col_type,
    single = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal" %in% names(table_data)) table_data$pVal
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q,
      phase  = if ("phase" %in% names(table_data)) table_data$phase
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    wt_primary = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal_wt" %in% names(table_data)) table_data$pVal_wt
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q_wt,
      phase  = if ("phase_wt" %in% names(table_data)) table_data$phase_wt
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    ko_primary = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal_ko" %in% names(table_data)) table_data$pVal_ko
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q_ko,
      phase  = if ("phase_ko" %in% names(table_data)) table_data$phase_ko
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    NULL
  )
}

# ============================================================
# HELPER: load only the display_dt table for a dataset
# ============================================================
load_display_table <- function(ds_id) {
  meta <- CIRC_DATASETS[[ds_id]]
  if (is.null(meta)) stop("Unknown ds_id: ", ds_id)

  candidate_paths <- unique(c(
    meta$display_file,
    file.path(meta$folder, paste0(meta$prefix, "_display_dt.rds")),
    file.path(meta$folder, paste0(basename(meta$folder), "_display_dt.rds"))
  ))

  candidate_paths <- unique(candidate_paths[!is.na(candidate_paths) & nzchar(candidate_paths)])
  path <- candidate_paths[file.exists(candidate_paths)][1]

  if (is.na(path)) {
    warning("Display table not found for ", ds_id, ". Tried: ",
            paste(candidate_paths, collapse = ", "))
    return(NULL)
  }
  readRDS(path)
}

# ============================================================
# HELPER: extract rhythmic gene set + full gene list for a dataset
# Returns list(all = <character>, rhythmic = <character>)
# Returns NULL if the data file is missing.
# ============================================================
get_gene_sets <- function(ds_id, bhq_cutoff) {
  tbl <- load_display_table(ds_id)
  if (is.null(tbl)) return(NULL)

  gd <- extract_gene_data(tbl, ds_id)
  if (is.null(gd)) return(NULL)

  all_genes      <- gd$symbol[!is.na(gd$symbol) & nzchar(gd$symbol)]
  rhythmic_genes <- gd$symbol[!is.na(gd$bhq) & gd$bhq < bhq_cutoff]

  list(all = all_genes, rhythmic = rhythmic_genes)
}

# ============================================================
# HELPER: significance stars from a (BH-corrected) p-value
# ============================================================
pval_to_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", ""))))
}

# ============================================================
# HELPER: make & save a pheatmap (PDF + PNG)
# ============================================================
make_heatmap <- function(mat,
                         number_mat,   # character matrix for cell labels
                         title,
                         base_filename,
                         color_palette,
                         breaks        = NULL,
                         legend_breaks = NULL,
                         legend_labels = NULL) {

  n         <- nrow(mat)
  cell_size <- if (n >= 30) CELL_SIZE_LARGE else CELL_SIZE_DEFAULT

  # Build pheatmap object once (silent=TRUE avoids drawing to current device)
  ph <- pheatmap(
    mat,
    cluster_rows             = TRUE,
    cluster_cols             = TRUE,
    clustering_method        = "ward.D2",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    display_numbers          = number_mat,
    number_color             = "black",
    fontsize                 = if (n >= 30) 9 else 11,
    fontsize_number          = if (n >= 30) 6 else 8,
    color                    = color_palette,
    breaks                   = breaks,
    legend_breaks            = legend_breaks,
    legend_labels            = legend_labels,
    na_col                   = "black",
    cellwidth                = cell_size,
    cellheight               = cell_size,
    border_color             = "grey80",
    main                     = title,
    silent                   = TRUE   # don't draw yet
  )

  # Estimate output dimensions (inches)
  margin_in <- 5   # labels + title + legend
  plot_in   <- n * cell_size / 72 + margin_in

  pdf_file <- file.path(OUTPUT_DIR, paste0(base_filename, ".pdf"))
  png_file <- file.path(OUTPUT_DIR, paste0(base_filename, ".png"))

  # Render the cached gtable to each device explicitly
  pdf(pdf_file, width = plot_in, height = plot_in)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  png(png_file, width = plot_in, height = plot_in, units = "in", res = PNG_DPI)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  message("  Saved: ", pdf_file, " + PNG")
}

# ============================================================
# HELPER: duplicate selected paper figures
# ============================================================
copy_to_paper_figures <- function(base_filename, paper_filename = base_filename) {
  dir.create(PAPER_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

  extensions <- c(".pdf", ".png")
  source_files <- file.path(OUTPUT_DIR, paste0(base_filename, extensions))
  destination_files <- file.path(PAPER_FIGURE_DIR, paste0(paper_filename, extensions))

  for (i in seq_along(source_files)) {
    source_file <- source_files[[i]]
    if (!file.exists(source_file)) {
      warning("Paper figure source file not found: ", source_file)
      next
    }

    file.copy(source_file, destination_files[[i]], overwrite = TRUE)
    message("  Paper copy: ", destination_files[[i]])
  }
}

# ============================================================
# MAIN: run analysis per species
# ============================================================
dir.create(OUTPUT_DIR, showWarnings = FALSE)

for (sp in SPECIES_TO_RUN) {
  message("\n========== Species: ", sp, " ==========")

  # --- 1. Collect datasets for this species (exclude KO duplicates) ----------
  sp_ids <- names(CIRC_DATASETS)[
    vapply(CIRC_DATASETS, function(m) {
      isTRUE(m$species == sp)
    }, logical(1))
  ]
  sp_ids <- setdiff(sp_ids, EXCLUDE_IDS)

  message("Datasets to analyse: ", length(sp_ids))

  # --- 2. Load gene sets for all datasets ------------------------------------
  gene_sets <- list()
  labels    <- character(0)   # human-readable label per ds_id

  sp_bhq <- BHQ_CUTOFF[sp]

  for (ds_id in sp_ids) {
    gs <- get_gene_sets(ds_id, bhq_cutoff = sp_bhq)
    if (is.null(gs)) {
      message("  SKIP (no data): ", ds_id)
      next
    }
    if (length(gs$rhythmic) == 0) {
      message("  SKIP (0 rhythmic genes at BH.Q < ", sp_bhq, "): ", ds_id)
      next
    }
    gene_sets[[ds_id]] <- gs
    labels[ds_id]      <- CIRC_DATASETS[[ds_id]]$dataset_label
    message("  Loaded ", ds_id, ": ", length(gs$rhythmic), " rhythmic / ",
            length(gs$all), " total")
  }

  ids <- names(gene_sets)
  n   <- length(ids)

  if (n < 2) {
    message("  Too few datasets with rhythmic genes - skipping species ", sp)
    next
  }

  # --- 3. Compute all-pairs metrics ------------------------------------------
  # Flat data frame of upper-triangle pairs
  pairs_df <- data.frame(
    i               = integer(0),
    j               = integer(0),
    id1             = character(0),
    id2             = character(0),
    jaccard         = numeric(0),
    fold_enrichment = numeric(0),
    overlap         = integer(0),
    size1           = integer(0),
    size2           = integer(0),
    p_value         = numeric(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      id1 <- ids[i]
      id2 <- ids[j]
      gs1 <- gene_sets[[id1]]
      gs2 <- gene_sets[[id2]]

      # Background = genes measured in BOTH datasets (intersection).
      # Only these genes could possibly be rhythmic in either set, so using
      # the union would inflate the universe and deflate the p-value.
      bg_genes <- intersect(gs1$all, gs2$all)
      bg_size  <- length(bg_genes)

      res <- hypergeometric_test(
        intersect(gs1$rhythmic, bg_genes),
        intersect(gs2$rhythmic, bg_genes),
        bg_size
      )

      pairs_df <- rbind(pairs_df, data.frame(
        i               = i,
        j               = j,
        id1             = id1,
        id2             = id2,
        jaccard         = res$jaccard,
        fold_enrichment = res$fold_enrichment,
        overlap         = res$overlap,
        size1           = res$size1,
        size2           = res$size2,
        p_value         = res$p_value,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- 4. BH-correct p-values across all pairs in this species ---------------
  pairs_df$padj <- p.adjust(pairs_df$p_value, method = "BH")

  # --- 5. Build square symmetric matrices ------------------------------------
  jaccard_mat <- matrix(NA_real_, nrow = n, ncol = n,
                        dimnames = list(labels[ids], labels[ids]))
  fe_mat      <- matrix(NA_real_, nrow = n, ncol = n,
                        dimnames = list(labels[ids], labels[ids]))
  pval_mat    <- matrix(NA_real_, nrow = n, ncol = n,
                        dimnames = list(labels[ids], labels[ids]))
  diag(jaccard_mat) <- 1
  diag(fe_mat)      <- NA   # fold enrichment on diagonal is undefined / Inf

  for (k in seq_len(nrow(pairs_df))) {
    i    <- pairs_df$i[k]
    j    <- pairs_df$j[k]
    li   <- labels[ids[i]]
    lj   <- labels[ids[j]]
    jaccard_mat[li, lj] <- jaccard_mat[lj, li] <- pairs_df$jaccard[k]
    fe_mat[li, lj]      <- fe_mat[lj, li]      <- pairs_df$fold_enrichment[k]
    pval_mat[li, lj]    <- pval_mat[lj, li]    <- pairs_df$padj[k]
  }

  # Asterisk character matrix (off-diagonal only)
  star_mat <- matrix("", nrow = n, ncol = n,
                     dimnames = list(labels[ids], labels[ids]))
  for (k in seq_len(nrow(pairs_df))) {
    li <- labels[ids[pairs_df$i[k]]]
    lj <- labels[ids[pairs_df$j[k]]]
    s  <- pval_to_stars(pairs_df$padj[k])
    star_mat[li, lj] <- star_mat[lj, li] <- s
  }

  # --- 6. Save raw matrices as CSV -------------------------------------------
  sp_csv_prefix <- file.path(OUTPUT_DIR, sp)
  write.csv(as.data.frame(jaccard_mat), paste0(sp_csv_prefix, "_jaccard.csv"))
  write.csv(as.data.frame(fe_mat),      paste0(sp_csv_prefix, "_fold_enrichment.csv"))
  write.csv(pairs_df,                   paste0(sp_csv_prefix, "_pairs_table.csv"),
            row.names = FALSE)
  message("  CSVs saved.")

  # --- 7. Build -log10(padj) matrix ------------------------------------------
  logp_mat <- matrix(0, nrow = n, ncol = n,
                     dimnames = list(labels[ids], labels[ids]))
  for (k in seq_len(nrow(pairs_df))) {
    li <- labels[ids[pairs_df$i[k]]]
    lj <- labels[ids[pairs_df$j[k]]]
    v  <- -log10(pmax(pairs_df$padj[k], 1e-300))
    logp_mat[li, lj] <- logp_mat[lj, li] <- v
  }

  # --- 8. Colour palettes ----------------------------------------------------
  # Jaccard: viridis plasma or YlOrRd; breaks computed dynamically in 9a
  if (has_viridis) {
    jac_colors <- viridis::plasma(100)
  } else {
    jac_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
  }

  # Fold enrichment: diverging blue→white→red with white FIXED at exactly 1.
  # Breaks are asymmetric: the number of colour steps on each side of 1 is
  # proportional to how far the observed range extends below / above 1.
  fe_vals  <- na.omit(as.vector(fe_mat))
  fe_min   <- max(0, min(fe_vals, na.rm = TRUE))
  fe_max   <- min(max(fe_vals, na.rm = TRUE), 10)
  # Ensure 1 is inside the range (cap pathological cases)
  fe_min   <- min(fe_min, 0.99)
  fe_max   <- max(fe_max, 1.01)

  n_total  <- 200   # total colour steps
  n_low    <- round(n_total * (1 - fe_min) / (fe_max - fe_min))
  n_high   <- n_total - n_low

  fe_colors <- c(
    colorRampPalette(c("#2166AC", "#F7F7F7"))(n_low),   # blue  → white (FE < 1)
    colorRampPalette(c("#F7F7F7", "#B2182B"))(n_high)   # white → red   (FE > 1)
  )
  # n_total colours require n_total + 1 breakpoints; 1 is pinned at the join
  fe_breaks <- c(
    seq(fe_min, 1,      length.out = n_low  + 1),
    seq(1,      fe_max, length.out = n_high + 1)[-1]   # drop duplicate 1
  )

  # -log10(padj): sequential purple
  lp_max    <- max(logp_mat, na.rm = TRUE)
  lp_colors <- colorRampPalette(brewer.pal(9, "Purples"))(100)
  lp_breaks <- seq(0, max(lp_max, 1), length.out = 101)

  # --- 9. Plot heatmaps ------------------------------------------------------
  message("  Plotting heatmaps...")

  # 9a. Jaccard
  # Diagonal is shown as black (NA) so the colour scale spans only actual
  # between-dataset comparisons, not the trivial self-score of 1.
  jac_display <- jaccard_mat
  diag(jac_display) <- NA
  jac_off    <- jac_display[upper.tri(jac_display)]
  jac_max    <- ceiling(max(jac_off, na.rm = TRUE) * 20) / 20  # ceil to nearest 0.05
  jac_breaks <- seq(0, jac_max, length.out = 101)
  make_heatmap(
    mat           = jac_display,
    number_mat    = star_mat,
    title         = paste0(tools::toTitleCase(sp), " - Jaccard similarity of rhythmic gene sets (BH.Q < ", sp_bhq, ")"),
    base_filename = paste0(sp, "_jaccard"),
    color_palette = jac_colors,
    breaks        = jac_breaks,
    legend_breaks = round(seq(0, jac_max, length.out = 6), 2),
    legend_labels = as.character(round(seq(0, jac_max, length.out = 6), 2))
  )

  # 9b. Fold enrichment (diagonal shown as black NA)
  fe_display <- fe_mat
  diag(fe_display) <- NA
  make_heatmap(
    mat           = fe_display,
    number_mat    = star_mat,
    title         = paste0(tools::toTitleCase(sp), " - Fold enrichment of rhythmic gene overlap (BH.Q < ", sp_bhq, ")"),
    base_filename = paste0(sp, "_fold_enrichment"),
    color_palette = fe_colors,
    breaks        = fe_breaks,
    legend_breaks = round(sort(unique(c(seq(fe_min, 1, length.out = 3),
                                        seq(1, fe_max, length.out = 4)))), 2),
    legend_labels = as.character(round(sort(unique(c(seq(fe_min, 1, length.out = 3),
                                                      seq(1, fe_max, length.out = 4)))), 2))
  )

  if (identical(sp, "mouse")) {
    copy_to_paper_figures("mouse_fold_enrichment", "Fig1F_mouse_fold_enricment")
  }

  # 9c. -log10(BH.Q) significance
  logp_display <- logp_mat
  diag(logp_display) <- NA
  make_heatmap(
    mat           = logp_display,
    number_mat    = star_mat,
    title         = paste0(tools::toTitleCase(sp), " - Overlap significance: -log10(BH.Q) [cutoff ", sp_bhq, "]"),
    base_filename = paste0(sp, "_neg_log_padj"),
    color_palette = lp_colors,
    breaks        = lp_breaks
  )

  message("  Done: ", sp)
}

message("\nAll done. Outputs in: ", file.path(PROJECT_ROOT, OUTPUT_DIR))

# ============================================================
# FOCUSED ANALYSIS: Chondrocyte synchronisation comparison
#
# Tests the hypothesis that different in vitro synchronisation methods
# (dex, hs, osmo) produce more divergent rhythmic gene sets than the
# natural variation between different tissues in vivo.
#
# Three groups:
#   CC - Chondrocyte vs Chondrocyte  (3 pairs: dex/hs/osmo cross)
#   CV - Chondrocyte vs In Vivo      (each chondrocyte vs all other mouse tissues)
#   VV - In Vivo vs In Vivo          (all non-chondrocyte mouse tissue pairs)
# ============================================================

message("\n========== Chondrocyte synchronisation comparison ==========")

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 is required for the chondrocyte comparison plot.")
library(ggplot2)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

CHONDRO_IDS <- c("osmo", "dex", "hs")

# Read from CSV so this section is runnable independently
cp <- read.csv(file.path(OUTPUT_DIR, "mouse_pairs_table.csv"),
               stringsAsFactors = FALSE)

# Short label: strip "Primary Chondrocytes - " prefix
make_pair_label <- function(id) {
  lbl <- CIRC_DATASETS[[id]]$dataset_label
  sub("Primary Chondrocytes - ", "", lbl)
}

# ---- Classify pairs -------------------------------------------------------
classify_pair <- function(a, b) {
  ac <- a %in% CHONDRO_IDS
  bc <- b %in% CHONDRO_IDS
  if (ac & bc)  "Chondrocyte\nvs Chondrocyte"
  else if (ac | bc) "Chondrocyte\nvs In Vivo"
  else          "In Vivo\nvs In Vivo"
}
cp$group <- factor(
  mapply(classify_pair, cp$id1, cp$id2),
  levels = c("Chondrocyte\nvs Chondrocyte",
             "Chondrocyte\nvs In Vivo",
             "In Vivo\nvs In Vivo")
)

# Label for the 3 CC points
cp$pair_label <- NA_character_
is_cc <- cp$group == "Chondrocyte\nvs Chondrocyte"
cp$pair_label[is_cc] <- paste0(
  vapply(cp$id1[is_cc], make_pair_label, character(1)), "\nvs\n",
  vapply(cp$id2[is_cc], make_pair_label, character(1))
)

# ---- Console summary -------------------------------------------------------
cc_j <- cp$jaccard[is_cc]
cv_j <- cp$jaccard[cp$group == "Chondrocyte\nvs In Vivo"]
vv_j <- cp$jaccard[cp$group == "In Vivo\nvs In Vivo"]
cc_fe <- cp$fold_enrichment[is_cc]
cv_fe <- cp$fold_enrichment[cp$group == "Chondrocyte\nvs In Vivo"]
vv_fe <- cp$fold_enrichment[cp$group == "In Vivo\nvs In Vivo"]

cp_export <- cp
cp_export$group <- gsub("\n", " ", as.character(cp_export$group))
cp_export$pair_label <- gsub("\n", " ", cp_export$pair_label)
write.csv(cp_export,
          file.path(OUTPUT_DIR, "chondrocyte_comparison_pairs.csv"),
          row.names = FALSE)

chondro_summary <- do.call(rbind, lapply(levels(cp$group), function(g) {
  rows <- cp[cp$group == g, ]
  data.frame(
    group                = gsub("\n", " ", g),
    n                    = nrow(rows),
    mean_jaccard         = mean(rows$jaccard),
    median_jaccard       = median(rows$jaccard),
    min_jaccard          = min(rows$jaccard),
    max_jaccard          = max(rows$jaccard),
    mean_fold_enrichment = mean(rows$fold_enrichment),
    median_fold_enrichment = median(rows$fold_enrichment),
    min_fold_enrichment  = min(rows$fold_enrichment),
    max_fold_enrichment  = max(rows$fold_enrichment),
    stringsAsFactors     = FALSE
  )
}))
# Pairwise rows share datasets, so a row-wise test would treat non-independent
# observations as independent. Instead, permute which three dataset labels are
# called "chondrocyte" and compare the observed group mean differences against
# that label-permutation null.
classify_pair_for_selection <- function(a, b, selected_ids) {
  ac <- a %in% selected_ids
  bc <- b %in% selected_ids
  if (ac & bc)  "Chondrocyte\nvs Chondrocyte"
  else if (ac | bc) "Chondrocyte\nvs In Vivo"
  else          "In Vivo\nvs In Vivo"
}

group_metric_means <- function(selected_ids, metric) {
  groups <- factor(
    mapply(classify_pair_for_selection, cp$id1, cp$id2,
           MoreArgs = list(selected_ids = selected_ids)),
    levels = levels(cp$group)
  )
  tapply(cp[[metric]], groups, mean)
}

all_pair_ids <- sort(unique(c(cp$id1, cp$id2)))
label_sets <- combn(all_pair_ids, length(CHONDRO_IDS), simplify = FALSE)
group_pairs <- combn(levels(cp$group), 2, simplify = FALSE)
metrics_to_test <- c("jaccard", "fold_enrichment")

permutation_tests <- do.call(rbind, lapply(metrics_to_test, function(metric) {
  observed_means <- group_metric_means(CHONDRO_IDS, metric)
  do.call(rbind, lapply(group_pairs, function(gp) {
    observed_difference <- unname(observed_means[gp[1]] - observed_means[gp[2]])
    null_differences <- vapply(label_sets, function(selected_ids) {
      null_means <- group_metric_means(selected_ids, metric)
      unname(null_means[gp[1]] - null_means[gp[2]])
    }, numeric(1))

    data.frame(
      metric              = metric,
      comparison          = paste(gsub("\n", " ", gp), collapse = " minus "),
      group_1             = gsub("\n", " ", gp[1]),
      group_2             = gsub("\n", " ", gp[2]),
      group_1_mean        = unname(observed_means[gp[1]]),
      group_2_mean        = unname(observed_means[gp[2]]),
      observed_difference = observed_difference,
      p_value             = mean(abs(null_differences) >= abs(observed_difference)),
      method              = paste0(
        "Exact dataset-label permutation over all ",
        length(label_sets), " choices of ", length(CHONDRO_IDS), " labels"
      ),
      stringsAsFactors    = FALSE
    )
  }))
}))
permutation_tests$padj_bh <- p.adjust(permutation_tests$p_value, method = "BH")

make_column_slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("(^_|_$)", "", x)
}

add_pairwise_test_columns <- function(summary_df, tests_df) {
  summary_groups <- summary_df$group
  for (metric in unique(tests_df$metric)) {
    metric_tests <- tests_df[tests_df$metric == metric, ]
    prefix <- if (metric == "fold_enrichment") "fold_enrichment" else metric

    for (reference_group in summary_groups) {
      suffix <- make_column_slug(reference_group)
      difference_col <- paste0(prefix, "_difference_vs_", suffix)
      p_col <- paste0(prefix, "_p_value_vs_", suffix)
      padj_col <- paste0(prefix, "_padj_bh_vs_", suffix)

      summary_df[[difference_col]] <- NA_real_
      summary_df[[p_col]] <- NA_real_
      summary_df[[padj_col]] <- NA_real_

      reference_mean <- summary_df[summary_df$group == reference_group,
                                   paste0("mean_", prefix)]
      summary_df[[difference_col]] <- summary_df[[paste0("mean_", prefix)]] - reference_mean
      summary_df[[difference_col]][summary_df$group == reference_group] <- 0

      for (i in seq_len(nrow(metric_tests))) {
        test <- metric_tests[i, ]
        if (!(reference_group %in% c(test$group_1, test$group_2))) next

        other_group <- if (test$group_1 == reference_group) test$group_2 else test$group_1
        other_idx <- summary_df$group == other_group
        if (!any(other_idx)) next

        summary_df[[p_col]][other_idx] <- test$p_value
        summary_df[[padj_col]][other_idx] <- test$padj_bh
      }
    }
  }
  summary_df
}

chondro_summary <- add_pairwise_test_columns(chondro_summary, permutation_tests)

write.csv(chondro_summary,
          file.path(OUTPUT_DIR, "chondrocyte_comparison_summary.csv"),
          row.names = FALSE)
write.csv(permutation_tests,
          file.path(OUTPUT_DIR, "chondrocyte_comparison_permutation_tests.csv"),
          row.names = FALSE)

message("  Saved: chondrocyte_comparison_pairs.csv + chondrocyte_comparison_summary.csv + ",
        "chondrocyte_comparison_permutation_tests.csv")

cat("\n--- Jaccard index by group ---\n")
for (g in levels(cp$group)) {
  v <- cp$jaccard[cp$group == g]
  cat(sprintf("  %-36s n=%3d  mean=%.4f  median=%.4f  [%.4f, %.4f]\n",
              gsub("\n", " ", g), length(v), mean(v), median(v), min(v), max(v)))
}

cat("\n--- Individual chondrocyte cross-pairs ---\n")
cc_rows <- cp[is_cc, ]
for (k in seq_len(nrow(cc_rows))) {
  pct <- mean(vv_j < cc_rows$jaccard[k]) * 100
  cat(sprintf("  %-28s  Jaccard=%.4f  FE=%.2fx  overlap=%d (%d+%d)  BH.Q=%.4f  [%4.1f%% of VV]\n",
              gsub("\n", " vs ", cc_rows$pair_label[k]),
              cc_rows$jaccard[k], cc_rows$fold_enrichment[k],
              cc_rows$overlap[k], cc_rows$size1[k], cc_rows$size2[k],
              cc_rows$padj[k], pct))
}

cat(sprintf("\n  Mean Jaccard:  CC=%.4f  CV=%.4f  VV=%.4f\n",
            mean(cc_j), mean(cv_j), mean(vv_j)))
cat(sprintf("  Pct of VV pairs BELOW mean CC Jaccard: %.1f%%\n",
            mean(vv_j < mean(cc_j)) * 100))
cat(sprintf("  Pct of VV pairs BELOW min  CC Jaccard: %.1f%%\n",
            mean(vv_j < min(cc_j)) * 100))

# ---- Plot ------------------------------------------------------------------
group_colours <- c(
  "In Vivo\nvs In Vivo"          = "#4E79A7",
  "Chondrocyte\nvs In Vivo"      = "#59A14F",
  "Chondrocyte\nvs Chondrocyte"  = "#E15759"
)

p_chondro <- ggplot(cp, aes(x = group, y = jaccard, fill = group, colour = group)) +
  geom_violin(
    data      = ~ subset(., group != "Chondrocyte\nvs Chondrocyte"),
    alpha     = 0.25, linewidth = 0.3, trim = TRUE
  ) +
  geom_boxplot(
    data      = ~ subset(., group != "Chondrocyte\nvs Chondrocyte"),
    width     = 0.18, alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data  = ~ subset(., group == "Chondrocyte\nvs Chondrocyte"),
    size  = 4, shape = 18
  ) +
  {if (has_ggrepel)
    ggrepel::geom_label_repel(
      data          = ~ subset(., !is.na(pair_label)),
      aes(label     = pair_label),
      size          = 2.8, fill = "white", alpha = 0.9, colour = "#E15759",
      label.size    = 0.2,
      nudge_x       = 0.5, direction = "y", segment.size = 0.3, min.segment.length = 0
    )
  else
    geom_text(
      data      = ~ subset(., !is.na(pair_label)),
      aes(label = pair_label),
      hjust = -0.12, size = 2.6, colour = "#E15759"
    )
  } +
  geom_hline(yintercept = mean(cc_j), linetype = "dashed",
             colour = "#E15759", linewidth = 0.5, alpha = 0.7) +
  annotate("text", x = 0.52, y = mean(cc_j),
           label = sprintf("mean CC = %.3f", mean(cc_j)),
           size = 2.6, colour = "#E15759", vjust = -0.5, hjust = 0) +
  scale_fill_manual(values   = group_colours, guide = "none") +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_y_continuous(
    name   = "Jaccard index (rhythmic gene overlap)",
    limits = c(0, NA), expand = expansion(mult = c(0, 0.06))
  ) +
  scale_x_discrete(name = NULL) +
  labs(
    title    = "Chondrocyte synchronisation methods vs in vivo tissue diversity",
    subtitle = sprintf(
      "CC: dex/hs/osmo cross-pairs (n=3)  |  CV: chondrocyte vs in-vivo (n=%d)  |  VV: all in-vivo pairs (n=%d)",
      sum(cp$group == "Chondrocyte\nvs In Vivo"),
      sum(cp$group == "In Vivo\nvs In Vivo")
    ),
    caption  = sprintf(
      "Mouse, BH.Q < %.2f.   Mean Jaccard: CC=%.4f vs VV=%.4f.   %.1f%% of in-vivo tissue pairs fall below the CC mean.",
      BHQ_CUTOFF["mouse"], mean(cc_j), mean(vv_j),
      mean(vv_j < mean(cc_j)) * 100
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x        = element_text(size = 12, colour = "black"),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, colour = "grey40"),
    plot.caption       = element_text(size = 8, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey92")
  )

p_chondro_fe <- ggplot(cp, aes(x = group, y = fold_enrichment, fill = group, colour = group)) +
  geom_violin(
    data      = ~ subset(., group != "Chondrocyte\nvs Chondrocyte"),
    alpha     = 0.25, linewidth = 0.3, trim = TRUE
  ) +
  geom_boxplot(
    data      = ~ subset(., group != "Chondrocyte\nvs Chondrocyte"),
    width     = 0.18, alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data  = ~ subset(., group == "Chondrocyte\nvs Chondrocyte"),
    size  = 4, shape = 18
  ) +
  {if (has_ggrepel)
    ggrepel::geom_label_repel(
      data          = ~ subset(., !is.na(pair_label)),
      aes(label     = pair_label),
      size          = 2.8, fill = "white", alpha = 0.9, colour = "#E15759",
      label.size    = 0.2,
      nudge_x       = 0.5, direction = "y", segment.size = 0.3, min.segment.length = 0
    )
  else
    geom_text(
      data      = ~ subset(., !is.na(pair_label)),
      aes(label = pair_label),
      hjust = -0.12, size = 2.6, colour = "#E15759"
    )
  } +
  geom_hline(yintercept = 1, linetype = "dotted",
             colour = "grey45", linewidth = 0.45, alpha = 0.8) +
  geom_hline(yintercept = mean(cc_fe), linetype = "dashed",
             colour = "#E15759", linewidth = 0.5, alpha = 0.7) +
  annotate("text", x = 0.52, y = mean(cc_fe),
           label = sprintf("mean CC = %.2fx", mean(cc_fe)),
           size = 2.6, colour = "#E15759", vjust = -0.5, hjust = 0) +
  scale_fill_manual(values   = group_colours, guide = "none") +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_y_continuous(
    name   = "Fold enrichment of rhythmic gene overlap",
    limits = c(0, NA), expand = expansion(mult = c(0, 0.06))
  ) +
  scale_x_discrete(name = NULL) +
  labs(
    title    = "Chondrocyte synchronisation methods vs in vivo tissue diversity",
    subtitle = sprintf(
      "CC: dex/hs/osmo cross-pairs (n=3)  |  CV: chondrocyte vs in-vivo (n=%d)  |  VV: all in-vivo pairs (n=%d)",
      sum(cp$group == "Chondrocyte\nvs In Vivo"),
      sum(cp$group == "In Vivo\nvs In Vivo")
    ),
    caption  = sprintf(
      "Mouse, BH.Q < %.2f.   Mean fold enrichment: CC=%.2fx vs VV=%.2fx.   FE=1 marks overlap expected by chance.",
      BHQ_CUTOFF["mouse"], mean(cc_fe), mean(vv_fe)
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x        = element_text(size = 12, colour = "black"),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, colour = "grey40"),
    plot.caption       = element_text(size = 8, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey92")
  )

ggsave(file.path(OUTPUT_DIR, "chondrocyte_comparison.pdf"),
       p_chondro, width = 6.5, height = 4.5)
ggsave(file.path(OUTPUT_DIR, "chondrocyte_comparison.png"),
       p_chondro, width = 6.5, height = 4.5, dpi = PNG_DPI)
ggsave(file.path(OUTPUT_DIR, "chondrocyte_comparison_fold_enrichment.pdf"),
       p_chondro_fe, width = 6.5, height = 4.5)
ggsave(file.path(OUTPUT_DIR, "chondrocyte_comparison_fold_enrichment.png"),
       p_chondro_fe, width = 6.5, height = 4.5, dpi = PNG_DPI)
copy_to_paper_figures("chondrocyte_comparison_fold_enrichment", "Fig1E_chondrocyte_comparison_fold_enrichment")
message("  Saved: chondrocyte_comparison.pdf/PNG + chondrocyte_comparison_fold_enrichment.pdf/PNG")

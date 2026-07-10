# Script Inventory

This inventory should be completed during the first pass through each analysis
script. It records inputs, outputs, external services, and custom preprocessing
that must be preserved.

| New script | Source script(s) | Inputs | Outputs | External services | Notes |
|---|---|---|---|---|---|
| `analysis/01_cartilage.Rmd` | `cartilage.Rmd`, `circ_cartilage_timecourse_bodyclocks.rmd` | `data/raw/wt_cartilage.csv`, `data/raw/ko_cartilage.csv` | cartilage RAIN/circ tables and RDS files | Ensembl, STRINGdb | Paired WT/KO dataset |
| `analysis/02_tendon.Rmd` | `tendon.Rmd`, `circ_tendon_timecourse_bodyclocks.rmd` | `data/raw/tendon.csv` | tendon RAIN/circ tables and RDS files | Ensembl, STRINGdb | Preserve 12 h rhythm handling |
| `analysis/03_chondrocytes_dexamethasone.Rmd` | `chondrocytes_dexamethasone.Rmd` | `data/raw/dex_raw.csv`, `data/meta/dex_meta.csv` | dex RDS files | Ensembl, STRINGdb | DESeq2 preprocessing |
| `analysis/04_chondrocytes_heat_shock.Rmd` | `chondrocytes_heat_shock.rmd` | `data/raw/hs_raw.csv`, `data/meta/hs_meta.csv` | heat shock RDS files | Ensembl, STRINGdb | DESeq2 preprocessing |
| `analysis/05_chondrocytes_osmotic_stress.Rmd` | `chondrocytes_osmotic_stress.rmd` | `data/raw/osmo_raw.csv`, `data/meta/osmo_meta.csv` | osmotic stress RDS files | Ensembl, STRINGdb | Custom contrast and RAIN setup |
| `analysis/06_glomeruli.Rmd` | `glomeruli_for_bodyclocks.rmd` | `data/raw/gloms_full.csv` | glomeruli RDS files | Ensembl, STRINGdb | Kidney dataset |
| `analysis/07_mammary_gland.Rmd` | `mammary_gland_processing_for_Bodyclocks.rmd` | `data/raw/mammary_gland.txt` | mammary gland RDS files | AnnotationDbi, STRINGdb | Array annotation |
| `analysis/08_nih3t3.Rmd` | `nih3t3_bodyclocks.rmd` | `data/raw/NIH3T3_Forskolin.csv` | NIH3T3 RDS files | Ensembl, STRINGdb | 1 h sampling |
| `analysis/09_podocytes.Rmd` | `podocytes_dexamethasone.rmd` | `data/raw/podocyte_data.csv`, `data/meta/podocyte_meta.csv` | podocyte RDS files | Ensembl, STRINGdb | Audit before deeper refactor |
| `analysis/10_xiphoid.Rmd` | `xipphoid_for_bodyclocks.rmd` | `data/raw/xiphoid_array.csv` | xiphoid RDS files | AnnotationDbi, STRINGdb | Filename spelling fixed |
| `analysis/11_mouse_atlas.Rmd` | `mouse_atlas_processing_for_BodyClocks.rmd` | `data/raw/GSE54650/*.csv` | mouse atlas RDS files | Ensembl, STRINGdb | Batch script |
| `analysis/12_baboon_atlas.Rmd` | `baboon_batch_processing_for_BodyClocks.rmd` | `data/raw/GSE98965_csvs/*.csv`, `data/meta/baboon_meta.csv` | baboon atlas RDS files | Ensembl REST/BioMart, STRINGdb | Batch script |
| `analysis/13_liver_rnaseq.Rmd` | `liver_RNAseq_batch.rmd` | `data/raw/GSE158600_genes_fpkm.csv` | liver RNA-seq RDS files | Ensembl, STRINGdb | Preserve 3-pass symbol mapping |
| `analysis/comparison_analysis_paper.rmd` | Consolidated paper comparison analysis | `results/tables/*_circ.csv`, `results/tables/*_enrichment.csv`, chondrocyte/cartilage Shiny RDS files | Manuscript figures and supplementary tables in `results/publication_figures/` | None | Paper-only comparison figures; run via `--scripts=paper` after scripts 01, 03, 04, and 05 |
| `analysis/pairwise_similarity_analysis.R` | New standalone pairwise similarity analysis | Display tables from `results/shiny_data/` via `config/output_manifest.csv` | Pairwise similarity CSVs/heatmaps in `results/pairwise_similarity/`; selected Fig1E/F copies in `results/publication_figures/` | None | Computes Jaccard, fold enrichment, overlap significance, and chondrocyte-vs-tissue comparison; run as part of `--scripts=paper` |


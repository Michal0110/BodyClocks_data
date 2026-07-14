# Script Inventory

This inventory summarises the active analysis scripts, their source scripts,
inputs, generated outputs, and external annotation services.

For every app-selectable dataset, the network outputs include an enrichment
bundle containing both dataset-tested and all-STRING backgrounds. The former is
defined by the unique genes actually supplied to rhythmicity testing after the
dataset-specific expression and probe filters. Both variants include GO
Biological Process, KEGG, Reactome and WikiPathways terms at FDR <= 0.05; legacy
unsuffixed enrichment RDS files contain the all-STRING variant.

| Analysis script | Source script(s) | Inputs | Outputs | External services |
|---|---|---|---|---|
| `analysis/01_cartilage.Rmd` | `cartilage.Rmd`, `circ_cartilage_timecourse_bodyclocks.rmd` | `data/raw/wt_cartilage.csv`, `data/raw/ko_cartilage.csv` | Circadian-result and enrichment CSVs in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/cartilage_circ/` | Ensembl, STRINGdb |
| `analysis/02_tendon.Rmd` | `tendon.Rmd`, `circ_tendon_timecourse_bodyclocks.rmd` | `data/raw/tendon.csv` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/tendon/` | Ensembl, STRINGdb |
| `analysis/03_chondrocytes_dexamethasone.Rmd` | `chondrocytes_dexamethasone.Rmd` | `data/raw/dex_raw.csv`, `data/meta/dex_meta.csv` | Circadian-result and enrichment CSVs in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/dex/`; heatmaps in `results/dexamethasone/plots/` and `results/publication_figures/` | Ensembl, STRINGdb |
| `analysis/04_chondrocytes_heat_shock.Rmd` | `chondrocytes_heat_shock.rmd` | `data/raw/hs_raw.csv`, `data/meta/hs_meta.csv` | Circadian-result and enrichment CSVs in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/hs/`; heatmaps in `results/heat_shock/plots/` and `results/publication_figures/` | Ensembl, STRINGdb |
| `analysis/05_chondrocytes_osmotic_stress.Rmd` | `chondrocytes_osmotic_stress.rmd` | `data/raw/osmo_raw.csv`, `data/meta/osmo_meta.csv` | Circadian-result and enrichment CSVs in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/osmo/`; heatmaps in `results/osmo/plots/` and `results/publication_figures/` | Ensembl, STRINGdb |
| `analysis/06_glomeruli.Rmd` | `glomeruli_for_bodyclocks.rmd` | `data/raw/gloms_full.csv` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/gloms/` | Ensembl, STRINGdb |
| `analysis/07_mammary_gland.Rmd` | `mammary_gland_processing_for_Bodyclocks.rmd` | `data/raw/mammary_gland.txt` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/mammary_gland/` | AnnotationDbi, STRINGdb |
| `analysis/08_nih3t3.Rmd` | `nih3t3_bodyclocks.rmd` | `data/raw/NIH3T3_Forskolin.csv` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/nih3t3/` | Ensembl, STRINGdb |
| `analysis/09_podocytes.Rmd` | `podocytes_dexamethasone.rmd` | `data/raw/podocyte_data.csv`, `data/meta/podocyte_meta.csv` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/podocytes_dexamethasone/` | Ensembl, STRINGdb |
| `analysis/10_xiphoid.Rmd` | `xipphoid_for_bodyclocks.rmd` | `data/raw/xiphoid_array.csv` | Circadian-result CSV in `results/tables/`; app-ready RDS files in `results/shiny_data/mouse/xiphoid_cartilage/` | AnnotationDbi, STRINGdb |
| `analysis/11_mouse_atlas.Rmd` | `mouse_atlas_processing_for_BodyClocks.rmd` | `data/raw/GSE54650/*.csv` | Per-tissue circadian-result CSVs in `results/tables/`; per-tissue app-ready RDS files in `results/shiny_data/mouse/` | Ensembl, STRINGdb |
| `analysis/12_baboon_atlas.Rmd` | `baboon_batch_processing_for_BodyClocks.rmd` | `data/raw/GSE98965_csvs/*.csv`, `data/meta/baboon_meta.csv` | Per-tissue circadian-result CSVs in `results/tables/`; per-tissue app-ready RDS files in `results/shiny_data/baboon/` | Ensembl REST/BioMart, STRINGdb |
| `analysis/13_liver_rnaseq.Rmd` | `liver_RNAseq_batch.rmd` | `data/raw/GSE158600_genes_fpkm.csv` | Per-condition app-ready RDS files in `results/shiny_data/mouse/liver_RNAseq/` | Ensembl, STRINGdb |
| `analysis/comparison_analysis_paper.rmd` | Consolidated paper comparison analysis | Circadian-result CSVs in `results/tables/`; chondrocyte/cartilage app-ready RDS files; STRING and KEGG pathway mappings downloaded to `results/cache/pathway_annotations/` | Manuscript figure PDF/PNG files and supplementary/summary CSVs in `results/publication_figures/` | STRING downloads, KEGG REST API |
| `analysis/pairwise_similarity_analysis.R` | Standalone pairwise similarity analysis | Display-table RDS files listed in `config/output_manifest.csv` | Pairwise similarity CSVs and heatmap PDF/PNG files in `results/pairwise_similarity/`; Figure 1E/F copies in `results/publication_figures/` | None |

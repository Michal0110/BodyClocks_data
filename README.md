# BodyClocks Data Processing

This repository contains the analysis pipeline used to process raw circadian
transcriptomic datasets into RDS files consumed by the BodyClocks Shiny app,
and to generate the manuscript figures and supplementary tables.

The active analysis scripts are in `analysis/`. Raw input data are tracked with
Git LFS so the repository can remain self-contained without embedding large
files directly in Git history.

## Repository Layout

- `analysis/` — numbered dataset processing scripts.
- `R/` — shared utilities for I/O, annotation, BioMart, STRINGdb, and validation.
- `data/raw/` — raw input data files tracked with Git LFS.
- `data/meta/` — metadata files.
- `config/output_manifest.csv` — expected Shiny app RDS outputs.
- `data_manifest.csv` — raw/meta data inventory with SHA-256 checksums.
- `results/` — generated outputs; ignored by Git except `results/publication_figures/`.
- `results/publication_figures/` — manuscript-ready figures and supplementary tables.

## Environment

Install Git LFS before cloning or before adding raw data, then create and
activate the conda environment and restore renv:

```sh
git lfs install
conda env create -f environment.yml
conda activate bodyclocks_data
Rscript -e "renv::restore()"
```

Conda pins R and system libraries. renv pins R package versions inside that
environment.

## Running

Validate the pipeline structure:

```sh
Rscript run_analysis.R --dry-run
```

Run all analyses:

```sh
Rscript run_analysis.R
```

Run selected scripts:

```sh
Rscript run_analysis.R --scripts=03,09 --fail-fast
```

Run the standalone pairwise similarity analysis only:

```sh
Rscript run_analysis.R --scripts=pairwise --fail-fast
```

Reproduce only the paper comparison figures and tables:

```sh
Rscript run_analysis.R --scripts=paper --fail-fast
```

The `paper` target automatically runs the required upstream scripts first: `01` cartilage, `03` dexamethasone chondrocytes, `04` heat-shock chondrocytes, and `05` osmotic-stress chondrocytes. It then renders `analysis/comparison_analysis_paper.rmd` and runs `analysis/pairwise_similarity_analysis.R`. Use this target to regenerate the manuscript figures and supplementary tables in `results/publication_figures/` without running the full app data-processing pipeline.

## Publication Notes

Before depositing on Zenodo:

1. Confirm `data_manifest.csv` source/accession fields are complete.
2. Run `git lfs pull` so raw input data are present locally.
3. Confirm the archive contains full data files rather than LFS pointer files.
4. Include source code, raw/meta data, environment files, and `results/publication_figures/`; do not include the rest of `results/` unless you intentionally want to archive generated intermediates.
5. Record the final release tag and DOI in `docs/zenodo_deposit.md`.

External annotation services used by this pipeline, especially Ensembl BioMart
and STRINGdb, can fail transiently. Shared wrappers in `R/` are designed to
retry and fail loudly rather than silently producing partial annotation.


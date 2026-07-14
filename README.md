# BodyClocks data and analysis pipeline

This repository contains the reproducible R workflow supporting the BodyClocks
study of circadian transcriptomes. It processes curated mouse and baboon
transcriptomic datasets into a common set of circadian statistics, functional
annotations and STRING association networks; prepares data files consumed by
the [BodyClocks](https://www.bodyclocks.org) web application; and generates the
figures and supplementary tables accompanying the manuscript.

The workflow starts from the gene-level count or expression matrices supplied
in `data/raw/`. FASTQ quality control, trimming and genome alignment for newly
generated sequencing data are described in the manuscript but are not executed
by this repository.

## Repository contents

- `analysis/` — dataset-processing notebooks and manuscript comparison analyses.
- `R/` — shared functions for input/output, validation, annotation, BioMart and
  STRINGdb processing.
- `data/raw/` — input count and expression matrices, with large files managed by
  Git LFS.
- `data/meta/` — sample metadata.
- `data_manifest.csv` — input-file provenance, accession information and SHA-256
  checksums.
- `config/output_manifest.csv` — expected data products for the BodyClocks app.
- `results/publication_figures/` — manuscript figures and supplementary tables.
- `docs/` — data-source, environment, script-inventory and archival documentation.

The [script inventory](docs/script_inventory.md) lists the inputs, outputs and
external services used by each analysis.

## Installation

The reproducible environment combines conda, which supplies R and system
libraries, with renv, which pins R package versions. Install conda and Git LFS,
then run:

```sh
git lfs install
git lfs pull
conda env create -f environment.yml
conda activate bodyclocks_data
Rscript -e "renv::restore()"
```

The environment is based on R 4.5.2. See [environment documentation](docs/environment.md)
and `renv.lock` for the complete software specification.

## Reproducing the analysis

Run a structural validation before starting the analyses:

```sh
Rscript run_analysis.R --dry-run
```

Run the focused offline test for the background-aware enrichment schema and
STRING request contract with:

```sh
Rscript tests/test_enrichment_backgrounds.R
```

To reproduce the manuscript comparison figures and supplementary tables:

```sh
Rscript run_analysis.R --scripts=paper --fail-fast
```

The `paper` target runs all mouse dataset-processing scripts required by the
cross-dataset analysis; it does not process the baboon atlas. It then renders
the focal chondrocyte/cartilage comparison and runs the mouse pairwise-similarity
analysis that generates Figures 1E and 1F. The pairwise analysis applies its
configured exclusions to knockout datasets and xiphoid cartilage. Results are
written to `results/publication_figures/`.

To process every configured dataset and generate all app and publication data
products:

```sh
Rscript run_analysis.R --fail-fast
```

Individual analyses can be selected by their identifiers. For example:

```sh
Rscript run_analysis.R --scripts=03,09 --fail-fast
```

The standalone cross-dataset similarity analysis can be run after its required
dataset outputs have been generated:

```sh
Rscript run_analysis.R --scripts=pairwise --fail-fast
```

Generated reports, intermediate tables, caches and app-ready RDS files are
written below `results/`. Only publication-ready outputs are retained in version
control.

Each app network dataset includes enrichment results for two alternative
backgrounds: genes retained for rhythmicity testing in that dataset and all
proteins available for the species in STRING. Both variants use STRING v12.0,
FDR <= 0.05 and the GO Biological Process, KEGG, Reactome and WikiPathways
categories. The background-aware enrichment bundle is stored alongside the
network RDS files; the unsuffixed enrichment files retain the all-STRING result
for compatibility with earlier versions of the app.

After regenerating every dataset, validate the complete app artifact set and,
when the application repository is checked out alongside this repository,
synchronize it with:

```sh
Rscript -e 'source("R/validation_utils.R"); validate_rds_exists()'
rsync -a results/shiny_data/ ../BodyClocks/data/
```

## External resources and reproducibility

Some analyses query Ensembl BioMart and STRING v12.0, so an internet connection
is required when the corresponding cached annotations are unavailable. The
dataset-specific STRING enrichment requests use the unique genes that passed
the dataset's preprocessing filters and were supplied to RAIN as their
background; the rhythmic foreground retains the thresholds documented in each
analysis script. Shared wrappers retry transient failures and stop rather than
silently returning partial annotations. For long-term reproducibility, release
archives should include the input matrices, metadata, manifests, environment
files and final publication outputs.

The machine-readable data inventory is `data_manifest.csv`; source information
is described further in the [data-source documentation](docs/data_sources.md).
The analysis uses stable Ensembl identifiers for cross-dataset comparisons where
available and supplies MGI or external gene symbols for display.

## Related resources

- Web application: [BodyClocks.org](https://www.bodyclocks.org)
- Application source code: [Michal0110/BodyClocks](https://github.com/Michal0110/BodyClocks)
- Analysis source code: [Michal0110/BodyClocks_data](https://github.com/Michal0110/BodyClocks_data)

## Citation

If you use this workflow or its processed data, please cite the associated
publication and the archived software release. Machine-readable citation
metadata are provided in [`CITATION.cff`](CITATION.cff).

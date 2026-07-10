# Environment

This project uses conda and renv together.

Conda provides R 4.5.2 and system libraries such as libcurl, OpenSSL,
libxml2, and pandoc. renv then locks the R packages used inside that conda
environment.

Recommended setup:

```sh
conda env create -f environment.yml
conda activate bodyclocks_data
Rscript -e "renv::restore()"
Rscript run_analysis.R --dry-run
```

The intended order is conda first, then renv. Do not create the renv library
with one R installation and later run it from a different R installation.

To refresh documentation after package changes:

```sh
Rscript -e "renv::snapshot()"
Rscript -e "writeLines(capture.output(sessionInfo()), 'docs/sessionInfo.txt')"
Rscript -e "writeLines(capture.output(renv::diagnostics()), 'docs/renv_diagnostics.txt')"
conda env export --from-history > docs/conda_from_history.yml
```


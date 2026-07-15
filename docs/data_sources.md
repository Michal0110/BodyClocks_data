# Data Sources

The canonical machine-readable inventory is `data_manifest.csv`.

Each row records the local file, dataset label, source accession or URL when
known, checksum, and whether the file should be included in the Zenodo archive.
Before publication, verify these fields and record any source-specific licensing
or redistribution restrictions in the notes.

The primary-chondrocyte dexamethasone and heat-shock datasets share the assigned
ArrayExpress accession `E-MTAB-17163`. The accession is recorded in the manifest,
but the ArrayExpress study is not yet publicly accessible.

## Downloaded pathway annotations

`analysis/comparison_analysis_paper.rmd` downloads the following reference
files when they are not already available locally and caches them under
`results/cache/pathway_annotations/`:

- STRING v12.0 mouse protein enrichment terms:
  `https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0/10090.protein.enrichment.terms.v12.0.txt.gz`
- STRING v12.0 mouse protein aliases:
  `https://stringdb-downloads.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz`
- KEGG mouse pathway names:
  `https://rest.kegg.jp/list/pathway/mmu`
- KEGG mouse pathway-to-gene mappings:
  `https://rest.kegg.jp/link/mmu/pathway`

The KEGG files are downloaded at analysis time and are not tracked in Git or
included in the Zenodo archive. Their recorded checksums identify the files
used for the publication analysis; KEGG content may change between retrievals.

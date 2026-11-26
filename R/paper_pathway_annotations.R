# Original paper notebook mapping of STRING and KEGG releases to Ensembl genes.
read_paper_pathway_annotations <- function() {
  string_terms_path <- resolve_reference_file(
    "10090.protein.enrichment.terms.v12.0.txt.gz",
    paste0(
      "https://stringdb-downloads.org/download/",
      "protein.enrichment.terms.v12.0/",
      "10090.protein.enrichment.terms.v12.0.txt.gz"
    )
  )
  string_aliases_path <- resolve_reference_file(
    "10090.protein.aliases.v12.0.txt.gz",
    paste0(
      "https://stringdb-downloads.org/download/protein.aliases.v12.0/",
      "10090.protein.aliases.v12.0.txt.gz"
    ),
    existing_paths = "results/cache/stringdb/10090.protein.aliases.v12.0.txt.gz"
  )
  kegg_pathway_names_path <- resolve_reference_file(
    "mmu_kegg_pathway_names_raw.tsv",
    "https://rest.kegg.jp/list/pathway/mmu"
  )
  kegg_pathway_genes_path <- resolve_reference_file(
    "mmu_kegg_pathway_genes_raw.tsv",
    "https://rest.kegg.jp/link/mmu/pathway"
  )

  string_terms <- read_filtered_tsv_gz(
    string_terms_path,
    pattern = "\t(Reactome Pathways|WikiPathways)\t",
    column_names = c(
      "string_protein_id", "string_category", "term", "description"
    )
  )

  string_aliases <- read_filtered_tsv_gz(
    string_aliases_path,
    pattern = "\t(Ensembl_gene|KEGG_GENEID)$",
    column_names = c("string_protein_id", "alias", "source")
  )

  # Collapse protein isoforms to unique Ensembl gene-to-term relationships.
  protein_gene_map <- string_aliases %>%
    filter(source == "Ensembl_gene", grepl("^ENSMUSG", alias)) %>%
    transmute(
      string_protein_id,
      gene_id = normalise_ensembl_gene_ids(alias)
    ) %>%
    filter(!is.na(gene_id)) %>%
    distinct()

  string_annotation_map <- string_terms %>%
    filter(string_category %in% c("Reactome Pathways", "WikiPathways")) %>%
    inner_join(protein_gene_map, by = "string_protein_id") %>%
    transmute(
      gene_id,
      category = recode(
        string_category,
        "Reactome Pathways" = "RCTM",
        "WikiPathways" = "WikiPathways"
      ),
      term,
      description
    ) %>%
    distinct()

  kegg_pathway_genes <- read.delim(
    kegg_pathway_genes_path,
    header = FALSE, sep = "\t", quote = "", comment.char = "",
    col.names = c("pathway_id", "kegg_gene_id"),
    stringsAsFactors = FALSE
  ) %>%
    transmute(
      term = sub("^path:", "", pathway_id),
      ncbi_gene_id = sub("^mmu:", "", kegg_gene_id)
    ) %>%
    distinct()

  kegg_pathway_names <- read.delim(
    kegg_pathway_names_path,
    header = FALSE, sep = "\t", quote = "", comment.char = "",
    col.names = c("term", "description"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(description = sub(
      " - Mus musculus \\(house mouse\\)$", "", description
    )) %>%
    distinct(term, .keep_all = TRUE)

  kegg_protein_map <- string_aliases %>%
    filter(source == "KEGG_GENEID") %>%
    transmute(string_protein_id, ncbi_gene_id = alias) %>%
    distinct()

  kegg_annotation_map <- kegg_pathway_genes %>%
    inner_join(kegg_protein_map, by = "ncbi_gene_id") %>%
    inner_join(protein_gene_map, by = "string_protein_id") %>%
    inner_join(kegg_pathway_names, by = "term") %>%
    transmute(gene_id, category = "KEGG", term, description) %>%
    distinct()

  pathway_annotation_map <- bind_rows(
    kegg_annotation_map,
    string_annotation_map
  ) %>%
    distinct(gene_id, category, term, .keep_all = TRUE)
  list(map = pathway_annotation_map, mapped = unique(na.omit(protein_gene_map$gene_id)))
}

"""Join corrected RAIN and CircaCompare on the exact manuscript hypotheses."""
from pathlib import Path

import numpy as np
import pandas as pd

package = Path(__file__).resolve().parent
root = package.parents[1]
rain_figs = root / "results/publication_figures"
cc_figs = package / "figures"
out = package / "results"
methods = ["RAIN", "CircaCompare"]


def join_tests(a, b, key):
    # Every hypothesis is eligible in both arms because the tested backgrounds match.
    wide = a.merge(b, on=key, how="outer", validate="one_to_one", indicator=True,
                   suffixes=("_RAIN", "_CircaCompare"))
    if not wide["_merge"].eq("both").all():
        raise ValueError(f"Unmatched hypotheses for {key}")
    return wide.drop(columns="_merge")


def transition(a, b, inclusive=False):
    # Individual pathway FDR uses ≤ 0.05. The other decision columns use < 0.05.
    passed_a = a.le(.05) if inclusive else a.lt(.05)
    passed_b = b.le(.05) if inclusive else b.lt(.05)
    return np.select([a.isna() | b.isna(), passed_a & passed_b, passed_a, passed_b],
                     ["unavailable", "both", "RAIN_only", "CircaCompare_only"], default="neither")


def compare_file(name, keys):
    return join_tests(pd.read_csv(root / "revision/rain" / name),
                      pd.read_csv(out / name), keys)


# Match comparisons on both the call threshold and the clock exclusion rule.
overlap = compare_file("overlap_tests.csv", ["pair", "threshold", "drop_clock"])
overlap["decision"] = transition(overlap.q_RAIN, overlap.q_CircaCompare)
overlap.to_csv(out / "gene_overlap_decisions.csv", index=False)


# Use the six-pair BH values for pairwise tests and the permutation p for the single
# three-culture test. Missing phase tests retain an explicit unavailable state.
phase = compare_file("phase_tests.csv", ["comparison", "rule"])
triple = phase.comparison.eq("Osmo_vs_Dex_vs_HS")
phase["decision"] = transition(phase.q_RAIN.where(~triple, phase.p_RAIN),
                               phase.q_CircaCompare.where(~triple, phase.p_CircaCompare))
phase.to_csv(out / "phase_decisions.csv", index=False)


# Keep a dataset's pathway tests separate for each pair-specific background. A pathway
# category overlap and an individual pathway enrichment have different units of analysis.
category_key = ["Background_Strategy", "Pair", "Category"]
term_key = ["Background_Strategy", "Pair_Context", "dataset", "category", "term"]
term_frames = []
category_frames = []
for directory in [rain_figs, cc_figs]:
    category = pd.read_csv(directory / "pathway_overlap_background_comparison.csv")
    assert len(category) == 36
    category_frames.append(category)
    term = pd.read_csv(directory / "all_eligible_pathway_tests.csv.gz")
    reported = pd.read_csv(directory / "pathway_enrichment_background_comparison.csv")
    # Compare even zero-hit terms, and check all significant IDs against the paper export.
    assert set(map(tuple, term.loc[term.fdr.le(.05), term_key].values)) == set(map(tuple, reported[term_key].values))

    # Recover each rhythmic foreground size in its own background. The all-mouse foreground
    # for a dataset is shared across its pairwise comparisons.
    sizes = []
    for i in (1, 2):
        f = category[["Background_Strategy", "Pair", f"Dataset{i}", f"Rhythmic_Genes_D{i}"]].copy()
        f.columns = ["Background_Strategy", "Pair_Context", "dataset", "foreground"]
        f.loc[f.Background_Strategy == "all_mouse", "Pair_Context"] = "All_pairs"
        sizes.append(f)
    contexts = ["Background_Strategy", "Pair_Context", "dataset"]
    sizes = pd.concat(sizes).drop_duplicates()
    term = term.merge(sizes, on=contexts, validate="many_to_one", how="left")
    assert term.foreground.notna().all()
    term["fold_enrichment"] = (term.number_of_genes / term.foreground /
                               (term.number_of_genes_in_background / term.Background_Genes))
    term_frames.append(term)


# Compare the two null models for category-level overlap, each with its own BH correction.
category = join_tests(*category_frames, category_key)
for name, column in [("permutation", "Perm_P_BH"), ("conditional", "Conditional_Hypergeom_P_BH")]:
    category[name + "_decision"] = transition(category[column + "_RAIN"], category[column + "_CircaCompare"])
category.to_csv(out / "pathway_overlap_decisions.csv", index=False)


# Require identical term definitions and eligible background sizes for individual pathways.
term = join_tests(*term_frames, term_key)
for column in ("description", "Background_Genes", "number_of_genes_in_background"):
    assert term[column + "_RAIN"].fillna("").equals(term[column + "_CircaCompare"].fillna(""))
    term[column] = term.pop(column + "_RAIN")
    term.drop(columns=column + "_CircaCompare", inplace=True)
term["decision"] = transition(term.fdr_RAIN, term.fdr_CircaCompare, inclusive=True)
term.to_csv(out / "individual_pathway_decisions.csv.gz", index=False)
term.groupby(term_key[:4] + ["decision"]).size().rename("count").reset_index().to_csv(
    out / "individual_pathway_counts.csv", index=False)


# Match gene-level decisions to the raw RAIN p values from the corrected inputs.
genes = []
for dataset, filename in [("Osmo","osmo"), ("Dex","dex"), ("HS","hs"), ("Cartilage","cartilage_wt")]:
    f = pd.read_csv(package / "calls" / f"single_{filename}.csv")
    reference = pd.read_csv(root / "results/tables" / (("wt_cartilage" if filename=="cartilage_wt" else filename) + "_circ.csv"))
    id_col = "CycID" if dataset == "Cartilage" else "Row.names"
    f = f.merge(reference[[id_col,"pVal"]].rename(columns={id_col:"ensembl_gene_id", "pVal":"rain_p"}),
                on="ensembl_gene_id", validate="one_to_one")
    f["dataset"] = dataset
    f["decision"] = transition(f.rain_q, f.rhythmic_q)
    genes.append(f[["dataset","ensembl_gene_id","symbol","rain_p","rain_q","rhythmic_p","rhythmic_q","decision"]])
pd.concat(genes).to_csv(out / "rhythmicity_decisions.csv.gz", index=False)

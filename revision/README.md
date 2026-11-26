# Corrected RAIN manuscript figures

The revision compares osmotic stress, dexamethasone and heat shock in cultured mouse chondrocytes with WT cartilage. The treatment protocols, study, cell state and sampling design can all contribute to observed differences.

## Reproduce from this worktree

Use the repository R environment and run from the repository root. The raw Dex count file requires Git LFS. Python's standard library is sufficient for input preparation.

```sh
python3 revision/prepare_inputs.py
Rscript revision/rebuild_dex.R
Rscript revision/render_paper.R
Rscript analysis/pairwise_similarity_analysis.R
Rscript revision/compare_focal.R RAIN . revision/rain
```

The four compressed expression/RAIN tables, sample times and SHA-256 manifest in [inputs](inputs) also provide a common starting point for alternative detectors. The Dex replay evaluates the original notebook against frozen annotation and checks the SHA-256 of the full table. Fig. 1E/F use the saved mouse displays, including the recomputed Dex display. The annotation snapshot contains the paper notebook's mapped STRING v12 and KEGG gene sets. Both revised analyses use this snapshot. Two KEGG source files differ from the submission manifest. [Source digests](inputs/pathway_references.tsv)

## Manuscript update checklist

| Location | Revision to use |
|---|---|
| Dex methods and gene counts | The 35 retained samples include three at 24 h and two at 28 h. Corrected RAIN detects 2,043 of 12,323 genes at BH < 0.05 (submitted 2,008). 326 gained and 291 lost the call. Of the 1,717 retained calls, 220 changed peak phase. At BH ≤ 0.01 the count changes from 596 to 624. |
| Fig. 1B–D | Replace Dex heatmap, four-way traces and overlap matrix. Four-way overlap changes from 15 to 16 with the addition of *Lgalsl*. Nine of 16 belong to the paper's clock panel. Recompute Venn labels from exclusive Ensembl-ID regions. |
| Fig. 1E/F | Replace the fold-enrichment figures and use the revised mouse pair/group tables. Dex–podocytes crosses the BH threshold (0.117 to 0.00224). The 276 pairs excluding Dex retain their overlaps and effects. Global adjustment can change their q values. These figures use the original symbol-based display inputs. |
| Fig. 2A–E | Replace phase distributions, pairwise clock-relative comparisons and clock-panel heatmap. Use the three-culture Venn counts for 2D. Triple rhythmic genes change from 64 to 62 and concordant genes from 23 to 20 (35.9% to 32.3%, permutation p from 0.00070 to 0.00310). |
| Fig. 3, Supplementary Fig. 1 and Tables 1–4 | Use the regenerated pathway panels and all four supplementary tables. Preserve the distinction between the pair-specific tested background and the all-mouse mapped background. |

[Manuscript figures and tables](../results/publication_figures/) · [Exclusive Venn counts](rain/venn_regions.csv) · [Gene membership](rain/rhythmic_membership.csv) · [Submitted-to-corrected differences](rain/dex_correction/) · [Mouse comparison statistics](rain/mouse_pairs/)

The correction tables compare original and corrected timing with matched references. The rebuilt full Dex table, time-series RDS and display RDS matched the earlier corrected replay byte for byte. All nine shared paper CSVs, including the four supplementary tables, matched the corrected replay cell for cell.

The shared comparison script calls the paper notebook's hypergeometric and clock-relative functions and uses its clock panel. Overlap uses the genes tested in both studies, a one-sided upper tail, and BH across the six pairs. Clock exclusion removes genes from both lists and the universe. The primary phase test requires rhythmic Bmal1 in each dataset. Pairwise concordance uses ≤4 h on a 24 h circle and 1,000 permutations (empirical p, six-pair BH). Three-culture concordance requires all three distances ≤4 h and uses 10,000 permutations with a plus-one correction. The additional anchor-exclusion, clock-panel-only and tested-anchor results are labelled in the phase tables.


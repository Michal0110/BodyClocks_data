# CircaCompare sensitivity analysis for the manuscript revision

Corrected RAIN and CircaCompare support the conclusion that the rhythmic gene lists differ substantially between synchronisation conditions. Both methods identify the same three modestly enriched pairwise overlaps. The proportion of four-way genes belonging to the clock panel, and the evidence for excess three-culture phase concordance after Bmal1 alignment, depend on the detector.

## Rhythmic genes and agreement with corrected RAIN

All calls used a fixed period of 24 h, the same normalised expression and sample-level times, all biological replicates, and BH q < 0.05 within each dataset. The primary CircaCompare analysis uses the amplitude-coefficient p value returned by the pinned package's single-dataset fit.

| Dataset | Tested | Corrected RAIN | CircaCompare | Both | % of CircaCompare in RAIN | % of RAIN in CircaCompare |
|---|---:|---:|---:|---:|---:|---:|
| Osmotic stress | 12,398 | 879 | 3,738 | 775 | 20.7 | 88.2 |
| Dexamethasone | 12,323 | 2,043 | 2,252 | 1,199 | 53.2 | 58.7 |
| Heat shock | 12,375 | 5,255 | 4,321 | 3,225 | 74.6 | 61.4 |
| WT cartilage | 12,559 | 1,678 | 1,036 | 721 | 69.6 | 43.0 |

Within each dataset, agreement between detectors exceeded chance, including after removing the clock panel. As an additional sensitivity check, jointly testing the sine and cosine coefficients against an intercept-only model produced 457, 62, 611 and 128 discoveries, respectively. This check uses the same 24-hour waveform with two numerator degrees of freedom. At zero amplitude the fitted peak time is unidentified. [All gene-level decisions](results/rhythmicity_decisions.csv.gz) · [Within-dataset overlap](results/within_dataset_overlap_with_RAIN.csv) · [Additional harmonic test](calls/dataset_summary.csv)

## Gene-list overlap across conditions

| Pair | CircaCompare observed / expected | Fold over chance | BH q | Clock genes in overlap | BH q after clock exclusion |
|---|---:|---:|---:|---:|---:|
| Osmo–Dex | 686 / 661.8 | 1.04 | 0.166 | 9 | 0.239 |
| Osmo–HS | 1,273 / 1,260.0 | 1.01 | 0.359 | 9 | 0.427 |
| Osmo–cartilage | 240 / 277.3 | 0.87 | 0.998 | 5 | 0.999 |
| Dex–HS | 837 / 774.9 | 1.08 | 0.00375 | 9 | 0.00698 |
| Dex–cartilage | 194 / 168.5 | 1.15 | 0.0277 | 5 | 0.0549 |
| HS–cartilage | 367 / 322.3 | 1.14 | 0.00375 | 4 | 0.00678 |

The six primary decisions agree with corrected RAIN. After clock exclusion, Dex–cartilage remains significant with RAIN alone. At the stricter rhythmicity cutoff of q ≤ 0.01, HS–cartilage passes only with RAIN when all genes are included. Dex–HS passes only with RAIN after removing clock genes at that cutoff. The tests use the genes assayed in both datasets and assess excess overlap. A high p value leaves the random-overlap hypothesis compatible with the data. [Full matched tests](results/gene_overlap_decisions.csv) · [Clock-gene identities](results/overlap_tests.csv)

Across the 10,764 genes tested in all four datasets, the four-way intersection contains 22 CircaCompare calls, including four clock-panel genes, versus 16 corrected RAIN calls, including nine clock-panel genes. Five genes survive both analyses: *Bmal1*, *Cry1*, *Dbp*, *Npas2* and *Sema4d*. A clock-rich subset therefore recurs across detectors. Calling the entire four-way intersection predominantly clock genes is supported by corrected RAIN (9/16) and would overstate the CircaCompare result (4/22). The [submitted manuscript](https://doi.org/10.64898/2026.07.20.739603) reported 15 four-way RAIN genes. [Four-way gene identities](calls/four_way_genes.csv) · [Exclusive Venn counts](results/venn_regions.csv)

## Clock-relative phase

Relative to the fitted Bmal1 peak in each culture, 14/246 (5.7%) shared CircaCompare genes had peaks within four hours in **every** culture pair. The shuffled mean was 6.7% (p = 0.767, 10,000 permutations). Corrected RAIN gave 20/62 (32.3%, shuffled mean 17.9%, p = 0.00310). Submitted RAIN gave 23/64 (35.9%, p = 0.00070). Removing Bmal1 itself gave 13/245 with CircaCompare (p = 0.848) and 19/61 with corrected RAIN (p = 0.00600). Among the eight clock-panel genes shared by the three cultures under CircaCompare, 4/8 met the criterion (p = 0.00860, exploratory subset).

For the six pairwise comparisons, Osmo–Dex gained significant concordance (286/686, 41.7% versus a 30.4% null). Dex–HS lost it (250/837, 29.9% versus 32.0% null). Dex–cartilage and HS–cartilage retained it. For Osmo–HS and Osmo–cartilage the observed concordance was compatible with the null under both methods. Zero empirical pairwise p values in the files mean zero exceedances in 1,000 permutations.

Gene selection and phase estimates both contribute. Of 23 genes rhythmic in all three cultures under both detectors, 9 were concordant using RAIN phases, 7 using CircaCompare phases and 5 using both. RAIN returns phases on the sampling grid. CircaCompare estimates continuous peak times. [Combined manuscript panels](figures/CircaCompare_figures.pdf) · [All phase tests](results/phase_decisions.csv) · [Matched-gene classifications](calls/clock_relative_matched_genes.csv) · [Phase figures and supplementary tables](figures/)

## Functional enrichment

The same manuscript code, mapped annotation snapshot, eligibility filters, list-size permutation null, seeds and BH families were used for each detector. With pair-specific tested backgrounds, **0/18 pair-by-category overlap tests passed BH under either detector**. With the all-mouse mapped background, the submitted paper reported 11/18. The original timing replay with the revised references gave 12/18, corrected RAIN gave 12/18 and CircaCompare gave 14/18. Ten passed with both and two passed with neither. The six changes in permutation BH q are listed below.

| Pair | Category | Corrected RAIN q | CircaCompare q | Significant with |
|---|---|---:|---:|---|
| Osmo–Dex | KEGG | 0.159 | 0.0245 | CircaCompare |
| Osmo–Dex | WikiPathways | 0.140 | 0.0160 | CircaCompare |
| Osmo–cartilage | KEGG | 0.318 | 0.0437 | CircaCompare |
| Osmo–cartilage | Reactome | 0.514 | 0.0245 | CircaCompare |
| Osmo–HS | KEGG | 0.00654 | 0.511 | RAIN |
| HS–cartilage | WikiPathways | 0.00200 | 0.198 | RAIN |

The secondary conditional hypergeometric test changed decisions in the same six comparisons. The background dependence is robust. Individual pathway and annotation-category claims require qualification. Two KEGG reference source files differ from the July manifest. Both revised arms used the same frozen map, so these detector comparisons hold the annotation fixed. [All 36 tests, both nulls and effect sizes](results/pathway_overlap_decisions.csv) · [Reference digests](../../revision/inputs/pathway_references.tsv)

At individual-term FDR ≤ 0.05, the all-mouse analysis had 513 dataset–term calls in both, 354 RAIN-only and 471 CircaCompare-only (6,606 in neither). With the pair-specific backgrounds there were 40 calls in both, 88 RAIN-only, 231 CircaCompare-only and 20,677 in neither. Here a call is a dataset–term test **within its background context**. A pathway tested for the same dataset against several pair backgrounds contributes several tests. The KEGG circadian rhythm pathway was enriched in all four studies with both methods under the all-mouse background. [All 28,980 term tests including zero hits, effect sizes and exact IDs](results/individual_pathway_decisions.csv.gz) · [Counts by dataset, context and category](results/individual_pathway_counts.csv)

## Native between-group CircaCompare tests

| Pair | Rhythmic in both | Mesor difference | Amplitude difference | Peak-time difference |
|---|---:|---:|---:|---:|
| Osmo–Dex | 686 | 603 | 7 | 161 |
| Osmo–HS | 1,273 | 1,008 | 0 | 289 |
| Osmo–cartilage | 240 | 200 | 59 | 145 |
| Dex–HS | 837 | 581 | 0 | 426 |
| Dex–cartilage | 194 | 180 | 146 | 149 |
| HS–cartilage | 367 | 308 | 145 | 243 |

The native fits use the package's nominal single-series rhythmicity gate in both datasets (14,528 eligible pair fits). The table uses the subset with BH q < 0.05 for rhythmicity in both groups. BH is applied across that subset for each endpoint and pair. Effects are group 2 minus group 1. Peak-time differences use hours after synchronisation for cultures and ZT for cartilage. Mesor and amplitude can reflect study-specific expression scales or tissue composition. Pairwise peak differences and Bmal1-relative concordance address different hypotheses. [All eligible fits, both adjustment families and effect sizes](calls/) · [Pair summaries](calls/pair_summary.csv) · [CircaCompare package](https://github.com/RWParsons/circacompare/tree/5df1144f5b5127893438ec0fc42e72cc3ac75ec1)

## Reproduction

From the repository root, use the repository R environment and Python with NumPy, pandas and SciPy. The figure workflow starts from the four bundled fit tables. It invokes the original comparison notebook and the original three Fig. 1B heatmap chunks. All 49,655 individual fits and all 14,528 native pairwise fits converged. A fresh fit reproduced all reported fit fields. Independent harmonic least-squares checks reproduced the single and pairwise significance classifications. [Fit verification](results/fit_verification.txt) · [Corrected RAIN inputs](../../revision/README.md)

```sh
python3 revision/prepare_inputs.py
p=alternative_processing/circacompare
Rscript "$p/summarize_fits.R" "$p/calls" results/tables
Rscript revision/compare_focal.R CircaCompare "$p/calls" "$p/results"
Rscript revision/render_paper.R "$p/calls" "$p/figures"
Rscript "$p/heatmaps.R"
Rscript revision/plot_venn.R "$p/results" "$p/figures"
Rscript revision/plot_sensitivity.R CircaCompare "$p/results" "$p/figures"
python3 "$p/compare_tests.py"
python3 "$p/verify_fits.py"
# Figure checks use Pillow and pdftotext. pdfunite combines manuscript panels.
python3 "$p/verify_figures.py"
pdfunite "$p"/figures/Fig*.pdf "$p"/figures/Supplementary_Figure1.pdf "$p"/figures/CircaCompare_figures.pdf

# To repeat the native fit in a fresh output directory
git clone https://github.com/RWParsons/circacompare.git /tmp/circacompare
git -C /tmp/circacompare checkout 5df1144f5b5127893438ec0fc42e72cc3ac75ec1
mkdir -p results/circacompare/library
Rscript -e 'install.packages("/tmp/circacompare", lib="results/circacompare/library", repos=NULL, type="source")'
Rscript "$p/fit.R" results/circacompare_repeat 4
```


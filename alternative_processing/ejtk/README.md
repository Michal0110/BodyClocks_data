# Fixed-24-hour eJTK sensitivity analysis

Across the three chondrocyte cultures, 92–99% of eJTK calls were also rhythmic under corrected RAIN. The six genes called rhythmic in all three cultures were *Cry1*, *Npas2*, *Nr1d1*, *Nr1d2*, *Per3* and *Tef*. All six belong to the manuscript's 12-gene clock panel. WT cartilage had zero eJTK discoveries among 12,559 tested genes. This limits inference about four-way overlap and phase alignment in cartilage.

| Dataset | Tested | eJTK | Corrected RAIN | Both | % of eJTK in RAIN | % of RAIN in eJTK |
|---|---:|---:|---:|---:|---:|---:|
| Osmotic stress | 12,398 | 354 | 879 | 325 | 91.8 | 37.0 |
| Dexamethasone | 12,323 | 65 | 2,043 | 63 | 96.9 | 3.1 |
| Heat shock | 12,375 | 886 | 5,255 | 874 | 98.6 | 16.6 |
| WT cartilage | 12,559 | 0 | 1,678 | 0 | undefined | 0 |

Within each culture, agreement with corrected RAIN exceeded random expectation (13.0×, 5.85× and 2.32× for Osmo, Dex and HS). After removing the clock panel from lists and background, 316, 53 and 866 calls still overlapped RAIN, respectively (all four-dataset BH q < 10⁻³⁸). Most of the individual eJTK calls in these cultures lie outside the clock panel. [Within-dataset overlap](results/within_dataset_overlap_with_RAIN.csv) · [Figure](figures/detection_and_overlap.pdf)

## Overlap between synchronisation conditions

| Pair | Observed / expected | Clock genes | Enrichment | BH q | Without clock, observed / expected | Without clock, BH q |
|---|---:|---:|---:|---:|---:|---:|
| Osmo–Dex | 11 / 1.83 | 7 | 6.03× | 0.00000977 | 4 / 1.49 | 0.185 |
| Osmo–HS | 27 / 24.29 | 7 | 1.11× | 0.622 | 20 / 23.46 | 1.00 |
| Dex–HS | 15 / 4.56 | 7 | 3.29× | 0.0000957 | 8 / 3.81 | 0.185 |
| Each of the three cartilage pairs | 0 / 0 | 0 | undefined | 1.00 | 0 / 0 | 1.00 |

Two culture pairs showed above-chance overlap before clock exclusion. None passed six-pair BH after exclusion. For Dex–HS the non-clock overlap was nominally enriched (raw p = 0.0346). The identities of the seven clock genes differ by pair. [Pairwise tests and clock-gene names](results/overlap_tests.csv) · [Gene-level overlaps](results/overlap_genes.csv)

## Bmal1-relative phase

Dex–HS was the only pair with shared calls and a rhythmic Bmal1 anchor in both studies. Nine of 15 genes peaked within four hours after alignment (60%, shuffled mean 37.1%, 55 of 1,000 permutations at least as concordant, p = 0.055, six-pair BH q = 0.330). Excluding Bmal1 gave 8/14 (q = 0.510). Within the clock panel, 5/7 were concordant (q = 0.720). Osmo Bmal1 had q = 0.0618. The primary three-culture phase statistic was therefore unavailable. In the sensitivity analysis using the measured Osmo Bmal1 phase, 0/6 shared clock genes met all three pairwise distance criteria (p = 1). [Anchor phases](results/bmal1.csv) · [Phase tests](results/phase_tests.csv) · [Per-gene classifications](results/phase_genes.csv)

## Reproduction

The analysis uses all retained replicates, the corrected Dex timestamps, a fixed 24 h period, a two-hour phase/asymmetry grid, 1,000 empirical null draws per dataset and within-dataset BH on GammaP. The paper threshold is q < 0.05. The tables also include the manuscript's q ≤ 0.01 overlap analysis. The archived call files and a fresh fit from the pinned upstream code matched byte for byte for all four datasets. The comparison uses the same downstream R functions as corrected RAIN. [Settings](calls/settings.json) · [Inputs and corrected RAIN workflow](../../revision/README.md)

Run from the repository root in the project R environment. Python 3.12 dependencies are recorded in [python_packages.txt](calls/python_packages.txt).

```sh
python3 revision/prepare_inputs.py
p=alternative_processing/ejtk
python3 "$p/verify_calls.py"
Rscript revision/compare_focal.R eJTK "$p/calls" "$p/results"
Rscript revision/plot_sensitivity.R eJTK "$p/results" "$p/figures"
Rscript revision/plot_venn.R "$p/results" "$p/figures"

# To repeat the fit into a fresh output directory
git clone https://github.com/alanlhutchison/empirical-JTK_CYCLE-with-asymmetry.git /tmp/ejtk-upstream
git -C /tmp/ejtk-upstream checkout 0159108946cab048ed0ec77810d793c6be0f4166
python3 "$p/run_ejtk.py" /tmp/ejtk-upstream results/tables results/ejtk_repeat --periods 24 --grid-step 2
```


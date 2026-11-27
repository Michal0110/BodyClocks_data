"""Verify the four archived eJTK exports against the corrected RAIN inputs."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


# Each dataset has its own tested genes, sample layout and BH adjustment.
package = Path(__file__).resolve().parent
root = package.parents[1]
inputs = root / "revision/inputs"
calls = package / "calls"
manifest = pd.read_csv(inputs / "inputs.csv").set_index("dataset")
layout = pd.read_csv(inputs / "sample_hours.csv")
summaries = pd.read_csv(calls / "summary.csv").set_index("dataset")
assert set(summaries.index) == set(manifest.index)


for dataset, item in manifest.iterrows():
    # Check file contents as well as names. Both Dex exports are called dex_circ.csv,
    # but only the corrected export has the checksum recorded in the revision manifest.
    source = root / "results/tables" / item.file
    summary = summaries.loc[dataset]
    checksum = hashlib.sha256(source.read_bytes()).hexdigest()
    assert checksum == item.sha256 == summary.rain_input_sha256
    id_column = "CycID" if dataset == "cartilage_wt" else "Row.names"
    rain = pd.read_csv(source).set_index(id_column)

    # Match by sample ID so sorting cannot hide a misassigned replicate. Compare
    # both the recorded collection hour and the hour on the folded 24 h clock.
    dataset_dir = calls / dataset
    samples = pd.read_csv(dataset_dir / "samples.csv").set_index("sample").sort_index()
    expected = layout.loc[layout.dataset.eq(dataset)].drop(columns="dataset")
    pd.testing.assert_frame_equal(samples, expected.set_index("sample").sort_index())


    # Compare every tested gene, including genes above threshold. GammaP already
    # accounts for the template search. Apply BH to that column within this dataset.
    genes = pd.read_csv(dataset_dir / "gene_comparison.csv").set_index("ensembl_gene_id")
    assert genes.index.is_unique and set(genes.index) == set(rain.index)
    assert genes.Period.eq(24).all()
    np.testing.assert_allclose(genes.RAIN_BH, rain.loc[genes.index, "BH.Q"], atol=1e-15, rtol=1e-10)
    expected_gamma_bh = false_discovery_control(genes.GammaP)
    np.testing.assert_allclose(genes.GammaBH, expected_gamma_bh, atol=1e-15, rtol=1e-10)

    # Reproduce the paper's strict q < 0.05 count from the numeric q values.
    rhythmic_count = genes.GammaBH.lt(0.05).sum()
    assert rhythmic_count == summary.eJTK_significant
    print(dataset, len(genes), rhythmic_count)

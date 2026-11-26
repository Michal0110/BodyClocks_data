"""Verify the four fixed-period exports against the corrected inputs and recomputed BH."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

pkg = Path(__file__).resolve().parent
root = pkg.parents[1]
manifest = pd.read_csv(root / "revision/inputs/inputs.csv").set_index("dataset")
layout = pd.read_csv(root / "revision/inputs/sample_hours.csv")
summaries = pd.read_csv(pkg / "calls/summary.csv").set_index("dataset")
assert set(summaries.index) == set(manifest.index)


# Check each fit against the same corrected expression file and replicate times used by RAIN.
for dataset, item in manifest.iterrows():
    source = root / "results/tables" / item.file
    assert hashlib.sha256(source.read_bytes()).hexdigest() == item.sha256 == summaries.loc[dataset,"rain_input_sha256"]
    x = pd.read_csv(source).set_index("CycID" if dataset == "cartilage_wt" else "Row.names")
    folder = pkg / "calls" / dataset
    samples = pd.read_csv(folder / "samples.csv").set_index("sample").sort_index()
    expected = layout[layout.dataset.eq(dataset)].drop(columns="dataset").set_index("sample").sort_index()
    pd.testing.assert_frame_equal(samples, expected)

    # Match tested genes, verify the fixed period and independently adjust GammaP across genes.
    f = pd.read_csv(folder / "gene_comparison.csv").set_index("ensembl_gene_id")
    assert f.index.is_unique and set(f.index) == set(x.index) and f.Period.eq(24).all()
    np.testing.assert_allclose(f.RAIN_BH, x.loc[f.index, "BH.Q"], atol=1e-15, rtol=1e-10)
    np.testing.assert_allclose(f.GammaBH, false_discovery_control(f.GammaP), atol=1e-15, rtol=1e-10)
    assert f.GammaBH.lt(.05).sum() == summaries.loc[dataset,"eJTK_significant"]
    print(dataset, len(f), f.GammaBH.lt(.05).sum())

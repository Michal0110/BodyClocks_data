#!/usr/bin/env python3
"""Detect rhythmic genes with empirical JTK_CYCLE (eJTK) on the four RAIN expression tables.

For each gene, eJTK compares the ranks of its observations with rhythmic templates.
The winning template supplies peak phase and trough time. The null distribution
accounts for searching over templates. Multiplicity across genes is handled within
each dataset. The paper's overlap and Bmal1-relative tests run in revision/compare_focal.R.

Requires Python 3.12, NumPy, SciPy, pandas, statsmodels, Cython and setuptools.
Reviewer run from the repository root:
    python alternative_processing/ejtk/run_ejtk.py UPSTREAM_CLONE RAIN_TABLE_DIR NEW_OUTPUT_DIR

The defaults reproduce the fixed 24 h run. See --help for the optional grid and sample layout.
"""

import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import comb
from pathlib import Path
import re
import runpy
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.special import erfc
from scipy.stats import false_discovery_control


REPO_ROOT = Path(__file__).resolve().parents[2]
UPSTREAM_COMMIT = "0159108946cab048ed0ec77810d793c6be0f4166"
SEED = 20261124  # Restart the upstream random generator at this seed for each dataset.
CUTOFF = 0.05    # GammaBH threshold for the native eJTK summary.
PERIOD = 24      # Reference period in hours for the fixed-period reviewer run.
TABLES = {
    "osmo": "osmo_circ.csv",
    "dex": "dex_circ.csv",
    "hs": "hs_circ.csv",
    "cartilage_wt": "wt_cartilage_circ.csv",
}
# Column order and ZT mapping from the cartilage notebook's RAIN input.
CARTILAGE_HOURS = {
    "WT_1AM_2": 18, "WT_5AM_2": 22, "WT_9AM_2": 2, "WT_1PM_1": 6,
    "WT_5PM_2": 10, "WT_9PM_2": 14, "WT_1AM_3": 18, "WT_5AM_3": 22,
    "WT_9AM_3": 2, "WT_1PM_2": 6, "WT_5PM_3": 10, "WT_9PM_3": 14,
}


def sample_hours(dataset):
    """Return elapsed hours since synchronisation for cultures and Zeitgeber hours for cartilage."""
    if dataset == "cartilage_wt":
        return CARTILAGE_HOURS

    # Use elapsed time by sample ID. This keeps the Dex 24 h and 28 h replicates
    # aligned with corrected RAIN. dex_24h_b and dex_24h_c both belong at 24 h.
    # The excluded Dex outlier is dex_28h_b.
    with (REPO_ROOT / "data/meta" / f"{dataset}_meta.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    sample_column, hour_column = ("symbol", "timepoint") if dataset == "osmo" else ("sample", "time")
    hours = {}
    for row in rows:
        hour = int(row[hour_column].removesuffix("h"))
        if dataset == "osmo":
            if hour == 0:
                continue  # RAIN excluded the three 0 h observations.
            # R prefixes column names that start with a digit when it exports Osmo samples.
            hours["X" + row[sample_column]] = hour
        else:
            hours[row[sample_column]] = hour
    return hours


def prepare_dataset(dataset, rain_file, dataset_dir, legacy_grid, sample_times=None):
    frame = pd.read_csv(rain_file)
    gene_column = "CycID" if dataset == "cartilage_wt" else "Row.names"
    if frame[gene_column].isna().any() or not frame[gene_column].is_unique:
        raise ValueError(f"Missing or duplicate Ensembl ID in {rain_file}")
    hours = sample_hours(dataset) if sample_times is None else sample_times
    samples = sorted(hours, key=hours.__getitem__)
    expression = frame.loc[:, samples]  # Normalised expression used by RAIN, one column per replicate.
    if not np.isfinite(expression.to_numpy()).all() or np.any(expression.nunique(axis=1) < 2):
        raise ValueError(f"Missing or constant expression series in {rain_file}")


    # Folding maps culture samples at, for example, 4 h and 28 h to phase 4 while
    # retaining separate expression values. Cartilage hours are already on the ZT clock.
    # Period searches need elapsed hours so a non-24 h template can drift across cycles.
    time_points = [hours[s] % PERIOD if legacy_grid else hours[s] for s in samples]
    # ZT is an upstream header convention. For cultures these numbers are still
    # hours since synchronisation modulo 24. Repeated headers retain the replicate values.
    input_path = dataset_dir / "expression.txt"
    with input_path.open("w") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["#", *(f"ZT{hour}" for hour in time_points)])
        writer.writerows(zip(frame[gene_column], *expression.to_numpy().T))
    pd.DataFrame({
        "sample": samples,
        "elapsed_or_ZT_hour": [hours[s] for s in samples],
        "eJTK_input_hour": time_points,
    }).to_csv(dataset_dir / "samples.csv", index=False)
    return frame.set_index(gene_column), input_path


def run_template_search(engine, dataset_dir, input_path, legacy_grid):
    """Run the pinned waveform search and its null calibration for one sample layout."""
    raw_path = dataset_dir / "ejtk_jtkout.txt"
    upstream_args = [
        str(engine / "eJTK-CalcP.py"),
        "-f", str(input_path),
        "-o", str(raw_path),
        "-w", str(engine / "waveform_cosine.txt"),
        "-p", str(engine / "periods.txt"),
        "-s", str(engine / "phases.txt"),
        "-a", str(engine / "asymmetries.txt"),
    ]

    # Each dataset gets its own interpreter, random seed and log. The reviewer
    # configuration calls the upstream matcher for both expression and null profiles.
    seed_and_run = (
        "import sys; sys.path.insert(0, sys.argv.pop(1)); "
        "from run_ejtk import run_upstream; run_upstream()"
    )
    command = [
        sys.executable, "-B", "-c", seed_and_run, str(Path(__file__).resolve().parent),
        str(SEED), "legacy" if legacy_grid else "fast", *upstream_args,
    ]
    with (dataset_dir / "run.log").open("w") as log:
        subprocess.run(command, cwd=engine, stdout=log, stderr=subprocess.STDOUT, check=True)


def run_dataset(dataset, rain_dir, output, engine, periods, legacy_grid, sample_times=None):
    """Prepare samples, fit rhythmic templates, then compare calls with the matching RAIN table."""
    dataset_dir = output / dataset
    dataset_dir.mkdir()
    rain_file = rain_dir / TABLES.get(dataset, f"{dataset}_circ.csv")
    rain, input_path = prepare_dataset(dataset, rain_file, dataset_dir, legacy_grid, sample_times)
    run_template_search(engine, dataset_dir, input_path, legacy_grid)


    # The pinned engine scores 1,000 random normal profiles on the same sample layout
    # and template grid as the genes. CalcP takes the minimum of the plus-one empirical
    # tail, the template-count Bonferroni p (BF), and the fitted Gamma tail as GammaP.
    # GammaBH adjusts those values across all tested genes in this dataset.
    calibrated = pd.read_table(dataset_dir / "ejtk_jtkout_GammaP.txt", index_col="ID")
    if not calibrated.index.is_unique or set(calibrated.index) != set(rain.index):
        raise ValueError(f"eJTK and RAIN gene IDs differ for {dataset}")
    if not calibrated["Period"].isin(periods).all() or not calibrated["GammaP"].between(0, 1).all():
        raise ValueError(f"Invalid eJTK period or p-values for {dataset}")
    np.testing.assert_allclose(calibrated["GammaBH"], false_discovery_control(calibrated["GammaP"]),
                               rtol=1e-10, atol=1e-15)

    # Join to the corrected RAIN table by Ensembl ID before counting overlap.
    rain_columns = rain[["mgi_symbol", "pVal", "BH.Q"]].rename(
        columns={"pVal": "RAIN_p", "BH.Q": "RAIN_BH"}
    )
    joined = calibrated.join(rain_columns)

    # The original export flags q <= 0.05. The manuscript comparison applies q < 0.05
    # to the numeric columns. Both conventions give the same calls on these four tables.
    ejtk_significant = joined["GammaBH"] <= CUTOFF
    rain_significant = joined["RAIN_BH"] <= CUTOFF
    joined.assign(ejtk_significant=ejtk_significant, rain_significant=rain_significant).to_csv(
        dataset_dir / "gene_comparison.csv", index_label="ensembl_gene_id")
    return {
        "dataset": dataset,
        "tested_genes": len(joined),
        "RAIN_significant": rain_significant.sum(),
        "eJTK_significant": ejtk_significant.sum(),
        "shared": (ejtk_significant & rain_significant).sum(),
        "RAIN_only": (rain_significant & ~ejtk_significant).sum(),
        "eJTK_only": (ejtk_significant & ~rain_significant).sum(),
        "rain_input_sha256": hashlib.sha256(rain_file.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("upstream", type=Path, help="Clone of the pinned upstream eJTK repository.")
    parser.add_argument("rain_tables", type=Path, help="Normalised expression and corrected RAIN exports.")
    parser.add_argument("output", type=Path, help="New directory for this run and its gene-level tables.")
    parser.add_argument("--periods", nargs="+", type=int, default=[PERIOD],
                        help="Candidate periods in hours. The reviewer analysis fixes the period at 24.")
    parser.add_argument("--grid-step", type=int, choices=(1, 2),
                        help="Spacing of peak phases and peak-to-trough widths in hours (reviewer run: 2).")
    parser.add_argument("--sample-hours", type=Path, help="CSV with dataset, sample and hour columns")
    args = parser.parse_args()
    samples = dict.fromkeys(TABLES)
    if args.sample_hours:
        layout = pd.read_csv(args.sample_hours)
        if layout[["dataset", "sample", "hour"]].isna().any().any() or layout.duplicated(["dataset", "sample"]).any():
            parser.error("Sample hours must uniquely identify each dataset and sample with a finite hour.")
        if not np.isfinite(layout["hour"]).all() or not layout["dataset"].str.fullmatch(r"[a-zA-Z0-9_]+").all():
            parser.error("Sample hours must be numeric and dataset names must be alphanumeric or underscores.")
        samples = {name: dict(zip(rows["sample"], rows["hour"]))
                   for name, rows in layout.groupby("dataset", sort=False)}
        if not samples:
            parser.error("The sample hours CSV is empty.")


    # The reviewer run uses the upstream two-hour phase and peak-to-trough width grid.
    # The optional one-hour grid uses the cached rank scan for the larger template search.
    periods = sorted(set(args.periods))
    if min(periods) < 4:
        parser.error("Periods must be at least 4 hours.")
    step = args.grid_step or (2 if periods == [PERIOD] else 1)
    legacy_grid = periods == [PERIOD] and step == 2
    if step == 2 and not legacy_grid:
        parser.error("Use a 1-hour grid for a search outside the original 24-hour configuration.")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    engine = prepare_upstream(args.upstream.resolve(), output, legacy_grid)


    # Phase is the template peak. Width is the forward time from peak to trough.
    # At 24 h, width 12 is symmetric. Shorter and longer widths allow asymmetric rhythms.
    # The two-hour grid has 12 possible peak positions and 11 widths. Reciprocal
    # peak/trough pairs share a template because the matcher considers both signs (66 templates).
    phases = list(range(0, max(periods), step))
    widths = list(range(step, max(periods), step))
    for name, grid in {"periods": periods, "phases": phases, "asymmetries": widths}.items():
        (engine / f"{name}.txt").write_text("\n".join(map(str, grid)) + "\n")
    settings = {
        "upstream_commit": UPSTREAM_COMMIT,
        "period_hours": periods,
        "phase_hours": phases,
        "width_hours": widths,
        "templates_by_period": {p: 66 if legacy_grid else comb(p, 2) for p in periods},
        "null_draws_per_dataset": 1000,  # Built into this version of upstream eJTK.
        "seed": SEED,
        "dataset_BH_cutoff": CUTOFF,
        "python": sys.version.split()[0],
    }
    if args.sample_hours:
        settings["sample_hours_sha256"] = hashlib.sha256(args.sample_hours.read_bytes()).hexdigest()
    (output / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")


    # Datasets can run independently. Save finished rows in input order as each fit returns.
    results = {}
    rain_tables = args.rain_tables.resolve()
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending = {pool.submit(run_dataset, dataset, rain_tables, output, engine, periods, legacy_grid,
                               sample_times): dataset for dataset, sample_times in samples.items()}
        for future in as_completed(pending):
            dataset = pending[future]
            results[dataset] = future.result()
            rows = [results[name] for name in samples if name in results]
            pd.DataFrame(rows).to_csv(output / "summary.csv", index=False)
            print(dataset, results[dataset]["eJTK_significant"], flush=True)


# Upstream compatibility
# The functions below adapt the pinned Python 2 engine to the project's Python 3
# environment. The fixed-24 h analysis uses its original Kendall template matcher.
def prepare_upstream(upstream, output, legacy_grid):
    """Build a run-local copy from the pinned Git objects so Python 3 edits stay reproducible."""
    commit = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if commit != UPSTREAM_COMMIT:
        raise ValueError(f"Expected upstream {UPSTREAM_COMMIT}, found {commit}")
    engine = output / "upstream"
    (engine / "bin").mkdir(parents=True)
    license_text = subprocess.check_output(["git", "-C", str(upstream), "show", f"{commit}:LICENSE.txt"])
    (engine / "LICENSE.txt").write_bytes(license_text)

    # Port the pinned Python 2 code in the output directory. In particular, the upstream
    # Cython extension expects map() to produce lists that can be indexed.
    for filename in ("eJTK-CalcP.py", "bin/eJTK.py", "bin/CalcP.py", "bin/accessories.pyx"):
        source = subprocess.check_output(["git", "-C", str(upstream), "show", f"{commit}:{filename}"], text=True)
        if filename.endswith(".py"):
            source = re.sub(r"^( +)print (.+)$", r"\1print(\2)", source, flags=re.MULTILINE)
            source = source.replace("xrange(", "range(")
            source = re.sub(r"^ +version=VERSION\s*\n", "", source, flags=re.MULTILINE)
        if filename == "bin/CalcP.py":
            # The optional --fit path depends on legacy mpfit. Use the default Gamma fit.
            for unused in ("import matplotlib.pylab as plt", "from mpfit import mpfit",
                           "import arbfit", "arbFit = arbfit.arbFit"):
                source = source.replace(unused + "\n", "")
        if filename == "bin/accessories.pyx":
            source = source.replace("np.array(map(lamb_get_matches,triples))", "np.array(list(map(lamb_get_matches,triples)))")
            source = source.replace("reference = map(float,reference)", "reference = list(map(float,reference))")
            source = source.replace("kkey = map(float,kkey)", "kkey = list(map(float,kkey))")
            source = source.replace("serie = map(float,list(kkey)    )", "serie = kkey")
            if not legacy_grid:
                # A phase at one period duplicates phase zero and breaks the reciprocal-pair count.
                source = source.replace("phase <=period", "phase < period").replace("phase <= period", "phase < period")
                source = source.replace("generate_base_reference(new_header,waveform,period,phase,width)",
                                        "generate_base_reference(np.remainder(new_header, period),waveform,period,phase,width)")
        (engine / filename).write_text(source)
    for name in ("waveform_cosine.txt", "period24.txt", "phases_00-22_by2.txt", "asymmetries_02-22_by2.txt"):
        source = subprocess.check_output(["git", "-C", str(upstream), "show", f"{commit}:ref_files/{name}"])
        (engine / name).write_bytes(source)
    subprocess.run([str(Path(sys.executable).with_name("cythonize")), "-2", "-i", "accessories.pyx"],
                   cwd=engine / "bin", check=True)
    return engine


class RankScan:
    """Optional one-hour grid matcher. Cache each template's ranks across genes and null draws."""

    def __init__(self, pick_best_match):
        self.pick_best_match = pick_best_match
        self.references = None

    def __call__(self, series, waveform, triples, references, hours):
        if references is not self.references:
            self.references = references
            self.left, self.right = np.triu_indices(len(hours), k=1)
            waves = np.array([references[tuple(triple)] for triple in triples])
            self.signs = np.sign(waves[:, self.left] - waves[:, self.right])
            self.untied_log = np.log(np.count_nonzero(self.signs, axis=1))
            periods, counts = np.unique(triples[:, 0].astype(int), return_counts=True)
            if list(counts) != [comb(period, 2) for period in periods]:
                raise ValueError("Incomplete eJTK phase/asymmetry template grid")

        # For the optional dense template grid, compare pairwise rank signs. The denominator
        # accounts for tied sample or template ranks.
        values = np.asarray(series[1:], dtype=float)
        signs = np.sign(values[self.left] - values[self.right])
        denom = np.exp(0.5 * (self.untied_log + np.log(np.count_nonzero(signs))))
        tau = (self.signs @ signs) / denom
        best_tau = np.max(np.abs(tau))
        selected = np.flatnonzero(np.abs(tau) == best_tau)
        periods, phases, widths = triples[selected].T
        nadirs = (phases + widths) % periods

        # A negative correlation reverses the template, exchanging its peak and trough.
        reverse = tau[selected] < 0
        n = len(hours)
        variance = (4.0 * n + 10.0) / (9.0 * n * (n - 1))
        p = erfc(best_tau / np.sqrt(variance) / 1.4142136) / 2.0
        matches = np.column_stack((np.full(len(selected), best_tau), np.full(len(selected), p), periods,
                                   np.where(reverse, nadirs, phases), np.where(reverse, phases, nadirs),
                                   np.full(len(selected), hours[values.argmax()]),
                                   np.full(len(selected), hours[values.argmin()])))
        # Delegate ties to the upstream rule, including its seeded random choice if ties remain.
        best_match = self.pick_best_match(matches)
        return [series[0], waveform, *best_match[2:], best_match[0], best_match[1]]


def run_upstream():
    """Child-process entry point: seed both random tie-breaking and simulated null profiles."""
    np.random.seed(int(sys.argv.pop(1)))
    mode = sys.argv.pop(1)
    sys.argv.pop(0)
    if mode == "fast":
        sys.path.insert(0, str(Path(sys.argv[0]).parent / "bin"))
        import eJTK
        from accessories import pick_best_match

        eJTK.a_get_best_match = RankScan(pick_best_match)
    runpy.run_path(sys.argv[0], run_name="__main__")


if __name__ == "__main__":
    main()

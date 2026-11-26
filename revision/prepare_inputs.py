"""Restore the four normalized paper tables from the revision's frozen inputs."""
import csv
import gzip
import hashlib
import shutil
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = root / "revision/inputs"
target = root / "results/tables"
target.mkdir(parents=True, exist_ok=True)
with (inputs / "inputs.csv").open() as stream:
    for row in csv.DictReader(stream):
        name = row["file"]
        data = gzip.decompress((inputs / "tables" / (name + ".gz")).read_bytes())
        if hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Input checksum mismatch for {name}")
        (target / name).write_bytes(data)


# The original RAIN plotting code reads manuscript time-series and display RDS files.
# Only the three unchanged time series and the 24 unchanged mouse displays are archived.
mouse = root / "results/shiny_data/mouse"
for source in (inputs / "plot").glob("*/*.rds"):
    dest = mouse / source.relative_to(inputs / "plot")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)
with (inputs / "display/sources.tsv").open() as stream:
    for row in csv.DictReader(stream, delimiter="\t"):
        if row["variant"] != "both":
            continue
        relative = row["source"].split(":data/mouse/")[1]
        source = inputs / "display/mouse" / relative
        if hashlib.sha256(source.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Display checksum mismatch for {relative}")
        dest = mouse / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, dest)

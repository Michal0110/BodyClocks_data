"""Check submitted-panel membership, timing, overlap and permutation results in the revision."""
import csv
from pathlib import Path
import subprocess

def rows(path):
    with path.open() as source:
        return list(csv.DictReader(source))

def page(stem):
    return subprocess.check_output(["pdftotext", "-layout", str(figures / (stem + ".pdf")), "-"], text=True)

figures = Path("results/publication_figures")
results = Path("revision/rain")
membership = rows(results / "rhythmic_membership.csv")
common = [g for g in membership if all(g[d] == "TRUE" for d in ("Osmo","Dex","HS","Cartilage"))]
assert len(common) == 16 and sum(g["clock_panel"] == "TRUE" for g in common) == 9
four_way_plot = page("Fig1C_overlapping_genes")
assert all(g["symbol"] in four_way_plot for g in common)
phase = [g for g in rows(results / "phase_genes.csv") if g["rule"] == "paper_anchor"
         and g["comparison"] == "Osmo_vs_Dex_vs_HS"]
supp = rows(figures / "Supplementary_table_1.csv")
assert len(phase) == len(supp) == 62
aligned = {g["symbol"] for g in phase if g["concordant"] == "TRUE"}
assert aligned == {g["gene"] for g in supp if g["three_way_concordant"] == "TRUE"}
assert len(aligned) == 20 and all(g in page("Fig2E_clock_relative_phase_heatmap") for g in aligned)
pairs = {g["pair"]: g for g in rows(figures / "crp_pairwise_summary.csv")}
for test in rows(results / "phase_tests.csv"):
    if test["rule"] == "paper_anchor" and test["comparison"] in pairs:
        orig = pairs[test["comparison"]]
        assert (test["shared"],test["concordant"],test["p"],test["q"]) == (
            orig["n_shared_genes"], orig["n_clock_controlled"], orig["perm_p_value"], orig["perm_p_BH"])
times = rows(results / "dex_correction/sample_timing.csv")
assert len(times) == 35 and all(t["corrected_hours"] == t["metadata_hours"] for t in times)
deltas = rows(results / "dex_correction/pairwise_mouse_changes.csv")
unchanged = [r for r in deltas if "dex" not in (r["id1"], r["id2"])]
assert len(unchanged) == 276 and all(float(r[f + "_change"]) == 0 for r in unchanged
                                      for f in ("overlap","jaccard","fold_enrichment","p_value"))
wider = rows(results / "mouse_pairs/mouse_pairs_table.csv")
assert len(wider) == 300
for f in figures.glob("*.pdf"):
    assert f.with_suffix(".png").stat().st_size > 1000 and len(page(f.stem).strip()) > 10
print("Verified 16 four-way genes, 20/62 three-way concordant genes, six phase tests and 300 mouse pairs.")


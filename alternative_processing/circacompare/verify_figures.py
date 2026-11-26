from pathlib import Path
import csv
import math
import subprocess
from PIL import Image, ImageStat

root = Path("alternative_processing/circacompare")
out, calls = root / "figures", root / "calls"
def read(p):
    with p.open() as f: return list(csv.DictReader(f))
def text_of(name):
    return subprocess.check_output(["pdftotext", "-layout", str(out/(name+".pdf")), "-"], text=True)
def stems(ext):
    return {x.stem for x in out.glob("*"+ext) if x.stem!="CircaCompare_figures"}
expected = {"Fig1C_overlapping_genes","Fig1D_jaccard_4way_comparison","Fig2A_rose_plots_combined",
            "Fig2B_pairwise_CRP_comparison","Fig2C_crp_pairwise_proportion","Fig2E_clock_relative_phase_heatmap",
            "Fig3A_enrichment_barchart","Supplementary_Figure1",*["Fig3B_enrichment_jaccard_"+k for k in ["KEGG","RCTM","WikiPathways"]]}
expected.update(["detection_and_overlap","Fig2D_three_culture_overlap", *["Fig1B_heatmap_"+d for d in ("osmo","dex","hs")]])
assert stems(".png")==stems(".pdf")==expected
for name in sorted(expected):
    with Image.open(out/(name+".png")) as image:
        image.load()
        assert min(image.size)>=1000
        image.thumbnail((300,300))
        assert max(ImageStat.Stat(image.convert("RGB")).stddev)>10
    page=text_of(name)
    assert len(page.strip())>15 and page.count("\f")==1
four = {r["symbol"] for r in read(calls/"four_way_genes.csv") if r["CircaCompare"]=="TRUE"}
assert len(four)==22 and all(sym in text_of("Fig1C_overlapping_genes") for sym in four)
counts=[int(r["circacompare"]) for r in read(calls/"dataset_summary.csv")]
rose=text_of("Fig2A_rose_plots_combined")
assert all(f"n = {n} rhythmic genes" in rose for n in counts)
a=read(out/"Supplementary_table_1.csv")
b=read(calls/"clock_relative_three_way_genes.csv")
concordance = lambda x, flag:{r["symbol" if "symbol" in r else "gene"] for r in x if r[flag]=="TRUE"}
assert len(a)==len(b)==246
assert concordance(a,"three_way_concordant")==concordance(b,"concordant")
assert all(g in text_of("Fig2E_clock_relative_phase_heatmap") for g in concordance(b,"concordant"))
three=read(out/"crp_3way_summary.csv")[0];baseline=read(calls/"clock_relative_three_way_summary.csv")[0]
assert int(three["n_concordant"])==int(baseline["concordant"])==14
assert math.isclose(float(three["perm_p_value"]), float(baseline["perm_p"]),abs_tol=1e-12),(three,baseline)
pairs=read(out/"crp_pairwise_summary.csv")
baseline=read(calls/"pair_summary.csv")
labels=dict(Osmo="osmo",Dex="dex",HS="hs",Cartilage="cartilage_wt")
for x in pairs:
    z=next(z for z in baseline if z["family"]=="both_rhythmic_set" and (z["first"],z["second"])==(labels[x["dataset1"]],labels[x["dataset2"]]))
    assert int(x["n_shared_genes"])==int(z["eligible"])
    assert f'n = {z["eligible"]} genes' in text_of("Fig2B_pairwise_CRP_comparison")
supp2=read(out/"Supplementary_table_2.csv")
overlap=read(calls/"overlaps.csv")
for x in (x for x in supp2 if x["BH.Q cut-off"]=="0.05"):
    z=next(z for z in overlap if z["method"]=="CircaCompare" and z["groups"]==labels[x["Dataset1"]]+"+"+labels[x["Dataset2"]])
    assert int(x["Overlap"])==int(z["shared_rhythmic"])
    assert math.isclose(float(x["Jaccard_Index"]),round(float(z["jaccard"]),3))
s4=read(out/"Supplementary_table_4.csv")
assert len(s4)==36 and len({(r["Background"],r["Dataset_1"],r["Dataset_2"],r["Annotation_category"]) for r in s4})==36
for k,v in [("KEGG","KEGG Pathways"),("RCTM","Reactome"),("WikiPathways","WikiPathways")]:
    page=text_of("Fig3B_enrichment_jaccard_"+k)
    assert (v in page) or (k == "KEGG" and "EGG Pathways" in page)
    for x in read(out/"pathway_overlap_all_mouse.csv"):
        if x["Category"]==v:
            assert f'{float(x["Jaccard_Index"]):.3f}' in page
assert (out/"CircaCompare.html").stat().st_size>100000
print(f"{len(expected)} PDF/PNG pairs verified. Fig 1C: {len(four)} genes. Fig 2E: {three['n_concordant']}/{three['n_genes']}, permutation p={three['perm_p_value']}. All six pairwise overlap counts and all three pathway matrices checked.")


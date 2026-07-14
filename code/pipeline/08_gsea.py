"""
08 Gsea

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 111  [python]
# ======================================================================
import gseapy, requests
print("gseapy", gseapy.__version__)
# probe Enrichr (gene-set source) reachability
for name,u in {
 "Enrichr geneset API":"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020",
 "Enrichr datasetStats":"https://maayanlab.cloud/Enrichr/datasetStatistics",
}.items():
    try:
        r=requests.get(u,timeout=25,stream=True); print(f"{name}: {r.status_code}")
    except Exception as e: print(f"{name}: ERR {type(e).__name__}")


# ======================================================================
# cell 120  [python]
# ======================================================================
import requests
def fetch_enrichr_library(libname):
    url="https://maayanlab.cloud/Enrichr/geneSetLibrary"
    r=requests.get(url, params={"mode":"text","libraryName":libname}, timeout=90)
    r.raise_for_status()
    d={}
    for line in r.text.strip().split("\n"):
        parts=line.split("\t")
        term=parts[0]
        genes=[g.split(",")[0].strip() for g in parts[2:] if g.strip()]
        d[term]=genes
    return d
hallmark=fetch_enrichr_library("MSigDB_Hallmark_2020")
gobp=fetch_enrichr_library("GO_Biological_Process_2021")
print("Hallmark sets:", len(hallmark), "| GO-BP sets:", len(gobp))
print("Hallmark IFN sets:", [k for k in hallmark if "INTERFERON" in k.upper()])
print("IFN-a set size:", len(hallmark.get("Interferon Alpha Response",[])),
      "| IFN-g:", len(hallmark.get("Interferon Gamma Response",[])))
import json
json.dump({"hallmark":hallmark,"gobp":gobp}, open("annot/genesets.json","w"))


# ======================================================================
# cell 121  [python]
# ======================================================================
import gseapy as gp
gsea_res={}
rank_by_cond={}
for cond in ["Rest","Stim8hr","Stim48hr"]:
    d=setdb1[setdb1.condition==cond].copy()
    # ranking metric: signed Wald stat (direction-aware). collapse dup symbols by max|stat|.
    d=d.dropna(subset=["stat"])
    d["absstat"]=d.stat.abs()
    d=d.sort_values("absstat",ascending=False).drop_duplicates("gene_name")
    rnk=d[["gene_name","stat"]].sort_values("stat",ascending=False).reset_index(drop=True)
    rank_by_cond[cond]=rnk
    pre=gp.prerank(rnk=rnk, gene_sets=hallmark, min_size=5, max_size=500,
                   permutation_num=1000, seed=42, threads=4, no_plot=True,
                   outdir=None)
    res=pre.res2d.copy()
    res["condition"]=cond
    gsea_res[cond]=res
    ifna=res[res.Term=="Interferon Alpha Response"]
    print(f"{cond}: top positive NES ->", res.sort_values('NES',ascending=False).Term.iloc[0],
          f"(NES {res.sort_values('NES',ascending=False).NES.iloc[0]:.2f})",
          f"| IFN-a NES={ifna.NES.iloc[0]:.2f} FDR={ifna['FDR q-val'].iloc[0]:.3g}")


# ======================================================================
# cell 125  [python]
# ======================================================================
# check hallmark gene format
ifna_genes=hallmark["Interferon Alpha Response"]
print("IFN-a hallmark genes sample:", ifna_genes[:8])
print("derep_erv sample:", derep_erv[:8])
overlap=set(derep_erv)&set(ifna_genes)
print("direct overlap derep_erv ∩ IFN-a:", len(overlap), sorted(overlap)[:10])
# genes are uppercase in both -> overlap exists. The gseapy 'enrich' background issue is likely the background arg format.
# Do ORA manually with hypergeometric (cleaner, no gseapy quirks)
from scipy.stats import hypergeom
import numpy as np
bg_set=set(bg); q=set(derep_erv)&bg_set; N=len(bg_set); n=len(q)
ora_rows=[]
for term,genes in hallmark.items():
    gs=set(genes)&bg_set; K=len(gs)
    if K<5: continue
    k=len(q&gs)
    if k==0: continue
    p=hypergeom.sf(k-1,N,K,n)
    # odds ratio
    a,b,c,d=k,n-k,K-k,N-n-(K-k)
    orr=(a*d)/(b*c) if b*c>0 else np.inf
    ora_rows.append(dict(Term=term,k=k,K=K,OR=orr,p=p,genes=";".join(sorted(q&gs))))
ora=pd.DataFrame(ora_rows)
from statsmodels.stats.multitest import multipletests
ora["padj"]=multipletests(ora.p,method="fdr_bh")[1]
ora=ora.sort_values("padj")
print("\nORA (hypergeometric) top Hallmark — ERV-proximal de-repressed, SETDB1 Stim48hr:")
print(ora.head(8)[["Term","k","K","OR","padj"]].to_string(index=False))


# ======================================================================
# cell 127  [python]
# ======================================================================
# assemble the GSEA deliverable table: all Hallmark terms x 3 conditions
gsea_all = pd.concat(gsea_res.values(), ignore_index=True)
gsea_all = gsea_all.rename(columns={"NOM p-val":"pval","FDR q-val":"fdr","Lead_genes":"leading_edge"})
keep=["condition","Term","ES","NES","pval","fdr","Tag %","Gene %","leading_edge"]
gsea_all=gsea_all[keep].copy()
gsea_all["NES"]=gsea_all.NES.astype(float); gsea_all["fdr"]=gsea_all.fdr.astype(float); gsea_all["pval"]=gsea_all.pval.astype(float)
gsea_all.to_csv("setdb1_gsea_hallmark.csv", index=False)
print("saved setdb1_gsea_hallmark.csv:", gsea_all.shape)

# IFN summary across conditions (the honest per-condition significance table)
ifn=gsea_all[gsea_all.Term.isin(["Interferon Alpha Response","Interferon Gamma Response"])]
print("\n=== Interferon programs per condition ===")
print(ifn.pivot_table(index="Term",columns="condition",values=["NES","fdr"]).round(3).to_string())
# rank of IFN-a among positive-NES terms per condition
print("\nIFN-a rank among all 50 Hallmark by NES (1 = top positive):")
for cond in ["Rest","Stim8hr","Stim48hr"]:
    g=gsea_res[cond].sort_values("NES",ascending=False).reset_index(drop=True)
    r=g.index[g.Term=="Interferon Alpha Response"][0]+1
    print(f"  {cond}: #{r}")


# ======================================================================
# cell 129  [python]
# ======================================================================
import matplotlib.pyplot as plt
import matplotlib as mpl
apply_figure_style()

# recompute running ES for IFN-a per condition for panel A (gseapy stores it)
# re-run prerank keeping plots data for IFN-a
es_curves={}
for cond in ["Rest","Stim8hr","Stim48hr"]:
    pre=gp.prerank(rnk=rank_by_cond[cond], gene_sets={"IFNa":hallmark["Interferon Alpha Response"]},
                   min_size=5,max_size=500,permutation_num=100,seed=42,threads=4,no_plot=True,outdir=None)
    t=pre.results["IFNa"]
    es_curves[cond]=dict(RES=np.array(t["RES"]), hits=t["hits"], nes=t["nes"], fdr=t["fdr"])
print("ES curves computed:", {k:round(v["nes"],2) for k,v in es_curves.items()})


# ======================================================================
# cell 160  [python]
# ======================================================================
# running-ES curves (shape only, for the line); label with CANONICAL NES from the table
es_curves={}
canon_nes={c: gsea_all[(gsea_all.condition==c)&(gsea_all.Term=='Interferon Alpha Response')].NES.iloc[0] for c in ["Rest","Stim8hr","Stim48hr"]}
canon_fdr={c: gsea_all[(gsea_all.condition==c)&(gsea_all.Term=='Interferon Alpha Response')].fdr.iloc[0] for c in ["Rest","Stim8hr","Stim48hr"]}
for cond in ["Rest","Stim8hr","Stim48hr"]:
    pre=gp.prerank(rnk=rank_by_cond[cond],gene_sets={"IFNa":hallmark["Interferon Alpha Response"]},
                   min_size=5,max_size=500,permutation_num=1000,seed=42,threads=4,no_plot=True,outdir=None)
    es_curves[cond]=np.array(pre.results["IFNa"]["RES"])
print("curve lengths:", {k:len(v) for k,v in es_curves.items()})
# rebuild heatmap matrix + ORA for panels C/D
conds=["Rest","Stim8hr","Stim48hr"]
genes_show=["IFITM1","IFITM2","IFITM3","OAS1","MX1","RSAD2","STAT2","IRF9","APOL1","APOL2","APOL3","APOL6","CD274"]
mat=np.array([[setdb1[(setdb1.gene_name==g)&(setdb1.condition==cx)].log2FC.iloc[0] for cx in conds] for g in genes_show])
padj_mat=np.array([[setdb1[(setdb1.gene_name==g)&(setdb1.condition==cx)].padj.iloc[0] for cx in conds] for g in genes_show])
ora_full=pd.read_csv("setdb1_ora_hallmark.csv")
ora_full=ora_full[ora_full.geneset=="all de-repressed"].sort_values("padj").head(5).copy()
print("state rebuilt")


# ======================================================================
# cell 229  [python]
# ======================================================================
import pandas as pd, numpy as np, anndata as ad
# Finding 2: genuine NES comparison (Panel A legend vs canonical table)
g=pd.read_csv("setdb1_gsea_hallmark.csv")
ifna=g[g.Term=="Interferon Alpha Response"].set_index("condition")["NES"]
legend={"Rest":1.53,"Stim8hr":2.28,"Stim48hr":2.30}  # values printed in Panel A legend
table={c:round(float(ifna[c]),2) for c in ["Rest","Stim8hr","Stim48hr"]}
nes_match={c:(legend[c]==table[c]) for c in legend}
print("table NES:",table,"| legend:",legend,"| all_match:",all(nes_match.values()))

# inspect isg embedding for figure rebuild
A=ad.read_h5ad("isg_umap_cells.h5ad")
print("obs cols:",list(A.obs.columns))
print("has X_umap:", "X_umap" in A.obsm)
print(A.obs[["group","condition"]].value_counts() if {"group","condition"}.issubset(A.obs.columns) else "no group/condition")

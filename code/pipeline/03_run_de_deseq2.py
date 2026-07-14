"""
03 Run De Deseq2

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 30  [python]
# ======================================================================
import os, numba
# patch njit to avoid disk-cache locator failure in sandbox
_o=numba.njit
numba.njit=lambda *a,**k:(k.__setitem__('cache',False) or _o(*a,**k))
numba.core.decorators.njit=numba.njit
import scanpy as sc, pertpy, pydeseq2, formulaic_contrasts
import sys; sys.path.insert(0,".")
from MultiStatePerturbSeqDataset import MultistatePerturbSeqDataset
from pertpy.tools._differential_gene_expression._pydeseq2 import PyDESeq2
print("OK | scanpy",sc.__version__,"pertpy",pertpy.__version__,"pydeseq2",pydeseq2.__version__,"fc",formulaic_contrasts.__version__)


# ======================================================================
# cell 49  [python]
# ======================================================================
import subprocess
print(subprocess.run(["tail","-6","de_run.log"],capture_output=True,text=True).stdout)
print("partial:", subprocess.run(["bash","-c","wc -l de_module_recomputed.csv.partial 2>/dev/null || echo none"],capture_output=True,text=True).stdout)


# ======================================================================
# cell 56  [bash]
# ======================================================================
cd "$(pwd)" && DE_NCPUS=6 python run_de_module.py > de_run.log 2>&1; echo "EXIT $?"


# ======================================================================
# cell 58  [python]
# ======================================================================
import pandas as pd, numpy as np
de = pd.read_csv("de_module_recomputed.csv")
print("shape:", de.shape)
print("targets:", sorted(de.target.unique()))
print("per condition rows:", de.condition.value_counts().to_dict())
# on-target sanity: each target's logFC on itself should be strongly negative
print("\non-target KD (recomputed log2FC of target on itself):")
for t in ["SETDB1","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2","SUV39H1","SUV39H2","MPHOSPH8"]:
    sub = de[(de.target==t)&(de.gene_name==t)]
    if len(sub): print(f"  {t:9s}:", {r.condition: round(r.log_fc,2) for r in sub.itertuples()})
# tidy columns for the deliverable
de_out = de[["target","condition","gene_name","gene_id","log_fc","lfcSE","stat","p_value","adj_p_value","baseMean"]].copy()
de_out = de_out.rename(columns={"log_fc":"log2FC","p_value":"pvalue","adj_p_value":"padj"})
de_out.to_csv("de_module_recomputed.csv", index=False)
print("\nfinal columns:", de_out.columns.tolist())


# ======================================================================
# cell 91  [python]
# ======================================================================
import gzip, re
# --- GENCODE v44: gene records -> canonical TSS (strand-aware) ---
# parse 'gene' features; TSS = start if + strand else end
rows=[]
with gzip.open("annot/gencode.v44.basic.annotation.gtf.gz","rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if f[2]!="gene": continue
        chrom,start,end,strand,attr=f[0],int(f[3]),int(f[4]),f[6],f[8]
        gid=re.search(r'gene_id "([^"]+)"',attr)
        gname=re.search(r'gene_name "([^"]+)"',attr)
        gtype=re.search(r'gene_type "([^"]+)"',attr)
        rows.append((chrom,start,end,strand,
                     gid.group(1) if gid else None,
                     gname.group(1) if gname else None,
                     gtype.group(1) if gtype else None))
genes=pd.DataFrame(rows,columns=["chrom","start","end","strand","gene_id_full","gene_name","gene_type"])
genes["gene_id"]=genes.gene_id_full.str.split(".").str[0]   # strip version
genes["tss"]=np.where(genes.strand=="+", genes.start, genes.end)
genes["gene_length"]=genes.end-genes.start
print("GENCODE gene records:", len(genes), "| unique gene_id:", genes.gene_id.nunique())

# background = authors' 10,282 tested genes (Ensembl IDs) from recomputed DE
de = pd.read_csv("de_module_recomputed.csv")
bg_ids = pd.Index(de.gene_id.unique())
print("background gene_ids (authors' tested universe):", len(bg_ids))
# map onto GENCODE by Ensembl ID
gmap = genes[genes.gene_id.isin(bg_ids)].drop_duplicates("gene_id").set_index("gene_id")
matched = bg_ids.isin(gmap.index)
print(f"matched to GENCODE TSS: {matched.sum()} / {len(bg_ids)} ({100*matched.mean():.1f}%)")
print("unmatched examples:", list(bg_ids[~matched][:8]))


# ======================================================================
# cell 94  [python]
# ======================================================================
from scipy.stats import fisher_exact
ann = pd.read_parquet("erv_proximity_annotation.parquet").set_index("gene_id")
de = pd.read_csv("de_module_recomputed.csv")

PADJ, LFC = 0.10, 0.5
# de-repressed = up-regulated past threshold in >=1 condition, on-target excluded
de["derep"] = (de.padj<PADJ) & (de.log2FC>LFC) & (de.target!=de.gene_name)
derep_by_target = {t: set(g[g.derep].gene_id.unique()) for t,g in de.groupby("target")}
print("de-repressed gene counts per target (pooled across conditions):")
for t in ["SETDB1","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2","SUV39H1","SUV39H2"]:
    print(f"  {t:8s}: {len(derep_by_target.get(t,set()))}")

def fisher_or(derep_ids, near_flag_col, universe_idx, exclude_ids=None):
    idx = universe_idx
    if exclude_ids is not None:
        idx = idx.difference(exclude_ids)
    near = ann.loc[idx, near_flag_col].values
    dr = idx.isin(derep_ids)
    a = int(( dr &  near).sum()); b = int(( dr & ~near).sum())
    c = int((~dr &  near).sum()); d = int((~dr & ~near).sum())
    orr, p = fisher_exact([[a,b],[c,d]])
    # Haldane-Anscombe CI
    from math import log, sqrt, exp
    a2,b2,c2,d2 = a+0.5,b+0.5,c+0.5,d+0.5
    se = sqrt(1/a2+1/b2+1/c2+1/d2); lo=exp(log(a2*d2/(b2*c2))-1.96*se); hi=exp(log(a2*d2/(b2*c2))+1.96*se)
    return dict(a=a,b=b,c=c,d=d,OR=orr,CI_low=lo,CI_high=hi,p=p,n_derep=int(dr.sum()),n_near_derep=a)

uni = ann.index
targets_mod = ["SETDB1","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2"]
raw_rows=[]
for t in targets_mod:
    for lab in ["50kb","100kb"]:
        r=fisher_or(derep_by_target[t], f"near_silenced_erv_{lab}", uni)
        raw_rows.append(dict(target=t,window=lab,analysis="raw",**r))
raw=pd.DataFrame(raw_rows)
print("\n=== RAW enrichment (uncorrected) — silenced-ERV proximity of de-repressed genes ===")
print(raw[["target","window","n_derep","OR","CI_low","CI_high","p"]].round(3).to_string(index=False))


# ======================================================================
# cell 116  [python]
# ======================================================================
import pandas as pd, numpy as np
de = pd.read_csv("de_module_recomputed.csv")
ann = pd.read_parquet("erv_proximity_annotation.parquet")[["gene_id","gene_name","near_silenced_erv_100kb","near_silenced_erv_50kb","is_KRABZNF"]]

# focus target = SETDB1 (writer, module core, strongest signal; deliverables are setdb1_*)
setdb1 = de[de.target=="SETDB1"].copy()
print("SETDB1 rows:", len(setdb1), "| conditions:", sorted(setdb1.condition.unique()))
print("genes per condition:", setdb1.groupby("condition").size().to_dict())
# merge ERV-proximity + KRAB flag
setdb1 = setdb1.merge(ann, on=["gene_id","gene_name"], how="left")
# use HGNC symbol for MSigDB/Enrichr; drop rows w/o symbol
setdb1 = setdb1[setdb1.gene_name.notna() & (setdb1.gene_name!="")].copy()
print("with symbol:", len(setdb1))
# sanity: key ISGs present in the tested universe?
isg_core=["IFITM1","IFITM2","IFITM3","OAS1","MX1","RSAD2","STAT2","IRF9","CD274","APOL1","APOL2","APOL3","APOL6"]
present={g: g in set(setdb1.gene_name) for g in isg_core}
print("ISG-core present in SETDB1-tested universe:", present)


# ======================================================================
# cell 134  [python]
# ======================================================================
import matplotlib.pyplot as plt, matplotlib as mpl
# regenerate fig_erv_enrichment.png with corrected suptitle (data shows raw up to ~14x)
# reuse the exact panel code; only the suptitle text changes.
exec(open("run_de_module.py").read()) if False else None
# rather than re-run everything, reload state needed
import pandas as pd, numpy as np
from math import exp
R = pd.read_csv("erv_enrichment_results.csv")
ann = pd.read_parquet("erv_proximity_annotation.parquet").set_index("gene_id")
# rebuild pooled set for panel C
de = pd.read_csv("de_module_recomputed.csv")
core_targets=["SETDB1","TASOR","PPHLN1","ATF7IP"]
de["derep"]=(de.padj<0.10)&(de.log2FC>0.5)&(de.target!=de.gene_name)
pooled=set().union(*[set(de[(de.target==t)&de.derep].gene_id) for t in core_targets])
def get(target,window,analysis):
    r=R[(R.target==target)&(R.window==window)&(R.analysis==analysis)]
    return (r.OR.iloc[0],r.CI_low.iloc[0],r.CI_high.iloc[0]) if len(r) else (np.nan,)*3
apply_figure_style()
craw,cadj="#c6702a","#08519c"
tord=["SETDB1","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2"]
groups=[("MODULE_CORE_POOLED",""),("SUV39_POOLED","")]

fig=plt.figure(figsize=(12,5.2)); gs=fig.add_gridspec(1,3,width_ratios=[1.15,1.0,1.0],wspace=0.42)
axA=fig.add_subplot(gs[0,0]); y=np.arange(len(tord))[::-1]
for i,t in enumerate(tord):
    orr,lo,hi=get(t,"50kb","raw"); orj,loj,hij=get(t,"50kb","logit_adj")
    axA.plot([lo,hi],[y[i]+0.16]*2,c=craw,lw=1.4); axA.scatter([orr],[y[i]+0.16],c=craw,s=26,zorder=3)
    axA.plot([loj,hij],[y[i]-0.16]*2,c=cadj,lw=1.4); axA.scatter([orj],[y[i]-0.16],c=cadj,s=26,zorder=3)
axA.axvline(1,ls=":",c="#555",lw=1); axA.set_yticks(y); axA.set_yticklabels(tord)
axA.set_xscale("log"); axA.set_xlim(0.8,32); axA.set_xticks([1,2,5,10,20]); axA.set_xticklabels(["1","2","5","10","20"])
axA.set_xlabel("Odds ratio (silenced-ERV proximity, ±50 kb)")
axA.set_title("A  Per-target enrichment: raw vs KRAB-ZNF-corrected",fontsize=9,loc="left")
axA.scatter([],[],c=craw,s=26,label="Raw (unco…[+2550 chars]


# ======================================================================
# cell 143  [python]
# ======================================================================
# Load ERV-proximity annotation + SETDB1 de-repressed set for a focused antagonism test
ann=pd.read_parquet("erv_proximity_annotation.parquet")
erv_prox_ids=set(ann[ann.near_silenced_erv_100kb==True].gene_id)
de=pd.read_csv("de_module_recomputed.csv")
setdb1_up = de[(de.target=="SETDB1")&(de.condition=="Stim48hr")&(de.padj<0.10)&(de.log2FC>0.5)]
setdb1_up_ids=set(setdb1_up.gene_id)
setdb1_up_erv = setdb1_up_ids & erv_prox_ids
print(f"SETDB1 de-repressed (Stim48hr): {len(setdb1_up_ids)}, of which ERV-proximal: {len(setdb1_up_erv)}")

# columns for those gene sets
def cols_for(ids): return [np.where(gene_ids_de==g)[0][0] for g in ids if (gene_ids_de==g).any()]
c_up=cols_for(setdb1_up_ids); c_up_erv=cols_for(setdb1_up_erv)

print("\n=== Mean log2FC on SETDB1-de-repressed genes (Stim48hr) — antagonism = negative ===")
print(f"{'node':8s} {'all_SETDB1up':>13s} {'ERVprox_SETDB1up':>17s}")
for n in ["SETDB1","KDM4A","KDM4B","KDM4C","CBX1","CBX3","CBX5"]:
    if (n,"Stim48hr") in node_rows:
        v=prof(n,"Stim48hr")
        print(f"{n:8s} {np.nanmean(v[c_up]):+13.2f} {np.nanmean(v[c_up_erv]):+17.2f}")


# ======================================================================
# cell 154  [python]
# ======================================================================
import json, pandas as pd, numpy as np
tg=json.load(open("handoff/chembl_targets.json"))
ba=json.load(open("handoff/chembl_bioact.json"))
cl=json.load(open("handoff/chembl_clinical.json"))
antag=pd.read_csv("annot/node_antagonism.csv")

# node metadata: axis role + whether directly measured as KD target + on-target KD efficacy (from our DE)
de=pd.read_csv("de_module_recomputed.csv")
# on-target KD for writer nodes we recomputed
axis_map={"SETDB1":"writer","TRIM28":"writer","ATF7IP":"writer","TASOR":"writer","PPHLN1":"writer",
          "KDM4A":"eraser","KDM4B":"eraser","KDM4C":"eraser","CBX1":"reader","CBX3":"reader","CBX5":"reader"}
role_detail={"SETDB1":"H3K9 methyltransferase (writer)","TRIM28":"KAP1 scaffold/corepressor (writer complex)",
 "ATF7IP":"SETDB1 stabiliser (writer complex)","TASOR":"HUSH core (writer complex)","PPHLN1":"HUSH subunit (writer complex)",
 "KDM4A":"H3K9me3 demethylase (eraser)","KDM4B":"H3K9me3 demethylase (eraser)","KDM4C":"H3K9me3 demethylase (eraser)",
 "CBX1":"HP1β reader","CBX3":"HP1γ reader","CBX5":"HP1α reader"}

rows=[]
for g,axis in axis_map.items():
    b=ba.get(g,{}); c=cl.get(g,{})
    # antagonism (Stim48hr) if eraser/reader; on-target KD if writer
    a48=antag[(antag.node==g)&(antag.condition=="Stim48hr")]
    rows.append(dict(
        node=g, axis=axis, role=role_detail[g],
        chembl_target_id=tg[g].get("target_chembl_id"),
        has_chembl_target=tg[g].get("target_chembl_id") is not None,
        n_compounds_uM_or_better=b.get("n_compounds_uM_or_better",0),
        n_compounds_sub100nM=b.get("n_compounds_sub100nM",0),
        best_pchembl=b.get("best_pchembl"),
        max_clinical_phase=c.get("max_phase_among_potent"),
        n_curated_mechanisms=c.get("n_mechanisms_curated"),
        directly_measured_KD_target=True,  # all screened (confirmed earlier)
        antag_median_log2FC_stim48=(a48.median_log2FC.iloc[0] if len(a48) else np.nan),
        antag_wilcoxon_p_stim48=(a48.wilcoxon_p.iloc[0] if len(a48) else n…[+266 chars]


# ======================================================================
# cell 168  [python]
# ======================================================================
import matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.patches import Patch
apply_figure_style()
axis_col={"writer":"#08519c","eraser":"#d94801","reader":"#6a51a3"}
c=pd.read_csv("part3_ranked_candidates.csv")
# need antag scatter data + prof(); reload DE-stats node profiles
import h5py, fsspec
url="s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
fs=fsspec.filesystem("s3",anon=True)
f=h5py.File(fs.open(url,"rb",block_size=4*1024*1024,cache_type="readahead"),"r")
tc=f["obs"]["target_contrast_gene_name"]["categories"][:].astype(str)[f["obs"]["target_contrast_gene_name"]["codes"][:]]
cc_=f["obs"]["culture_condition"]["categories"][:].astype(str)[f["obs"]["culture_condition"]["codes"][:]]
gid=f["var"]["_index"][:].astype(str)
r_setdb1=np.where((tc=="SETDB1")&(cc_=="Stim48hr"))[0][0]
r_kdm4c=np.where((tc=="KDM4C")&(cc_=="Stim48hr"))[0][0]
rows=sorted([r_setdb1,r_kdm4c])
lfc=f["layers"]["log_fc"][rows,:]; loc={r:i for i,r in enumerate(rows)}
de=pd.read_csv("de_module_recomputed.csv")
ann2=pd.read_parquet("erv_proximity_annotation.parquet")
erv_ids=set(ann2[ann2.near_silenced_erv_100kb==True].gene_id)
su=de[(de.target=="SETDB1")&(de.condition=="Stim48hr")&(de.padj<0.10)&(de.log2FC>0.5)]
up_erv=set(su.gene_id)&erv_ids
cols=[np.where(gid==g)[0][0] for g in up_erv if (gid==g).any()]
sv=lfc[loc[r_setdb1],cols]; kv=lfc[loc[r_kdm4c],cols]
m=np.isfinite(sv)&np.isfinite(kv)
print("scatter points:",m.sum())


# ======================================================================
# cell 178  [python]
# ======================================================================
import gseapy as gp, os, requests
os.environ["HOME"]=os.getcwd()
gsea_all=pd.read_csv("setdb1_gsea_hallmark.csv")
# rebuild ranked sigs + IFN-a running ES (1000 perms) for the curves
de=pd.read_csv("de_module_recomputed.csv")
s1=de[de.target=="SETDB1"].dropna(subset=["stat"]).copy()
def lib(n):
    r=requests.get("https://maayanlab.cloud/Enrichr/geneSetLibrary",params={"mode":"text","libraryName":n},timeout=90);r.raise_for_status()
    d={};
    for line in r.text.strip().split("\n"):
        p=line.split("\t"); d[p[0]]=[g.split(",")[0].strip() for g in p[2:] if g.strip()]
    return d
hallmark=lib("MSigDB_Hallmark_2020")
conds=["Rest","Stim8hr","Stim48hr"]
es={}
for c in conds:
    d=s1[s1.condition==c].copy(); d["a"]=d.stat.abs()
    rnk=d.sort_values("a",ascending=False).drop_duplicates("gene_name")[["gene_name","stat"]].sort_values("stat",ascending=False)
    pre=gp.prerank(rnk=rnk,gene_sets={"IFNa":hallmark["Interferon Alpha Response"]},min_size=5,max_size=500,permutation_num=1000,seed=42,threads=4,no_plot=True,outdir=None)
    es[c]=np.array(pre.results["IFNa"]["RES"])
canon_nes={c:gsea_all[(gsea_all.condition==c)&(gsea_all.Term=="Interferon Alpha Response")].NES.iloc[0] for c in conds}
print("IFN-a NES:",{c:round(canon_nes[c],2) for c in conds},"| ES curve lens:",{c:len(es[c]) for c in conds})


# ======================================================================
# cell 181  [python]
# ======================================================================
de=pd.read_csv("de_module_recomputed.csv")
v=de[(de.target=="SETDB1")&(de.condition=="Stim48hr")].copy()
v["neglog10padj"]=-np.log10(v.padj.clip(lower=1e-300))
isg=["IFITM1","IFITM2","IFITM3","OAS1","MX1","RSAD2"]
apol=["APOL1","APOL2","APOL3","APOL6"]
padj_thr=0.10; lfc_thr=0.5
v["sig"]=(v.padj<padj_thr)&(v.log2FC.abs()>lfc_thr)
print("SETDB1 Stim48hr:", len(v), "genes |", "up:",int(((v.log2FC>lfc_thr)&(v.padj<padj_thr)).sum()),
      "down:",int(((v.log2FC<-lfc_thr)&(v.padj<padj_thr)).sum()))
print("ISG present:",{g:g in set(v.gene_name) for g in isg})
print("APOL present:",{g:g in set(v.gene_name) for g in apol})
print("SETDB1 on-target:", v[v.gene_name=="SETDB1"][["log2FC","padj"]].values)
# cap neglog10 for display
print("max neglog10padj (non-ontarget):", v[v.gene_name!="SETDB1"].neglog10padj.max().round(1))


# ======================================================================
# cell 222  [python]
# ======================================================================
# rebuild heatmap + ORA state for the 4-panel deliverable figure
de=pd.read_csv("de_module_recomputed.csv")
s1=de[de.target=="SETDB1"].copy()
genes_show=["IFITM1","IFITM2","IFITM3","OAS1","MX1","RSAD2","STAT2","IRF9","APOL1","APOL2","APOL3","APOL6","CD274"]
mat=np.array([[s1[(s1.gene_name==gg)&(s1.condition==cx)].log2FC.iloc[0] for cx in conds] for gg in genes_show])
padj_mat=np.array([[s1[(s1.gene_name==gg)&(s1.condition==cx)].padj.iloc[0] for cx in conds] for gg in genes_show])
ora_full=pd.read_csv("setdb1_ora_hallmark.csv"); ora_full=ora_full[ora_full.geneset=="all de-repressed"].sort_values("padj").head(5).copy()
print("state ok; heatmap",mat.shape,"ora",ora_full.shape)

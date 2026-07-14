"""
09 Ora

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 126  [python]
# ======================================================================
def ora_hyper(query, bg, libraries):
    bg_set=set(bg); q=set(query)&bg_set; N=len(bg_set); n=len(q); rows=[]
    for libname,lib in libraries.items():
        for term,genes in lib.items():
            gs=set(genes)&bg_set; K=len(gs)
            if K<5: continue
            k=len(q&gs)
            if k<2: continue
            p=hypergeom.sf(k-1,N,K,n)
            a,b,c,d=k,n-k,K-k,N-n-(K-k)
            orr=(a*d)/(b*c) if b*c>0 else np.inf
            rows.append(dict(library=libname,Term=term,k=k,K=K,OR=orr,p=p,genes=";".join(sorted(q&gs))))
    df=pd.DataFrame(rows)
    df["padj"]=multipletests(df.p,method="fdr_bh")[1]
    return df.sort_values("padj")

# full de-repressed set (all up, padj<0.10, |lfc|>0.5), Stim48hr
derep_full = setdb1[(setdb1.condition=="Stim48hr")&(setdb1.padj<0.10)&(setdb1.log2FC>0.5)].gene_name.dropna().unique().tolist()
print(f"Full de-repressed set (SETDB1 Stim48hr): {len(derep_full)} genes")
ora_full=ora_hyper(derep_full, bg, {"Hallmark":hallmark})
print("\nORA top Hallmark (FULL de-repressed set, Stim48hr):")
print(ora_full.head(6)[["Term","k","K","OR","padj"]].to_string(index=False))
# GO-BP ORA
ora_go=ora_hyper(derep_full, bg, {"GO_BP":gobp})
print("\nORA top GO-BP (FULL de-repressed, Stim48hr):")
print(ora_go.head(6)[["Term","k","K","OR","padj"]].to_string(index=False)[:900])


# ======================================================================
# cell 128  [python]
# ======================================================================
# save ORA tables too (both the ERV-proximal subset and full set, both libraries)
ora_erv = ora_hyper(derep_erv, bg, {"Hallmark":hallmark})
ora_erv["geneset"]="ERV-proximal de-repressed (100kb)"
ora_full2=ora_full.copy(); ora_full2["geneset"]="all de-repressed"
ora_out=pd.concat([ora_full2.assign(library="Hallmark"), ora_erv.assign(library="Hallmark")], ignore_index=True)
ora_out.to_csv("setdb1_ora_hallmark.csv", index=False)
print("saved setdb1_ora_hallmark.csv")
# also save the ranked signatures for reproducibility
for cond,rnk in rank_by_cond.items():
    rnk.to_csv(f"annot/rank_setdb1_{cond}.csv", index=False)
print("saved ranked signatures")


# ======================================================================
# cell 223  [python]
# ======================================================================
ora=pd.read_csv("setdb1_ora_hallmark.csv")
print("genesets:",ora.geneset.unique())
# use the ERV-proximal de-repressed set for ORA panel (matches the ORA-vs-GSEA contrast)
sets=ora.geneset.unique()
ora_show=ora[ora.geneset==sets[0]].sort_values("p").head(6)
print(ora_show[["Term","OR","padj"]].to_string())

"""
13 Fig Iii Setdb1 Volcano

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 184  [python]
# ======================================================================
from adjustText import adjust_text
CAP=90  # display cap for on-target
vd=v.copy(); vd["y"]=vd.neglog10padj.clip(upper=CAP)
isg_set,apol_set=set(isg),set(apol)
def build_volcano(theme):
    grid=("#888" if theme=="light" else "#aaa"); base=("#c9c9c9" if theme=="light" else "#555")
    fig,ax=plt.subplots(figsize=(7.4,6.2))
    ns=vd[~vd.sig]; ax.scatter(ns.log2FC,ns.y,s=6,c=base,alpha=0.45,edgecolors="none",zorder=1,rasterized=True)
    up=vd[vd.sig&(vd.log2FC>0)]; dn=vd[vd.sig&(vd.log2FC<0)]
    ax.scatter(up.log2FC,up.y,s=9,c="#d94801",alpha=0.7,edgecolors="none",zorder=2,rasterized=True,label=f"Up (de-repressed), n={len(up)}")
    ax.scatter(dn.log2FC,dn.y,s=9,c="#3182bd",alpha=0.7,edgecolors="none",zorder=2,rasterized=True,label=f"Down, n={len(dn)}")
    ax.axvline(lfc_thr,ls=":",c=grid,lw=0.8); ax.axvline(-lfc_thr,ls=":",c=grid,lw=0.8)
    ax.axhline(-np.log10(padj_thr),ls=":",c=grid,lw=0.8)
    texts=[]
    # ISG core (orange), APOL (purple), on-target SETDB1 (dark)
    for _,row in vd.iterrows():
        g=row.gene_name
        if g in isg_set:
            ax.scatter(row.log2FC,row.y,s=42,facecolors="none",edgecolors="#7f2704",lw=1.4,zorder=4)
            texts.append(ax.text(row.log2FC,row.y,g,fontsize=7.5,color="#7f2704",fontweight="bold",zorder=5))
        elif g in apol_set:
            ax.scatter(row.log2FC,row.y,s=42,facecolors="none",edgecolors="#6a51a3",lw=1.4,zorder=4)
            texts.append(ax.text(row.log2FC,row.y,g,fontsize=7.5,color="#6a51a3",fontweight="bold",zorder=5))
        elif g=="SETDB1":
            ax.scatter(row.log2FC,row.y,s=55,marker="D",facecolors="none",edgecolors=("#111" if theme=="light" else "#eee"),lw=1.5,zorder=4)
            texts.append(ax.text(row.log2FC,row.y,"SETDB1\n(on-target)",fontsize=7,color=("#111" if theme=="light" else "#eee"),zorder=5))
    adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle="-",color=grid,lw=0.6),
                expand=(1.4,1.8),force_text=(0.4,0.6))
    ax.set_xlabel("log₂ fold-change (SETDB…[+864 chars]


# ======================================================================
# cell 185  [python]
# ======================================================================
def build_volcano(theme):
    grid=("#888" if theme=="light" else "#aaa"); base=("#c9c9c9" if theme=="light" else "#555")
    fgc=("#111" if theme=="light" else "#eee")
    fig,ax=plt.subplots(figsize=(7.4,6.2))
    ns=vd[~vd.sig]; ax.scatter(ns.log2FC,ns.y,s=6,c=base,alpha=0.45,edgecolors="none",zorder=1,rasterized=True)
    up=vd[vd.sig&(vd.log2FC>0)]; dn=vd[vd.sig&(vd.log2FC<0)]
    ax.scatter(up.log2FC,up.y,s=9,c="#d94801",alpha=0.7,edgecolors="none",zorder=2,rasterized=True,label=f"Up (de-repressed), n={len(up)}")
    ax.scatter(dn.log2FC,dn.y,s=9,c="#3182bd",alpha=0.7,edgecolors="none",zorder=2,rasterized=True,label=f"Down, n={len(dn)}")
    ax.axvline(lfc_thr,ls=":",c=grid,lw=0.8); ax.axvline(-lfc_thr,ls=":",c=grid,lw=0.8)
    ax.axhline(-np.log10(padj_thr),ls=":",c=grid,lw=0.8)
    # class handles
    ax.scatter([],[],s=42,facecolors="none",edgecolors="#7f2704",lw=1.4,label="ISG core")
    ax.scatter([],[],s=42,facecolors="none",edgecolors="#6a51a3",lw=1.4,label="APOL family")
    ax.scatter([],[],s=55,marker="D",facecolors="none",edgecolors=fgc,lw=1.5,label="SETDB1 (on-target)")
    texts=[]
    for _,row in vd.iterrows():
        g=row.gene_name
        if g in isg_set:
            ax.scatter(row.log2FC,row.y,s=42,facecolors="none",edgecolors="#7f2704",lw=1.4,zorder=4)
            texts.append(ax.text(row.log2FC,row.y,g,fontsize=7.5,color="#7f2704",fontweight="bold",zorder=5))
        elif g in apol_set:
            ax.scatter(row.log2FC,row.y,s=42,facecolors="none",edgecolors="#6a51a3",lw=1.4,zorder=4)
            texts.append(ax.text(row.log2FC,row.y,g,fontsize=7.5,color="#6a51a3",fontweight="bold",zorder=5))
    adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle="-",color=grid,lw=0.6),expand=(1.4,1.8),force_text=(0.4,0.6))
    # on-target: manual diamond + leader into empty mid-left
    ot=vd[vd.gene_name=="SETDB1"].iloc[0]
    ax.scatter(ot.log2FC,ot.y,s=60,marker="D",facecolors="none",edgecolors=fgc,lw=1.6,zorder=6)
    ax.annotate("SETDB1 (on-target)…[+818 chars]

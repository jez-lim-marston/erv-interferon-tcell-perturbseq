"""
16 Fig Setdb1 Gsea 4Panel

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 224  [python]
# ======================================================================
from matplotlib.colors import TwoSlopeNorm
cond_colors={"Rest":"#9ecae1","Stim8hr":"#4292c6","Stim48hr":"#084594"}
def build_setdb1_gsea(theme):
    fgc=("#111" if theme=="light" else "#eee"); grid=("#888" if theme=="light" else "#aaa")
    fig=plt.figure(figsize=(13,9))
    gs=fig.add_gridspec(2,2,hspace=0.34,wspace=0.28,height_ratios=[1,1.05])
    # A running ES
    axA=fig.add_subplot(gs[0,0])
    for c in conds: axA.plot(np.arange(len(es[c])),es[c],c=cond_colors[c],lw=1.8,label=f"{c} (NES {canon_nes[c]:.2f})")
    axA.axhline(0,c=grid,lw=0.8); axA.set_xlabel("Gene rank (de-repression → repression)"); axA.set_ylabel("Running ES")
    axA.set_title("A  Hallmark IFN-α enrichment (pre-ranked GSEA)",fontsize=9.5,loc="left")
    axA.legend(frameon=False,fontsize=7.5,loc="upper right",title="Activation",title_fontsize=7.5)
    # B NES bars
    axB=fig.add_subplot(gs[0,1]); x=np.arange(len(conds)); w=0.36
    for i,(term,pl) in enumerate(zip(prog,plabs)):
        nes=[g[(g.condition==c)&(g.Term==term)].NES.iloc[0] for c in conds]; fdrv=[g[(g.condition==c)&(g.Term==term)].fdr.iloc[0] for c in conds]
        col="#d94801" if i==0 else "#6a51a3"; axB.bar(x+(i-0.5)*w,nes,w,color=col,label=pl)
        for xi,(n,fd) in enumerate(zip(nes,fdrv)):
            st="***" if fd<0.001 else "**" if fd<0.01 else "*" if fd<0.05 else "ns"
            axB.text(x[xi]+(i-0.5)*w,n+0.05,st,ha="center",va="bottom",fontsize=8,color=(fgc if fd<0.05 else grid))
    axB.set_xticks(x); axB.set_xticklabels(conds); axB.set_ylabel("NES"); axB.set_ylim(0,3.0)
    axB.set_title("B  IFN-α top-ranked all 3; FDR-sig in stimulated; IFN-γ stim-dependent",fontsize=8.2,loc="left")
    axB.legend(frameon=False,fontsize=8,loc="upper left")
    axB.text(0.52,0.9,"*** FDR<0.001  ** <0.01  * <0.05  ns",transform=axB.transAxes,ha="left",va="top",fontsize=6.3,color=grid)
    # C heatmap
    axC=fig.add_subplot(gs[1,0]); norm=TwoSlopeNorm(vmin=-0.5,vcenter=0,vmax=1.5)
    im=axC.imshow(mat,aspect="auto",cmap="RdBu_r",…[+1495 chars]


# ======================================================================
# cell 225  [python]
# ======================================================================
short={"Interferon Alpha Response":"IFN-α Response","Epithelial Mesenchymal Transition":"Epithelial-Mes. Transition",
       "Wnt-beta Catenin Signaling":"Wnt/β-catenin","Inflammatory Response":"Inflammatory Resp.","KRAS Signaling Dn":"KRAS Signaling Dn"}
def build_setdb1_gsea(theme):
    fgc=("#111" if theme=="light" else "#eee"); grid=("#888" if theme=="light" else "#aaa")
    fig=plt.figure(figsize=(13.5,9))
    gs=fig.add_gridspec(2,2,hspace=0.34,wspace=0.42,height_ratios=[1,1.05])
    axA=fig.add_subplot(gs[0,0])
    for c in conds: axA.plot(np.arange(len(es[c])),es[c],c=cond_colors[c],lw=1.8,label=f"{c} (NES {canon_nes[c]:.2f})")
    axA.axhline(0,c=grid,lw=0.8); axA.set_xlabel("Gene rank (de-repression → repression)"); axA.set_ylabel("Running ES")
    axA.set_title("A  Hallmark IFN-α enrichment (pre-ranked GSEA)",fontsize=9.5,loc="left")
    axA.legend(frameon=False,fontsize=7.5,loc="upper right",title="Activation",title_fontsize=7.5)
    axB=fig.add_subplot(gs[0,1]); x=np.arange(len(conds)); w=0.36
    for i,(term,pl) in enumerate(zip(prog,plabs)):
        nes=[g[(g.condition==c)&(g.Term==term)].NES.iloc[0] for c in conds]; fdrv=[g[(g.condition==c)&(g.Term==term)].fdr.iloc[0] for c in conds]
        col="#d94801" if i==0 else "#6a51a3"; axB.bar(x+(i-0.5)*w,nes,w,color=col,label=pl)
        for xi,(n,fd) in enumerate(zip(nes,fdrv)):
            st="***" if fd<0.001 else "**" if fd<0.01 else "*" if fd<0.05 else "ns"
            axB.text(x[xi]+(i-0.5)*w,n+0.05,st,ha="center",va="bottom",fontsize=8,color=(fgc if fd<0.05 else grid))
    axB.set_xticks(x); axB.set_xticklabels(conds); axB.set_ylabel("NES"); axB.set_ylim(0,3.0)
    axB.set_title("B  IFN-α top-ranked all 3; FDR-sig in stimulated; IFN-γ stim-dependent",fontsize=8.2,loc="left")
    axB.legend(frameon=False,fontsize=8,loc="upper left")
    axB.text(0.52,0.9,"*** FDR<0.001  ** <0.01  * <0.05  ns",transform=axB.transAxes,ha="left",va="top",fontsize=6.3,color=grid)
    axC=fig.add_subplot(gs[1,0]); norm=Tw…[+1977 chars]

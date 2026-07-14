"""
12 Fig Ii Interferon Gsea

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 180  [python]
# ======================================================================
def build_gsea(theme):
    fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={"width_ratios":[1.25,1],"wspace":0.32})
    grid=("#888" if theme=="light" else "#aaa")
    for c in conds:
        axA.plot(np.arange(len(es[c])),es[c],c=cond_colors[c],lw=1.9,label=f"{c} (NES {canon_nes[c]:.2f})")
    axA.axhline(0,c=grid,lw=0.8)
    axA.set_xlabel("Gene rank (de-repression → repression)"); axA.set_ylabel("Running enrichment score")
    axA.set_title("A  Hallmark IFN-α enrichment (pre-ranked GSEA, SETDB1 KD)",fontsize=9.5,loc="left")
    axA.legend(frameon=False,fontsize=8,loc="upper right",title="Activation",title_fontsize=8)
    x=np.arange(len(conds)); w=0.36
    for i,(term,pl) in enumerate(zip(prog,plabs)):
        nes=[gsea_all[(gsea_all.condition==c)&(gsea_all.Term==term)].NES.iloc[0] for c in conds]
        fdr=[gsea_all[(gsea_all.condition==c)&(gsea_all.Term==term)].fdr.iloc[0] for c in conds]
        col="#d94801" if i==0 else "#6a51a3"
        axB.bar(x+(i-0.5)*w,nes,w,color=col,label=pl)
        for xi,(n,fd) in enumerate(zip(nes,fdr)):
            star="***" if fd<0.001 else "**" if fd<0.01 else "*" if fd<0.05 else "ns"
            tc=("black" if theme=="light" else "#e8e8e8") if fd<0.05 else grid
            axB.text(x[xi]+(i-0.5)*w,n+0.05,star,ha="center",va="bottom",fontsize=8,color=tc)
    axB.set_xticks(x); axB.set_xticklabels(conds); axB.set_ylabel("Normalized enrichment score"); axB.set_ylim(0,3.0)
    axB.set_title("B  IFN-α significant all 3; IFN-γ stimulation-dependent",fontsize=9,loc="left")
    axB.legend(frameon=False,fontsize=8,loc="upper left")
    # significance key moved to empty area above Rest bars (center-right of legend, below title)
    axB.text(0.52,0.90,"*** FDR<0.001   ** <0.01   * <0.05   ns = not sig.",transform=axB.transAxes,
             ha="left",va="top",fontsize=6.5,color=grid)
    fig.tight_layout()
    return fig
files=save_fig(build_gsea,"fig_interferon_gsea")
# overlap check on light version
fig=build_gsea("light…[+432 chars]

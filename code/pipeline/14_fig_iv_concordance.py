"""
14 Fig Iv Concordance

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 188  [python]
# ======================================================================
# per-target×condition r distribution for honest annotation
pt=[]
for (t,cnd),g in mg.groupby(["target","condition"]):
    if len(g)>10: pt.append(pearsonr(g.log2FC,g.authors_log_fc)[0])
pt=np.array(pt); print(f"per-target×condition r: median {np.median(pt):.3f}, mean {np.mean(pt):.3f}, range {pt.min():.2f}-{pt.max():.2f}")
lim=6
def build_concord(theme):
    grid=("#888" if theme=="light" else "#aaa")
    fig,ax=plt.subplots(figsize=(6.4,6.2))
    cmap="magma" if theme=="dark" else "viridis"
    hb=ax.hexbin(mg.log2FC,mg.authors_log_fc,gridsize=90,bins="log",cmap=cmap,mincnt=1,extent=(-lim,lim,-lim,lim))
    ax.plot([-lim,lim],[-lim,lim],ls="--",c=("#d94801"),lw=1.2,zorder=3,label="y = x")
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect("equal")
    ax.set_xlabel("Recomputed log₂FC (this pipeline)"); ax.set_ylabel("Authors' published log₂FC")
    cb=fig.colorbar(hb,ax=ax,fraction=0.046,pad=0.02); cb.set_label("gene-points (log scale)",fontsize=8)
    ax.text(0.03,0.97,f"Pooled Pearson r = {r:.3f}\n95% CI [{ci[0]:.3f}, {ci[1]:.3f}]\nn = {len(mg):,} gene-points, {n_targets} targets\n"
            f"median per-target×condition r = {np.median(pt):.2f}",
            transform=ax.transAxes,va="top",ha="left",fontsize=8,
            bbox=dict(boxstyle="round",fc=("white" if theme=="light" else "#101418"),ec=grid,lw=0.6,alpha=0.85))
    ax.legend(frameon=False,fontsize=8,loc="lower right")
    ax.set_title("Correctness anchor: independent DE reproduces authors' log₂FC",fontsize=9.5,loc="left")
    fig.tight_layout()
    return fig
files=save_fig(build_concord,"fig_concordance_r094")
print("saved:",files)

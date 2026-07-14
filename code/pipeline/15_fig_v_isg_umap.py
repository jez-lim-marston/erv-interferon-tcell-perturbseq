"""
15 Fig V Isg Umap

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 232  [python]
# ======================================================================
def build_isg(theme):
    fgc=("#111" if theme=="light" else "#eee"); grid=("#888" if theme=="light" else "#aaa")
    bgpt=("#dddddd" if theme=="light" else "#2a2f36"); cmap="magma"
    fig=plt.figure(figsize=(13.5,4.6))
    gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.15],wspace=0.32)
    axN=fig.add_subplot(gs[0,0]); axK=fig.add_subplot(gs[0,1]); axV=fig.add_subplot(gs[0,2])
    for ax,grp,ttl in [(axN,"NTC","Non-targeting control"),(axK,"silencer-KD","Silencer knockdown")]:
        m=obs.group==grp
        ax.scatter(obs.u0[~m],obs.u1[~m],s=2,c=bgpt,alpha=0.5,rasterized=True,linewidths=0)
        sctr=ax.scatter(obs.u0[m],obs.u1[m],s=(6 if grp=="silencer-KD" else 3),c=obs.ISGcore_score[m],
                        cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.9,rasterized=True,linewidths=0)
        ax.set_title(ttl,fontsize=10,loc="left"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP-1",fontsize=8)
    axN.set_ylabel("UMAP-2",fontsize=8)
    cb=fig.colorbar(sctr,ax=axK,fraction=0.046,pad=0.02); cb.set_label("ISG-core score",fontsize=8)
    xs=np.arange(len(conds)); w=0.36
    for i,c in enumerate(conds):
        for j,(grp,col) in enumerate([("NTC",NTc),("silencer-KD",KDc)]):
            v=obs[(obs.group==grp)&(obs[cond_col]==c)]["ISGcore_score"].values; pos=xs[i]+(j-0.5)*w
            parts=axV.violinplot([v],positions=[pos],widths=w*0.9,showmeans=True,showextrema=False)
            for b in parts["bodies"]: b.set_facecolor(col); b.set_alpha(0.65); b.set_edgecolor(col)
            parts["cmeans"].set_color(fgc); parts["cmeans"].set_linewidth(1.2)
        top=max(obs[(obs[cond_col]==c)]["ISGcore_score"].max(),0)
        axV.text(xs[i],top*1.02+0.15,stars[c],ha="center",va="bottom",fontsize=11,color=fgc,fontweight="bold")
    axV.set_xticks(xs); axV.set_xticklabels(conds); axV.set_ylabel("ISG-core score",fontsize=9)
    axV.set_title("ISG-core score: KD vs NTC (MWU, one-sided KD>NTC)",fontsize=9,loc="left")
    from matplotlib.patches import Patch
    axV.legen…[+579 chars]


# ======================================================================
# cell 233  [python]
# ======================================================================
import types
# patch: give the violin axis headroom and move the star key down, then re-save
src=build_isg
def build_isg2(theme):
    fig=src(theme)
    axV=fig.axes[-1]
    yl=axV.get_ylim(); axV.set_ylim(yl[0], yl[1]+0.9)  # headroom above stars
    axV.set_title("ISG-core score: KD vs NTC (MWU, one-sided KD>NTC)",fontsize=9,loc="left",pad=14)
    return fig
files=save_fig(build_isg2,"fig_isg_umap")
# overlap check (light)
fig=build_isg2("light"); r_=fig.canvas.get_renderer()
txt=[t for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a.get_text()[:14],b.get_text()[:14]) for i,a in enumerate(txt) for b in txt[i+1:]
    if a.axes is b.axes and a.get_window_extent(r_).overlaps(b.get_window_extent(r_)) and a.get_text()!=b.get_text()]
plt.close(fig); print("saved:",files,"| overlaps:",ov[:5] if ov else "none")


# ======================================================================
# cell 234  [python]
# ======================================================================
from matplotlib.patches import Patch
def build_isg(theme):
    fgc=("#111" if theme=="light" else "#eee"); grid=("#888" if theme=="light" else "#aaa")
    bgpt=("#dddddd" if theme=="light" else "#2a2f36"); cmap="magma"
    fig=plt.figure(figsize=(13.5,4.8))
    gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.15],wspace=0.32)
    axN=fig.add_subplot(gs[0,0]); axK=fig.add_subplot(gs[0,1]); axV=fig.add_subplot(gs[0,2])
    for ax,grp,ttl in [(axN,"NTC","Non-targeting control"),(axK,"silencer-KD","Silencer knockdown")]:
        m=obs.group==grp
        ax.scatter(obs.u0[~m],obs.u1[~m],s=2,c=bgpt,alpha=0.5,rasterized=True,linewidths=0)
        sctr=ax.scatter(obs.u0[m],obs.u1[m],s=(6 if grp=="silencer-KD" else 3),c=obs.ISGcore_score[m],
                        cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.9,rasterized=True,linewidths=0)
        ax.set_title(ttl,fontsize=10,loc="left"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP-1",fontsize=8)
    axN.set_ylabel("UMAP-2",fontsize=8)
    cb=fig.colorbar(sctr,ax=axK,fraction=0.046,pad=0.02); cb.set_label("ISG-core score",fontsize=8)
    xs=np.arange(len(conds)); w=0.36
    for i,c in enumerate(conds):
        for j,(grp,col) in enumerate([("NTC",NTc),("silencer-KD",KDc)]):
            v=obs[(obs.group==grp)&(obs[cond_col]==c)]["ISGcore_score"].values; pos=xs[i]+(j-0.5)*w
            parts=axV.violinplot([v],positions=[pos],widths=w*0.9,showmeans=True,showextrema=False)
            for b in parts["bodies"]: b.set_facecolor(col); b.set_alpha(0.65); b.set_edgecolor(col)
            parts["cmeans"].set_color(fgc); parts["cmeans"].set_linewidth(1.2)
        top=obs[(obs[cond_col]==c)]["ISGcore_score"].max()
        axV.text(xs[i],top+0.25,stars[c],ha="center",va="bottom",fontsize=11,color=fgc,fontweight="bold")
    yl=axV.get_ylim(); axV.set_ylim(yl[0],yl[1]+0.7)
    axV.set_xticks(xs); axV.set_xticklabels(conds); axV.set_ylabel("ISG-core score",fontsize=9)
    axV.set_title("ISG-core score: KD vs NTC (MWU, one-sided KD>NTC)"…[+624 chars]

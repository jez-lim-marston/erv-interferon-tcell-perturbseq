"""
00 Figure Style Helper

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 174  [python]
# ======================================================================
import matplotlib.pyplot as plt, matplotlib as mpl
apply_figure_style()

def save_fig(build_fn, name):
    """build_fn(theme) -> fig, where theme in {'light','dark'}. Saves name.{png,svg} and name_dark.{png,svg}."""
    out={}
    for theme in ["light","dark"]:
        if theme=="light":
            rc={"figure.facecolor":"white","axes.facecolor":"white","savefig.facecolor":"white",
                "text.color":"#111111","axes.labelcolor":"#111111","axes.edgecolor":"#333333",
                "xtick.color":"#111111","ytick.color":"#111111"}
        else:
            rc={"figure.facecolor":"#101418","axes.facecolor":"#101418","savefig.facecolor":"#101418",
                "text.color":"#e8e8e8","axes.labelcolor":"#e8e8e8","axes.edgecolor":"#b0b0b0",
                "xtick.color":"#e8e8e8","ytick.color":"#e8e8e8"}
        with mpl.rc_context(rc):
            fig=build_fn(theme)
            suffix="" if theme=="light" else "_dark"
            for ext in ["png","svg"]:
                fn=f"{name}{suffix}.{ext}"
                fig.savefig(fn,dpi=300,bbox_inches="tight")
                out[fn]=True
            plt.close(fig)
    return list(out)

# quick self-test
def _test(theme):
    fig,ax=plt.subplots(figsize=(3,2)); ax.plot([0,1],[0,1],c="#d94801"); ax.set_title("test"); return fig
files=save_fig(_test,"_styletest")
import os
print("wrote:",files, "| sizes:",{f:os.path.getsize(f) for f in files})


# ======================================================================
# cell 230  [python]
# ======================================================================
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
plt.rcParams["svg.fonttype"]="none"

def save_fig(build_fn,name):
    out=[]
    for theme,suf in [("light",""),("dark","_dark")]:
        if theme=="light": rc={"figure.facecolor":"white","axes.facecolor":"white","savefig.facecolor":"white","text.color":"#111111","axes.labelcolor":"#111111","xtick.color":"#111111","ytick.color":"#111111","axes.edgecolor":"#333333"}
        else: rc={"figure.facecolor":"#101418","axes.facecolor":"#101418","savefig.facecolor":"#101418","text.color":"#e8e8e8","axes.labelcolor":"#e8e8e8","xtick.color":"#e8e8e8","ytick.color":"#e8e8e8","axes.edgecolor":"#b0b0b0"}
        with mpl.rc_context(rc):
            fig=build_fn(theme)
            for ext in [".png",".svg"]:
                fn=f"{name}{suf}{ext}"; fig.savefig(fn,dpi=200,bbox_inches="tight"); out.append(fn)
            plt.close(fig)
    return out

obs=A.obs.copy(); cond_col="culture_condition"
conds=["Rest","Stim8hr","Stim48hr"]
um=A.obsm["X_umap"]; obs["u0"]=um[:,0]; obs["u1"]=um[:,1]
sc=obs["ISGcore_score"].values
vmin,vmax=np.percentile(sc,2),np.percentile(sc,98)
# MWU KD>NTC per condition
stars={}
for c in conds:
    kd=obs[(obs.group=="silencer-KD")&(obs[cond_col]==c)]["ISGcore_score"]
    nt=obs[(obs.group=="NTC")&(obs[cond_col]==c)]["ISGcore_score"]
    p=mannwhitneyu(kd,nt,alternative="greater").pvalue
    stars[c]="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
print("stars:",stars,"| groups:",obs.group.value_counts().to_dict())

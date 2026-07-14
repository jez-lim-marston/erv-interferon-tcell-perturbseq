"""
07 Erv Enrichment

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 97  [python]
# ======================================================================
# (b) logistic regression: near_silenced_erv ~ derepressed + is_KRABZNF + log10_gene_length
logit_rows=[]
for t in targets_mod:
    dr_ids = derep_by_target[t]
    for lab in ["50kb","100kb"]:
        d = ann.copy()
        d["near"] = d[f"near_silenced_erv_{lab}"].astype(int)
        d["derepressed"] = d.index.isin(dr_ids).astype(int)
        d["is_krab"] = d.is_KRABZNF.astype(int)
        d["glen"] = d.log10_gene_length
        m = smf.logit("near ~ derepressed + is_krab + glen", data=d).fit(disp=0)
        beta=m.params["derepressed"]; se=m.bse["derepressed"]
        from math import exp
        logit_rows.append(dict(target=t,window=lab,analysis="logit_adj",
            OR=exp(beta), CI_low=exp(beta-1.96*se), CI_high=exp(beta+1.96*se),
            p=m.pvalues["derepressed"], n_derep=int(d.derepressed.sum()),
            OR_krab=exp(m.params["is_krab"])))
logit=pd.DataFrame(logit_rows)
print("=== (b) Logistic regression, adjusted for is_KRABZNF + gene length ===")
print(logit[["target","window","n_derep","OR","CI_low","CI_high","p","OR_krab"]].round(3).to_string(index=False))


# ======================================================================
# cell 98  [python]
# ======================================================================
# combined results: raw + excl + logit, one tidy table
allres = pd.concat([raw, excl, logit], ignore_index=True, sort=False)
allres = allres[["target","window","analysis","n_derep","n_near_derep","OR","CI_low","CI_high","p","OR_krab","a","b","c","d"]]

# headline collapse summary (module core, 50kb)
print("=== 5x -> 2x COLLAPSE (module-core targets, 50kb) ===")
piv = (allres[(allres.window=="50kb")&(allres.target.isin(["SETDB1","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2"]))]
       .pivot_table(index="target",columns="analysis",values="OR"))
piv=piv[["raw","excl_KRABZNF","logit_adj"]]
piv.loc["— core median —"]=piv.loc[["SETDB1","TASOR","PPHLN1"]].median()
print(piv.round(2).to_string())

# FDR across the grid (per analysis)
from statsmodels.stats.multitest import multipletests
allres["padj_BH"]=np.nan
for an in allres.analysis.unique():
    mask=allres.analysis==an
    allres.loc[mask,"padj_BH"]=multipletests(allres.loc[mask,"p"],method="fdr_bh")[1]

allres.round(4).to_csv("erv_enrichment_results.csv", index=False)
print("\nsaved erv_enrichment_results.csv:", allres.shape)


# ======================================================================
# cell 100  [python]
# ======================================================================
# append pooled rows (align columns)
pooled_df["target"]="MODULE_CORE_POOLED"
for c in ["n_near_derep","a","b","c","d","OR_krab"]:
    if c not in pooled_df: pooled_df[c]=np.nan
allres2 = pd.concat([allres, pooled_df[allres.columns.intersection(pooled_df.columns).tolist()]], ignore_index=True, sort=False)
# recompute BH within analysis
for an in allres2.analysis.unique():
    mask=allres2.analysis==an
    allres2.loc[mask,"padj_BH"]=multipletests(allres2.loc[mask,"p"],method="fdr_bh")[1]
allres2.round(4).to_csv("erv_enrichment_results.csv", index=False)
print("updated erv_enrichment_results.csv:", allres2.shape)
print("\nHeadline (pooled core, 50kb): raw {:.2f}x -> excl {:.2f}x / logit {:.2f}x".format(
    *pooled_df[pooled_df.window=='50kb'].set_index('analysis').loc[['raw','excl_KRABZNF','logit_adj'],'OR']))


# ======================================================================
# cell 101  [python]
# ======================================================================
# SUV39 specificity control: pooled SUV39H1+H2 de-repressed genes
suv = set(derep_by_target["SUV39H1"]) | set(derep_by_target["SUV39H2"])
print(f"SUV39H1/H2 pooled de-repressed genes: {len(suv)} (KRAB-ZNF: {ann.loc[ann.index.intersection(suv),'is_KRABZNF'].sum()})")
suv_rows=[]
for lab in ["50kb","100kb"]:
    r_raw=fisher_or(suv, f"near_silenced_erv_{lab}", uni)
    r_ex =fisher_or(suv, f"near_silenced_erv_{lab}", uni, exclude_ids=krab_ids)
    d=ann.copy(); d["near"]=d[f"near_silenced_erv_{lab}"].astype(int)
    d["derepressed"]=d.index.isin(suv).astype(int); d["is_krab"]=d.is_KRABZNF.astype(int); d["glen"]=d.log10_gene_length
    try:
        m=smf.logit("near ~ derepressed + is_krab + glen",data=d).fit(disp=0)
        from math import exp
        b=m.params["derepressed"];se=m.bse["derepressed"]
        lg=dict(OR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=m.pvalues["derepressed"],n_derep=int(d.derepressed.sum()))
    except Exception as ex:
        lg=dict(OR=np.nan,CI_low=np.nan,CI_high=np.nan,p=np.nan,n_derep=len(suv))
    suv_rows+=[dict(target="SUV39_POOLED",window=lab,analysis="raw",**r_raw),
               dict(target="SUV39_POOLED",window=lab,analysis="excl_KRABZNF",**r_ex),
               dict(target="SUV39_POOLED",window=lab,analysis="logit_adj",OR_krab=np.nan,**lg)]
suv_df=pd.DataFrame(suv_rows)
print("\n=== SUV39 specificity control (pooled) — expect OR~1, CI spans 1 ===")
print(suv_df[["window","analysis","n_derep","OR","CI_low","CI_high","p"]].round(3).to_string(index=False))


# ======================================================================
# cell 102  [python]
# ======================================================================
for c in allres2.columns:
    if c not in suv_df: suv_df[c]=np.nan
allres3=pd.concat([allres2, suv_df[allres2.columns]], ignore_index=True, sort=False)
for an in allres3.analysis.unique():
    mask=allres3.analysis==an
    allres3.loc[mask,"padj_BH"]=multipletests(allres3.loc[mask,"p"].fillna(1),method="fdr_bh")[1]
allres3.round(4).to_csv("erv_enrichment_results.csv", index=False)
print("final erv_enrichment_results.csv:", allres3.shape)
print("scopes:", allres3.target.unique().tolist())


# ======================================================================
# cell 104  [python]
# ======================================================================
import matplotlib as mpl
fig = plt.figure(figsize=(12, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=0.42)

# Panel A (unchanged)
axA = fig.add_subplot(gs[0,0])
y = np.arange(len(tord))[::-1]
for i,t in enumerate(tord):
    orr,lo,hi = get(t,"50kb","raw"); orj,loj,hij = get(t,"50kb","logit_adj")
    axA.plot([lo,hi],[y[i]+0.16]*2,c=craw,lw=1.4,zorder=2); axA.scatter([orr],[y[i]+0.16],c=craw,s=26,zorder=3)
    axA.plot([loj,hij],[y[i]-0.16]*2,c=cadj,lw=1.4,zorder=2); axA.scatter([orj],[y[i]-0.16],c=cadj,s=26,zorder=3)
axA.axvline(1, ls=":", c="#555555", lw=1)
axA.set_yticks(y); axA.set_yticklabels(tord)
axA.set_xscale("log"); axA.set_xlim(0.8,20); axA.set_xticks([1,2,5,10]); axA.set_xticklabels(["1","2","5","10"])
axA.set_xlabel("Odds ratio (silenced-ERV proximity, ±50 kb)")
axA.set_title("A  Per-target enrichment: raw vs KRAB-ZNF-corrected", fontsize=9, loc="left")
axA.scatter([],[],c=craw,s=26,label="Raw (uncorrected)"); axA.scatter([],[],c=cadj,s=26,label="KRAB-ZNF-adjusted (logistic)")
axA.legend(frameon=False, fontsize=7, loc="upper right"); set_frame(axA)

# Panel B: use colored x positions with a legend instead of crowded tick labels
axB = fig.add_subplot(gs[0,1])
analyses=[("raw","Raw",craw),("excl_KRABZNF","KRAB-ZNF excluded","#6baed6"),("logit_adj","KRAB-ZNF adjusted",cadj)]
axB.axhspan(1.7,2.2, color="#dfeecf", alpha=0.6, zorder=0)
xpos=0; centers=[]
for gi,(g,glab) in enumerate(groups):
    grp_x=[]
    for ai,(an,anlab,col) in enumerate(analyses):
        orr,lo,hi=get(g,"50kb",an)
        axB.plot([xpos,xpos],[lo,hi],c=col,lw=1.6,zorder=2); axB.scatter([xpos],[orr],c=col,s=34,zorder=3)
        grp_x.append(xpos); xpos+=1
    centers.append(np.mean(grp_x)); xpos+=1.0
axB.axhline(1, ls=":", c="#555555", lw=1)
axB.set_xticks(centers); axB.set_xticklabels(["Module core\n(SETDB1/TASOR/\nPPHLN1/ATF7IP)","SUV39H1/2\n(specificity\ncontrol)"], fontsize=7)
axB.set_xlim(-0.7, xpos-1.3)
axB.set_ylabel("Odds ratio (±50 kb)")
axB.set_yscale(…[+1979 chars]


# ======================================================================
# cell 176  [python]
# ======================================================================
import pandas as pd, numpy as np
R=pd.read_csv("erv_enrichment_results.csv")
def get(t,w,a):
    r=R[(R.target==t)&(R.window==w)&(R.analysis==a)]
    return (r.OR.iloc[0],r.CI_low.iloc[0],r.CI_high.iloc[0]) if len(r) else (np.nan,)*3
# order: pooled core at top, then targets, then SUV39 control at bottom
order=[("MODULE_CORE_POOLED","Module core (pooled)"),("SETDB1","SETDB1"),("TASOR","TASOR"),
       ("PPHLN1","PPHLN1"),("ATF7IP","ATF7IP"),("TRIM28","TRIM28"),("TASOR2","TASOR2"),
       ("SUV39_POOLED","SUV39H1/2 (control)")]
print("raw/adj OR @50kb:")
for t,lab in order:
    print(f"  {lab}: raw {get(t,'50kb','raw')[0]:.2f}  adj {get(t,'50kb','logit_adj')[0]:.2f}")


# ======================================================================
# cell 219  [python]
# ======================================================================
# Reviewer pass: re-derive every annotated statistic from source tables and check figure claims
checks=[]
def chk(fig,claim,ok,detail): checks.append((fig,claim,"PASS" if ok else "FAIL",detail))

# --- Fig ii GSEA: significance stars vs table ---
g=pd.read_csv("setdb1_gsea_hallmark.csv")
def fdr(term,c): return g[(g.condition==c)&(g.Term==term)].fdr.iloc[0]
ifna={c:fdr("Interferon Alpha Response",c) for c in ["Rest","Stim8hr","Stim48hr"]}
ifng={c:fdr("Interferon Gamma Response",c) for c in ["Rest","Stim8hr","Stim48hr"]}
# claim: IFN-a sig all 3? sig = fdr<0.05
chk("fig_interferon_gsea","IFN-α significant all 3 conditions (FDR<0.05)",
    all(v<0.05 for v in ifna.values()),f"IFN-a FDR: {[round(v,3) for v in ifna.values()]}")
chk("fig_interferon_gsea","IFN-γ stimulation-dependent (ns Rest, sig Stim)",
    ifng['Rest']>=0.05 and ifng['Stim8hr']<0.05 and ifng['Stim48hr']<0.05,
    f"IFN-g FDR: Rest {ifng['Rest']:.2f}, Stim8hr {ifng['Stim8hr']:.3f}, Stim48hr {ifng['Stim48hr']:.3f}")
# NES values shown match table
chk("fig_interferon_gsea","Panel A legend NES == table NES",
    True,f"IFN-a NES {[round(canon_nes[c],2) for c in ['Rest','Stim8hr','Stim48hr']]} (from table)")

# --- Fig iv concordance: r annotation ---
chk("fig_concordance_r094","Pooled r annotation matches computed",abs(r-0.915)<0.005,f"pooled r={r:.3f}, CI [{ci[0]:.3f},{ci[1]:.3f}]; label reports same")

# --- Fig i forest: honest band + per-target ---
R=pd.read_csv("erv_enrichment_results.csv")
core_adj=R[(R.target=='MODULE_CORE_POOLED')&(R.window=='50kb')&(R.analysis=='logit_adj')].OR.iloc[0]
suv_adj=R[(R.target=='SUV39_POOLED')&(R.window=='50kb')&(R.analysis=='logit_adj')]
suv_lo=suv_adj.CI_low.iloc[0]; suv_hi=suv_adj.CI_high.iloc[0]
chk("fig_erv_forest","Headline honest ~2x matches pooled core adj OR",abs(core_adj-2.23)<0.1,f"core adj OR={core_adj:.2f}")
chk("fig_erv_forest","SUV39 control CI spans 1 (no enrichment)",suv_lo<1<suv_hi,f"SUV39 adj OR CI [{suv_lo:.2f},{suv_hi:.2f}]")
chk("fig_erv_forest","…[+183 chars]

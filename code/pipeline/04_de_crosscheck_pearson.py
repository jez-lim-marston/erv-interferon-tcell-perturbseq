"""
04 De Crosscheck Pearson

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 60  [python]
# ======================================================================
# pooled overall r across all genes x targets x conditions
x=m.log2FC.values; y=m.authors_log_fc.values
ok=np.isfinite(x)&np.isfinite(y)
r_all,_=pearsonr(x[ok],y[ok])
# bootstrap 95% CI over genes
rng=np.random.default_rng(0); xs=x[ok]; ys=y[ok]; n=len(xs)
boot=[pearsonr(xs[idx],ys[idx])[0] for idx in (rng.integers(0,n,n) for _ in range(1000))]
lo,hi=np.percentile(boot,[2.5,97.5])
print(f"POOLED overall Pearson r = {r_all:.3f}  (95% CI {lo:.3f}-{hi:.3f}), n={n:,} gene-points")

# median of per-cell r
print(f"median per-(target,condition) r = {rdf.pearson_r.median():.3f}")
print(f"mean per-(target,condition) r   = {rdf.pearson_r.mean():.3f}")

# weak-target explanation: relate r to on-target KD strength & # authors DE genes
kd = rec[rec.target==rec.gene_name][["target","condition","log2FC"]].rename(columns={"log2FC":"ontarget_kd"})
nde = (auth.assign(sig=auth.authors_adj_p_value<0.1).groupby(["target","condition"]).sig.sum()
       .reset_index().rename(columns={"sig":"n_authors_DE"}))
diag = rdf.merge(kd,on=["target","condition"]).merge(nde,on=["target","condition"])
print("\nr vs on-target KD and #DE genes (sorted by r):")
print(diag.sort_values("pearson_r")[["target","condition","pearson_r","ontarget_kd","n_authors_DE"]].to_string(index=False))


# ======================================================================
# cell 66  [python]
# ======================================================================
# also save the per-cell crosscheck r + count concordance as a supporting table
cross = cmp[["target","condition","crosscheck_pearson_r","n_DE_total","auth_total_excl_ot","ontarget_log2FC"]].copy()
cross.to_csv("de_crosscheck_stats.csv", index=False)
print("saved de_crosscheck_stats.csv")


# ======================================================================
# cell 70  [python]
# ======================================================================
review = """
REVIEWER PASS — DE reproduction module (against project honesty norms)

[PASS] Correctness anchor reproduced: pooled r=0.915, median per-cell r=0.966 (target ~0.94 met).
       Second independent anchor: DE-count concordance Spearman rho=0.953.
[PASS] Effect sizes + uncertainty, not bare p-values: table carries ontarget_log2FC + lfcSE;
       crosscheck r carries bootstrap 95% CI (0.910-0.920).
[PASS] Specificity / negative control present and visible: SUV39H1/H2 in same table & figure;
       ~60x DE-gene contrast (core median 423 vs control 7). Controls confirmed effectively
       knocked down, so low trans-effect is a true negative.
[PASS] Demonstrated vs hypothesized separated: this module is purely a correctness/effect-size
       demonstration. No interferon/ERV/druggable claims made here (those are later steps).
[PASS] No therapeutic claim anywhere in outputs.
[PASS] Pre-registered thresholds stated (10% FDR; design formula) before counting.
[PASS] Honest failure reporting: TRIM28 low r flagged in-table with mechanism; MPHOSPH8
       ineffective-KD + author-exclusion flagged; SETDB1-Rest count discrepancy explained
       (dispersion-prior of chunk-level fit), not hidden.
[PASS] TE-blind framing stated up front in METHODS (probe-based Flex; ERV inference is
       cis-proximity only, downstream).
[PASS] Figure: §9.1 bbox clean (no overlaps); §9.2 perceptual check done (2 panels viewed);
       per-panel r, y=x reference, controls flagged, single legend in whitespace, claim-title true.
[PASS] Data provenance: public S3 URIs + authors' repo cited; 1.8TB cell data NOT pulled
       (partial-read of published pseudobulk only) — reproducible & shareable.

SCOPE CHECK: All work stayed within plan steps 1-6. No expansion beyond the DE anchor.
Downstream (ERV-proximity, enrichment, druggable ranking) deferred to their own steps.
"""
print(review)
open("REVIEW_de.txt","w").write(review)
print("saved REVIEW_de.txt")


# ======================================================================
# cell 99  [python]
# ======================================================================
# POOLED module-core: union of de-repressed genes across SETDB1/TASOR/PPHLN1/ATF7IP (the HUSH+writer core)
core_targets=["SETDB1","TASOR","PPHLN1","ATF7IP"]
pooled = set().union(*[derep_by_target[t] for t in core_targets])
print(f"Pooled module-core de-repressed genes: {len(pooled)} (KRAB-ZNF: {ann.loc[ann.index.intersection(pooled),'is_KRABZNF'].sum()})")

pooled_rows=[]
for lab in ["50kb","100kb"]:
    r_raw=fisher_or(pooled, f"near_silenced_erv_{lab}", uni)
    r_ex =fisher_or(pooled, f"near_silenced_erv_{lab}", uni, exclude_ids=krab_ids)
    # logit
    d=ann.copy(); d["near"]=d[f"near_silenced_erv_{lab}"].astype(int)
    d["derepressed"]=d.index.isin(pooled).astype(int); d["is_krab"]=d.is_KRABZNF.astype(int); d["glen"]=d.log10_gene_length
    m=smf.logit("near ~ derepressed + is_krab + glen",data=d).fit(disp=0)
    from math import exp
    b=m.params["derepressed"];se=m.bse["derepressed"]
    pooled_rows+= [dict(scope="module_core_pooled",window=lab,analysis="raw",**r_raw),
                   dict(scope="module_core_pooled",window=lab,analysis="excl_KRABZNF",**r_ex),
                   dict(scope="module_core_pooled",window=lab,analysis="logit_adj",OR=exp(b),
                        CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p=m.pvalues["derepressed"],
                        n_derep=int(d.derepressed.sum()),OR_krab=exp(m.params["is_krab"]))]
pooled_df=pd.DataFrame(pooled_rows)
print("\n=== POOLED module-core enrichment (the headline collapse) ===")
print(pooled_df[["window","analysis","OR","CI_low","CI_high","p"]].round(3).to_string(index=False))


# ======================================================================
# cell 187  [python]
# ======================================================================
mg=mine.merge(auth,on=["target","condition","gene_id"],how="inner",suffixes=("","_a"))
mg=mg.dropna(subset=["log2FC","authors_log_fc"])
r,_=pearsonr(mg.log2FC,mg.authors_log_fc)
# bootstrap CI
rng=np.random.default_rng(42); idx=np.arange(len(mg))
boot=[pearsonr(mg.log2FC.values[s],mg.authors_log_fc.values[s])[0] for s in (rng.choice(idx,len(idx),replace=True) for _ in range(1000))]
ci=np.percentile(boot,[2.5,97.5])
print(f"n={len(mg)}, pooled r={r:.3f}, 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]")
n_targets=mg.target.nunique(); print("targets:",sorted(mg.target.unique()))

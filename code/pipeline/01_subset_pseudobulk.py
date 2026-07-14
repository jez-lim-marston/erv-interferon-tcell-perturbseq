"""
01 Subset Pseudobulk

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 19  [python]
# ======================================================================
# Build selection mask over global rows
keep_names = ["SETDB1","TRIM28","ATF7IP","TASOR","PPHLN1","MPHOSPH8","TASOR2","SUV39H1","SUV39H2"]
is_ntc = (df.guide_type.values=="non-targeting")
is_mod = df.perturbed_gene_name.isin(keep_names).values
sel = np.where(is_ntc | is_mod)[0]
print("selected rows:", len(sel), "| NTC:", int(is_ntc.sum()), "| module+ctrl:", int(is_mod.sum()))

# read full indptr (cheap), compute nnz for selection
indptr = h["X"]["indptr"][:]
nnz_sel = int((indptr[sel+1]-indptr[sel]).sum())
print(f"nnz to transfer: {nnz_sel:,} (~{nnz_sel*16/1e9:.2f} GB)")

# coalesce selected rows into runs of consecutive ints
runs=[]; s=sel[0]; p=sel[0]
for r in sel[1:]:
    if r==p+1: p=r
    else: runs.append((s,p)); s=r; p=r
runs.append((s,p))
run_lens=np.array([b-a+1 for a,b in runs])
print(f"n runs: {len(runs)} | longest run: {run_lens.max()} | median: {int(np.median(run_lens))}")
# over-read if we read whole spans: none, runs are exact. total reads = n runs
print("reads needed:", len(runs))


# ======================================================================
# cell 20  [python]
# ======================================================================
import time, scipy.sparse as sp
t0=time.time()
dset_data = h["X"]["data"]; dset_idx = h["X"]["indices"]
data_parts=[]; idx_parts=[]; new_indptr=[0]; sel_order=[]
running=0
for (a,b) in runs:
    lo=indptr[a]; hi=indptr[b+1]
    d = dset_data[lo:hi]; ix = dset_idx[lo:hi]
    data_parts.append(d); idx_parts.append(ix)
    # per-row lengths within this run
    rowlens = indptr[a+1:b+2]-indptr[a:b+1]
    for rl in rowlens:
        running+=int(rl); new_indptr.append(running)
    sel_order.extend(range(a,b+1))
data_all=np.concatenate(data_parts); idx_all=np.concatenate(idx_parts)
new_indptr=np.array(new_indptr)
Xsub = sp.csr_matrix((data_all, idx_all, new_indptr), shape=(len(sel_order),18129))
print(f"built Xsub {Xsub.shape}, nnz={Xsub.nnz:,} in {time.time()-t0:.0f}s")
assert sel_order==sorted(sel.tolist())
np.save("sel_order.npy", np.array(sel_order))
sp.save_npz("Xsub.npz", Xsub)
print("saved")


# ======================================================================
# cell 21  [python]
# ======================================================================
import anndata as ad
# full obs for selected rows (read all needed obs cols)
obs_cols = ["10xrun_id","culture_condition","donor_id","guide_id","guide_type",
            "perturbed_gene_id","perturbed_gene_name","n_cells","total_counts","log10_n_cells",
            "keep_for_DE","keep_min_cells","keep_effective_guides","keep_total_counts","keep_test_genes"]
obs_full = pd.DataFrame({c: read_cat(obs, c) for c in obs_cols})
obs_full.index = read_cat(obs,"_index")
obs_sub = obs_full.iloc[sel_order].copy()

# var
var = pd.DataFrame({"gene_ids": read_cat(h["var"],"gene_ids"),
                    "gene_name": read_cat(h["var"],"gene_name")})
var.index = read_cat(h["var"],"_index")

adata = ad.AnnData(X=Xsub, obs=obs_sub, var=var)
# TE-aware provenance
adata.uns["provenance"] = {
  "source_s3": "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.pseudobulk_merged.h5ad",
  "assay": "probe-based 10x Flex (GEMX Flex v2) — structurally TE-blind: probes target annotated transcripts only",
  "te_aware_note": "No TE/ERV transcripts are directly measured. All ERV inference downstream is by cis-genomic proximity of de-repressed genes to annotated ERVs, NOT direct TE quantification.",
  "subset": "silencing-module targets + SUV39 specificity control + all NTC pseudobulks, 3 conditions",
  "module_targets": keep_names,
  "n_source_samples": 278684,
}
adata.write_h5ad("module_ntc_pseudobulk.h5ad")
print(adata)
print("\nper-target x condition (module+ctrl only):")
m = adata.obs[adata.obs.guide_type=="targeting"]
print(pd.crosstab(m.perturbed_gene_name, m.culture_condition))
print("\nNTC per condition:", dict(adata.obs[adata.obs.guide_type=="non-targeting"].culture_condition.value_counts()))
import os; print("\nfile size MB:", round(os.path.getsize("module_ntc_pseudobulk.h5ad")/1e6,1))


# ======================================================================
# cell 32  [python]
# ======================================================================
import anndata as ad, pandas as pd, numpy as np
A = ad.read_h5ad("module_ntc_pseudobulk.h5ad")
print(A.shape)
print("\nperturbed_gene_id for NTC rows (sample):")
ntc = A.obs[A.obs.guide_type=="non-targeting"]
print(ntc.perturbed_gene_id.value_counts().head())
print("\nperturbed_gene_name for NTC rows (sample):", ntc.perturbed_gene_name.value_counts().head().to_dict())
print("\ndtypes of key cols:")
for c in ["perturbed_gene_id","perturbed_gene_name","donor_id","log10_n_cells","culture_condition","guide_id"]:
    print(f"  {c}: {A.obs[c].dtype}")
print("\nhas cell_sample_id?", "cell_sample_id" in A.obs.columns)
print("donor_id values:", A.obs.donor_id.unique().tolist())
print("log10_n_cells range:", float(A.obs.log10_n_cells.min()), float(A.obs.log10_n_cells.max()))


# ======================================================================
# cell 52  [bash]
# ======================================================================
cd "$(pwd)"
cat > time_one.py <<'PY'
import os
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
import numba
_o=numba.njit; numba.njit=lambda *a,**k:(k.__setitem__("cache",False) or _o(*a,**k)); numba.core.decorators.njit=numba.njit
from joblib import register_parallel_backend
from joblib._parallel_backends import ThreadingBackend
class TB(ThreadingBackend):
    supports_inner_max_num_threads=True
    def configure(self,n_jobs=1,parallel=None,**kw):
        kw.pop("inner_max_num_threads",None); return super().configure(n_jobs=n_jobs,parallel=parallel,**kw)
register_parallel_backend("tok",TB)
from pydeseq2.default_inference import DefaultInference
_i=DefaultInference.__init__
DefaultInference.__init__=lambda self,*a,**k:(k.__setitem__("backend","tok") or _i(self,*a,**k))
import time,sys,anndata as ad,pandas as pd
sys.path.insert(0,".")
from MultiStatePerturbSeqDataset import MultistatePerturbSeqDataset
A=ad.read_h5ad("module_ntc_pseudobulk.h5ad")
au=set(pd.read_parquet("authors_de_module.parquet").gene_id.unique())
A.var["gene_ids"]=A.var["gene_ids"].astype(str)
A=A[:,A.var["gene_ids"].isin(au).values].copy()
ncpu=int(os.environ.get("NC","6"))
ms=MultistatePerturbSeqDataset(A,sample_cols=["donor_id"],perturbation_type="CRISPRi",
    target_col="perturbed_gene_id",sgrna_col="guide_id",state_col="culture_condition",control_level="NTC")
t0=time.time()
res=ms.run_target_DE(design_formula="~ log10_n_cells + donor_id + target",test_state=["Rest"],
    test_targets=None,min_counts_per_gene=10,return_model=False,n_cpus=ncpu)
print(f"REST_TIME={time.time()-t0:.0f}s rows={len(res)}",flush=True)
res.to_csv("de_rest_only.csv",index=False)
PY
NC=6 timeout 700 python time_one.py 2>&1 | grep -vE "cpu_info|InitGoogle|sysctl" | tail -15; echo "EXIT ${PIPESTATUS[0]}"


# ======================================================================
# cell 63  [python]
# ======================================================================
import anndata as ad
A0 = ad.read_h5ad("module_ntc_pseudobulk.h5ad")  # for n samples / n guides per target x cond
FDR=0.10

rec = pd.read_csv("de_module_recomputed.csv")
rec["sig"] = rec.padj < FDR
# exclude on-target gene from trans-DE counts (like authors' n_downstream)
rec["is_ontarget"] = rec.target == rec.gene_name

rows=[]
for (t,c), g in rec.groupby(["target","condition"]):
    trans = g[~g.is_ontarget]
    sig = trans[trans.sig]
    n_up=int((sig.log2FC>0).sum()); n_dn=int((sig.log2FC<0).sum())
    ot = g[g.is_ontarget]
    ot_lfc = float(ot.log2FC.iloc[0]) if len(ot) else np.nan
    ot_padj= float(ot.padj.iloc[0]) if len(ot) else np.nan
    ot_se  = float(ot.lfcSE.iloc[0]) if len(ot) else np.nan
    # sample & guide counts from checkpoint (keep_for_DE only, matching what DE used effectively)
    sub = A0.obs[(A0.obs.perturbed_gene_name==t)&(A0.obs.culture_condition==c)]
    n_samp=int(len(sub)); n_guide=int(sub.guide_id.nunique())
    rows.append(dict(target=t, condition=c,
        n_DE_up=n_up, n_DE_down=n_dn, n_DE_total=n_up+n_dn,
        ontarget_log2FC=round(ot_lfc,3), ontarget_lfcSE=round(ot_se,3), ontarget_padj=ot_padj,
        ontarget_significant=bool(ot_padj<FDR) if np.isfinite(ot_padj) else False,
        n_pseudobulk_samples=n_samp, n_guides=n_guide,
        crosscheck_pearson_r=round(rlu.get((t,c),np.nan),3)))
summ = pd.DataFrame(rows)
# role & note columns
role = {"SETDB1":"writer (module core)","TASOR":"HUSH","PPHLN1":"HUSH","TASOR2":"HUSH-adjacent",
        "ATF7IP":"module core","TRIM28":"module core","SUV39H1":"specificity control","SUV39H2":"specificity control"}
summ["role"]=summ.target.map(role)
def note(r):
    n=[]
    if not r.ontarget_significant: n.append("on-target KD not significant at 10% FDR")
    elif abs(r.ontarget_log2FC)<0.5: n.append("weak on-target KD (|log2FC|<0.5)")
    if r.crosscheck_pearson_r<0.9: n.append("crosscheck r<0.9 (few real DE genes; logFC dominated by noise)")
    if r.target in ("SUV39H1","SUV39H2"): n.ap…[+716 chars]


# ======================================================================
# cell 209  [python]
# ======================================================================
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
A=ad.read_h5ad("module_ntc_pseudobulk.h5ad")
print("pseudobulk:",A.shape)
# label KD vs NTC + module target
obs=A.obs
print("perturbed_gene_id sample:",obs["perturbed_gene_id"].value_counts().head(12).to_dict())
print("conditions:",obs["culture_condition"].value_counts().to_dict() if "culture_condition" in obs else obs.columns.tolist())

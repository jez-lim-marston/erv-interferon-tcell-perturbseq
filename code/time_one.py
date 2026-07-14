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

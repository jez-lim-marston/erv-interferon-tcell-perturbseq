#!/usr/bin/env python
"""Reproduce authors' pseudobulk DESeq2 for the silencing-module targets + NTC,
across Rest/Stim8hr/Stim48hr, using their exact wrapper and design formula.
Writes de_module_recomputed.csv incrementally."""
import os, sys, time, warnings, json
# BLAS single-threaded so joblib threading parallelizes cleanly (else oversubscription thrash)
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
warnings.filterwarnings("ignore")

# --- sandbox numba fix: disable disk cache (no writable locator) ---
import numba
_o = numba.njit
numba.njit = lambda *a, **k: (k.__setitem__("cache", False) or _o(*a, **k))
numba.core.decorators.njit = numba.njit

# --- sandbox joblib fix: loky can't spawn (SC_SEM_NSEMS_MAX blocked) ---
from joblib import register_parallel_backend
from joblib._parallel_backends import ThreadingBackend
class ThreadingInnerOK(ThreadingBackend):
    supports_inner_max_num_threads = True
    def configure(self, n_jobs=1, parallel=None, **backend_args):
        backend_args.pop("inner_max_num_threads", None)
        return super().configure(n_jobs=n_jobs, parallel=parallel, **backend_args)
register_parallel_backend("threading_inner_ok", ThreadingInnerOK)

# Patch DefaultInference.__init__ ON THE CLASS so every instantiation
# (pertpy PyDESeq2.fit AND pydeseq2.ds.DeseqStats) uses the threading backend.
from pydeseq2.default_inference import DefaultInference
_orig_init = DefaultInference.__init__
def _init(self, *a, **k):
    k["backend"] = "threading_inner_ok"
    _orig_init(self, *a, **k)
DefaultInference.__init__ = _init

import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, ".")
from MultiStatePerturbSeqDataset import MultistatePerturbSeqDataset

N_CPUS = int(os.environ.get("DE_NCPUS", "8"))
DESIGN = "~ log10_n_cells + donor_id + target"
CONDS  = ["Rest", "Stim8hr", "Stim48hr"]

print("Loading checkpoint + author gene set...", flush=True)
A = ad.read_h5ad("module_ntc_pseudobulk.h5ad")
authors = pd.read_parquet("authors_de_module.parquet")
author_gene_ids = set(pd.unique(authors.gene_id))
A.var["gene_ids"] = A.var["gene_ids"].astype(str)
A = A[:, A.var["gene_ids"].isin(author_gene_ids).values].copy()   # 10,282 tested genes
print(f"  subset to author test-gene universe: {A.shape}", flush=True)

id2name = dict(zip(A.var["gene_ids"], A.var["gene_name"]))
all_res = []
done_conds = []
outcsv = "de_module_recomputed.csv"

for cond in CONDS:
    t0 = time.time()
    print(f"[{cond}] building wrapper + fitting DESeq2 ...", flush=True)
    ms = MultistatePerturbSeqDataset(
        A, sample_cols=["donor_id"], perturbation_type="CRISPRi",
        target_col="perturbed_gene_id", sgrna_col="guide_id",
        state_col="culture_condition", control_level="NTC",
    )
    res = ms.run_target_DE(
        design_formula=DESIGN, test_state=[cond], test_targets=None,
        min_counts_per_gene=10, return_model=False, n_cpus=N_CPUS,
    )
    res["condition"] = cond
    all_res.append(res)
    done_conds.append(cond)
    # incremental save (partial ok)
    pd.concat(all_res, ignore_index=True).to_csv(outcsv + ".partial", index=False)
    print(f"[{cond}] done in {time.time()-t0:.0f}s | rows={len(res):,}", flush=True)

de = pd.concat(all_res, ignore_index=True)
# 'contrast' column carries the target gene_id; map to gene symbol; add gene symbol for measured genes
de = de.rename(columns={"contrast": "target_gene_id"})
de["target"] = de["target_gene_id"].map(id2name)
if "variable" in de.columns:
    de["gene_id"] = de["variable"]
    de["gene_name"] = de["variable"].map(id2name)
de.to_csv(outcsv, index=False)
print("WROTE", outcsv, "shape", de.shape, flush=True)
print("columns:", de.columns.tolist(), flush=True)
print(de.head(3).to_string(), flush=True)

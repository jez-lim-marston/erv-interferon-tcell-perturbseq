"""
02 Authors Anchor

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 22  [python]
# ======================================================================
# h2 = authors' DE_stats.h5ad (already open). obs = target x condition rows.
# Build obs frame with the columns we need
o2cols = ["target_contrast_gene_name","culture_condition","target_contrast","n_cells_target",
          "n_up_genes","n_down_genes","n_total_de_genes","ontarget_effect_size","ontarget_significant",
          "target_baseMean","n_downstream","n_guides","donor_correlation_all_mean",
          "guide_correlation_all","low_target_gex"]
def read_any(grp,col):
    node=grp[col]
    if isinstance(node,h5py.Group) and "categories" in node:
        return read_cat(grp,col)
    arr=node[:]
    if arr.dtype.kind=="S": arr=np.array([x.decode() for x in arr])
    return arr
o2 = pd.DataFrame({c: read_any(h2["obs"],c) for c in o2cols})
mods=["SETDB1","TRIM28","ATF7IP","TASOR","PPHLN1","MPHOSPH8","TASOR2","SUV39H1","SUV39H2"]
mask2 = o2.target_contrast_gene_name.isin(mods).values
rows2 = np.where(mask2)[0]
print("author DE-stats module rows:", len(rows2))
print(o2.loc[mask2, ["target_contrast_gene_name","culture_condition","n_total_de_genes","ontarget_effect_size","ontarget_significant","target_baseMean"]].to_string())


# ======================================================================
# cell 23  [python]
# ======================================================================
# Extract layer values for the 24 module rows, all 10282 genes
gene_names = read_any(h2["var"],"gene_name")
gene_ids = read_any(h2["var"],"gene_ids")
layers = {}
for L in ["log_fc","lfcSE","p_value","adj_p_value","baseMean","zscore"]:
    layers[L] = h2["layers"][L][rows2, :]   # (24, 10282)
    print(L, layers[L].shape)

# long-form dataframe
recs=[]
tgt = o2.target_contrast_gene_name.values[rows2]
cnd = o2.culture_condition.values[rows2]
for i in range(len(rows2)):
    recs.append(pd.DataFrame({
        "target": tgt[i], "condition": cnd[i],
        "gene_name": gene_names, "gene_id": gene_ids,
        "authors_log_fc": layers["log_fc"][i],
        "authors_lfcSE": layers["lfcSE"][i],
        "authors_p_value": layers["p_value"][i],
        "authors_adj_p_value": layers["adj_p_value"][i],
        "authors_baseMean": layers["baseMean"][i],
        "authors_zscore": layers["zscore"][i],
    }))
authors_de = pd.concat(recs, ignore_index=True)
print("\nlong-form authors DE:", authors_de.shape)
print("NaN log_fc fraction:", authors_de.authors_log_fc.isna().mean().round(3))
authors_de.to_parquet("authors_de_module.parquet", index=False)
print("saved authors_de_module.parquet")
print(authors_de.head(3).to_string())


# ======================================================================
# cell 24  [python]
# ======================================================================
for col in ["gene_name","gene_id"]:
    authors_de[col] = authors_de[col].map(lambda x: x.decode() if isinstance(x,(bytes,bytes)) else x)
authors_de.to_parquet("authors_de_module.parquet", index=False)
# per-target measured-gene count (non-NaN)
print("author DE genes per target/condition (all should be 10282):",
      authors_de.groupby(["target","condition"]).size().unique())
print("\ndecoded head:")
print(authors_de.head(3)[["target","condition","gene_name","gene_id","authors_log_fc"]].to_string())
# sanity: on-target self logFC should be strongly negative
for t in ["SETDB1","TASOR","SUV39H1"]:
    r=authors_de[(authors_de.target==t)&(authors_de.gene_name==t)]
    print(f"  {t} on-target logFC by cond:", dict(zip(r.condition, r.authors_log_fc.round(2))))


# ======================================================================
# cell 33  [python]
# ======================================================================
import time
# map authors' gene set (by gene_id) onto pseudobulk var
authors_de = pd.read_parquet("authors_de_module.parquet")
author_gene_ids = pd.unique(authors_de.gene_id)   # 10282
print("author test genes:", len(author_gene_ids))
# pseudobulk var index is gene symbol? check
print("A.var head:\n", A.var.head(2))
print("A.var_names sample:", list(A.var_names[:3]))
# match by gene_ids column
A.var["gene_ids"]=A.var["gene_ids"].astype(str)
in_author = A.var["gene_ids"].isin(set(author_gene_ids))
print("pseudobulk genes matching author set:", int(in_author.sum()), "/", len(author_gene_ids))


# ======================================================================
# cell 59  [python]
# ======================================================================
import pandas as pd, numpy as np
from scipy.stats import pearsonr

rec = pd.read_csv("de_module_recomputed.csv")   # target, condition, gene_id, log2FC ...
auth = pd.read_parquet("authors_de_module.parquet")  # target, condition, gene_id, authors_log_fc ...

# merge on target + condition + gene_id
m = rec.merge(auth[["target","condition","gene_id","authors_log_fc","authors_lfcSE","authors_adj_p_value"]],
              on=["target","condition","gene_id"], how="inner")
print("merged rows:", len(m), "| targets:", sorted(m.target.unique()))
# authors excluded MPHOSPH8 -> should have 8 targets in merge
print("targets in merge (authors tested):", sorted(m.target.unique()))
print("MPHOSPH8 in authors anchor?", "MPHOSPH8" in auth.target.unique())

# per (target, condition) Pearson r over genes (both logFC finite)
rows=[]
for (t,c), g in m.groupby(["target","condition"]):
    x=g.log2FC.values; y=g.authors_log_fc.values
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()>10:
        r,p = pearsonr(x[ok],y[ok])
        rows.append({"target":t,"condition":c,"pearson_r":r,"n_genes":int(ok.sum())})
rdf=pd.DataFrame(rows)
print("\nper-target x condition Pearson r:")
print(rdf.pivot(index="target",columns="condition",values="pearson_r").round(3).to_string())


# ======================================================================
# cell 135  [python]
# ======================================================================
import pandas as pd
# is KDM4C (or any KDM4/eraser, HP1/CBX reader) in the recomputed DE? The screen targeted the silencing module;
# check which druggable-axis nodes were actually perturbed & measured.
de = pd.read_csv("de_module_recomputed.csv")
print("targets in recomputed DE:", sorted(de.target.unique()))
# check the AUTHORS' full DE-stats for eraser/reader nodes (they screened genome-scale)
auth = pd.read_parquet("authors_de_module.parquet")
print("targets in authors_de_module subset:", sorted(auth.target.unique()))
# these are only the module subset we pulled. The genome-scale screen has ALL targets.
# check node genes as MEASURED genes (they're in the 10,282 tested universe regardless)
nodes = ["SETDB1","TRIM28","ATF7IP","TASOR","PPHLN1","MPHOSPH8","TASOR2",  # writer/HUSH
         "KDM4A","KDM4B","KDM4C","KDM4D","KDM4E",  # eraser
         "CBX1","CBX3","CBX5"]  # reader (HP1 beta/gamma/alpha)
present_as_gene = {n: n in set(de.gene_name) for n in nodes}
print("\nnode present as MEASURED gene in universe:", present_as_gene)


# ======================================================================
# cell 186  [python]
# ======================================================================
from scipy.stats import pearsonr
auth=pd.read_parquet("authors_de_module.parquet")
mine=pd.read_csv("de_module_recomputed.csv")
print("authors cols:",list(auth.columns)[:12])
print("mine cols:",list(mine.columns))

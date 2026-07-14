"""
10 Druggable Nodes

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 145  [python]
# ======================================================================
# save the antagonism evidence + node ISG/de-rep direction table for the ranking
antag_rows=[]
for c in ["Rest","Stim8hr","Stim48hr"]:
    su=de[(de.target=="SETDB1")&(de.condition==c)&(de.padj<0.10)&(de.log2FC>0.5)]
    ids=set(su.gene_id)&erv_prox_ids; cc=cols_for(ids)
    s=prof("SETDB1",c)
    for n in ["KDM4A","KDM4B","KDM4C","CBX1","CBX3","CBX5","TASOR","PPHLN1","ATF7IP","TRIM28","TASOR2"]:
        if (n,c) in node_rows and len(cc)>10:
            v=prof(n,c)[cc]; vf=v[np.isfinite(v)]
            sm=np.isfinite(s[cc])&np.isfinite(v)
            r,pr=pearsonr(s[cc][sm], v[sm])
            antag_rows.append(dict(node=n,condition=c,n_genes=len(cc),
                median_log2FC=float(np.median(vf)), frac_down=float(np.mean(vf<0)),
                wilcoxon_p=float(wilcoxon(vf)[1]), anticorr_vs_SETDB1_r=float(r), anticorr_p=float(pr)))
antag=pd.DataFrame(antag_rows)
antag.to_csv("annot/node_antagonism.csv", index=False)
print("KDM4C rows:")
print(antag[antag.node=="KDM4C"].round(4).to_string(index=False))
print("\nAll erasers, Stim48hr (antagonism ranking, most-negative = strongest):")
print(antag[(antag.condition=="Stim48hr")&(antag.node.str.startswith("KDM"))].sort_values("median_log2FC")[["node","median_log2FC","frac_down","wilcoxon_p","anticorr_vs_SETDB1_r"]].round(4).to_string(index=False))


# ======================================================================
# cell 147  [python]
# ======================================================================
import json
nodes = json.load(open("handoff/nodes_query.json"))
allnodes = [g for lst in nodes.values() for g in lst]

targets = {}
for g in allnodes:
    try:
        r = host.mcp("chembl","target_search", gene_symbol=g, organism="Homo sapiens", target_type="SINGLE PROTEIN", limit=10)
        ts = r.get("targets",[])
        # pick the human single-protein whose components match the gene
        pick=None
        for t in ts:
            if t.get("organism")=="Homo sapiens":
                pick=t; break
        if pick is None and ts: pick=ts[0]
        targets[g] = {"target_chembl_id": pick.get("target_chembl_id") if pick else None,
                      "pref_name": pick.get("pref_name") if pick else None,
                      "n_hits": len(ts)}
    except Exception as e:
        targets[g] = {"error": str(e)[:120]}
    print(g, "->", targets[g])
json.dump(targets, open("handoff/chembl_targets.json","w"))


# ======================================================================
# cell 148  [python]
# ======================================================================
import json
targets = json.load(open("handoff/chembl_targets.json"))
bioact = {}
for g,info in targets.items():
    tid = info.get("target_chembl_id")
    if not tid:
        bioact[g] = {"target_chembl_id": None, "n_activities": 0, "n_potent_100nM": 0, "best_pchembl": None, "n_compounds": 0}
        print(g, "-> no target"); continue
    try:
        # potent activities (sub-100nM = pchembl>=7), binding+functional
        r = host.mcp("chembl","get_bioactivity", target_chembl_id=tid, min_pchembl=7, limit=1000)
        acts = r.get("activities",[])
        total = r.get("total", r.get("count", len(acts)))
        mols = set(a.get("molecule_chembl_id") for a in acts if a.get("molecule_chembl_id"))
        pchembls = [a.get("pchembl_value") for a in acts if a.get("pchembl_value") is not None]
        pchembls = [float(x) for x in pchembls]
        # also grab overall activity count (any potency)
        r_all = host.mcp("chembl","get_bioactivity", target_chembl_id=tid, min_pchembl=5, limit=1000)
        allacts=r_all.get("activities",[]); allmols=set(a.get("molecule_chembl_id") for a in allacts if a.get("molecule_chembl_id"))
        bioact[g] = {"target_chembl_id": tid, "pref_name": info.get("pref_name"),
                     "n_potent_activities_sub100nM": len(acts),
                     "n_compounds_sub100nM": len(mols),
                     "best_pchembl": max(pchembls) if pchembls else None,
                     "n_compounds_uM_or_better": len(allmols)}
    except Exception as e:
        bioact[g] = {"target_chembl_id": tid, "error": str(e)[:150]}
    print(g, "->", bioact[g])
json.dump(bioact, open("handoff/chembl_bioact.json","w"))


# ======================================================================
# cell 149  [python]
# ======================================================================
import json
targets = json.load(open("handoff/chembl_targets.json"))
# For each node with chemistry, find the max clinical phase across its potent compounds + any mechanism/approved drug.
clinical = {}
for g in ["SETDB1","TRIM28","KDM4A","KDM4B","KDM4C"]:
    tid = targets[g]["target_chembl_id"]
    maxphase=None; approved=[]
    try:
        # mechanism = curated drug-target (clinical-stage indicator)
        mech = host.mcp("chembl","get_mechanism", target_chembl_id=tid, limit=50)
        mechs = mech.get("mechanisms",[])
        # bioactivity molecules, then check their max_phase via compound_search
        r = host.mcp("chembl","get_bioactivity", target_chembl_id=tid, min_pchembl=6, limit=1000)
        mols = list({a.get("molecule_chembl_id") for a in r.get("activities",[]) if a.get("molecule_chembl_id")})
        phases=[]
        # sample up to 40 molecules for phase (avoid excessive calls)
        for mid in mols[:40]:
            cs = host.mcp("chembl","compound_search", name=mid, chembl_id=mid, limit=1)
            for c in cs.get("compounds",[]):
                mp=c.get("max_phase")
                if mp is not None: phases.append(mp)
        clinical[g]={"n_mechanisms_curated":len(mechs),
                     "max_phase_among_potent": max(phases) if phases else None,
                     "n_molecules_checked": len(mols[:40]),
                     "any_approved_or_clinical": any((p or 0)>=1 for p in phases) or len(mechs)>0}
    except Exception as e:
        clinical[g]={"error":str(e)[:150]}
    print(g, "->", clinical[g])
json.dump(clinical, open("handoff/chembl_clinical.json","w"))


# ======================================================================
# cell 155  [python]
# ======================================================================
# Transparent composite scoring (documented, not black-box)
# tractability_score (0-1): normalized log10(compounds sub-100nM +1) scaled + potency component
c=cand.copy()
c["log_potent"]=np.log10(c.n_compounds_sub100nM.fillna(0)+1)
c["tract_chem"]=c.log_potent/c.log_potent.max() if c.log_potent.max()>0 else 0
c["tract_potency"]=(c.best_pchembl.fillna(0)/9.5).clip(0,1)  # 9.5 ~ ceiling
c["tractability_score"]=(0.6*c["tract_chem"]+0.4*c["tract_potency"]).round(3)

# mechanistic evidence score: for erasers, antagonism strength (negative anticorr = good); writers = on-target de-repression role
# Use |anticorr| when negative (antagonistic) for erasers/readers; writers scored by module-core role (SETDB1 anchor)
def mech_score(r):
    if r.axis=="writer":
        # writer nodes ARE the de-repression drivers; SETDB1 strongest. Score by demonstrated DE (from de_summary) role
        return {"SETDB1":1.0,"TASOR":0.8,"PPHLN1":0.75,"ATF7IP":0.5,"TRIM28":0.4}.get(r.node,0.3)
    else:
        # eraser/reader: antagonism = negative shift on SETDB1-derepressed ERV-prox genes
        rr=r.antag_anticorr_r_stim48
        return float(np.clip(-rr/0.25,0,1)) if rr<0 else 0.0
c["mechanism_score"]=c.apply(mech_score,axis=1).round(3)

# directly-measured bonus: all screened, so uniform; note it
c["directly_measured"]=True
# composite = tractability x mechanism (both must be present to be a lead)
c["composite_score"]=(0.5*c.tractability_score+0.5*c.mechanism_score).round(3)
c=c.sort_values("composite_score",ascending=False).reset_index(drop=True)
c["rank"]=c.index+1
print(c[["rank","node","axis","tractability_score","mechanism_score","composite_score",
         "n_compounds_sub100nM","best_pchembl","antag_anticorr_r_stim48","max_clinical_phase"]].to_string(index=False))


# ======================================================================
# cell 156  [python]
# ======================================================================
# part3_ranked_candidates.csv — the ranked node table with evidence + honesty flags
ranked=c[["rank","node","axis","role","directly_measured","composite_score","tractability_score","mechanism_score",
          "has_chembl_target","n_compounds_uM_or_better","n_compounds_sub100nM","best_pchembl",
          "max_clinical_phase","n_curated_mechanisms",
          "antag_median_log2FC_stim48","antag_wilcoxon_p_stim48","antag_anticorr_r_stim48"]].copy()
ranked["clinical_stage_compound"]="none (preclinical only)"
ranked["is_lead_node"]=ranked.node=="KDM4C"
ranked.to_csv("part3_ranked_candidates.csv", index=False)
print("saved part3_ranked_candidates.csv:", ranked.shape)

# part3_chembl_tractability.csv — the raw ChEMBL evidence per node
tract=c[["node","axis","role","chembl_target_id","has_chembl_target",
         "n_compounds_uM_or_better","n_compounds_sub100nM","best_pchembl",
         "max_clinical_phase","n_curated_mechanisms"]].copy()
tract["best_potency_nM"]=(10**(9-tract.best_pchembl)).round(2)  # pchembl->nM
tract["clinical_status"]="preclinical (max_phase = null/none)"
tract["tractability_call"]=np.where(tract.n_compounds_sub100nM>=50,"high (rich sub-100nM chemistry)",
                            np.where(tract.n_compounds_sub100nM>=1,"low (isolated tool compounds)",
                            np.where(tract.has_chembl_target,"minimal (few weak binders)","none (no ChEMBL target)")))
tract.to_csv("part3_chembl_tractability.csv", index=False)
print("saved part3_chembl_tractability.csv:", tract.shape)
print("\n=== Tractability by axis ===")
print(tract[["node","axis","n_compounds_sub100nM","best_potency_nM","tractability_call","clinical_status"]].to_string(index=False))


# ======================================================================
# cell 165  [python]
# ======================================================================
import json
targets=json.load(open("handoff/chembl_targets.json"))
# For each node target, pull the FULL potent compound set and check max_phase on every molecule.
clinical_full={}
for g in ["SETDB1","TRIM28","PPHLN1","KDM4A","KDM4B","KDM4C","CBX1","CBX3","CBX5"]:
    tid=targets[g].get("target_chembl_id")
    if not tid:
        clinical_full[g]={"max_phase_seen":None,"n_checked":0}; continue
    try:
        # all activities uM-or-better -> unique molecules
        r=host.mcp("chembl","get_bioactivity",target_chembl_id=tid,min_pchembl=5,limit=1000)
        mols=sorted({a.get("molecule_chembl_id") for a in r.get("activities",[]) if a.get("molecule_chembl_id")})
        phases=[]
        # check every molecule's max_phase
        for mid in mols:
            cs=host.mcp("chembl","compound_search",name=mid,chembl_id=mid,limit=1)
            for c in cs.get("compounds",[]):
                mp=c.get("max_phase")
                phases.append(mp if mp is not None else 0)
        mx=max(phases) if phases else None
        n_clin=sum(1 for p in phases if (p or 0)>=1)
        clinical_full[g]={"n_checked":len(mols),"max_phase_seen":mx,"n_clinical_stage":n_clin}
    except Exception as e:
        clinical_full[g]={"error":str(e)[:120]}
    print(g,"->",clinical_full[g])
json.dump(clinical_full,open("handoff/chembl_clinical_full.json","w"))


# ======================================================================
# cell 166  [python]
# ======================================================================
import json
targets=json.load(open("handoff/chembl_targets.json"))
detail={}
for g in ["KDM4A","KDM4B","KDM4C","CBX3"]:
    tid=targets[g]["target_chembl_id"]
    r=host.mcp("chembl","get_bioactivity",target_chembl_id=tid,min_pchembl=5,limit=1000)
    # map molecule -> best pchembl against this target
    from collections import defaultdict
    molpc=defaultdict(list)
    for a in r.get("activities",[]):
        mid=a.get("molecule_chembl_id"); pc=a.get("pchembl_value")
        if mid and pc is not None: molpc[mid].append(float(pc))
    mols=sorted(molpc)
    clin=[]
    for mid in mols:
        cs=host.mcp("chembl","compound_search",name=mid,chembl_id=mid,limit=1)
        for c in cs.get("compounds",[]):
            mp=c.get("max_phase")
            if mp is not None and mp>=1:
                clin.append({"id":mid,"name":c.get("pref_name"),"max_phase":mp,
                             "best_pchembl_vs_target":max(molpc[mid])})
    # sort by potency
    clin=sorted(clin,key=lambda x:-x["best_pchembl_vs_target"])
    detail[g]=clin
    print(f"\n{g}: {len(clin)} clinical-stage compounds w/ recorded activity")
    for c in clin[:6]:
        print(f"   {c['name'] or c['id']} (phase {c['max_phase']}): pChEMBL vs {g} = {c['best_pchembl_vs_target']:.1f} ({10**(9-c['best_pchembl_vs_target']):.0f} nM)")
json.dump(detail,open("handoff/chembl_clinical_detail.json","w"))


# ======================================================================
# cell 167  [python]
# ======================================================================
import json, pandas as pd, numpy as np
cf=json.load(open("handoff/chembl_clinical_full.json"))
cd=json.load(open("handoff/chembl_clinical_detail.json"))

# corrected clinical annotation per node
# mechanism-directed clinical compound = zavondemstat (KDM4 inhibitor). Others are off-target/promiscuous.
mechanism_directed={"KDM4A":"zavondemstat (Ph2, KDM4 inhibitor)","KDM4B":"zavondemstat (Ph2, KDM4 inhibitor)",
                    "KDM4C":"zavondemstat (Ph2, KDM4 inhibitor)"}
offtarget_only={"CBX3":"molibresib (BET inhibitor, off-target 6 uM)"}

ranked=pd.read_csv("part3_ranked_candidates.csv")
tract=pd.read_csv("part3_chembl_tractability.csv")

def clin_status(node):
    mx=cf.get(node,{}).get("max_phase_seen")
    if node in mechanism_directed:
        return f"clinical-stage mechanism inhibitor exists: {mechanism_directed[node]}"
    if mx is not None and mx>=1:
        # only off-target/promiscuous approved compounds
        return f"no mechanism-directed clinical compound; only off-target hits (max_phase={mx:g}, uM/promiscuous)"
    return "preclinical only (no clinical-stage compound with recorded activity)"

for df in (ranked,tract):
    df["max_clinical_phase"]=df.node.map(lambda n: cf.get(n,{}).get("max_phase_seen"))
    df["clinical_status_corrected"]=df.node.map(clin_status)
# drop stale columns that asserted 'none'
if "clinical_stage_compound" in ranked: ranked=ranked.drop(columns=["clinical_stage_compound"])
if "clinical_status" in tract: tract=tract.drop(columns=["clinical_status"])
ranked.to_csv("part3_ranked_candidates.csv",index=False)
tract.to_csv("part3_chembl_tractability.csv",index=False)
print(tract[["node","axis","n_compounds_sub100nM","max_clinical_phase","clinical_status_corrected"]].to_string(index=False))

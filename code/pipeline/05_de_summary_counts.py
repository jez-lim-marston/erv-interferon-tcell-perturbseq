"""
05 De Summary Counts

As-run analysis code, exported verbatim (per-cell) from the Claude Science execution
log for this project. Executed interactively in conda env `scte` — python 3.11.15;
scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3, formulaic-contrasts
1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs. Public data are read by partial-range
requests from s3://genome-scale-tcell-perturb-seq/marson2025_data/ (no full download).
This is a faithful record of what produced the named artifacts under results/, with
cell boundaries preserved — not a refactored package.
"""


# ======================================================================
# cell 65  [python]
# ======================================================================
summ.loc[summ.target=="MPHOSPH8","role"]="HUSH-adjacent (excluded by authors)"
summ.loc[summ.target=="MPHOSPH8","notes"]=(
    "on-target KD not significant at 10% FDR; absent from authors' published DE-stats "
    "(failed guide-effectiveness filter) — no crosscheck possible")
summ.to_csv("de_summary_counts.csv", index=False)
print("Saved de_summary_counts.csv:", summ.shape)
# headline contrast for the arc
core = summ[summ.target.isin(["SETDB1","TASOR","PPHLN1"])].n_DE_total
ctrl = summ[summ.target.isin(["SUV39H1","SUV39H2"])].n_DE_total
print(f"\nModule core (SETDB1/TASOR/PPHLN1) DE genes: median={int(core.median())}, range {core.min()}-{core.max()}")
print(f"SUV39 specificity control DE genes:        median={int(ctrl.median())}, range {ctrl.min()}-{ctrl.max()}")
print(f"Contrast ratio (median core / median control): {core.median()/max(ctrl.median(),1):.0f}x")

# Correction: Th1/Th2 module figure

**What changed.** The figure `results/fig_module_th1th2_shift.png` was removed and replaced by
`results/figures/fig_module_th1th2_null.{png,svg}` (light + dark). Only the framing/title
changed — the underlying data, heatmap cells, module-score bars, and significance stars are
identical (recomputed from the same sources, `results/de_crosscheck_geneset.csv` and
`results/genesetscore_summary.csv`; module scores reproduce exactly).

**Why.** The old title read *"Core H3K9me3/ERV-silencing knockdowns shift the Th1/Th2 module
in the predicted direction under stimulation; SUV39H2 shifts oppositely."* That claim is
contradicted by the figure's own significance annotation and by the project's central finding
(the Th1/Th2 hypothesis scores null; the real consequence is an ERV-gated interferon program).

**What the data actually support** (two-sided permutation p on the module score = mean(Th1
log2FC) − mean(Th2 log2FC)):

| Target | Role | Rest | Stim-8h | Stim-48h | Significant? |
|---|---|---|---|---|---|
| SETDB1 | core | 0.119 | 0.053 | 0.213 | no |
| ATF7IP | core | 0.754 | 0.520 | 0.354 | no |
| PPHLN1 | core | 0.157 | 0.541 | 0.182 | no |
| TASOR | core | 0.490 | 0.298 | 0.302 | no |
| TRIM28 | core (pleiotropic) | 0.302 | **0.003** | **0.006** | yes |
| SUV39H1 | off-ERV comparator | 0.655 | 0.118 | 0.646 | no |
| SUV39H2 | off-ERV comparator | 0.679 | **0.003** | **0.006** | yes (opposite direction) |

The four core ERV-silencing enzymes with a clean mechanistic identity (SETDB1, ATF7IP, PPHLN1,
TASOR) are **n.s. in every condition**. The only significant movers are (i) **TRIM28**, which is
pleiotropic (KAP1 has many H3K9me3-independent roles) and cannot be read as clean support for
the ERV→Th1/Th2 hypothesis, and (ii) **SUV39H2**, the off-ERV *negative control*, which moves in
the *opposite* direction. The honest read: **the Th1/Th2 signature scores null for the core
module.** This is the null result that motivates the data-led pivot to the ERV-gated interferon
program (see `README.md`, `docs/SUBMISSION_SUMMARY.md`).

The new title states this directly: *"Th1/Th2 signature scores NULL for the core ERV-silencing
enzymes … only pleiotropic TRIM28 and the off-ERV control SUV39H2 move Th2, in opposite
directions."*

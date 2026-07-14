# Pre-registration — ERV-cis-proximity enrichment

**Written before any enrichment statistic is computed.** Part of the TE-aware
reanalysis of the Marson/Pritchard CD4⁺ T-cell CRISPRi Perturb-seq. All
thresholds below are fixed here and not revised after seeing enrichment results.

## Hypothesis (pre-declared)

Genes **de-repressed** (up-regulated) by knockdown of the SETDB1/HUSH silencing
module are **enriched near silenced (H3K9me3-marked) ERVs** relative to the
tested-gene background. The raw enrichment is expected to be **inflated (~5×)**
because KRAB-ZNF genes are simultaneously (a) H3K9me3/ERV-proximal, (b)
genomically clustered (chr19 KZNF clusters), and (c) disproportionately
de-repressed. After accounting for the KRAB-ZNF confounder the enrichment is
expected to **collapse to an honest ~2× (odds ratio ≈ 1.7–2.2)** and remain
significant. The **SUV39H1/SUV39H2** methyltransferase knockdowns (a
SETDB1-independent, pericentromeric/non-ERV H3K9 axis) are the **specificity
control** and are expected to show **no enrichment (OR ≈ 1, CI spanning 1)**.

This is a *demonstration of an axis*, not a discovery claim; the assay is
probe-based 10x Flex (TE-blind), so enrichment is measured as **cis-genomic
proximity of protein-coding genes to ERVs**, never as direct TE quantification.

## Frozen definitions & thresholds

| Choice | Pre-registered value |
|--------|---------------------|
| **De-repressed gene** | `padj < 0.10` **AND** `log2FC > +0.5` (up-regulated) in that target×condition, from `de_module_recomputed.csv`. On-target gene excluded. |
| **Direction** | Up-regulated only (silencing removed ⇒ de-repression). Down-regulated genes are not "de-repressed". |
| **Background universe** | The **10,282 authors'-tested genes** (the assay's actual measurable universe), NOT the whole genome. Enrichment is de-repressed vs the rest of this universe. |
| **Silenced ERV** | An ERV/LTR RepeatMasker element (`repClass == "LTR"`) whose interval **overlaps an H3K9me3 ChIP-seq peak** (≥1 bp). |
| **ERV classes** | LTR `repClass`; reported overall and by `repFamily` (ERV1, ERVK, ERVL, ERVL-MaLR, and other LTR). |
| **H3K9me3 peaks** | Union/consensus of 3 ENCODE primary-CD4⁺-T-cell H3K9me3 ChIP-seq experiments (ENCSR787WLV, ENCSR453GNY, ENCSR692ICP), GRCh38 narrowPeak. Consensus = peak present in ≥2 of 3 experiments (primary); union reported as sensitivity. |
| **Proximity windows** | **TSS ± 50 kb** and **TSS ± 100 kb** — BOTH pre-registered, reported side by side. TSS is strand-aware, canonical transcript, GENCODE v44. |
| **Genome build** | GRCh38 / hg38 throughout. |
| **Primary enrichment statistic** | Fisher exact **odds ratio** + **95% CI** (two-sided), on the 2×2 [de-repressed × near-silenced-ERV] table over the 10,282-gene universe. |
| **Confounder-adjusted statistic** | Logistic regression `near_silenced_ERV ~ derepressed + is_KRABZNF + log10_gene_length`; adjusted OR = `exp(β_derepressed)` + 95% CI. Also report Fisher OR with KRAB-ZNF genes **excluded** as a complementary correction. |
| **KRAB-ZNF definition** | Curated KRAB-domain zinc-finger gene set mapped to Ensembl IDs (HGNC ZNF family filtered to KRAB-domain via InterPro IPR001909, cross-checked against the Imbeault/Trono KZNF catalogue). Frozen before enrichment. |
| **Targets tested** | Silencing module: SETDB1, TASOR, PPHLN1, ATF7IP, TRIM28, TASOR2. Specificity control: SUV39H1, SUV39H2. (MPHOSPH8 excluded — ineffective KD, per DE step.) |
| **Multiple testing** | OR + CI is the primary readout; p-values reported with BH-FDR across the target×window grid, flagged but not gating the effect-size interpretation. |
| **Per-condition handling** | Primary analysis pools de-repressed calls across the 3 activation conditions (a gene is "de-repressed" for a target if it passes threshold in ≥1 condition); per-condition results reported in the table as a robustness check. |

## Pre-declared success/interpretation rules

1. **Report BOTH raw and corrected OR** in the same table — the ~5×→~2× collapse
   is the headline, not something to bury. Failing to show the raw inflation
   would itself be a reporting failure.
2. The corrected enrichment is considered **supported** if the module-core
   targets (SETDB1, TASOR, PPHLN1) show adjusted OR CIs that **exclude 1** and
   lie in the ~1.5–2.5 range at ≥1 window.
3. The specificity control is considered **passed** if SUV39H1/H2 adjusted OR
   CIs **include 1** (no ERV enrichment), OR are clearly below the module-core
   effect. A SUV39 enrichment as strong as the module would falsify specificity.
4. Effect sizes + CIs are reported throughout; no bare p-values.
5. Window sensitivity (50 vs 100 kb) and peak-set choice (consensus vs union)
   are reported as robustness, not cherry-picked.

## Honest limitations (declared up front)

- H3K9me3 peaks come from **reference** primary CD4⁺ T cells (ENCODE), not the
  exact screen donors; H3K9me3 landscapes are largely stable but this is an
  approximation stated openly.
- "Silenced ERV" conflates the H3K9me3 mark with functional silencing; we do not
  measure ERV transcription (TE-blind assay).
- Proximity ≠ causation: cis-proximity is consistent with, but does not prove,
  direct ERV-mediated regulation of the nearby gene.

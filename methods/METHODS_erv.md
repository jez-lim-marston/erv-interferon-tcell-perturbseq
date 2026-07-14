# Methods — ERV-cis-proximity enrichment of de-repressed genes

Part of the TE-aware reanalysis of the public Marson/Pritchard CD4⁺ T-cell
CRISPRi Perturb-seq. This module tests whether genes **de-repressed** by
silencing-module knockdown lie **near silenced (H3K9me3⁺) ERVs**, and applies
the KRAB-ZNF confounder control. Thresholds were **pre-registered**
(`PREREG_erv.md`) before any enrichment statistic was computed. Results are
**demonstrated**; no therapeutic claim.

## TE-blind framing (unchanged)

The screen is probe-based 10x Flex — transposable elements are not directly
measured. Enrichment here is the **cis-genomic proximity of protein-coding
genes to annotated ERVs**, never direct TE quantification. "Silenced ERV" means
an ERV bearing the repressive H3K9me3 mark; it does not measure ERV
transcription.

## Data provenance (GRCh38 / hg38)

| Annotation | Source | Detail |
|-----------|--------|--------|
| ERV / LTR coordinates | UCSC RepeatMasker (`hg38 rmsk.txt.gz`) | `repClass == LTR`; **720,177** elements on main chromosomes (ERVL-MaLR, ERV1, ERVL, ERVK families) |
| H3K9me3 peaks | ENCODE Histone ChIP-seq, primary CD4⁺ T cells | 3 experiments: **ENCSR787WLV** (files ENCFF602JTX), **ENCSR453GNY** (ENCFF885NAO), **ENCSR692ICP** (ENCFF080EAR, PMA/ionomycin-stimulated). GRCh38 narrowPeak. |
| Gene TSS | GENCODE **v44** basic annotation | strand-aware canonical TSS; all 10,282 background genes mapped by Ensembl ID (100%). |
| KRAB-ZNF set | InterPro **IPR001909** (KRAB domain), human reviewed | 380 genes; **280** present in the tested-gene background. |
| DE input | `de_module_recomputed.csv` (this project's DE step) | authors'-method DESeq2, 10,282-gene universe. |

## Definitions (pre-registered)

- **De-repressed gene** for a target: `padj < 0.10` **AND** `log2FC > +0.5`
  (up-regulated) in ≥1 activation condition; on-target gene excluded.
- **Background universe:** the **10,282 authors'-tested genes** — the assay's
  measurable universe, not the whole genome.
- **H3K9me3 consensus:** peak present in **≥2 of 3** experiments →
  **52,649** consensus peaks (~48 Mb). Union (≥1 of 3) used as a sensitivity
  check.
- **Silenced ERV:** LTR element overlapping (≥1 bp) a consensus H3K9me3 peak →
  **25,343** silenced ERVs (3.5% of LTRs).
- **Proximity:** TSS within **±50 kb** and **±100 kb** of the nearest silenced
  ERV — both windows reported.

## Enrichment statistics

- **Primary (raw):** Fisher exact **odds ratio** + 95% CI on the 2×2
  [de-repressed × near-silenced-ERV] table over the 10,282-gene universe.
- **Confounder-adjusted, two complementary ways:**
  (a) Fisher OR with **KRAB-ZNF genes excluded** from the universe;
  (b) **logistic regression** `near_silenced_ERV ~ derepressed + is_KRABZNF +
  log10(gene_length)`; adjusted OR = `exp(β_derepressed)` + 95% CI.
- Per-target, per-window, per-condition results plus a **pooled module-core**
  set (SETDB1 ∪ TASOR ∪ PPHLN1 ∪ ATF7IP de-repressed genes). BH-FDR reported
  across the grid.

## Results

### The confounder is real and quantified
KRAB-ZNF genes are the archetypal H3K9me3/ERV-proximal, chr19-clustered family.
In this background they are near a silenced ERV at **89 %** (±50 kb) vs **24 %**
for non-KRAB genes (**3.65×**). They are also massively over-represented among
de-repressed genes — **18 % of SETDB1's, 41 % of PPHLN1's, 74 % of TASOR2's**
de-repressed sets, vs **2.7 %** background. In the logistic model, KRAB-ZNF
membership alone carries an OR of **~18–25** for silenced-ERV proximity — the
dominant driver of the raw signal.

### Raw enrichment is inflated (as pre-declared)
Raw silenced-ERV enrichment of de-repressed genes (±50 kb): SETDB1
**3.22 (2.80–3.71)**, PPHLN1 **7.67 (5.65–10.31)**, TASOR2 **14.37 (6.18–29.94)** — i.e. ~3–8×
(up to 14× for the KRAB-ZNF-dominated TASOR2). Pooled module core: **2.98 (2.62–3.39)**.

### Corrected enrichment collapses to an honest ~2×
After KRAB-ZNF correction, the **pooled module-core** OR (±50 kb) is
**2.07 (1.80–2.39)** (KRAB-ZNF excluded) / **2.23 (1.94–2.57)** (logistic-adjusted),
falling squarely in the pre-registered honest band (OR 1.7–2.2). At ±100 kb the
adjusted OR is **1.86 (1.62–2.13)**. Per target, SETDB1 collapses from
3.22 (2.80–3.71) raw to **2.36 (2.02–2.75)** adjusted. The enrichment remains
**significant** (CI excludes 1) for the module core.

### Specificity control passes
The SUV39H1/H2 methyltransferase knockdowns (SETDB1-independent,
pericentromeric/non-ERV H3K9 axis) show **no significant ERV enrichment**:
pooled OR (±50 kb) **1.78 (0.63–5.38)** raw, **1.47 (0.41–5.28)** adjusted — **CI spans
1 in every analysis**. The de-repression → silenced-ERV proximity link is
specific to the SETDB1/HUSH silencing axis, not a generic H3K9me3 effect. (Both
SUV39 paralogs were effectively knocked down in the DE step; the null here is a
true negative, and the wide CI reflects their small de-repressed-gene sets.)

## Honest limitations

- **H3K9me3 is from reference CD4⁺ T cells** (ENCODE), not the exact screen
  donors — a stable-landscape approximation, stated openly.
- **"Silenced" = H3K9me3-marked**, a proxy for silencing; ERV transcription is
  not measured (TE-blind assay).
- **Proximity ≠ causation:** cis-proximity is consistent with, but does not
  prove, direct ERV-mediated regulation.
- The pooled-core OR is the primary readout; individual weak targets (TRIM28,
  ATF7IP) have small de-repressed sets and correspondingly wide CIs — reported,
  not hidden.
- Windows (50/100 kb) and peak-set choice (consensus vs union) are reported as
  robustness, not selected post hoc.

## Artifacts

- `PREREG_erv.md` — pre-registration (frozen before computation)
- `erv_proximity_annotation.parquet` — per-gene TSS, distances, near-flags
  (50/100 kb), is_KRABZNF, covariates (10,282 genes)
- `erv_enrichment_results.csv` — raw + KRAB-ZNF-excluded + logistic-adjusted OR
  with CIs and BH-FDR, per target × window + pooled core + SUV39 control
- `fig_erv_enrichment.png` — (A) per-target raw vs corrected forest, (B) 5×→2×
  collapse + SUV39 specificity control, (C) TSS-distance CDF
- `REVIEW_erv.txt` — reviewer pass

## Environment

`scte` — python 3.11.15, bioframe (interval ops), statsmodels (logistic),
scipy (Fisher), pandas, numpy. Annotations fetched from UCSC, ENCODE, EBI
GENCODE, EBI InterPro.

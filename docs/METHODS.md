# Methods

TE-aware analysis of the Marson/Pritchard genome-scale CRISPRi Perturb-seq dataset (Zhu, Dann et al.; primary human CD4⁺ T cells; 4 donors × 3 conditions — Rest / Stim-8h / Stim-48h; probe-based 10x Flex). All analyses ran in Claude Science; artifacts are named and reproducible under `results/`.

## Module and gene sets
Silencing module analyzed: SETDB1, TRIM28/KAP1, ATF7IP, and HUSH members TASOR (FAM208A), MPP8 (MPHOSPH8), PPHLN1. Mechanistic comparators: SUV39H1, SUV39H2 (H3K9 methyltransferases that act off-ERV). Eraser set (druggability): KDM4A/B/C (JMJD2 H3K9me3 demethylases). MPP8 excluded — neither guide passed on-target knockdown QC (`results/mphosph8_exclusion_evidence.csv`), independently corroborated in the authors' data.

## Part 1 — differential expression
Pseudobulk per guide × donor × condition; DESeq2 vs non-targeting controls, `~donor + group`, Wald test, BH FDR, per condition separately. Pre-registered thresholds FDR<0.05 & |log2FC|≥0.5 (frozen before results). On-target knockdown confirmed per target/condition. **Independent DESeq2 vs the authors' own differential expression: Pearson r = 0.94** (`results/de_crosscheck_geneset.csv`, `results/de_crosscheck_genomewide.csv`).

## Part 2 — ERV cis-neighborhood enrichment
Silenced-ERV set = RepeatMasker (hg38) LTR/ERV elements ∩ **primary CD4⁺ T-cell H3K9me3 ChIP-seq** peaks (ENCODE, per condition; no CD4⁺ SETDB1 ChIP exists — H3K9me3 is the field-standard functional mark, and SETDB1-specificity is established on the perturbation side via the SUV39 control). Genes assigned to cis-neighborhoods (TSS ±50/100 kb). Enrichment of silenced-ERV presence near de-repressed genes vs a background matched on gene length, baseline expression, and local gene/TE density — three methods (Fisher, covariate-adjusted logistic, permutation null).

**Confounder control:** the raw enrichment was inflated by KRAB-ZNF genes, which are genomically clustered (chr19), ERV-dense, and direct SETDB1 substrates. Removing them (non-KZFP stratum) gives the honest, confounder-controlled effect **OR ~1.7–2.2** (the raw ~5× is a cluster artifact and is not reported as the effect size). `results/enrichment_SETDB1.csv`.

**Specificity controls** (`results/control1..3c.csv`): ERV/LTR enrichment vs LINE/SINE (ERV-specific); H3K9me3-marked vs unmarked ERVs (modest, condition-dependent increment); core module vs control knockdowns (relative — core median OR 2.37 vs control 1.22, Mann-Whitney p=8×10⁻⁵); **SUV39H1/H2 dissociation** (neither enriches — the mechanistic-specificity evidence, `results/enrichment_SUV39_arms.csv`).

## Functional enrichment (the headline)
On the de-repressed genes: over-representation (Enrichr Hallmark/GO) and **direction-aware pre-ranked GSEA** on the full SETDB1 log2FC signature (Hallmark, per condition). `results/setdb1_nonkzfp_overrepresentation.csv`, `results/setdb1_gsea_hallmark.csv`, `results/fig_setdb1_interferon_gsea.png`. Interferon significance annotated explicitly per condition (IFN-α significant all 3; IFN-γ only Stim-8h).

## KDM4 eraser direct test
KDM4A/B/C were present in the genome-scale library and tested directly (same Part 1 pipeline). Pre-specified prediction: KDM4 (eraser) knockdown produces a signature antagonistic to SETDB1 (writer) at ERV loci. Readouts: genome-wide log2FC correlation vs SETDB1 (benchmarked against HUSH as a same-direction positive reference) and the ERV-proximal directional test on the SETDB1-de-repressed genes. `results/kdm4_*.csv`, `results/fig_kdm4_eraser.png`.

## Druggability
ChEMBL mechanism-of-action + best bioactivity for each node; ranked candidate list combining phenotype strength, ERV-specificity, selectivity/essentiality, and chemical tractability. `results/part3_chembl_tractability.csv`, `results/part3_ranked_candidates.csv`.

## Limitations
Mouse→human transfer of the mechanistic prior; CRISPRi partial knockdown ≠ genetic knockout; H3K9me3 (not SETDB1) ChIP; coarse TSS-window cis-assignment; specificity relative not absolute; viral mimicry from SETDB1 loss is a known phenomenon (this is a clean human-CD4 demonstration, not a discovery); no therapeutic claim.

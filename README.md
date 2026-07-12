# TE-aware target discovery in a T-cell Perturb-seq dataset: the SETDB1/HUSH machinery gates an ERV-encoded interferon regulatory program

**Built with Claude: Life Sciences Hackathon — Research track**
Dataset: Marson & Pritchard genome-scale CRISPRi Perturb-seq in primary human CD4⁺ T cells (Zhu, Dann et al., public).

## One-paragraph summary

The Marson/Pritchard Perturb-seq readout is **probe-based 10x Flex**, which is structurally blind to transposable elements (TEs) / endogenous retroviruses (ERVs) — the exact regulatory layer the H3K9me3 silencing machinery acts on. We built a **TE-aware analysis** of this TE-blind dataset to ask what the ERV-silencing module (SETDB1, TRIM28/KAP1, ATF7IP, and the HUSH complex) does when knocked down, and which nodes are druggable. Genes de-repressed by SETDB1/HUSH knockdown are **modestly but robustly enriched near silenced (H3K9me3-marked) ERVs** (confounder-controlled odds ratio ~1.7–2.2, three methods × three conditions), with a clean **mechanistic-specificity control** (the SUV39H1/2 H3K9 methyltransferases, which act off-ERV, do not enrich). Functionally, the ERV-proximal de-repressed genes are **dominated by an interferon program** — i.e., releasing the H3K9me3 lock de-represses interferon/ISG genes that sit near silenced ERVs, consistent with ERVs acting as **cis-regulatory elements of interferon genes** (the MER41/AIM2 paradigm) rather than as mere pseudo-viral triggers. We then nominate and rank **druggable nodes** in the axis, and directly test the KDM4/JMJD2 "eraser."

## Headline result — SETDB1/HUSH gates an ERV-encoded interferon regulatory program

ERVs are not silenced junk: specific ERV/LTR sequences have been co-opted as **cis-regulatory elements (enhancers/promoters) for innate-immune and interferon genes** — the canonical example is the MER41 LTR family supplying STAT1-binding, interferon-inducible enhancers to genes such as AIM2 and the APOL family (Chuong, Elde & Feschotte, *Science* 2016). In the resting state these ERV-derived regulatory elements are held **OFF** by SETDB1/HUSH-deposited H3K9me3.

In this dataset, knocking down the SETDB1/HUSH machinery **de-represses interferon/ISG genes that are specifically enriched near silenced ERVs** — releasing the H3K9me3 lock licenses an ERV-associated interferon regulatory program *in cis*:

- The de-repressed genes near silenced ERVs are **massively interferon-enriched** (over-representation: Interferon-γ Response FDR 3×10⁻¹³, Interferon-α Response FDR 3.6×10⁻¹²), including known ERV-enhancer-regulated genes (**APOL family**) and checkpoint-relevant **CD274/PD-L1**.
- **Direction-aware GSEA** on the full SETDB1 signature confirms a broader interferon program (Interferon-α top positive-NES in all 3 conditions: Rest +1.66, Stim-8h +1.80, Stim-48h +1.74, q<0.03; Interferon-γ positive but stimulation-dependent, significant only at Stim-8h; leading edge = the classic ISG core IFITM1/2/3, OAS1, OASL, RSAD2, BST2, STAT2, IRF9) — a *trans* interferon layer on top of the cis de-repression.

**Read precisely: ERVs constitute a SETDB1/HUSH-gated cis-regulatory layer of the interferon program in human CD4⁺ T cells** — a repressive/inducible regulatory axis (MER41/AIM2-type), gated by H3K9me3. The spatial linkage (interferon genes enriched *near* silenced ERVs) is the signature of the cis-regulatory-element mechanism, not pure trans sensing. Whether the broader signature is additionally amplified by classic "viral mimicry" (trans sensing of ERV transcripts) **cannot be separated in this probe-based, TE-blind assay and is not claimed**.

Figure: `results/fig_setdb1_interferon_gsea.png`. Data: `results/setdb1_nonkzfp_overrepresentation.csv`, `results/setdb1_gsea_hallmark.csv`.

## The honest analytical arc (the journey is the finding)

1. **Hypothesis** (from Adoue et al., *Immunity* 2019, mouse): SETDB1/ERV silencing tunes Th1/Th2 balance.
2. **Confirmed** a modest, **confounder-controlled** ERV-cis-enrichment (OR ~1.7–2.2; the raw ~5× was a KRAB-ZNF genomic-clustering artifact and is **not** the effect size), with a clean **SUV39 dissociation** establishing mechanistic specificity.
3. **But the Th1/Th2 module scored null.** Rather than force the hypothesis, we ran an unbiased functional enrichment of the de-repressed genes.
4. **The real, robust immune consequence is the ERV-gated interferon program above — not a Th1/Th2 lineage shift.** This explains the null Th1/Th2: the knockdown releases the ERV-silencing layer that gates interferon-gene regulatory elements, not the T-helper lineage machinery.

## Drug targets (the deliverable)

Ranked druggable nodes in the axis (`results/part3_ranked_candidates.csv`, `results/fig_candidate_ranking.png`), with chemical tractability from ChEMBL (`results/part3_chembl_tractability.csv`):

- **KDM4C** (JMJD2 H3K9me3 "eraser") — the lead **directly-measured** node. Its knockdown produces an ERV-proximal transcriptional signature **antagonistic** to SETDB1 (writer↔eraser), confirmed by direct perturbation (`results/kdm4_ervprox_directional.csv`, `results/fig_kdm4_eraser.png`); selective sub-100 nM chemistry exists. Effect is a modest, locus-specific rheostat (does not reverse phenotype).
- **KDM4A** — shares the eraser mechanism (directly measured) but weak chemistry; supporting evidence, not the lead.
- **KDM4B** — distinct: no ERV antagonism here; rationale rests on a separate (published) cGAS-axis mechanism, flagged as **not** ERV-validated by this analysis.
- **SETDB1 / HUSH (PPHLN1, TASOR)** — strongest mechanistic hits but pan-essential / undrugged; the writer-side of the axis.

Two directions of one axis: **writer** (SETDB1/HUSH) inhibition → viral mimicry (pro-immune); **eraser** (KDM4) inhibition → reinforce silencing → dampen interferon. Which direction is therapeutic is disease-dependent and **not demonstrated here** — these are motivated hypotheses, not claims.

## Rigor and reproducibility

- **Independent DESeq2 reproduces the authors' own differential expression** (Pearson r = 0.94) — pipeline-correctness anchor (`results/de_crosscheck_*.csv`, `results/fig_crosscheck_vs_authors.png`).
- **Pre-registered thresholds**, confounder control (KRAB-ZNF clustering: raw ~5× → honest ~2×), and specificity controls (ERV vs LINE/SINE; H3K9me3-marked vs unmarked; core vs control knockdowns; SUV39 dissociation). See `docs/METHODS.md`; all numeric claims trace to named artifacts in `results/`.

## How Claude Science contributed

- Autonomous multi-stage pipeline (remote S3 data slicing on a laptop; pseudobulk DESeq2; ERV-cis enrichment; ORA + direction-aware GSEA; ChEMBL tractability) with named, reproducible artifacts and provenance.
- **Plan-mode checkpoints** before each genome-scale run; **Reviewer** agent caught a figure-significance-annotation error (interferon-γ marked significant where it wasn't) and confounder issues, forcing corrections.
- MCP/skill-driven orchestration of ChEMBL and Enrichr/GSEA.

## Honest boundaries

- Two distinct ERV→interferon mechanisms exist: (i) **ERVs as cis-regulatory elements** of interferon genes (MER41/AIM2-type co-option; Chuong 2016) and (ii) **"viral mimicry"** — trans sensing of ERV transcripts as non-self nucleic acid. Both are **published**. Our cis-proximity result (interferon genes enriched near silenced ERVs) specifically supports mechanism (i); the probe-based Flex assay is TE-blind, so ERV transcripts are **not measured** and the sensing step of mechanism (ii) is **not demonstrated** here.
- The contribution is the **clean demonstration in primary human CD4⁺ T cells via genome-scale Perturb-seq** that SETDB1/HUSH loss de-represses an ERV-proximal interferon program, unifying the cis-ERV enrichment and the interferon signature — **not** discovery of the ERV→interferon relationship.
- Interferon-α is the robust all-conditions signal; interferon-γ is stimulation-dependent. ERV-cis-enrichment is modest (OR ~2), and is partly ERV family/location, not solely the H3K9me3 mark. Cis-proximity is not proof of enhancer function.
- **No therapeutic claim is made.** Cross-species transfer (mouse→human), CRISPRi partial knockdown vs. genetic knockout, and coarse cis-window assignment are limitations.

## License

Apache 2.0 — see `LICENSE`.

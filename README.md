# TE-aware target discovery in a T-cell Perturb-seq dataset: the SETDB1/ERV silencing machinery drives a viral-mimicry interferon program

**Built with Claude: Life Sciences Hackathon — Research track**
Dataset: Marson & Pritchard genome-scale CRISPRi Perturb-seq in primary human CD4⁺ T cells (Zhu, Dann et al., public).

## One-paragraph summary

The Marson/Pritchard Perturb-seq readout is **probe-based 10x Flex**, which is structurally blind to transposable elements (TEs) / endogenous retroviruses (ERVs) — the exact regulatory layer the H3K9me3 silencing machinery acts on. We built a **TE-aware analysis** of this TE-blind dataset to ask what the ERV-silencing module (SETDB1, TRIM28/KAP1, ATF7IP, and the HUSH complex) does when knocked down, and which nodes are druggable. Genes de-repressed by SETDB1/HUSH knockdown are **modestly but robustly enriched near silenced (H3K9me3-marked) ERVs** (confounder-controlled odds ratio ~1.7–2.2, three methods × three conditions), with a clean **mechanistic-specificity control** (the SUV39H1/2 H3K9 methyltransferases, which act off-ERV, do not enrich). Functionally, the dominant, direction-correct consequence is the **canonical viral-mimicry interferon program**: interferon-α response is the top positive-enrichment hallmark in all three conditions, the leading edge is the classic ISG core (IFITM1/2/3, OAS1, OASL, RSAD2, BST2, STAT2, IRF9), and **CD274/PD-L1** sits in the interferon-γ leading edge. We then nominate and rank **druggable nodes** in the axis, and directly test the KDM4/JMJD2 "eraser."

## Headline result — SETDB1/HUSH knockdown → viral-mimicry interferon program

De-repressing the ERV-silencing machinery in primary human CD4⁺ T cells triggers the textbook **ERV-de-silencing → viral-mimicry → interferon** response — the same axis exploited by DNMT/HDAC-inhibitor "viral mimicry" in cancer immunotherapy — shown cleanly here in T cells:

- **Over-representation** (de-repressed non-KZFP genes): Interferon-γ Response FDR 3×10⁻¹³, Interferon-α Response FDR 3.6×10⁻¹² — dominant by orders of magnitude.
- **Direction-aware GSEA** (full SETDB1 signature): Interferon-α is the **top positive-NES hallmark in all 3 conditions** (Rest +1.66, Stim-8h +1.80, Stim-48h +1.74; q<0.03). Interferon-γ is positive but **stimulation-dependent** (significant only at Stim-8h).
- Leading edge = classic ISG core; **CD274/PD-L1** in the IFN-γ leading edge (a checkpoint-relevant readout).

Figure: `results/fig_setdb1_interferon_gsea.png`. Data: `results/setdb1_nonkzfp_overrepresentation.csv`, `results/setdb1_gsea_hallmark.csv`.

## The honest analytical arc (the journey is the finding)

1. **Hypothesis** (from Adoue et al., *Immunity* 2019, mouse): SETDB1/ERV silencing tunes Th1/Th2 balance.
2. **Confirmed** a modest, **confounder-controlled** ERV-cis-enrichment (OR ~1.7–2.2; the raw ~5× was a KRAB-ZNF genomic-clustering artifact and is **not** the effect size), with a clean **SUV39 dissociation** establishing mechanistic specificity.
3. **But the Th1/Th2 module scored null.** Rather than force the hypothesis, we ran an unbiased functional enrichment of the de-repressed genes.
4. **The real, robust immune consequence is the viral-mimicry interferon program above — not a Th1/Th2 lineage shift.** This is a *bystander innate-immune / trans* program (ISG induction downstream of ERV sensing), and it explains the null Th1/Th2: the knockdown hits the ERV-silencing / interferon layer, not the T-helper lineage machinery.

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

- Viral mimicry from SETDB1 loss is an **established** phenomenon (cancer immunotherapy). The contribution here is the **clean demonstration in primary human CD4⁺ T cells via genome-scale Perturb-seq**, the integrated cis-ERV / trans-ISG decomposition, and the PD-L1 T-cell readout — **not** discovery of viral mimicry.
- Interferon-α is the robust all-conditions signal; interferon-γ is stimulation-dependent.
- ERV-cis-enrichment is modest (OR ~2), and the enrichment is partly ERV family/location, not solely the H3K9me3 mark.
- **No therapeutic claim is made.** Cross-species transfer (mouse→human), CRISPRi partial knockdown vs. genetic knockout, and coarse cis-window assignment are limitations.

## License

Apache 2.0 — see `LICENSE`.

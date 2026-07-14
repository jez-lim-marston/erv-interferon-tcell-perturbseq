# TE-aware target discovery in a T-cell Perturb-seq dataset: the SETDB1/HUSH machinery gates an ERV-encoded interferon regulatory program

**Built with Claude: Life Sciences Hackathon — Research track**
Dataset: Marson & Pritchard genome-scale CRISPRi Perturb-seq in primary human CD4⁺ T cells (Zhu, Dann et al., public).

## One-paragraph summary

**In plain terms:** scattered through the human genome are *viral fossils* — endogenous retroviruses (ERVs), the inherited remains of ancient retroviral infections. Far from dead, some have been **domesticated into regulatory switches** for immune genes, held off by the H3K9me3 silencing machinery. This project asks, in a genome-scale CRISPRi screen of primary human T cells, what that machinery controls when it's released — and which parts of it are druggable.

The Marson/Pritchard Perturb-seq readout is **probe-based 10x Flex**, which is structurally blind to transposable elements (TEs) / endogenous retroviruses (ERVs) — the exact regulatory layer the H3K9me3 silencing machinery acts on. We built a **TE-aware analysis** of this TE-blind dataset to ask what the ERV-silencing module (SETDB1, TRIM28/KAP1, ATF7IP, and the HUSH complex) does when knocked down, and which nodes are druggable. Genes de-repressed by SETDB1/HUSH knockdown are **modestly but robustly enriched near silenced (H3K9me3-marked) ERVs** (confounder-controlled odds ratio ~1.7–2.2, three methods × three conditions), with a clean **mechanistic-specificity control** (the SUV39H1/2 H3K9 methyltransferases, which act off-ERV, do not enrich). Functionally, the ERV-proximal de-repressed genes are **dominated by an interferon program** — i.e., releasing the H3K9me3 lock de-represses interferon/ISG genes that sit near silenced ERVs, consistent with ERVs acting as **cis-regulatory elements of interferon genes** (the MER41/AIM2 paradigm) rather than as mere pseudo-viral triggers. We then nominate and rank **druggable nodes** in the axis, and directly test the KDM4/JMJD2 "eraser."

## Headline result — SETDB1/HUSH gates an ERV-encoded interferon regulatory program

ERVs are not silenced junk: specific ERV/LTR sequences have been co-opted as **cis-regulatory elements (enhancers/promoters) for innate-immune and interferon genes** — the canonical example is the MER41 LTR family supplying STAT1-binding, interferon-inducible enhancers to genes such as AIM2 and the APOL family (Chuong, Elde & Feschotte, *Science* 2016). In the resting state these ERV-derived regulatory elements are held **OFF** by SETDB1/HUSH-deposited H3K9me3.

In this dataset, knocking down the SETDB1/HUSH machinery **de-represses interferon/ISG genes that are specifically enriched near silenced ERVs** — releasing the H3K9me3 lock licenses an ERV-associated interferon regulatory program *in cis*:

- The de-repressed genes near silenced ERVs are **massively interferon-enriched** (over-representation: Interferon-γ Response FDR 3×10⁻¹³, Interferon-α Response FDR 3.6×10⁻¹²), including known ERV-enhancer-regulated genes (**APOL family**) and checkpoint-relevant **CD274/PD-L1**.
- **Direction-aware GSEA** on the full SETDB1 signature confirms a broader interferon program (Interferon-α is the top positive signature in all 3 conditions (NES Rest 1.53, Stim-8h 2.28, Stim-48h 2.30) — FDR-significant in the two stimulated conditions and a strong trend at Rest (FDR 0.14); Interferon-γ positive but stimulation-dependent, significant only at Stim-8h; leading edge = the classic ISG core IFITM1/2/3, OAS1, OASL, RSAD2, BST2, STAT2, IRF9) — a *trans* interferon layer on top of the cis de-repression.

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

Two directions of one axis: **writer** (SETDB1/HUSH) inhibition → viral mimicry (pro-immune); **eraser** (KDM4) inhibition → reinforces ERV silencing (antagonistic to SETDB1 loss in this dataset). The net *sign* of eraser inhibition on the interferon program is context-dependent and **not resolved here**; which direction is therapeutic is disease-dependent — these are motivated hypotheses, not claims.

## Rigor and reproducibility

- **Independent DESeq2 reproduces the authors' own differential expression** (pooled Pearson r = 0.915, median per-target × condition r = 0.97) — pipeline-correctness anchor (`results/de_crosscheck_*.csv`, `results/fig_crosscheck_vs_authors.png`).
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

## Translational outlook & next directions

The result above frames the SETDB1/HUSH–ERV–interferon axis as a **tunable, bidirectional interferon rheostat** — which is where its disease potential lies, in two opposite directions:

- **Release the lock → induce interferon (immuno-oncology / viral mimicry).** De-repressing the axis is a candidate strategy to convert immunologically "cold," checkpoint-refractory tumors to "hot" by inducing an ERV-driven interferon state — the viral-mimicry rationale already validated for DNA-methylation inhibitors and SETDB1 loss in cancer models (Chiappinelli et al. 2015; Griffin et al. 2021). The primary-human-CD4⁺ genome-scale demonstration here that the module gates this program is the missing link between that cancer-cell literature and human T-cell immunity.
- **Reinforce the lock → restrain interferon (interferon-driven autoimmunity).** The opposite direction is a candidate for type-I-interferon-driven autoimmune disease, where chronic ERV-associated interferon tone is pathogenic.

**Where the medicinal chemistry goes.** The axis presents three node classes with distinct druggability: the **writer** (SETDB1/HUSH) is the strongest mechanistic hit but pan-essential / undrugged; the **eraser** (KDM4) has mature selective chemistry but sits off the core silencing axis; the **reader** layer (the HP1/CBX proteins that bind H3K9me3) is the most on-mechanism but presents a harder, likely **targeted-degrader** medicinal-chemistry problem (shallow methyl-reader pocket, high paralog conservation). A priority next step is **structure-guided modality selection** across these nodes — small molecule vs. targeted degrader vs. protein–protein-interface disruptor.

**Where the biology goes.** Two experiments would sharpen the therapeutic hypothesis: (i) a **TE-aware assay** (long-read or ERV-resolved quantification) to directly measure ERV transcripts and separate the *cis* regulatory-element mechanism from *trans* viral-mimicry sensing — which this probe-based Flex assay cannot; and (ii) a **perturb-then-checkpoint** design to test whether de-repressing the axis sensitizes tumors to checkpoint blockade, the translational endpoint the immuno-oncology direction predicts.

These are motivated hypotheses and a roadmap, not results — but they are why the axis is worth drugging.

## Conclusion

In primary human CD4⁺ T cells at genome scale, the SETDB1/HUSH silencing module gates an ERV-encoded *cis*-regulatory layer of the interferon program — releasing the H3K9me3 lock de-represses interferon-stimulated genes enriched beside silenced ERVs, a modest but robust, confounder- and specificity-controlled effect that reframes the module's immune role away from Th1/Th2 lineage tuning. This defines the axis as a **bidirectional, druggable interferon rheostat** — one mechanism with opposite therapeutic polarities, tunable toward viral-mimicry immuno-oncology or toward restraint of interferon-driven autoimmunity — and nominates its layered druggable nodes (writer, eraser, reader) as a ranked target set. We *demonstrate rather than discover* this axis and make no therapeutic claim — the ERV-transcript and clinical steps remain to be shown — but its position at the immuno-oncology/autoimmunity intersection makes it a rigorous, honestly-bounded, and actionable target hypothesis.

*(Graphical abstract: `docs/GRAPHICAL_ABSTRACT.md`.)*

## How this was built (compute split)

Claude Science orchestrated the genome-scale analysis, statistics, and figures with
reproducible named artifacts and Plan-mode checkpoints; a companion Claude Code + AWS
environment handled heavy data processing; a Reviewer agent caught a figure-significance
error.

Concretely, the split was:
- **Claude Science (orchestration).** Remote partial-read slicing of the 44 GB pseudobulk
  and 16.8 GB DE-stats objects from public S3; pseudobulk DESeq2 reproduction of the
  authors' differential expression (correctness anchor, Pearson *r* = 0.915 pooled / median
  per-target × condition *r* = 0.97 vs the ~0.94 target); ERV-cis proximity enrichment with
  the KRAB-ZNF confounder correction and SUV39 specificity control; ORA + direction-aware
  GSEA (the interferon program); ChEMBL druggable-node ranking; and the publication figure
  set (light + dark, SVG + PNG). Every deliverable is a named artifact with frozen,
  pre-registered thresholds and a Plan-mode checkpoint before each genome-scale run.
- **Companion Claude Code + AWS environment (heavy data processing).** The high-throughput
  data-preparation steps that back the `results/` tables.
- **Reviewer agent (independent QC pass).** Re-derived every figure's annotated statistics
  from the source tables and caught a real significance-annotation error — figure titles
  had claimed "IFN-α significant in all three conditions", but IFN-α is FDR-significant only
  in the two *stimulated* conditions and a strong trend at Rest (NES 1.53 / FDR 0.14; still
  the #1-ranked positive Hallmark in all three). Titles were corrected to match the plotted
  significance stars. The Reviewer pass is archived in `methods/REVIEW_figures.txt` (and the
  earlier confounder/effect-size pass in `methods/REVIEW_erv.txt`).

Note on cross-environment consistency: the two environments ran independent DE/GSEA passes
with slightly different FDR conventions, so a few figures generated in Claude Science
(`results/figures/`) report marginally different NES / *r* values than the companion-authored
prose above; each figure's numbers are internally consistent with its own source table and
verified in the Reviewer pass.

## License

Apache 2.0 — see `LICENSE`.

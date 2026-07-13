# 3-minute demo video — script / storyboard

Target 3:00. Times are cumulative. Visual = the figure/screen to show.

---

**(a) The question — 0:00–0:30**
Visual: title slide → the SETDB1/HUSH + ERV/H3K9me3 cartoon.

"Scattered through our DNA are **viral fossils** — the burned-out remains of viruses that
infected our ancestors millions of years ago. But they're not dead: the genome has 
domesticated some of them into switches that control the immune system's interferon 
alarm, and the SETDB1/HUSH machinery keeps those switches locked off with a molecular
seal called H3K9me3. So we asked, in the Marson–Pritchard genome-scale CRISPRi 
Perturb-seq of primary human CD4 T cells: what happens when you break that seal — and is
any part of the machinery druggable?"

**(b) The unique angle + rigor — 0:30–1:00**
Visual: 'probe-based Flex = TE-blind' graphic → the author-cross-check figure (r=0.94) → 
the confounder 5×→2× note.

"The catch: this dataset is read on a probe panel that's blind to transposable elements —
the exact layer we care about. So we built a TE-aware analysis of a TE-blind dataset. 
We pre-registered thresholds, reproduced the authors' own differential expression at 
r=0.94, and — importantly — controlled a confounder: a raw five-fold ERV enrichment 
collapsed to an honest two-fold once we removed clustered KRAB-zinc-finger genes."

**(c) The pivot + headline — 1:00–2:10 (the core, most time)**
Visual: SETDB1 ERV-enrichment + SUV39 dissociation figure → then the interferon GSEA/ORA 
figure (IFN-α starred all 3) → highlight APOL / CD274-PD-L1.

"We started from a mouse hypothesis that this axis tunes Th1 versus Th2 balance. The ERV 
enrichment held — modest but robust, with a clean specificity control: the off-ERV SUV39 
enzymes don't enrich. But the Th1/Th2 signature came back null. So instead of forcing it,
we asked, unbiased: what ARE these de-repressed genes?

The answer is an interferon program — and here's the precise part: the interferon genes 
are enriched right next to the silenced ERVs. That's the signature of ERVs acting as 
cis-regulatory elements of interferon genes — the MER41–AIM2 mechanism — gated by 
SETDB1/HUSH. Interferon-alpha is the top hallmark in every condition; the leading edge is 
the classic ISG core; and it includes checkpoint-relevant PD-L1 and known 
ERV-enhancer-controlled genes. So this isn't Th1/Th2 tuning — it's an ERV-interferon 
regulatory axis, and that reframing explains why the lineage signature was null."

**(d) Drug targets, honestly — 2:10–2:45**
Visual: ranked-candidate figure; SETDB1-writer / KDM4-eraser two-arrow schematic.

"For targets: the writer, SETDB1/HUSH, is the strongest hit but pan-essential. The 
tractable node is the eraser — KDM4C — which we tested directly: its knockdown antagonizes
SETDB1 at exactly these ERV-proximal genes, and it has selective sub-100-nanomolar 
chemistry. Writer versus eraser are two directions of one axis. Which direction is 
therapeutic is disease-dependent — that's a hypothesis this work makes testable, not a 
claim we're making.

*[Outlook teaser — first to cut if over 3:00]* And this axis is a two-way interferon dial 
with real disease stakes: release it to warm up cold, checkpoint-refractory tumors, or 
reinforce it to cool interferon-driven autoimmunity — which is why the next chapter is 
drugging it across its writer, eraser, and reader nodes."

**(e) How Claude Science ran it — 2:45–3:00**
Visual: Claude Science pipeline / artifact tree; Reviewer 'finding' chip.

"Claude Science ran the whole pipeline autonomously — slicing a 45-gigabyte matrix on a 
laptop, DESeq2, enrichment, druggability — with reproducible artifacts and Plan-mode 
checkpoints. Its Reviewer even caught us marking an interferon bar significant when it 
wasn't. The result: a rigorous, honest, TE-aware framework that followed the data to a 
better answer than its hypothesis."

---

**Cut priority if over time:** first drop the (d) outlook-teaser sentence; then trim (b) 
to one sentence (keep the r=0.94 and 5×→2×); then trim (d)'s writer/eraser detail. Never 
cut the (c) pivot or the honesty caveats.

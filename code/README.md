# Analysis code

As-run analysis code for the TE-aware reanalysis of the Marson/Pritchard genome-scale CRISPRi
Perturb-seq dataset (primary human CD4+ T cells). This is the **orchestration-side** code that
ran in Claude Science and produced the named artifacts under `../results/`.

## What this is (and isn't)

- **Faithful, not refactored.** `pipeline/*.py` and `as_run_analysis_log.py` are exported
  verbatim from the execution log, cell by cell, in execution order. They record exactly what
  was run — including interactive iteration — rather than a cleaned-up importable package.
- **What produced this code.** The Claude Science session ran the entire orchestration-side
  analysis itself — remote partial-read S3 slicing, pseudobulk DESeq2 reproduction, ERV-cis
  enrichment, GSEA/ORA, ChEMBL ranking, and the figure set — with reproducible named artifacts
  and Plan-mode checkpoints, and a Reviewer pass that caught a figure-significance error. The
  pre-existing `../results/` tables and `../docs/` narrative in this repo predate this code and
  were prepared separately; their generating scripts are not part of this export.
- **Credentials scrubbed.** All cells that touched credentials/tokens or git authentication,
  and the code-export machinery itself, are excluded. No secrets are present.

## Layout

- `pipeline/` — analysis grouped by stage, in order:
  - `00_figure_style_helper.py` — light/dark SVG+PNG `save_fig` helper
  - `01_subset_pseudobulk.py` — partial-range S3 read; module+SUV39+NTC pseudobulk subset
  - `02_authors_anchor.py` — pull the authors' published DE (correctness anchor)
  - `03_run_de_deseq2.py` — pseudobulk DESeq2 (pertpy/pydeseq2), sandbox-parallelism fix
  - `04_de_crosscheck_pearson.py` — reproduce authors' DE, pooled Pearson r + bootstrap CI
  - `05_de_summary_counts.py` — per-target x condition DE summary table
  - `06_erv_proximity_annotation.py` — RepeatMasker LTRs INT H3K9me3 consensus; TSS distances; KRAB-ZNF flag
  - `07_erv_enrichment.py` — Fisher OR (raw + KRAB-ZNF-excluded) + logistic adjustment; SUV39 control
  - `08_gsea.py` — direction-aware pre-ranked GSEA (Hallmark, per condition)
  - `09_ora.py` — hypergeometric over-representation (contrast to GSEA)
  - `10_druggable_nodes.py` — ChEMBL tractability + KDM4/writer/reader ranking
  - `11..18_fig_*.py` — publication figures (light+dark, SVG+PNG)
- `as_run_analysis_log.py` — complete chronological log of every analysis cell
- `run_de_module.py` — standalone DESeq2 driver (BLAS pinned to 1 thread; threading joblib
  backend; numba disk-cache disabled — the fixes that let pydeseq2 run in the sandbox)
- `MultiStatePerturbSeqDataset.py` — the authors' dataset wrapper (from their public repo),
  used unchanged for methodological identity with their DE

## Environment

conda env `scte`: python 3.11.15; scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, pertpy 1.0.3,
formulaic-contrasts 1.0.0, gseapy 1.3.0, bioframe, boto3, s3fs, statsmodels.

Public data: `s3://genome-scale-tcell-perturb-seq/marson2025_data/` (AWS no-sign-request; read
by partial-range requests, never fully downloaded). GEO GSE314342 / SRA SRP643211.

Licensed under Apache-2.0 (see `../LICENSE`).

# Meningioma molecular aggressiveness and recurrence program

This repository contains the analysis code and lightweight release materials for a retrospective cross-cohort meningioma transcriptomics study.

## Manuscript

Working title:

> Cross-cohort transcriptomics defines a grade-associated meningioma aggressiveness program enriched in recurrent tumours

Primary target: Scientific Reports.  
Fallback target: BMC Medical Genomics.

## Analysis Overview

The study defines a grade-associated molecular aggressiveness program from GSE136661, validates recurrence and grade associations in independent bulk cohorts, and resolves the program in single-cell meningioma data. Geneformer zero-shot embeddings are retained as a secondary, analysis-dependent benchmark against classical expression features.

Main public datasets:

- GSE136661
- GSE16581
- GSE74385
- GSE183655
- GSE206647

Raw GEO downloads and large model/intermediate files are intentionally not included in the release repository.

## Reproducibility

Run commands from the repository root. The canonical run order is documented in:

- `docs/RUN_ORDER.md`

Key entry points:

- `scripts/05_aggressiveness_program.R` — bulk program definition
- `scripts/06_crosscohort_validation.R` — GSE16581 validation
- `scripts/16_GSE74385_recurrence.R` — GSE74385 recurrence validation
- `scripts/24_P0_robustness.R` — grade-adjusted Firth and batch checks
- `scripts/15_GSE206647_grade.R` and `scripts/31_pseudobulk_cellcycle.R` — single-cell grade validation and pseudobulk controls
- `scripts/17a_export_tumor_for_geneformer.R`, `scripts/17b_geneformer_embed.py`, `scripts/51_repair_geneformer_patient_mapping.R`, `scripts/52_geneformer_patient_level_audit.R` — secondary Geneformer benchmark and patient-level sensitivity analysis
- `scripts/40_independent_headline_audit.R` — independent raw/assembled-input headline audit
- `scripts/42_finalize_corrected_target_matrix.R` — apply the independently verified full-library pseudobulk correction

## Environment

Environment snapshots are stored in:

- `docs/sessionInfo_R.txt`
- `docs/geneformer_requirements.txt`

The completed Geneformer analysis used local model assets and CPU inference. The dependency snapshot contains the pinned Python stack; install Geneformer separately from its official model repository as described in `docs/geneformer_requirements.txt`. The release repository does not include model weights.

## Release Scope

Include:

- analysis scripts
- manuscript-supporting markdown documents
- supplementary CSV tables
- publication-ready figures
- environment manifests
- repository metadata files (`CITATION.cff`, `.zenodo.json`, `LICENSE`)

Exclude:

- raw GEO downloads
- large `.rds`, `.h5ad`, `.mtx`, model-weight, and local virtual-environment files
- rendered manuscript QA caches
- personal machine logs

## Citation

The Zenodo concept DOI should be inserted after the first public archive is created. ORCID records and the DOI remain portal/release items; do not reuse the DOI from another project.
